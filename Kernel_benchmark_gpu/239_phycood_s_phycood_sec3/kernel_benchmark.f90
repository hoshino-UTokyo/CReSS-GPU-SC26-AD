program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: ni, nj, nk
  integer :: nkm1, nkm2
  real :: zsfc, zflat0, htuiv, htuivz

  ! Arrays
  real, allocatable :: ht(:,:)
  real, allocatable :: zsth(:)
  real, allocatable :: zsth_in(:)
  real, allocatable :: zph(:,:,:)
  real, allocatable :: zph_ref(:,:,:)
  real, allocatable :: zsth_ref(:)

  ! Benchmark variables
  integer :: num_iterations, warmup_iterations
  real(8) :: total_time, avg_time
  real(8) :: start_time, end_time
  character(len=256) :: data_dir
  real :: tolerance
  integer :: error_count
  real :: max_error
  integer :: iter, i, j, k

  ! Read benchmark configuration
  open(10, file='benchmark.conf', status='old')
  read(10, '(A)') data_dir
  read(10, *) num_iterations
  read(10, *) warmup_iterations
  read(10, *) tolerance
  close(10)

  ! Read parameters
  call read_params(trim(data_dir)//'/params.txt')

  ! Allocate arrays
  allocate(ht(0:ni+1, 0:nj+1))
  allocate(zsth(1:nk), zsth_in(1:nk))
  allocate(zph(0:ni+1, 0:nj+1, 1:nk))
  allocate(zph_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(zsth_ref(1:nk))

  ! Read input data
  call read_binary_2d(trim(data_dir)//'/ht.bin', ht, 0, ni+1, 0, nj+1)
  call read_binary_1d(trim(data_dir)//'/zsth_in.bin', zsth_in, 1, nk)

  ! Read reference output
  call read_binary_3d(trim(data_dir)//'/zph_ref.bin', zph_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_binary_1d(trim(data_dir)//'/zsth_ref.bin', zsth_ref, 1, nk)

  ! Warmup iterations
  do iter = 1, warmup_iterations
    zph = 0.0
    zsth = zsth_in

    call kernel_phycood(ni, nj, nk, nkm1, nkm2, zsfc, zflat0, htuiv, htuivz, &
                        ht, zsth, zph)
    !$acc wait
  end do

  ! Benchmark iterations
  total_time = 0.0d0
  do iter = 1, num_iterations
    zph = 0.0
    zsth = zsth_in

    !$acc wait
    start_time = omp_get_wtime()

    call kernel_phycood(ni, nj, nk, nkm1, nkm2, zsfc, zflat0, htuiv, htuivz, &
                        ht, zsth, zph)

    !$acc wait
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate zph
  error_count = 0
  max_error = 0.0
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        if (abs(zph_ref(i,j,k)) > 1.0e-30) then
          if (abs(zph(i,j,k) - zph_ref(i,j,k)) / abs(zph_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(zph(i,j,k) - zph_ref(i,j,k)) / abs(zph_ref(i,j,k)))
          end if
        else
          if (abs(zph(i,j,k) - zph_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(zph(i,j,k) - zph_ref(i,j,k)))
          end if
        end if
      end do
    end do
  end do

  ! Validate zsth
  do k = 1, nk
    if (abs(zsth_ref(k)) > 1.0e-30) then
      if (abs(zsth(k) - zsth_ref(k)) / abs(zsth_ref(k)) > tolerance) then
        error_count = error_count + 1
        max_error = max(max_error, abs(zsth(k) - zsth_ref(k)) / abs(zsth_ref(k)))
      end if
    else
      if (abs(zsth(k) - zsth_ref(k)) > tolerance) then
        error_count = error_count + 1
        max_error = max(max_error, abs(zsth(k) - zsth_ref(k)))
      end if
    end if
  end do

  ! Output results
  print '(A)', '=== Benchmark Results ==='
  print '(A,I0)', 'Kernel: phycood_s_phycood_sec3'
  print '(A,I0)', 'Iterations: ', num_iterations
  print '(A,F12.6,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,I0)', 'Errors: ', error_count
  print '(A,E12.5)', 'Max relative error: ', max_error
  if (error_count == 0) then
    print '(A)', 'PASSED'
  else
    print '(A)', 'FAILED'
  end if

  ! Cleanup
  deallocate(ht, zsth, zsth_in, zph, zph_ref, zsth_ref)

contains

  subroutine kernel_phycood(ni, nj, nk, nkm1, nkm2, zsfc, zflat0, htuiv, htuivz, &
                            ht, zsth, zph)
    implicit none
    integer, intent(in) :: ni, nj, nk, nkm1, nkm2
    real, intent(in) :: zsfc, zflat0, htuiv, htuivz
    real, intent(in) :: ht(0:ni+1, 0:nj+1)
    real, intent(inout) :: zsth(1:nk)
    real, intent(inout) :: zph(0:ni+1, 0:nj+1, 1:nk)
    integer :: i, j, k

    ! First kernel: compute zph for k=2 to nk-1
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      !$acc loop independent
      do j = 0, nj
        !$acc loop independent
        do i = 0, ni
          if (zsth(k) > zflat0) then
            zph(i,j,k) = zsth(k)
          else
            zph(i,j,k) = htuiv * (zflat0 - ht(i,j)) * (zsth(k) - zsfc) + ht(i,j)
          end if
        end do
      end do
    end do
    !$acc end kernels

    ! Second kernel: set boundary conditions for zph
    !$acc kernels
    !$acc loop independent
    do j = 0, nj
      !$acc loop independent
      do i = 0, ni
        zph(i,j,2) = ht(i,j)
        zph(i,j,1) = 2.0e0 * zph(i,j,2) - zph(i,j,3)
        zph(i,j,nk) = 2.0e0 * zph(i,j,nkm1) - zph(i,j,nkm2)
      end do
    end do
    !$acc end kernels

    ! Third kernel: update zsth for k=2 to nk-1
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      if (zsth(k) <= zflat0) then
        zsth(k) = htuivz * (zsth(k) - zsfc)
      end if
    end do
    !$acc end kernels

    ! Sequential updates for zsth boundaries (outside acc region)
    zsth(2) = 0.0e0
    zsth(1) = -zsth(3)
    zsth(nk) = 2.0e0 * zsth(nk-1) - zsth(nk-2)

  end subroutine kernel_phycood

  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    character(len=256) :: line
    character(len=64) :: key
    integer :: ios, eq_pos

    open(10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening params file: ', trim(filename)
      stop 1
    end if

    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        select case (trim(key))
        case ('ni')
          read(line(eq_pos+1:), *) ni
        case ('nj')
          read(line(eq_pos+1:), *) nj
        case ('nk')
          read(line(eq_pos+1:), *) nk
        case ('nkm1')
          read(line(eq_pos+1:), *) nkm1
        case ('nkm2')
          read(line(eq_pos+1:), *) nkm2
        case ('zsfc')
          read(line(eq_pos+1:), *) zsfc
        case ('zflat0')
          read(line(eq_pos+1:), *) zflat0
        case ('htuiv')
          read(line(eq_pos+1:), *) htuiv
        case ('htuivz')
          read(line(eq_pos+1:), *) htuivz
        end select
      end if
    end do
    close(10)
  end subroutine read_params

  subroutine read_binary_1d(filename, array, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: k1, k2
    real, intent(out) :: array(k1:k2)
    integer :: ios

    open(10, file=filename, form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(10) array
    close(10)
  end subroutine read_binary_1d

  subroutine read_binary_2d(filename, array, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: array(i1:i2, j1:j2)
    integer :: ios

    open(10, file=filename, form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(10) array
    close(10)
  end subroutine read_binary_2d

  subroutine read_binary_3d(filename, array, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: array(i1:i2, j1:j2, k1:k2)
    integer :: ios

    open(10, file=filename, form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(10) array
    close(10)
  end subroutine read_binary_3d

end program kernel_benchmark

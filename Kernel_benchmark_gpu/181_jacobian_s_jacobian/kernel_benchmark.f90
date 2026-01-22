!> GPU Kernel Benchmark: jacobian (s_jacobian)
!> Compute Jacobian transformation coefficients j31, j32, jcb
!> GPU Difficulty: Medium - Multiple 3D parallel loops
program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: ni, nj, nk

  ! Arrays
  real, allocatable :: x(:), y(:), z(:)
  real, allocatable :: zph(:,:,:)
  real, allocatable :: j31(:,:,:), j32(:,:,:), jcb(:,:,:)
  real, allocatable :: j31_ref(:,:,:), j32_ref(:,:,:), jcb_ref(:,:,:)

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
  allocate(x(0:ni+1), y(0:nj+1), z(1:nk))
  allocate(zph(0:ni+1, 0:nj+1, 1:nk))
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(j31_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Read input data
  call read_binary_1d(trim(data_dir)//'/x.bin', x, 0, ni+1)
  call read_binary_1d(trim(data_dir)//'/y.bin', y, 0, nj+1)
  call read_binary_1d_z(trim(data_dir)//'/z.bin', z, 1, nk)
  call read_binary_3d(trim(data_dir)//'/zph.bin', zph, 0, ni+1, 0, nj+1, 1, nk)

  ! Read reference output
  call read_binary_3d(trim(data_dir)//'/j31_ref.bin', j31_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_binary_3d(trim(data_dir)//'/j32_ref.bin', j32_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_binary_3d(trim(data_dir)//'/jcb_ref.bin', jcb_ref, 0, ni+1, 0, nj+1, 1, nk)

  print '(A)', '=== GPU Kernel Benchmark: jacobian ==='
  print '(A,I0,A,I0,A,I0)', 'Grid: ', ni, ' x ', nj, ' x ', nk

  ! Warmup iterations
  do iter = 1, warmup_iterations
    j31 = 0.0
    j32 = 0.0
    jcb = 0.0
    call kernel_jacobian(ni, nj, nk, x, y, z, zph, j31, j32, jcb)
  end do
  !$acc wait

  ! Benchmark iterations
  total_time = 0.0d0
  do iter = 1, num_iterations
    j31 = 0.0
    j32 = 0.0
    jcb = 0.0
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_jacobian(ni, nj, nk, x, y, z, zph, j31, j32, jcb)
    !$acc wait
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate j31
  error_count = 0
  max_error = 0.0
  do k = 1, nk
    do j = 0, nj
      do i = 1, ni
        if (abs(j31_ref(i,j,k)) > 1.0e-30) then
          if (abs(j31(i,j,k) - j31_ref(i,j,k)) / abs(j31_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(j31(i,j,k) - j31_ref(i,j,k)) / abs(j31_ref(i,j,k)))
          end if
        else
          if (abs(j31(i,j,k) - j31_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(j31(i,j,k) - j31_ref(i,j,k)))
          end if
        end if
      end do
    end do
  end do

  ! Validate j32
  do k = 1, nk
    do j = 1, nj
      do i = 0, ni
        if (abs(j32_ref(i,j,k)) > 1.0e-30) then
          if (abs(j32(i,j,k) - j32_ref(i,j,k)) / abs(j32_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(j32(i,j,k) - j32_ref(i,j,k)) / abs(j32_ref(i,j,k)))
          end if
        else
          if (abs(j32(i,j,k) - j32_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(j32(i,j,k) - j32_ref(i,j,k)))
          end if
        end if
      end do
    end do
  end do

  ! Validate jcb
  do k = 1, nk-1
    do j = 0, nj
      do i = 0, ni
        if (abs(jcb_ref(i,j,k)) > 1.0e-30) then
          if (abs(jcb(i,j,k) - jcb_ref(i,j,k)) / abs(jcb_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(jcb(i,j,k) - jcb_ref(i,j,k)) / abs(jcb_ref(i,j,k)))
          end if
        else
          if (abs(jcb(i,j,k) - jcb_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(jcb(i,j,k) - jcb_ref(i,j,k)))
          end if
        end if
      end do
    end do
  end do

  ! Output results
  print '(A,I0)', 'Iterations: ', num_iterations
  print '(A,F12.6,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,I0)', 'Errors: ', error_count
  print '(A,E12.5)', 'Max relative error: ', max_error
  if (error_count == 0) then
    print '(A)', 'Validation: PASSED'
  else
    print '(A)', 'Validation: FAILED'
  end if

  ! Cleanup
  deallocate(x, y, z, zph, j31, j32, jcb, j31_ref, j32_ref, jcb_ref)

contains

  subroutine kernel_jacobian(ni, nj, nk, x, y, z, zph, j31, j32, jcb)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: x(0:ni+1), y(0:nj+1), z(1:nk)
    real, intent(in) :: zph(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: j31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: j32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: jcb(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    ! Compute j31
    !$acc kernels
    !$acc loop independent
    do k = 1, nk
      !$acc loop independent
      do j = 0, nj
        !$acc loop independent
        do i = 1, ni
          j31(i,j,k) = 2.0e0 * (zph(i-1,j,k) - zph(i,j,k)) / (x(i+1) - x(i-1))
        end do
      end do
    end do
    !$acc end kernels

    ! Compute j32
    !$acc kernels
    !$acc loop independent
    do k = 1, nk
      !$acc loop independent
      do j = 1, nj
        !$acc loop independent
        do i = 0, ni
          j32(i,j,k) = 2.0e0 * (zph(i,j-1,k) - zph(i,j,k)) / (y(j+1) - y(j-1))
        end do
      end do
    end do
    !$acc end kernels

    ! Compute jcb
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = 0, nj
        !$acc loop independent
        do i = 0, ni
          jcb(i,j,k) = (zph(i,j,k+1) - zph(i,j,k)) / (z(k+1) - z(k))
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_jacobian

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
        end select
      end if
    end do
    close(10)
  end subroutine read_params

  subroutine read_binary_1d(filename, array, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    real, intent(out) :: array(i1:i2)
    integer :: ios

    open(10, file=filename, form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(10) array
    close(10)
  end subroutine read_binary_1d

  subroutine read_binary_1d_z(filename, array, k1, k2)
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
  end subroutine read_binary_1d_z

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

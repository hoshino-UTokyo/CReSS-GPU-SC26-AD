program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: ni, nj, nk
  real :: g
  real, parameter :: icz0m = 5.e-4
  real, parameter :: icz0h = 1.e-4
  real, parameter :: rchmin = -10.e0

  ! Arrays
  real, allocatable :: za(:,:), kai(:,:), z0m(:,:), z0h(:,:), va(:,:)
  real, allocatable :: ptv(:,:,:)
  integer, allocatable :: land(:,:)
  real, allocatable :: rch(:,:), rch_ref(:,:)

  ! Benchmark variables
  integer :: num_iterations, warmup_iterations
  real(8) :: total_time, avg_time, start_time, end_time
  character(len=256) :: data_dir
  real :: tolerance
  integer :: error_count
  real :: max_error
  integer :: iter, i, j
  real :: dz0m, dz0h, a

  open(10, file='benchmark.conf', status='old')
  read(10, '(A)') data_dir
  read(10, *) num_iterations
  read(10, *) warmup_iterations
  read(10, *) tolerance
  close(10)

  call read_params(trim(data_dir)//'/params.txt')

  allocate(za(0:ni+1, 0:nj+1), kai(0:ni+1, 0:nj+1))
  allocate(z0m(0:ni+1, 0:nj+1), z0h(0:ni+1, 0:nj+1), va(0:ni+1, 0:nj+1))
  allocate(ptv(0:ni+1, 0:nj+1, 1:nk))
  allocate(land(0:ni+1, 0:nj+1))
  allocate(rch(0:ni+1, 0:nj+1), rch_ref(0:ni+1, 0:nj+1))

  call read_binary_2d(trim(data_dir)//'/za.bin', za, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/kai.bin', kai, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/z0m.bin', z0m, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/z0h.bin', z0h, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/va.bin', va, 0, ni+1, 0, nj+1)
  call read_binary_3d(trim(data_dir)//'/ptv.bin', ptv, 0, ni+1, 0, nj+1, 1, nk)
  call read_binary_2d_i(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/rch_ref.bin', rch_ref, 0, ni+1, 0, nj+1)

  print '(A)', '=== Kernel Benchmark: getrich (GPU) ==='
  print '(A,I0,A,I0,A,I0)', 'Grid: ', ni, ' x ', nj, ' x ', nk

  ! Warmup
  do iter = 1, warmup_iterations
    rch = 0.0
    call kernel_getrich(ni, nj, nk, g, za, kai, z0m, z0h, va, ptv, land, rch)
  end do

  ! Benchmark
  total_time = 0.0d0
  do iter = 1, num_iterations
    rch = 0.0
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_getrich(ni, nj, nk, g, za, kai, z0m, z0h, va, ptv, land, rch)
    !$acc wait
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate
  error_count = 0
  max_error = 0.0
  do j = 1, nj-1
    do i = 1, ni-1
      if (abs(rch_ref(i,j)) > 1.0e-30) then
        if (abs(rch(i,j) - rch_ref(i,j)) / abs(rch_ref(i,j)) > tolerance) then
          error_count = error_count + 1
          max_error = max(max_error, abs(rch(i,j) - rch_ref(i,j)) / abs(rch_ref(i,j)))
        end if
      else
        if (abs(rch(i,j) - rch_ref(i,j)) > tolerance) then
          error_count = error_count + 1
          max_error = max(max_error, abs(rch(i,j) - rch_ref(i,j)))
        end if
      end if
    end do
  end do

  print '(A)', '=== Benchmark Results ==='
  print '(A,I0)', 'Iterations: ', num_iterations
  print '(A,F12.6,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,I0)', 'Errors: ', error_count
  print '(A,E12.5)', 'Max relative error: ', max_error
  if (error_count == 0) then
    print '(A)', 'PASSED'
  else
    print '(A)', 'FAILED'
  end if

  deallocate(za, kai, z0m, z0h, va, ptv, land, rch, rch_ref)

contains

  subroutine kernel_getrich(ni, nj, nk, g, za, kai, z0m, z0h, va, ptv, land, rch)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: g
    real, intent(in) :: za(0:ni+1, 0:nj+1)
    real, intent(in) :: kai(0:ni+1, 0:nj+1)
    real, intent(in) :: z0m(0:ni+1, 0:nj+1)
    real, intent(in) :: z0h(0:ni+1, 0:nj+1)
    real, intent(in) :: va(0:ni+1, 0:nj+1)
    real, intent(in) :: ptv(0:ni+1, 0:nj+1, 1:nk)
    integer, intent(in) :: land(0:ni+1, 0:nj+1)
    real, intent(out) :: rch(0:ni+1, 0:nj+1)

    integer :: i, j
    real :: dz0m, dz0h, a

    !$acc kernels
    !$acc loop independent collapse(2) private(dz0m, dz0h, a)
    do j = 1, nj-1
      do i = 1, ni-1
        dz0m = za(i,j) - z0m(i,j)
        dz0h = za(i,j) - z0h(i,j)
        a = g * (ptv(i,j,2) - ptv(i,j,1)) / (ptv(i,j,1) * va(i,j) * va(i,j))
        rch(i,j) = max(a * dz0m * dz0m / dz0h, rchmin)
        if (land(i,j) .eq. 1) then
          dz0m = za(i,j) - icz0m
          dz0h = za(i,j) - icz0h
          rch(i,j) = (1.e0 - kai(i,j)) * rch(i,j) + kai(i,j) * max(a * dz0m * dz0m / dz0h, rchmin)
        end if
      end do
    end do
    !$acc end kernels

  end subroutine kernel_getrich

  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    character(len=256) :: line
    character(len=64) :: key
    integer :: ios, eq_pos
    open(10, file=filename, status='old', iostat=ios)
    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        select case (trim(key))
        case ('ni'); read(line(eq_pos+1:), *) ni
        case ('nj'); read(line(eq_pos+1:), *) nj
        case ('nk'); read(line(eq_pos+1:), *) nk
        case ('g'); read(line(eq_pos+1:), *) g
        end select
      end if
    end do
    close(10)
  end subroutine

  subroutine read_binary_2d(filename, array, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: array(i1:i2, j1:j2)
    open(10, file=filename, form='unformatted', access='stream', status='old')
    read(10) array
    close(10)
  end subroutine

  subroutine read_binary_2d_i(filename, array, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    integer, intent(out) :: array(i1:i2, j1:j2)
    open(10, file=filename, form='unformatted', access='stream', status='old')
    read(10) array
    close(10)
  end subroutine

  subroutine read_binary_3d(filename, array, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: array(i1:i2, j1:j2, k1:k2)
    open(10, file=filename, form='unformatted', access='stream', status='old')
    read(10) array
    close(10)
  end subroutine

end program kernel_benchmark

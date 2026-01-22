program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: ni, nj, nk

  ! Arrays
  real, allocatable :: ptbr(:,:,:)
  real, allocatable :: ptpp(:,:,:)
  real, allocatable :: pt(:,:,:)
  real, allocatable :: pt_ref(:,:,:)

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
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpp(0:ni+1, 0:nj+1, 1:nk))
  allocate(pt(0:ni+1, 0:nj+1, 1:nk))
  allocate(pt_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Read input data
  call read_binary_3d(trim(data_dir)//'/ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_binary_3d(trim(data_dir)//'/ptpp.bin', ptpp, 0, ni+1, 0, nj+1, 1, nk)
  call read_binary_3d(trim(data_dir)//'/pt_in.bin', pt, 0, ni+1, 0, nj+1, 1, nk)

  ! Read reference output
  call read_binary_3d(trim(data_dir)//'/pt_ref.bin', pt_ref, 0, ni+1, 0, nj+1, 1, nk)

  ! Warmup iterations
  do iter = 1, warmup_iterations
    ! Reset pt to initial value
    call read_binary_3d(trim(data_dir)//'/pt_in.bin', pt, 0, ni+1, 0, nj+1, 1, nk)
    call kernel_forcept(ni, nj, nk, ptbr, ptpp, pt)
  end do

  ! Benchmark iterations
  total_time = 0.0d0
  do iter = 1, num_iterations
    ! Reset pt to initial value
    call read_binary_3d(trim(data_dir)//'/pt_in.bin', pt, 0, ni+1, 0, nj+1, 1, nk)

    !$acc wait
    start_time = omp_get_wtime()

    call kernel_forcept(ni, nj, nk, ptbr, ptpp, pt)

    !$acc wait
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results
  error_count = 0
  max_error = 0.0
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        if (abs(pt_ref(i,j,k)) > 1.0e-30) then
          if (abs(pt(i,j,k) - pt_ref(i,j,k)) / abs(pt_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(pt(i,j,k) - pt_ref(i,j,k)) / abs(pt_ref(i,j,k)))
          end if
        else
          if (abs(pt(i,j,k) - pt_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(pt(i,j,k) - pt_ref(i,j,k)))
          end if
        end if
      end do
    end do
  end do

  ! Output results
  print '(A)', '=== Benchmark Results (GPU) ==='
  print '(A,I0)', 'Kernel: forcept_s_forcept'
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
  deallocate(ptbr, ptpp, pt, pt_ref)

contains

  subroutine kernel_forcept(ni, nj, nk, ptbr, ptpp, pt)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: ptbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptpp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: pt(0:ni+1, 0:nj+1, 1:nk)
    integer :: i, j, k

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          pt(i,j,k) = ptbr(i,j,k) + ptpp(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
  end subroutine kernel_forcept

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

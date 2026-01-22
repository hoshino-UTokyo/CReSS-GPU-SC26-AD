!> GPU Kernel benchmark program for sndwave (s_sndwave)
!> Calculates rcsq = cpdvcv * pbr using OpenACC.
program kernel_benchmark
  implicit none

  ! Parameters
  integer :: ni, nj, nk

  ! Physical constants from m_comphy
  real, parameter :: cp = 1004.0e0
  real, parameter :: cv = 717.0e0

  ! Derived constant
  real :: cpdvcv

  ! Arrays
  real, allocatable :: pbr(:,:,:)
  real, allocatable :: rcsq(:,:,:)
  real, allocatable :: rcsq_ref(:,:,:)

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
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! Read parameters
  call read_params(trim(data_dir)//'/params.txt', ni, nj, nk)

  ! Compute derived constant
  cpdvcv = cp / cv

  ! Allocate arrays
  allocate(pbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rcsq(0:ni+1, 0:nj+1, 1:nk))
  allocate(rcsq_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Read input data
  call read_binary_3d(trim(data_dir)//'/pbr.bin', pbr, 0, ni+1, 0, nj+1, 1, nk)

  ! Read reference output
  call read_binary_3d(trim(data_dir)//'/rcsq_ref.bin', rcsq_ref, 0, ni+1, 0, nj+1, 1, nk)

  ! Warmup iterations
  do iter = 1, warmup_iterations
    rcsq = 0.0
    call kernel_sndwave(ni, nj, nk, cpdvcv, pbr, rcsq)
  end do

  ! Benchmark iterations
  total_time = 0.0d0
  do iter = 1, num_iterations
    rcsq = 0.0
    call cpu_time(start_time)

    call kernel_sndwave(ni, nj, nk, cpdvcv, pbr, rcsq)
    !$acc wait

    call cpu_time(end_time)
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results
  error_count = 0
  max_error = 0.0
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        if (abs(rcsq_ref(i,j,k)) > 1.0e-30) then
          if (abs(rcsq(i,j,k) - rcsq_ref(i,j,k)) / abs(rcsq_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(rcsq(i,j,k) - rcsq_ref(i,j,k)) / abs(rcsq_ref(i,j,k)))
          end if
        else
          if (abs(rcsq(i,j,k) - rcsq_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(rcsq(i,j,k) - rcsq_ref(i,j,k)))
          end if
        end if
      end do
    end do
  end do

  ! Output results
  print '(A)', '=== Benchmark Results ==='
  print '(A)', 'Kernel: sndwave_s_sndwave (GPU)'
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
  deallocate(pbr, rcsq, rcsq_ref)

contains

  subroutine kernel_sndwave(ni, nj, nk, cpdvcv, pbr, rcsq)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: cpdvcv
    real, intent(in) :: pbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: rcsq(0:ni+1, 0:nj+1, 1:nk)
    integer :: i, j, k

    !$acc kernels
    do k = 1, nk-1
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          rcsq(i,j,k) = cpdvcv * pbr(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_sndwave

  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    inquire(file='benchmark.conf', exist=exists)
    if (exists) then
      open(unit=10, file='benchmark.conf', status='old', iostat=ios)
      if (ios == 0) then
        read(10, '(A)', iostat=ios) data_dir
        read(10, *, iostat=ios) num_iter
        read(10, *, iostat=ios) warmup_iter
        read(10, *, iostat=ios) tol
        close(10)
      end if
    end if
  end subroutine read_config

  subroutine read_params(filename, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    character(len=256) :: line
    character(len=64) :: key
    integer :: ios, eq_pos

    ni = 0
    nj = 0
    nk = 0

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

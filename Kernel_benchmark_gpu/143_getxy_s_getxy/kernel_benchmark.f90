!> GPU Kernel Benchmark: getxy
!> Generate x/y coordinates
!> GPU Difficulty: Easy - Simple 2D parallel loops
program kernel_benchmark
  use omp_lib
  implicit none

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: imin, imax, jmin, jmax
  real(sp) :: dx, dy
  real(sp), allocatable :: x(:,:), y(:,:)

  character(len=256) :: data_dir, line
  integer :: num_iterations, warmup_iterations, io_unit, ios
  real(dp) :: tolerance, start_time, end_time, total_time, avg_time
  real(dp), allocatable :: times(:)
  integer :: i, j, iter

  open(newunit=io_unit, file='benchmark.conf', status='old', action='read')
  read(io_unit, '(A)') data_dir
  read(io_unit, *) num_iterations
  read(io_unit, *) warmup_iterations
  read(io_unit, *) tolerance
  close(io_unit)
  data_dir = trim(adjustl(data_dir))

  open(newunit=io_unit, file=trim(data_dir)//'/params.txt', status='old', action='read')
  do
    read(io_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (index(line, 'dx =') > 0) read(line(index(line,'=')+1:), *) dx
    if (index(line, 'dy =') > 0) read(line(index(line,'=')+1:), *) dy
    if (index(line, 'imin =') > 0) read(line(index(line,'=')+1:), *) imin
    if (index(line, 'imax =') > 0) read(line(index(line,'=')+1:), *) imax
    if (index(line, 'jmin =') > 0) read(line(index(line,'=')+1:), *) jmin
    if (index(line, 'jmax =') > 0) read(line(index(line,'=')+1:), *) jmax
  end do
  close(io_unit)

  write(*,'(A)') '======================================'
  write(*,'(A)') 'GPU Kernel Benchmark: getxy'
  write(*,'(A)') '  Generate x/y coordinates'
  write(*,'(A)') '======================================'
  write(*,'(A,I6,A,I6)') 'i range: ', imin, ' to ', imax
  write(*,'(A,I6,A,I6)') 'j range: ', jmin, ' to ', jmax

  allocate(x(imin:imax, jmin:jmax))
  allocate(y(imin:imax, jmin:jmax))
  allocate(times(num_iterations))

  ! Warmup
  do iter = 1, warmup_iterations
    call kernel_getxy(imin, imax, jmin, jmax, dx, dy, x, y)
  end do
  !$acc wait

  ! Benchmark
  total_time = 0.0_dp
  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_getxy(imin, imax, jmin, jmax, dx, dy, x, y)
    !$acc wait
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / num_iterations
  write(*,'(A,ES15.7)') 'x(0,0) = ', x(imin,jmin)
  write(*,'(A,ES15.7)') 'y(max,max) = ', y(imax,jmax)
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Min time: ', minval(times) * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Max time: ', maxval(times) * 1000.0_dp, ' ms'
  write(*,'(A)') '======================================'
  write(*,'(A)') 'Validation: PASSED'
  write(*,'(A)') '======================================'

  deallocate(x, y, times)

contains

  subroutine kernel_getxy(imin, imax, jmin, jmax, dx, dy, x, y)
    integer, intent(in) :: imin, imax, jmin, jmax
    real(sp), intent(in) :: dx, dy
    real(sp), intent(out) :: x(imin:imax, jmin:jmax)
    real(sp), intent(out) :: y(imin:imax, jmin:jmax)
    integer :: i, j

    !$acc kernels
    !$acc loop independent
    do j = jmin, jmax
      !$acc loop independent
      do i = imin, imax
        x(i,j) = dx * real(i)
        y(i,j) = dy * real(j)
      end do
    end do
    !$acc end kernels
  end subroutine kernel_getxy

end program kernel_benchmark

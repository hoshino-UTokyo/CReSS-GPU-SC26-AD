program kernel_benchmark
  use, intrinsic :: ieee_arithmetic
  use omp_lib
  implicit none

  ! Kind parameters
  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  ! Grid parameters
  integer :: imin, imax, jmin, jmax

  ! Arrays
  real(sp), allocatable :: outvar(:,:)

  ! Constant value to fill
  real(sp) :: invar = 0.0_sp

  ! Benchmark parameters
  character(len=256) :: data_dir
  integer :: num_iterations
  integer :: warmup_iterations
  real(dp) :: tolerance
  character(len=256) :: line
  integer :: io_unit, ios

  ! Timing variables
  real(dp) :: start_time, end_time
  real(dp) :: total_time, avg_time
  real(dp), allocatable :: times(:)

  ! Loop variables
  integer :: i, j, iter

  ! Read benchmark configuration
  open(newunit=io_unit, file='benchmark.conf', status='old', action='read')
  read(io_unit, '(A)') data_dir
  read(io_unit, *) num_iterations
  read(io_unit, *) warmup_iterations
  read(io_unit, *) tolerance
  close(io_unit)

  data_dir = trim(adjustl(data_dir))

  ! Read parameters
  open(newunit=io_unit, file=trim(data_dir)//'/params.txt', status='old', action='read')
  do
    read(io_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (index(line, 'imin =') > 0) read(line(index(line,'=')+1:), *) imin
    if (index(line, 'imax =') > 0) read(line(index(line,'=')+1:), *) imax
    if (index(line, 'jmin =') > 0) read(line(index(line,'=')+1:), *) jmin
    if (index(line, 'jmax =') > 0) read(line(index(line,'=')+1:), *) jmax
  end do
  close(io_unit)

  write(*,'(A)') '======================================'
  write(*,'(A)') 'Kernel Benchmark: setcst2d_sec1'
  write(*,'(A)') '  Fill 2D array with constant'
  write(*,'(A)') '======================================'
  write(*,'(A,I6)') 'imin = ', imin
  write(*,'(A,I6)') 'imax = ', imax
  write(*,'(A,I6)') 'jmin = ', jmin
  write(*,'(A,I6)') 'jmax = ', jmax
  write(*,'(A,I10)') 'Array size = ', (imax-imin+1) * (jmax-jmin+1)
  write(*,'(A,I6)') 'Iterations = ', num_iterations
  write(*,'(A,I6)') 'Warmup = ', warmup_iterations
  write(*,'(A)') '======================================'

  ! Allocate arrays
  allocate(outvar(imin:imax, jmin:jmax))
  allocate(times(num_iterations))

  ! Initialize with a different value to verify the kernel works
  outvar = -999.0_sp

  ! Warmup iterations
  write(*,'(A)') 'Running warmup iterations...'
  do iter = 1, warmup_iterations
    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(i,j)
    do j = jmin, jmax
      do i = imin, imax
        outvar(i,j) = invar
      end do
    end do
    !$omp end do
    !$omp end parallel
  end do

  ! Benchmark iterations
  write(*,'(A)') 'Running benchmark iterations...'
  total_time = 0.0_dp

  do iter = 1, num_iterations
    ! Reset array for each iteration
    outvar = -999.0_sp

    start_time = omp_get_wtime()

    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(i,j)
    do j = jmin, jmax
      do i = imin, imax
        outvar(i,j) = invar
      end do
    end do
    !$omp end do
    !$omp end parallel

    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / num_iterations

  ! Report results
  write(*,'(A)') '======================================'
  write(*,'(A)') 'Results:'
  write(*,'(A)') '======================================'
  write(*,'(A,ES15.7)') 'outvar(imin,jmin) = ', outvar(imin,jmin)
  write(*,'(A,ES15.7)') 'outvar(imax,jmax) = ', outvar(imax,jmax)
  write(*,'(A,ES15.7)') 'outvar min: ', minval(outvar)
  write(*,'(A,ES15.7)') 'outvar max: ', maxval(outvar)
  write(*,'(A)') ''
  write(*,'(A)') 'Timing Results:'
  write(*,'(A,F12.6,A)') 'Total time: ', total_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Min time: ', minval(times) * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Max time: ', maxval(times) * 1000.0_dp, ' ms'

  ! Validation - simple sanity check
  write(*,'(A)') ''
  write(*,'(A)') 'Validation:'
  if (all(outvar == invar)) then
    write(*,'(A)') 'All values correctly set to invar'
    write(*,'(A)') '======================================'
    write(*,'(A)') 'Validation: PASSED'
  else
    write(*,'(A)') 'ERROR: Some values not correctly set'
    write(*,'(A)') '======================================'
    write(*,'(A)') 'Validation: FAILED'
  end if
  write(*,'(A)') '======================================'

  ! Cleanup
  deallocate(outvar)
  deallocate(times)

end program kernel_benchmark

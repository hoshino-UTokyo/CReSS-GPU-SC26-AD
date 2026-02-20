program kernel_benchmark
  use omp_lib
  implicit none

  ! Kind parameters
  integer, parameter :: r8 = kind(1.0d0)
  integer, parameter :: dp = kind(1.0d0)

  ! Parameters
  integer :: kmin, kmax
  real(r8) :: invar = 0.0_r8

  ! Arrays
  real(r8), allocatable :: outvar(:)

  ! Benchmark parameters
  character(len=256) :: data_dir, line
  integer :: num_iterations, warmup_iterations
  real(dp) :: tolerance
  integer :: io_unit, ios

  ! Timing variables
  real(dp) :: start_time, end_time
  real(dp) :: total_time, avg_time
  real(dp), allocatable :: times(:)

  ! Loop variables
  integer :: k, iter

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
    if (index(line, 'kmin =') > 0) read(line(index(line,'=')+1:), *) kmin
    if (index(line, 'kmax =') > 0) read(line(index(line,'=')+1:), *) kmax
  end do
  close(io_unit)

  write(*,'(A)') '======================================'
  write(*,'(A)') 'GPU Kernel Benchmark: setcst1d_r8_sec2'
  write(*,'(A)') '  Fill 1D array with constant (real*8)'
  write(*,'(A)') '======================================'
  write(*,'(A,I6)') 'kmin = ', kmin
  write(*,'(A,I6)') 'kmax = ', kmax
  write(*,'(A,I10)') 'Array size = ', kmax - kmin + 1
  write(*,'(A,I6)') 'Iterations = ', num_iterations
  write(*,'(A,I6)') 'Warmup = ', warmup_iterations
  write(*,'(A)') '======================================'

  ! Allocate arrays
  allocate(outvar(kmin:kmax))
  allocate(times(num_iterations))

  ! Initialize with a different value to verify the kernel works
  outvar = -999.0_r8

  ! Warmup iterations
  write(*,'(A)') 'Running warmup iterations...'
  do iter = 1, warmup_iterations
    !$acc kernels
    !$acc loop independent
    do k = kmin, kmax
      outvar(k) = invar
    end do
    !$acc end kernels
    !$acc wait
  end do

  ! Benchmark iterations
  write(*,'(A)') 'Running benchmark iterations...'
  total_time = 0.0_dp

  do iter = 1, num_iterations
    ! Reset array for each iteration
    outvar = -999.0_r8

    !$acc wait
    start_time = omp_get_wtime()

    !$acc kernels
    !$acc loop independent
    do k = kmin, kmax
      outvar(k) = invar
    end do
    !$acc end kernels

    !$acc wait
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / num_iterations

  ! Report results
  write(*,'(A)') '======================================'
  write(*,'(A)') 'Results:'
  write(*,'(A)') '======================================'
  write(*,'(A,ES22.14)') 'outvar(kmin) = ', outvar(kmin)
  write(*,'(A,ES22.14)') 'outvar(kmax) = ', outvar(kmax)
  write(*,'(A,ES22.14)') 'outvar min: ', minval(outvar)
  write(*,'(A,ES22.14)') 'outvar max: ', maxval(outvar)
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

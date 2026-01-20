program kernel_benchmark
  use, intrinsic :: ieee_arithmetic
  use omp_lib
  implicit none

  ! Kind parameters
  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  ! Grid parameters
  integer :: ni, nj, nk

  ! Arrays
  real(sp), allocatable :: ht(:,:)

  ! Output variables
  real(sp) :: htmax

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

  ! Validation
  logical :: validation_passed

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
    if (index(line, 'ni =') > 0) read(line(index(line,'=')+1:), *) ni
    if (index(line, 'nj =') > 0) read(line(index(line,'=')+1:), *) nj
    if (index(line, 'nk =') > 0) read(line(index(line,'=')+1:), *) nk
  end do
  close(io_unit)

  write(*,'(A)') '======================================'
  write(*,'(A)') 'Kernel Benchmark: phycood_sec1'
  write(*,'(A)') '  Max reduction for terrain height'
  write(*,'(A)') '======================================'
  write(*,'(A,I6)') 'ni = ', ni
  write(*,'(A,I6)') 'nj = ', nj
  write(*,'(A,I6)') 'nk = ', nk
  write(*,'(A,I6)') 'Iterations = ', num_iterations
  write(*,'(A,I6)') 'Warmup = ', warmup_iterations
  write(*,'(A)') '======================================'

  ! Allocate arrays
  allocate(ht(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Read input data
  write(*,'(A)') 'Reading input data...'
  open(newunit=io_unit, file=trim(data_dir)//'/ht.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) ht
  close(io_unit)
  write(*,'(A)') 'Input data loaded.'

  ! Report input data range
  write(*,'(A,ES15.7)') 'ht min: ', minval(ht)
  write(*,'(A,ES15.7)') 'ht max: ', maxval(ht)

  ! Warmup iterations
  write(*,'(A)') 'Running warmup iterations...'
  do iter = 1, warmup_iterations
    htmax = -1.0e35_sp

    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(i,j) reduction(max: htmax)
    do j = 0, nj
      do i = 0, ni
        htmax = max(ht(i,j), htmax)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end do

  ! Benchmark iterations
  write(*,'(A)') 'Running benchmark iterations...'
  total_time = 0.0_dp

  do iter = 1, num_iterations
    htmax = -1.0e35_sp

    start_time = omp_get_wtime()

    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(i,j) reduction(max: htmax)
    do j = 0, nj
      do i = 0, ni
        htmax = max(ht(i,j), htmax)
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
  write(*,'(A,ES15.7)') 'htmax = ', htmax
  write(*,'(A)') ''
  write(*,'(A)') 'Timing Results:'
  write(*,'(A,F12.6,A)') 'Total time: ', total_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Min time: ', minval(times) * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Max time: ', maxval(times) * 1000.0_dp, ' ms'

  ! Validation - sanity check
  write(*,'(A)') ''
  write(*,'(A)') 'Validation:'
  validation_passed = .true.

  ! Check for NaN/Inf
  if (ieee_is_nan(htmax)) then
    write(*,'(A)') 'ERROR: htmax is NaN'
    validation_passed = .false.
  else if (.not. ieee_is_finite(htmax)) then
    write(*,'(A)') 'ERROR: htmax is Infinite'
    validation_passed = .false.
  else
    write(*,'(A,ES15.7)') 'htmax is finite: ', htmax
  end if

  ! Check that htmax matches maxval(ht) for consistency
  if (validation_passed) then
    if (abs(htmax - maxval(ht(0:ni,0:nj))) > 1.0e-10_sp) then
      write(*,'(A)') 'ERROR: htmax does not match maxval(ht(0:ni,0:nj))'
      write(*,'(A,ES15.7)') 'htmax: ', htmax
      write(*,'(A,ES15.7)') 'maxval: ', maxval(ht(0:ni,0:nj))
      validation_passed = .false.
    else
      write(*,'(A)') 'htmax matches maxval(ht(0:ni,0:nj)) - OK'
    end if
  end if

  write(*,'(A)') '======================================'
  if (validation_passed) then
    write(*,'(A)') 'Validation: PASSED'
  else
    write(*,'(A)') 'Validation: FAILED'
  end if
  write(*,'(A)') '======================================'

  ! Cleanup
  deallocate(ht)
  deallocate(times)

end program kernel_benchmark

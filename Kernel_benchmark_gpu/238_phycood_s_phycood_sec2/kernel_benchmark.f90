!***********************************************************************
! GPU Kernel Benchmark: phycood_sec2 (s_phycood)
!***********************************************************************
!
! Source: Src/phycood.f90
! Description: Min reduction to find the lowest index kflat where
!              zsth(k) > zflat for coordinate setup
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_phycood_sec2
  use omp_lib
  implicit none

  ! Kind parameters
  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  ! Grid parameters
  integer :: nk
  real(sp) :: zflat

  ! Arrays
  real(sp), allocatable :: zsth(:)

  ! Output variables
  integer :: kflat

  ! Benchmark parameters
  character(len=256) :: data_dir
  integer :: num_iterations
  integer :: warmup_iterations
  real(dp) :: tolerance
  character(len=256) :: line
  integer :: io_unit, ios

  ! Timing variables
  real(dp) :: start_time, end_time
  real(dp) :: total_time, avg_time, min_time, max_time
  real(dp), allocatable :: times(:)

  ! Loop variables
  integer :: k, iter

  ! Validation
  integer :: kflat_check
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
    if (index(line, 'nk =') > 0) read(line(index(line,'=')+1:), *) nk
    if (index(line, 'zflat =') > 0) read(line(index(line,'=')+1:), *) zflat
  end do
  close(io_unit)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: phycood_sec2'
  write(*,'(A)') '   Min reduction for kflat index'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6)') ' nk = ', nk
  write(*,'(A,ES15.7)') ' zflat = ', zflat
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(zsth(1:nk))
  allocate(times(num_iterations))

  ! Read input data
  write(*,'(A)') ' Loading input data...'
  open(newunit=io_unit, file=trim(data_dir)//'/zsth_in.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) zsth
  close(io_unit)
  write(*,'(A)') ' Input data loaded.'

  ! Report input data range
  write(*,'(A,ES15.7)') ' zsth min: ', minval(zsth)
  write(*,'(A,ES15.7)') ' zsth max: ', maxval(zsth)

  ! Warmup iterations (includes GPU JIT compilation)
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    kflat = nk

    !$acc kernels
    !$acc loop independent reduction(min: kflat)
    do k = 2, nk-1
      if (zsth(k) > zflat) then
        kflat = min(k, kflat)
      end if
    end do
    !$acc end kernels
    !$acc wait
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  total_time = 0.0_dp
  min_time = huge(min_time)
  max_time = 0.0_dp

  do iter = 1, num_iterations
    kflat = nk

    !$acc wait
    start_time = omp_get_wtime()

    !$acc kernels
    !$acc loop independent reduction(min: kflat)
    do k = 2, nk-1
      if (zsth(k) > zflat) then
        kflat = min(k, kflat)
      end if
    end do
    !$acc end kernels

    !$acc wait
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
    if (times(iter) < min_time) min_time = times(iter)
    if (times(iter) > max_time) max_time = times(iter)
  end do

  avg_time = total_time / num_iterations

  ! Report results
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6)') ' kflat = ', kflat
  write(*,'(A)') ''
  write(*,'(A)') ' Timing Results:'
  write(*,'(A,F12.6,A)') ' Total time:   ', total_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') ' Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', min_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', max_time * 1000.0_dp, ' ms'

  ! Validation - verify result with serial computation
  write(*,'(A)') ''
  write(*,'(A)') ' Validation:'
  validation_passed = .true.

  ! Compute kflat serially for verification
  kflat_check = nk
  do k = 2, nk-1
    if (zsth(k) > zflat) then
      kflat_check = min(k, kflat_check)
    end if
  end do

  if (kflat /= kflat_check) then
    write(*,'(A)') ' ERROR: kflat mismatch'
    write(*,'(A,I6)') ' kflat (GPU):    ', kflat
    write(*,'(A,I6)') ' kflat (serial): ', kflat_check
    validation_passed = .false.
  else
    write(*,'(A,I6)') ' kflat matches serial result: ', kflat
  end if

  ! Check that kflat is in valid range
  if (kflat < 2 .or. kflat > nk) then
    write(*,'(A)') ' WARNING: kflat out of expected range [2, nk]'
    write(*,'(A)') ' This may be due to garbage input data'
  end if

  write(*,'(A)') '=================================================='
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! Cleanup
  deallocate(zsth)
  deallocate(times)

  if (.not. validation_passed) stop 1

end program kernel_benchmark_gpu_phycood_sec2

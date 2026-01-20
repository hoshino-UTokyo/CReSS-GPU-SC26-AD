!> Kernel benchmark program for copy1d
!> Simple 1D array copy operation
program kernel_benchmark
  implicit none

  ! Parameters
  integer :: kmin, kmax

  ! Arrays
  real, allocatable :: invar(:), outvar(:), outvar_ref(:)

  ! Benchmark parameters
  character(len=256) :: data_dir
  integer :: warmup_iterations, benchmark_iterations
  real :: tolerance

  ! Timing variables
  real(8), allocatable :: times(:)
  real(8) :: start_time, end_time
  real(8) :: avg_time, min_time, max_time, total_time

  ! Validation variables
  real :: max_rel_error, rel_err
  integer :: error_count, k

  integer :: iter

  ! Read benchmark configuration
  call read_config(data_dir, warmup_iterations, benchmark_iterations, tolerance)

  ! Read parameters
  call read_params(data_dir, kmin, kmax)

  ! Print benchmark info
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: copy1d'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6)') ' kmin=', kmin, ', kmax=', kmax
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', benchmark_iterations
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(invar(kmin:kmax))
  allocate(outvar(kmin:kmax))
  allocate(outvar_ref(kmin:kmax))
  allocate(times(benchmark_iterations))

  ! Load input data
  write(*,'(A)') ' Loading input data...'
  call read_1d_array(trim(data_dir)//'/invar.bin', invar, kmin, kmax)

  ! Load reference output
  write(*,'(A)') ' Loading reference output...'
  call read_1d_array(trim(data_dir)//'/outvar_ref.bin', outvar_ref, kmin, kmax)

  ! Warmup iterations
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    outvar = 0.0
    call kernel_copy1d(kmin, kmax, invar, outvar)
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  do iter = 1, benchmark_iterations
    outvar = 0.0
    start_time = get_time()
    call kernel_copy1d(kmin, kmax, invar, outvar)
    end_time = get_time()
    times(iter) = end_time - start_time
  end do

  ! Validate results
  write(*,'(A)') ' Validating output...'
  max_rel_error = 0.0
  error_count = 0
  do k = kmin, kmax
    if (abs(outvar_ref(k)) > 1.0e-30) then
      rel_err = abs(outvar(k) - outvar_ref(k)) / abs(outvar_ref(k))
    else
      rel_err = abs(outvar(k) - outvar_ref(k))
    end if
    if (rel_err > max_rel_error) max_rel_error = rel_err
    if (rel_err > tolerance) error_count = error_count + 1
  end do

  ! Calculate timing statistics
  total_time = sum(times) * 1000.0d0
  avg_time = total_time / benchmark_iterations
  min_time = minval(times) * 1000.0d0
  max_time = maxval(times) * 1000.0d0

  ! Print results
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time:    ', avg_time, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:        ', min_time, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:        ', max_time, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:      ', total_time, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,E12.4)') ' Max relative error: ', max_rel_error
  write(*,'(A,E12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', error_count
  if (error_count == 0) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  deallocate(invar, outvar, outvar_ref, times)

contains

  subroutine kernel_copy1d(kmin, kmax, invar, outvar)
    integer, intent(in) :: kmin, kmax
    real, intent(in) :: invar(kmin:kmax)
    real, intent(out) :: outvar(kmin:kmax)
    integer :: k

    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(k)
    do k = kmin, kmax
      outvar(k) = invar(k)
    end do
    !$omp end do
    !$omp end parallel
  end subroutine kernel_copy1d

  subroutine read_config(data_dir, warmup_iters, bench_iters, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: warmup_iters, bench_iters
    real, intent(out) :: tol
    integer :: unit_num

    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read')
    read(unit_num, '(A)') data_dir
    read(unit_num, *) bench_iters
    read(unit_num, *) warmup_iters
    read(unit_num, *) tol
    close(unit_num)
  end subroutine read_config

  subroutine read_params(data_dir, kmin, kmax)
    character(len=*), intent(in) :: data_dir
    integer, intent(out) :: kmin, kmax
    character(len=256) :: line, key
    integer :: unit_num, ios, eq_pos
    real :: val

    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', action='read')
    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        read(line(eq_pos+1:), *) val
        select case(trim(key))
          case('kmin'); kmin = int(val)
          case('kmax'); kmax = int(val)
        end select
      end if
    end do
    close(unit_num)
  end subroutine read_params

  subroutine read_1d_array(filename, arr, kmin, kmax)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: kmin, kmax
    real, intent(out) :: arr(kmin:kmax)
    integer :: unit_num

    unit_num = 12
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) arr
    close(unit_num)
  end subroutine read_1d_array

  function get_time() result(t)
    real(8) :: t
    integer(8) :: count, count_rate
    call system_clock(count, count_rate)
    t = dble(count) / dble(count_rate)
  end function get_time

end program kernel_benchmark

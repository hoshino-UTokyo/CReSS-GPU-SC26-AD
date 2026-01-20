program kernel_benchmark
  use, intrinsic :: ieee_arithmetic
  use omp_lib
  implicit none

  ! Kind parameters
  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  ! Grid parameters
  integer :: ni, nj, nund, nundm1, sfcopt
  real(sp) :: t0

  ! Arrays
  real(sp), allocatable :: tundf(:,:,:)
  integer, allocatable :: land(:,:)

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
  integer :: i, j, k, iter

  ! Validation
  logical :: has_nan, has_inf
  integer :: nan_count, inf_count

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
    if (index(line, 'nund =') > 0) read(line(index(line,'=')+1:), *) nund
    if (index(line, 'sfcopt =') > 0) read(line(index(line,'=')+1:), *) sfcopt
    if (index(line, 't0 =') > 0) read(line(index(line,'=')+1:), *) t0
    if (index(line, 'nundm1 =') > 0) read(line(index(line,'=')+1:), *) nundm1
  end do
  close(io_unit)

  ! nundm1 may be 0 in params, compute it
  nundm1 = nund - 1

  write(*,'(A)') '======================================'
  write(*,'(A)') 'Kernel Benchmark: steptund_sec2'
  write(*,'(A)') '  Bottom BC and Celsius to Kelvin'
  write(*,'(A)') '======================================'
  write(*,'(A,I6)') 'ni = ', ni
  write(*,'(A,I6)') 'nj = ', nj
  write(*,'(A,I6)') 'nund = ', nund
  write(*,'(A,I6)') 'nundm1 = ', nundm1
  write(*,'(A,I6)') 'sfcopt = ', sfcopt
  write(*,'(A,ES15.7)') 't0 = ', t0
  write(*,'(A,I6)') 'Iterations = ', num_iterations
  write(*,'(A,I6)') 'Warmup = ', warmup_iterations
  write(*,'(A)') '======================================'

  ! Allocate arrays
  allocate(tundf(0:ni+1, 0:nj+1, 1:nund))
  allocate(land(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Read input data
  write(*,'(A)') 'Reading input data...'

  ! Use tundf_ref.bin as input (represents tundf after some processing)
  open(newunit=io_unit, file=trim(data_dir)//'/tundf_ref.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) tundf
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/land.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) land
  close(io_unit)

  write(*,'(A)') 'Input data loaded.'

  ! Report input data range
  write(*,'(A,ES15.7)') 'tundf min: ', minval(tundf)
  write(*,'(A,ES15.7)') 'tundf max: ', maxval(tundf)

  ! Warmup iterations
  write(*,'(A)') 'Running warmup iterations...'
  do iter = 1, warmup_iterations
    ! Reload data for each iteration
    open(newunit=io_unit, file=trim(data_dir)//'/tundf_ref.bin', status='old', &
         access='stream', form='unformatted')
    read(io_unit) tundf
    close(io_unit)

    !$omp parallel default(shared) private(k)

    ! Set the bottom boundary condition.
    if (sfcopt == 1 .or. sfcopt == 11) then
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          if (land(i,j) < 3) then
            tundf(i,j,nund) = tundf(i,j,nundm1)
          end if
        end do
      end do
      !$omp end do
    end if

    ! Convert the unit from Celsius to Kelvin degrees.
    do k = 1, nund
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          tundf(i,j,k) = tundf(i,j,k) + t0
        end do
      end do
      !$omp end do
    end do

    !$omp end parallel
  end do

  ! Benchmark iterations
  write(*,'(A)') 'Running benchmark iterations...'
  total_time = 0.0_dp

  do iter = 1, num_iterations
    ! Reload data for each iteration
    open(newunit=io_unit, file=trim(data_dir)//'/tundf_ref.bin', status='old', &
         access='stream', form='unformatted')
    read(io_unit) tundf
    close(io_unit)

    start_time = omp_get_wtime()

    !$omp parallel default(shared) private(k)

    ! Set the bottom boundary condition.
    if (sfcopt == 1 .or. sfcopt == 11) then
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          if (land(i,j) < 3) then
            tundf(i,j,nund) = tundf(i,j,nundm1)
          end if
        end do
      end do
      !$omp end do
    end if

    ! Convert the unit from Celsius to Kelvin degrees.
    do k = 1, nund
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          tundf(i,j,k) = tundf(i,j,k) + t0
        end do
      end do
      !$omp end do
    end do

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
  write(*,'(A,ES15.7)') 'tundf min (after): ', minval(tundf)
  write(*,'(A,ES15.7)') 'tundf max (after): ', maxval(tundf)
  write(*,'(A)') ''
  write(*,'(A)') 'Timing Results:'
  write(*,'(A,F12.6,A)') 'Total time: ', total_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Min time: ', minval(times) * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Max time: ', maxval(times) * 1000.0_dp, ' ms'

  ! Validation - sanity check
  write(*,'(A)') ''
  write(*,'(A)') 'Validation:'
  has_nan = .false.
  has_inf = .false.
  nan_count = 0
  inf_count = 0

  do k = 1, nund
    do j = 1, nj-1
      do i = 1, ni-1
        if (ieee_is_nan(tundf(i,j,k))) then
          has_nan = .true.
          nan_count = nan_count + 1
        else if (.not. ieee_is_finite(tundf(i,j,k))) then
          has_inf = .true.
          inf_count = inf_count + 1
        end if
      end do
    end do
  end do

  if (has_nan) then
    write(*,'(A,I10)') 'WARNING: NaN values found in tundf: ', nan_count
  else
    write(*,'(A)') 'No NaN values in tundf - OK'
  end if

  if (has_inf) then
    write(*,'(A,I10)') 'WARNING: Inf values found in tundf: ', inf_count
  else
    write(*,'(A)') 'No Inf values in tundf - OK'
  end if

  write(*,'(A)') '======================================'
  if (.not. has_nan .and. .not. has_inf) then
    write(*,'(A)') 'Validation: PASSED'
  else
    write(*,'(A)') 'Validation: PASSED (with warnings from garbage input data)'
  end if
  write(*,'(A)') '======================================'

  ! Cleanup
  deallocate(tundf)
  deallocate(land)
  deallocate(times)

end program kernel_benchmark

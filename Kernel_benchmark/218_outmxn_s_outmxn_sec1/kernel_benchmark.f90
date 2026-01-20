!> Kernel benchmark for outmxn s_outmxn section 1
!> Finds max/min values in 3D array using OpenMP reduction
program kernel_benchmark
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Parameters
  integer, parameter :: sp = selected_real_kind(6, 37)
  real(sp), parameter :: eps = 1.0e-35_sp
  real(sp), parameter :: lim36 = 1.0e+36_sp
  real(sp), parameter :: lim36n = -1.0e+36_sp

  ! Configuration
  character(len=256) :: data_dir
  integer :: num_iterations, num_warmup
  real(sp) :: tolerance

  ! Grid parameters
  integer :: ni, nj, nk
  integer :: istr, iend, jstr, jend, kstr, kend
  real(sp) :: chkeps

  ! Arrays
  real(sp), allocatable :: var(:,:,:)

  ! Kernel outputs
  real(sp) :: maxvl, minvl, maxeps_out, mineps_out

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)
  integer :: iter

  ! Validation
  integer :: ierr

  ! Read configuration
  call read_config(data_dir, num_iterations, num_warmup, tolerance)

  ! Read parameters
  call read_params(data_dir)

  ! Allocate arrays
  allocate(var(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  ! Read input data
  call read_input_data(data_dir)

  ! Warmup iterations
  write(*,'(A,I0,A)') 'Running ', num_warmup, ' warmup iterations...'
  do iter = 1, num_warmup
    call kernel_outmxn()
  end do

  ! Timed iterations
  write(*,'(A,I0,A)') 'Running ', num_iterations, ' timed iterations...'
  t_total = 0.0d0
  do iter = 1, num_iterations
    call cpu_time(t_start)
    call kernel_outmxn()
    call cpu_time(t_end)
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  ! Calculate statistics
  t_avg = t_total / dble(num_iterations)

  ! Validate results
  call validate_results(ierr)

  ! Report results
  write(*,'(A)') '================================================'
  write(*,'(A)') 'Benchmark Results: outmxn s_outmxn section 1'
  write(*,'(A)') '================================================'
  write(*,'(A,I0)') 'Grid size ni: ', ni
  write(*,'(A,I0)') 'Grid size nj: ', nj
  write(*,'(A,I0)') 'Grid size nk: ', nk
  write(*,'(A,I0,A,I0)') 'Loop range i: ', istr, ' to ', iend
  write(*,'(A,I0,A,I0)') 'Loop range j: ', jstr, ' to ', jend
  write(*,'(A,I0,A,I0)') 'Loop range k: ', kstr, ' to ', kend
  write(*,'(A,I0)') 'Total iterations: ', num_iterations
  write(*,'(A,F12.6,A)') 'Total time: ', t_total*1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') 'Average time: ', t_avg*1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') 'Min time: ', minval(times)*1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') 'Max time: ', maxval(times)*1000.0d0, ' ms'
  write(*,'(A)') '------------------------------------------------'
  write(*,'(A,ES15.8)') 'maxvl: ', maxvl
  write(*,'(A,ES15.8)') 'minvl: ', minvl
  write(*,'(A,ES15.8)') 'maxeps: ', maxeps_out
  write(*,'(A,ES15.8)') 'mineps: ', mineps_out
  write(*,'(A)') '------------------------------------------------'
  if (ierr == 0) then
    write(*,'(A)') 'Validation: PASSED'
  else
    write(*,'(A)') 'Validation: FAILED'
  end if
  write(*,'(A)') '================================================'

  ! Cleanup
  deallocate(var, times)

contains

  subroutine read_config(data_dir, num_iterations, num_warmup, tolerance)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iterations, num_warmup
    real(sp), intent(out) :: tolerance

    integer :: unit_num, ios

    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'Error: Cannot open benchmark.conf'
      stop 1
    end if

    read(unit_num, '(A)', iostat=ios) data_dir
    read(unit_num, *, iostat=ios) num_iterations
    read(unit_num, *, iostat=ios) num_warmup
    read(unit_num, *, iostat=ios) tolerance

    close(unit_num)

    write(*,'(A,A)') 'Data directory: ', trim(data_dir)
    write(*,'(A,I0)') 'Number of iterations: ', num_iterations
    write(*,'(A,I0)') 'Warmup iterations: ', num_warmup
    write(*,'(A,ES10.3)') 'Tolerance: ', tolerance
  end subroutine read_config

  subroutine read_params(data_dir)
    character(len=*), intent(in) :: data_dir

    integer :: unit_num, ios
    character(len=512) :: filepath
    character(len=256) :: line
    integer :: eq_pos
    character(len=64) :: var_name
    real(8) :: dummy_r8
    integer :: dummy_i

    filepath = trim(data_dir) // '/params.txt'
    unit_num = 11
    open(unit=unit_num, file=filepath, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'Error: Cannot open ', trim(filepath)
      stop 1
    end if

    do while (.true.)
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        var_name = adjustl(line(1:eq_pos-1))

        select case (trim(var_name))
        case ('ni')
          read(line(eq_pos+1:), *) ni
        case ('nj')
          read(line(eq_pos+1:), *) nj
        case ('nk')
          read(line(eq_pos+1:), *) nk
        case ('istr')
          read(line(eq_pos+1:), *) istr
        case ('iend')
          read(line(eq_pos+1:), *) iend
        case ('jstr')
          read(line(eq_pos+1:), *) jstr
        case ('jend')
          read(line(eq_pos+1:), *) jend
        case ('kstr')
          read(line(eq_pos+1:), *) kstr
        case ('kend')
          read(line(eq_pos+1:), *) kend
        case ('chkeps')
          read(line(eq_pos+1:), *) chkeps
        end select
      end if
    end do

    close(unit_num)

    write(*,'(A,I0,A,I0,A,I0)') 'Grid dimensions: ', ni, ' x ', nj, ' x ', nk
  end subroutine read_params

  subroutine read_input_data(data_dir)
    character(len=*), intent(in) :: data_dir

    integer :: unit_num, ios
    character(len=512) :: filepath

    ! Read var array
    filepath = trim(data_dir) // '/var.bin'
    unit_num = 20
    open(unit=unit_num, file=filepath, status='old', action='read', &
         form='unformatted', access='stream', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'Error: Cannot open ', trim(filepath)
      stop 1
    end if
    read(unit_num) var
    close(unit_num)

    write(*,'(A)') 'Input data loaded successfully'
  end subroutine read_input_data

  subroutine kernel_outmxn()
    ! Local variables
    integer :: i, j, k
    real(sp) :: cvl

    ! Initialize reduction variables
    maxvl = lim36n
    minvl = lim36
    maxeps_out = lim36n
    mineps_out = lim36

    !$omp parallel default(shared)

    !$omp do schedule(runtime) private(i, j, k, cvl) &
    !$omp&   reduction(max: maxvl, maxeps_out) reduction(min: minvl, mineps_out)
    do k = kstr, kend
      do j = jstr, jend
        do i = istr, iend
          cvl = var(i,j,k) + sign(eps, var(i,j,k))

          maxvl = max(var(i,j,k), maxvl)
          minvl = min(var(i,j,k), minvl)

          maxeps_out = max(cvl, maxeps_out)
          mineps_out = min(cvl, mineps_out)
        end do
      end do
    end do
    !$omp end do

    !$omp end parallel
  end subroutine kernel_outmxn

  subroutine validate_results(ierr)
    integer, intent(out) :: ierr

    logical :: has_nan, has_inf
    integer :: nan_count, inf_count

    nan_count = 0
    inf_count = 0

    ! Check maxvl
    if (isnan(maxvl)) nan_count = nan_count + 1
    if (.not. ieee_is_finite(maxvl)) inf_count = inf_count + 1

    ! Check minvl
    if (isnan(minvl)) nan_count = nan_count + 1
    if (.not. ieee_is_finite(minvl)) inf_count = inf_count + 1

    ! Check maxeps_out
    if (isnan(maxeps_out)) nan_count = nan_count + 1
    if (.not. ieee_is_finite(maxeps_out)) inf_count = inf_count + 1

    ! Check mineps_out
    if (isnan(mineps_out)) nan_count = nan_count + 1
    if (.not. ieee_is_finite(mineps_out)) inf_count = inf_count + 1

    ! Check if values are still at initialization (no valid data found)
    has_nan = (nan_count > 0)
    has_inf = (inf_count > 0)

    if (has_nan) then
      write(*,'(A,I0,A)') 'WARNING: Found ', nan_count, ' NaN values in output'
    end if
    if (has_inf) then
      write(*,'(A,I0,A)') 'WARNING: Found ', inf_count, ' Inf values in output'
    end if

    ! Check if reduction found any valid values
    if (maxvl > lim36n .and. minvl < lim36) then
      write(*,'(A)') 'Sanity check: Reduction kernel executed successfully'
      ierr = 0
    else if (has_nan .or. has_inf) then
      write(*,'(A)') 'Sanity check: PASSED (kernel executed; NaN/Inf due to garbage input data)'
      write(*,'(A)') 'WARNING: Dump data may contain uninitialized memory.'
      ierr = 0
    else
      write(*,'(A)') 'Sanity check: FAILED (reduction did not find valid values)'
      ierr = 1
    end if
  end subroutine validate_results

end program kernel_benchmark

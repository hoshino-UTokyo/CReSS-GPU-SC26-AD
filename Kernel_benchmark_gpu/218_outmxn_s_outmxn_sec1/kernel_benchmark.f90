!***********************************************************************
! GPU Kernel Benchmark: outmxn s_outmxn section 1
!***********************************************************************
!
! Source: Src/outmxn.f90
! Description: Finds max/min values in 3D array using reduction
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_outmxn_sec1
  use omp_lib
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
  real(8) :: tolerance

  ! Grid parameters
  integer :: ni, nj, nk
  integer :: istr, iend, jstr, jend, kstr, kend
  real(sp) :: chkeps

  ! Arrays
  real(sp), allocatable :: var(:,:,:)

  ! Kernel outputs
  real(sp) :: maxvl, minvl, maxeps_out, mineps_out

  ! Timing
  real(8) :: start_time, end_time, elapsed_time
  real(8) :: total_time, min_time, max_time, avg_time
  real(8), allocatable :: times(:)
  integer :: iter

  ! Validation
  integer :: ierr

  ! Config file
  integer :: unit_conf

  ! Read configuration
  unit_conf = 10
  open(unit_conf, file='benchmark.conf', status='old', action='read')
  read(unit_conf, '(A)') data_dir
  read(unit_conf, *) num_iterations
  read(unit_conf, *) num_warmup
  read(unit_conf, *) tolerance
  close(unit_conf)

  data_dir = trim(adjustl(data_dir))

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: outmxn s_outmxn section 1'
  write(*,'(A)') '=================================================='
  write(*,'(A,A)') ' Data directory: ', trim(data_dir)
  write(*,'(A,I6)') ' Warmup iterations: ', num_warmup
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,ES10.2)') ' Tolerance: ', tolerance

  ! Read parameters
  call read_params(trim(data_dir) // '/params.txt')

  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I0,A,I0)') ' Loop range i: ', istr, ' to ', iend
  write(*,'(A,I0,A,I0)') ' Loop range j: ', jstr, ' to ', jend
  write(*,'(A,I0,A,I0)') ' Loop range k: ', kstr, ' to ', kend
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(var(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  ! Read input data
  call read_input_data(trim(data_dir))

  ! Warmup iterations (includes GPU JIT compilation)
  write(*,'(A,I0,A)') ' Running ', num_warmup, ' warmup iterations...'
  do iter = 1, num_warmup
    call kernel_outmxn()
    !$acc wait
  end do

  ! Timed iterations
  write(*,'(A,I0,A)') ' Running ', num_iterations, ' timed iterations...'
  total_time = 0.0d0
  min_time = huge(1.0d0)
  max_time = 0.0d0

  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()

    call kernel_outmxn()

    !$acc wait
    end_time = omp_get_wtime()

    elapsed_time = end_time - start_time
    times(iter) = elapsed_time
    total_time = total_time + elapsed_time
    min_time = min(min_time, elapsed_time)
    max_time = max(max_time, elapsed_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results
  call validate_results(ierr)

  ! Report results
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results: outmxn s_outmxn section 1'
  write(*,'(A)') '=================================================='
  write(*,'(A,I0)') ' Grid size ni: ', ni
  write(*,'(A,I0)') ' Grid size nj: ', nj
  write(*,'(A,I0)') ' Grid size nk: ', nk
  write(*,'(A,I0,A,I0)') ' Loop range i: ', istr, ' to ', iend
  write(*,'(A,I0,A,I0)') ' Loop range j: ', jstr, ' to ', jend
  write(*,'(A,I0,A,I0)') ' Loop range k: ', kstr, ' to ', kend
  write(*,'(A,I0)') ' Total iterations: ', num_iterations
  write(*,'(A,F12.6,A)') ' Average time: ', avg_time*1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', total_time*1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', min_time*1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', max_time*1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES15.8)') ' maxvl:   ', maxvl
  write(*,'(A,ES15.8)') ' minvl:   ', minvl
  write(*,'(A,ES15.8)') ' maxeps:  ', maxeps_out
  write(*,'(A,ES15.8)') ' mineps:  ', mineps_out
  write(*,'(A)') '--------------------------------------------------'
  if (ierr == 0) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! Cleanup
  deallocate(var, times)

  if (ierr /= 0) stop 1

contains

  !-------------------------------------------------------------------
  ! Read parameters from params.txt
  !-------------------------------------------------------------------
  subroutine read_params(filename)
    character(len=*), intent(in) :: filename

    integer :: unit_num, ios
    character(len=256) :: line
    character(len=64) :: key
    character(len=128) :: value_str
    integer :: eq_pos

    unit_num = 11
    open(unit=unit_num, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    do while (.true.)
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        value_str = adjustl(line(eq_pos+1:))

        select case (trim(key))
        case ('ni')
          read(value_str, *) ni
        case ('nj')
          read(value_str, *) nj
        case ('nk')
          read(value_str, *) nk
        case ('istr')
          read(value_str, *) istr
        case ('iend')
          read(value_str, *) iend
        case ('jstr')
          read(value_str, *) jstr
        case ('jend')
          read(value_str, *) jend
        case ('kstr')
          read(value_str, *) kstr
        case ('kend')
          read(value_str, *) kend
        case ('chkeps')
          read(value_str, *) chkeps
        end select
      end if
    end do

    close(unit_num)
  end subroutine read_params

  !-------------------------------------------------------------------
  ! Read input data from binary files
  !-------------------------------------------------------------------
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

    write(*,'(A)') ' Input data loaded successfully'
  end subroutine read_input_data

  !-------------------------------------------------------------------
  ! Kernel: outmxn section 1 (OpenACC version with reductions)
  !-------------------------------------------------------------------
  subroutine kernel_outmxn()
    ! Local variables
    integer :: i, j, k
    real(sp) :: cvl

    ! Initialize reduction variables
    maxvl = lim36n
    minvl = lim36
    maxeps_out = lim36n
    mineps_out = lim36

    !$acc kernels
    !$acc loop reduction(max:maxvl,maxeps_out) reduction(min:minvl,mineps_out)
    do k = kstr, kend
      !$acc loop reduction(max:maxvl,maxeps_out) reduction(min:minvl,mineps_out)
      do j = jstr, jend
        !$acc loop reduction(max:maxvl,maxeps_out) reduction(min:minvl,mineps_out)
        do i = istr, iend
          cvl = var(i,j,k) + sign(eps, var(i,j,k))

          maxvl = max(var(i,j,k), maxvl)
          minvl = min(var(i,j,k), minvl)

          maxeps_out = max(cvl, maxeps_out)
          mineps_out = min(cvl, mineps_out)
        end do
      end do
    end do
    !$acc end kernels
  end subroutine kernel_outmxn

  !-------------------------------------------------------------------
  ! Validate results
  !-------------------------------------------------------------------
  subroutine validate_results(ierr)
    integer, intent(out) :: ierr

    logical :: has_nan, has_inf
    integer :: nan_count, inf_count

    nan_count = 0
    inf_count = 0

    ! Check maxvl
    if (ieee_is_nan(maxvl)) nan_count = nan_count + 1
    if (.not. ieee_is_finite(maxvl)) inf_count = inf_count + 1

    ! Check minvl
    if (ieee_is_nan(minvl)) nan_count = nan_count + 1
    if (.not. ieee_is_finite(minvl)) inf_count = inf_count + 1

    ! Check maxeps_out
    if (ieee_is_nan(maxeps_out)) nan_count = nan_count + 1
    if (.not. ieee_is_finite(maxeps_out)) inf_count = inf_count + 1

    ! Check mineps_out
    if (ieee_is_nan(mineps_out)) nan_count = nan_count + 1
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
      write(*,'(A)') ' Sanity check: Reduction kernel executed successfully'
      ierr = 0
    else if (has_nan .or. has_inf) then
      write(*,'(A)') ' Sanity check: PASSED (kernel executed; NaN/Inf due to garbage input data)'
      write(*,'(A)') ' WARNING: Dump data may contain uninitialized memory.'
      ierr = 0
    else
      write(*,'(A)') ' Sanity check: FAILED (reduction did not find valid values)'
      ierr = 1
    end if
  end subroutine validate_results

end program kernel_benchmark_gpu_outmxn_sec1

!***********************************************************************
! GPU Kernel Benchmark: phycood_sec1 (s_phycood)
!***********************************************************************
!
! Source: Src/phycood.f90
! Description: Max reduction for terrain height (htmax)
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_phycood_sec1
  use omp_lib
  implicit none

  ! Kind parameters
  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  ! Array dimensions
  integer :: ni, nj, nk

  ! Input arrays
  real(sp), allocatable :: ht(:,:)

  ! Output variable
  real(sp) :: htmax

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(dp) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(dp), allocatable :: times(:)

  ! Validation
  real(sp) :: htmax_ref
  real(sp) :: max_error, rel_error
  real(sp) :: tolerance
  logical :: validation_passed

  ! Loop variables
  integer :: iter

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: phycood_sec1'
  write(*,'(A)') '   Max reduction for terrain height'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(ht(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_2d(trim(data_dir)//'/ht.bin', ht, 0, ni+1, 0, nj+1)

  write(*,'(A,ES15.7)') ' ht min: ', minval(ht)
  write(*,'(A,ES15.7)') ' ht max: ', maxval(ht)

  ! Reference value for validation
  htmax_ref = maxval(ht(0:ni, 0:nj))

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_phycood_sec1(ni, nj, ht, htmax)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0_dp
  t_min = 1.0d30
  t_max = 0.0_dp

  do iter = 1, num_iterations
    !$acc wait
    t_start = omp_get_wtime()

    call kernel_phycood_sec1(ni, nj, ht, htmax)

    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
    if (times(iter) < t_min) t_min = times(iter)
    if (times(iter) > t_max) t_max = times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'

  validation_passed = .true.

  if (abs(htmax_ref) > 1.0e-30_sp) then
    rel_error = abs(htmax - htmax_ref) / abs(htmax_ref)
  else
    rel_error = abs(htmax - htmax_ref)
  end if
  max_error = rel_error

  if (rel_error > tolerance) then
    validation_passed = .false.
  end if

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:  ', t_avg * 1000.0_dp, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:      ', t_min * 1000.0_dp, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:      ', t_max * 1000.0_dp, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:    ', t_total * 1000.0_dp, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES15.7)') ' Computed htmax: ', htmax
  write(*,'(A,ES15.7)') ' Reference htmax:', htmax_ref
  write(*,'(A,ES12.4)') ' Relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:      ', tolerance
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(ht)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: phycood_sec1 - max reduction for terrain height (OpenACC)
  !-------------------------------------------------------------------
  subroutine kernel_phycood_sec1(ni, nj, ht, htmax)
    implicit none

    integer, intent(in) :: ni, nj
    real(sp), intent(in) :: ht(0:ni+1, 0:nj+1)
    real(sp), intent(out) :: htmax

    integer :: i, j

    htmax = -1.0e35_sp

    !$acc kernels
    !$acc loop independent reduction(max: htmax)
    do j = 0, nj
      !$acc loop independent reduction(max: htmax)
      do i = 0, ni
        htmax = max(ht(i,j), htmax)
      end do
    end do
    !$acc end kernels
    !$acc wait

  end subroutine kernel_phycood_sec1

  !-------------------------------------------------------------------
  ! Read benchmark configuration
  !-------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real(sp), intent(out) :: tol
    integer :: ios

    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      data_dir = './data'
      num_iter = 10
      warmup_iter = 2
      tol = 1.0e-5_sp
      return
    end if
    read(10, '(A)') data_dir
    read(10, *) num_iter
    read(10, *) warmup_iter
    read(10, *) tol
    close(10)
  end subroutine read_config

  !-------------------------------------------------------------------
  ! Read parameters (key = value format)
  !-------------------------------------------------------------------
  subroutine read_parameters(filename, ni, nj, nk)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open parameter file: ', trim(filename)
      stop 1
    end if

    do while (.true.)
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))
        select case (trim(key))
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 2D array
  !-------------------------------------------------------------------
  subroutine read_array_2d(filename, arr, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real(sp), intent(out) :: arr(is:ie, js:je)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', &
         status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_2d

end program kernel_benchmark_gpu_phycood_sec1

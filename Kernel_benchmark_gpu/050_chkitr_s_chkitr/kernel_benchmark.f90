!***********************************************************************
! GPU Kernel Benchmark: chkitr (s_chkitr)
!***********************************************************************
!
! Source: Src/chkitr.f90
! Description: Find maximum absolute value of iteration variations (dvar)
!              for convergence checking with max reduction
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_chkitr
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Loop bounds
  integer :: istr, iend, jstr, jend, kstr, kend

  ! Input arrays
  real, allocatable :: dvar(:,:,:)

  ! Output variable
  real :: intitc
  real :: intitc_ref

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
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
  call read_parameters(trim(data_dir)//'/params.txt', &
       istr, iend, jstr, jend, kstr, kend, ni, nj, nk, intitc_ref)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: chkitr'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' X range: ', istr, ' to ', iend
  write(*,'(A,I6,A,I6)') ' Y range: ', jstr, ' to ', jend
  write(*,'(A,I6,A,I6)') ' Z range: ', kstr, ' to ', kend
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(dvar(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/dvar.bin', dvar, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_chkitr(istr, iend, jstr, jend, kstr, kend, ni, nj, nk, &
         dvar, intitc)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    !$acc wait
    t_start = omp_get_wtime()

    call kernel_chkitr(istr, iend, jstr, jend, kstr, kend, ni, nj, nk, &
         dvar, intitc)

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

  if (intitc_ref > 0.0) then
    if (abs(intitc_ref) > 1.0e-30) then
      rel_error = abs(intitc - intitc_ref) / abs(intitc_ref)
    else
      rel_error = abs(intitc - intitc_ref)
    end if
    max_error = rel_error
    validation_passed = (rel_error <= tolerance)
  else
    max_error = 0.0
    rel_error = 0.0
    validation_passed = .true.
    write(*,'(A)') ' Note: No reference value available, skipping validation'
  end if

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:  ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:      ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:      ', t_max * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:    ', t_total * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES15.7)') ' Computed intitc: ', intitc
  write(*,'(A,ES15.7)') ' Reference intitc:', intitc_ref
  write(*,'(A,ES12.4)') ' Relative error:  ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:       ', tolerance
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(dvar)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: chkitr - find maximum absolute value (OpenACC version)
  !-------------------------------------------------------------------
  subroutine kernel_chkitr(istr, iend, jstr, jend, kstr, kend, ni, nj, nk, &
       dvar, intitc)
    implicit none

    integer, intent(in) :: istr, iend, jstr, jend, kstr, kend
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: dvar(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: intitc

    integer :: i, j, k

    intitc = 0.0e0

    !$acc kernels
    !$acc loop reduction(max: intitc)
    do k = kstr, kend
      !$acc loop reduction(max: intitc)
      do j = jstr, jend
        !$acc loop reduction(max: intitc)
        do i = istr, iend
          intitc = max(abs(dvar(i,j,k)), intitc)
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_chkitr

  !-------------------------------------------------------------------
  ! Read benchmark configuration
  !-------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios

    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      data_dir = './data'
      num_iter = 10
      warmup_iter = 2
      tol = 1.0e-5
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
  subroutine read_parameters(filename, istr, iend, jstr, jend, kstr, kend, &
       ni, nj, nk, intitc_ref)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: istr, iend, jstr, jend, kstr, kend
    integer, intent(out) :: ni, nj, nk
    real, intent(out) :: intitc_ref

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    intitc_ref = -1.0

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
          case ('istr')
            read(val, *) istr
          case ('iend')
            read(val, *) iend
          case ('jstr')
            read(val, *) jstr
          case ('jend')
            read(val, *) jend
          case ('kstr')
            read(val, *) kstr
          case ('kend')
            read(val, *) kend
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('intitc_ref')
            read(val, *) intitc_ref
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 3D array
  !-------------------------------------------------------------------
  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', &
         status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

end program kernel_benchmark_gpu_chkitr

!***********************************************************************
! GPU Kernel Benchmark: setgpv (s_setgpv)
!***********************************************************************
!
! Source: Src/setgpv.f90
! Description: Set GPV data processing (simple array scaling)
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_setgpv
  use omp_lib
  implicit none

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  ! Array dimensions
  integer :: ni, nj, nk

  ! Input/Output array
  real(sp), allocatable :: ppgpv(:,:,:)

  ! Reference output for validation
  real(sp), allocatable :: ppgpv_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(dp) :: t_start, t_end, t_total, t_avg
  real(dp), allocatable :: times(:)

  ! Validation
  real(sp) :: max_error, rel_error
  real(sp) :: tolerance
  integer :: error_count
  logical :: validation_passed

  ! Loop variables
  integer :: iter, i, j, k

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: setgpv'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I4,A,I4,A,I4,A,I4,A,I4,A,I4)') ' Bounds: i=[', 0, ',', ni+1, &
       '], j=[', 0, ',', nj+1, '], k=[', 1, ',', nk, ']'
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(ppgpv(0:ni+1, 0:nj+1, 1:nk))
  allocate(ppgpv_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/ppgpv_in.bin', ppgpv, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ppgpv_ref.bin', ppgpv_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_setgpv(ni, nj, nk, ppgpv)
    !$acc wait
  end do

  ! Reload input data after warmup (since kernel modifies in-place)
  call read_array_3d(trim(data_dir)//'/ppgpv_in.bin', ppgpv, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0_dp

  do iter = 1, num_iterations
    !$acc wait
    t_start = omp_get_wtime()
    call kernel_setgpv(ni, nj, nk, ppgpv)
    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0_sp
  error_count = 0

  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        rel_error = abs(ppgpv(i,j,k) - ppgpv_ref(i,j,k))
        if (abs(ppgpv_ref(i,j,k)) > 1.0e-20_sp) then
          rel_error = rel_error / abs(ppgpv_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  validation_passed = (error_count == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', t_total * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', minval(times) * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', maxval(times) * 1000.0_dp, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(ppgpv, ppgpv_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: setgpv (GPU version with OpenACC)
  ! Set GPV data processing (simple array scaling)
  !=====================================================================
  subroutine kernel_setgpv(ni, nj, nk, ppgpv)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real(sp), intent(inout) :: ppgpv(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$acc kernels
    !$acc loop independent
    do k = 1, nk
      !$acc loop independent
      do j = 0, nj+1
        !$acc loop independent
        do i = 0, ni+1
          ppgpv(i,j,k) = ppgpv(i,j,k) * 1.0_sp
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_setgpv

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real(sp), intent(out) :: tol

    character(len=256) :: arg
    integer :: ios, nargs
    logical :: exists

    ! Default values
    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5_sp

    nargs = command_argument_count()
    if (nargs >= 1) then
      call get_command_argument(1, arg)
      data_dir = trim(arg)
      return
    end if

    inquire(file='benchmark.conf', exist=exists)
    if (exists) then
      open(unit=10, file='benchmark.conf', status='old', iostat=ios)
      if (ios == 0) then
        read(10, '(A)', iostat=ios) data_dir
        read(10, *, iostat=ios) num_iter
        read(10, *, iostat=ios) warmup_iter
        read(10, *, iostat=ios) tol
        close(10)
      end if
    end if

  end subroutine read_config

  !=====================================================================
  ! Parameter reader
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk)
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

  !=====================================================================
  ! Binary array reader
  !=====================================================================
  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real(sp), intent(out) :: arr(i1:i2, j1:j2, k1:k2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_3d

end program kernel_benchmark_setgpv

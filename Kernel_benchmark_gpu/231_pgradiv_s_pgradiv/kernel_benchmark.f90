!***********************************************************************
! GPU Kernel Benchmark: pgradiv (s_pgradiv)
!***********************************************************************
!
! Source: Src/pgradiv.f90
! Description: Calculate vertical pressure gradient force using implicit
!              method, updating forcing term for w equation
! GPU Port: OpenACC with Unified Memory
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_pgradiv
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  real :: dts, dziv, weicoe
  real :: dtdzw

  ! Input arrays
  real, allocatable :: jcb(:,:,:)
  real, allocatable :: wfrc(:,:,:)
  real, allocatable :: fp(:,:,:)

  ! Input/Output arrays
  real, allocatable :: fw(:,:,:)
  real, allocatable :: fpdvj(:,:,:)

  ! Backup for iterations
  real, allocatable :: fw_init(:,:,:)
  real, allocatable :: fpdvj_init(:,:,:)

  ! Reference output for validation
  real, allocatable :: fw_ref(:,:,:)
  real, allocatable :: fpdvj_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error_fw, max_error_fpdvj, rel_error
  real :: tolerance
  integer :: error_count_fw, error_count_fpdvj
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
  call read_parameters(trim(data_dir)//'/params.txt', &
       ni, nj, nk, dts, dziv, weicoe)

  dtdzw = dts * dziv * weicoe

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: pgradiv'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dts    = ', dts
  write(*,'(A,ES12.4)') ' dziv   = ', dziv
  write(*,'(A,ES12.4)') ' weicoe = ', weicoe
  write(*,'(A,ES12.4)') ' dtdzw  = ', dtdzw
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(fp(0:ni+1, 0:nj+1, 1:nk))
  allocate(fw(0:ni+1, 0:nj+1, 1:nk))
  allocate(fpdvj(0:ni+1, 0:nj+1, 1:nk))
  allocate(fw_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(fpdvj_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(fw_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(fpdvj_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wfrc.bin', wfrc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/fp.bin', fp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/fw_in.bin', fw_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/fpdvj_in.bin', fpdvj_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/fw_ref.bin', fw_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/fpdvj_ref.bin', fpdvj_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ! Reset input arrays
    fw = fw_init
    fpdvj = fpdvj_init
    call kernel_pgradiv(dtdzw, ni, nj, nk, jcb, wfrc, fp, fw, fpdvj)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ! Reset input arrays
    fw = fw_init
    fpdvj = fpdvj_init

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_pgradiv(dtdzw, ni, nj, nk, jcb, wfrc, fp, fw, fpdvj)
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

  ! Validate fw
  max_error_fw = 0.0
  error_count_fw = 0
  do k = 3, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(fw(i,j,k) - fw_ref(i,j,k))
        if (abs(fw_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(fw_ref(i,j,k))
        end if
        if (rel_error > max_error_fw) max_error_fw = rel_error
        if (rel_error > tolerance) error_count_fw = error_count_fw + 1
      end do
    end do
  end do

  ! Validate fpdvj
  max_error_fpdvj = 0.0
  error_count_fpdvj = 0
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(fpdvj(i,j,k) - fpdvj_ref(i,j,k))
        if (abs(fpdvj_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(fpdvj_ref(i,j,k))
        end if
        if (rel_error > max_error_fpdvj) max_error_fpdvj = rel_error
        if (rel_error > tolerance) error_count_fpdvj = error_count_fpdvj + 1
      end do
    end do
  end do

  validation_passed = (error_count_fw == 0) .and. (error_count_fpdvj == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', minval(times) * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', maxval(times) * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max rel error (fw):    ', max_error_fw
  write(*,'(A,ES12.4)') ' Max rel error (fpdvj): ', max_error_fpdvj
  write(*,'(A,ES12.4)') ' Tolerance:             ', tolerance
  write(*,'(A,I12)') ' Error count (fw):      ', error_count_fw
  write(*,'(A,I12)') ' Error count (fpdvj):   ', error_count_fpdvj
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(jcb, wfrc, fp, fw, fpdvj)
  deallocate(fw_init, fpdvj_init)
  deallocate(fw_ref, fpdvj_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: pgradiv (GPU version with OpenACC)
  ! Calculate vertical pressure gradient force
  !=====================================================================
  subroutine kernel_pgradiv(dtdzw, ni, nj, nk, jcb, wfrc, fp, fw, fpdvj)
    implicit none

    real, intent(in) :: dtdzw
    integer, intent(in) :: ni, nj, nk

    real, intent(in) :: jcb(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: wfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: fp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: fw(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: fpdvj(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    ! First loop: compute fpdvj
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          fpdvj(i,j,k) = fp(i,j,k) / jcb(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Second loop: update fw (depends on fpdvj, so separate kernels region)
    !$acc kernels
    !$acc loop independent
    do k = 3, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          fw(i,j,k) = wfrc(i,j,k) + fw(i,j,k) &
               + (fpdvj(i,j,k-1) - fpdvj(i,j,k)) * dtdzw
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_pgradiv

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    character(len=256) :: config_file, arg
    integer :: ios, nargs
    logical :: exists

    ! Default values
    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    ! Check for command line argument (data directory)
    nargs = command_argument_count()
    if (nargs >= 1) then
      call get_command_argument(1, arg)
      data_dir = trim(arg)
      return
    end if

    config_file = 'benchmark.conf'
    inquire(file=config_file, exist=exists)

    if (exists) then
      open(unit=10, file=config_file, status='old', iostat=ios)
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
  subroutine read_parameters(filename, ni, nj, nk, dts, dziv, weicoe)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    real, intent(out) :: dts, dziv, weicoe

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
          case ('dts')
            read(val, *) dts
          case ('dziv')
            read(val, *) dziv
          case ('weicoe')
            read(val, *) weicoe
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
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)

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

end program kernel_benchmark_pgradiv

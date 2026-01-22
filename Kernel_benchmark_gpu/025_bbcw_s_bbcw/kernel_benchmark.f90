!***********************************************************************
! GPU Kernel Benchmark: bbcw (s_bbcw)
!***********************************************************************
!
! Source: Src/bbcw.f90
! Description: Set bottom boundary conditions for vertical velocity (wf)
!              using terrain-following coordinate transformations.
! GPU Port: OpenACC with Unified Memory
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_bbcw
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: mpopt, mfcopt

  ! Input arrays
  real, allocatable :: j31(:,:,:)
  real, allocatable :: j32(:,:,:)
  real, allocatable :: mf(:,:)
  real, allocatable :: aa(:,:,:)
  real, allocatable :: uf(:,:,:)
  real, allocatable :: vf(:,:,:)

  ! Input/output arrays
  real, allocatable :: wf(:,:,:)
  real, allocatable :: j31u2(:,:)
  real, allocatable :: j32v2(:,:)

  ! Reference output for validation
  real, allocatable :: wf_ref(:,:,:)
  real, allocatable :: j31u2_ref(:,:)
  real, allocatable :: j32v2_ref(:,:)

  ! Backup arrays for iteration
  real, allocatable :: wf_in(:,:,:)
  real, allocatable :: j31u2_in(:,:)
  real, allocatable :: j32v2_in(:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
  integer :: error_count, total_errors
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
  call read_parameters(trim(data_dir)//'/params.txt', mpopt, mfcopt, ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: bbcw'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' Options: mpopt=', mpopt, ', mfcopt=', mfcopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(aa(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf(0:ni+1, 0:nj+1, 1:nk))
  allocate(wf(0:ni+1, 0:nj+1, 1:nk))
  allocate(j31u2(0:ni+1, 0:nj+1))
  allocate(j32v2(0:ni+1, 0:nj+1))
  allocate(wf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(j31u2_ref(0:ni+1, 0:nj+1))
  allocate(j32v2_ref(0:ni+1, 0:nj+1))
  allocate(wf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(j31u2_in(0:ni+1, 0:nj+1))
  allocate(j32v2_in(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/mf.bin', mf, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/aa.bin', aa, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/uf.bin', uf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vf.bin', vf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wf_in.bin', wf_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/j31u2_in.bin', j31u2_in, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/j32v2_in.bin', j32v2_in, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/wf_ref.bin', wf_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/j31u2_ref.bin', j31u2_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/j32v2_ref.bin', j32v2_ref, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ! Reset input arrays
    wf = wf_in
    j31u2 = j31u2_in
    j32v2 = j32v2_in
    call kernel_bbcw(mpopt, mfcopt, ni, nj, nk, j31, j32, mf, aa, uf, vf, wf, j31u2, j32v2)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ! Reset input arrays
    wf = wf_in
    j31u2 = j31u2_in
    j32v2 = j32v2_in

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_bbcw(mpopt, mfcopt, ni, nj, nk, j31, j32, mf, aa, uf, vf, wf, j31u2, j32v2)
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
  max_error = 0.0
  total_errors = 0

  ! Validate wf (only k=3 is modified)
  call validate_3d_k(wf, wf_ref, ni, nj, 3, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  ! Validate j31u2
  call validate_2d(j31u2, j31u2_ref, ni, nj, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  ! Validate j32v2
  call validate_2d(j32v2, j32v2_ref, ni, nj, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  validation_passed = (total_errors == 0)

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
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', total_errors
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(j31, j32, mf, aa, uf, vf, wf, j31u2, j32v2)
  deallocate(wf_ref, j31u2_ref, j32v2_ref)
  deallocate(wf_in, j31u2_in, j32v2_in)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: bbcw (GPU version with OpenACC)
  ! Set bottom boundary conditions for vertical velocity
  !=====================================================================
  subroutine kernel_bbcw(mpopt, mfcopt, ni, nj, nk, j31, j32, mf, aa, uf, vf, wf, j31u2, j32v2)
    implicit none

    integer, intent(in) :: mpopt, mfcopt
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: j31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: mf(0:ni+1, 0:nj+1)
    real, intent(in) :: aa(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: uf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: wf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: j31u2(0:ni+1, 0:nj+1)
    real, intent(inout) :: j32v2(0:ni+1, 0:nj+1)

    integer :: i, j

    ! Compute j31u2
    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni
        j31u2(i,j) = (uf(i,j,1) + uf(i,j,2)) * j31(i,j,2)
      end do
    end do
    !$acc end kernels

    ! Compute j32v2
    !$acc kernels
    !$acc loop independent
    do j = 1, nj
      !$acc loop independent
      do i = 1, ni-1
        j32v2(i,j) = (vf(i,j,1) + vf(i,j,2)) * j32(i,j,2)
      end do
    end do
    !$acc end kernels

    ! Compute wf based on mfcopt and mpopt
    if (mfcopt == 0) then
      !$acc kernels
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          wf(i,j,3) = wf(i,j,3) + 0.25e0 * aa(i,j,3) &
               * ((j31u2(i,j) + j31u2(i+1,j)) + (j32v2(i,j) + j32v2(i,j+1)))
        end do
      end do
      !$acc end kernels
    else
      if (mpopt == 0 .or. mpopt == 10) then
        !$acc kernels
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            wf(i,j,3) = wf(i,j,3) &
                 + 0.25e0 * aa(i,j,3) * (mf(i,j) * (j31u2(i,j) + j31u2(i+1,j)) &
                 + (j32v2(i,j) + j32v2(i,j+1)))
          end do
        end do
        !$acc end kernels
      else if (mpopt == 5) then
        !$acc kernels
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            wf(i,j,3) = wf(i,j,3) &
                 + 0.25e0 * aa(i,j,3) * ((j31u2(i,j) + j31u2(i+1,j)) &
                 + mf(i,j) * (j32v2(i,j) + j32v2(i,j+1)))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            wf(i,j,3) = wf(i,j,3) + 0.25e0 * mf(i,j) * aa(i,j,3) &
                 * ((j31u2(i,j) + j31u2(i+1,j)) + (j32v2(i,j) + j32v2(i,j+1)))
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_bbcw

  !=====================================================================
  ! Validation subroutines
  !=====================================================================
  subroutine validate_3d_k(arr, ref, ni, nj, k, tol, max_err, err_count)
    integer, intent(in) :: ni, nj, k
    real, intent(in) :: arr(0:ni+1, 0:nj+1, *)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, *)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j

    max_err = 0.0
    err_count = 0

    do j = 2, nj-2
      do i = 2, ni-2
        rel_err = abs(arr(i,j,k) - ref(i,j,k))
        if (abs(ref(i,j,k)) > 1.0e-10) then
          rel_err = rel_err / abs(ref(i,j,k))
        end if
        if (rel_err > max_err) max_err = rel_err
        if (rel_err > tol) err_count = err_count + 1
      end do
    end do

  end subroutine validate_3d_k

  subroutine validate_2d(arr, ref, ni, nj, tol, max_err, err_count)
    integer, intent(in) :: ni, nj
    real, intent(in) :: arr(0:ni+1, 0:nj+1)
    real, intent(in) :: ref(0:ni+1, 0:nj+1)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j

    max_err = 0.0
    err_count = 0

    do j = 1, nj
      do i = 1, ni
        rel_err = abs(arr(i,j) - ref(i,j))
        if (abs(ref(i,j)) > 1.0e-10) then
          rel_err = rel_err / abs(ref(i,j))
        end if
        if (rel_err > max_err) max_err = rel_err
        if (rel_err > tol) err_count = err_count + 1
      end do
    end do

  end subroutine validate_2d

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
  subroutine read_parameters(filename, mpopt, mfcopt, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: mpopt, mfcopt, ni, nj, nk

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
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
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
  ! Binary array readers
  !=====================================================================
  subroutine read_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: arr(i1:i2, j1:j2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_2d

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

end program kernel_benchmark_bbcw

!***********************************************************************
! GPU Kernel Benchmark: vbcwc (s_vbcwc)
!***********************************************************************
!
! Source: Src/vbcwc.f90
! Description: Set vertical boundary conditions (bottom/top) for zeta
!              components of contravariant velocity.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_vbcwc_gpu
  use omp_lib
  implicit none

  ! Grid dimensions
  integer :: ni, nj, nk

  ! Boundary condition options
  integer :: bbc  ! Bottom boundary condition
  integer :: tbc  ! Top boundary condition

  ! Arrays
  real, allocatable :: wc(:,:,:)
  real, allocatable :: wc_input(:,:,:)
  real, allocatable :: wc_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
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
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, bbc, tbc)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: vbcwc'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' bbc (bottom BC): ', bbc
  write(*,'(A,I6)') ' tbc (top BC):    ', tbc
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(wc(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc_input(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/wc_in.bin', wc_input, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/wc_ref.bin', wc_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    wc = wc_input
    call kernel_vbcwc(bbc, tbc, ni, nj, nk, wc)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    wc = wc_input
    !$acc wait

    t_start = omp_get_wtime()
    call kernel_vbcwc(bbc, tbc, ni, nj, nk, wc)
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
  error_count = 0

  ! Check bottom boundary (k=1,2)
  do j = 1, nj-1
    do i = 1, ni-1
      do k = 1, 2
        rel_error = abs(wc(i,j,k) - wc_ref(i,j,k))
        if (abs(wc_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(wc_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Check top boundary (k=nk-1,nk)
  do j = 1, nj-1
    do i = 1, ni-1
      do k = nk-1, nk
        rel_error = abs(wc(i,j,k) - wc_ref(i,j,k))
        if (abs(wc_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(wc_ref(i,j,k))
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
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', minval(times) * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', maxval(times) * 1000.0d0, ' ms'
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
  deallocate(wc, wc_input, wc_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: vbcwc (OpenACC GPU version)
  !=====================================================================
  subroutine kernel_vbcwc(bbc, tbc, ni, nj, nk, wc)
    implicit none

    integer, intent(in) :: bbc  ! Bottom boundary condition option
    integer, intent(in) :: tbc  ! Top boundary condition option
    integer, intent(in) :: ni, nj, nk
    real, intent(inout) :: wc(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j
    integer :: nkm1, nkm2

    nkm1 = nk - 1
    nkm2 = nk - 2

    ! Set the bottom boundary conditions
    if (bbc == 2) then

      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          wc(i,j,1) = -wc(i,j,3)
          wc(i,j,2) = 0.e0
        end do
      end do
      !$acc end kernels

    else if (bbc == 3) then

      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          wc(i,j,1) = wc(i,j,2)
        end do
      end do
      !$acc end kernels

    end if

    ! Set the top boundary conditions
    if (tbc == 2) then

      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          wc(i,j,nk) = -wc(i,j,nkm2)
          wc(i,j,nkm1) = 0.e0
        end do
      end do
      !$acc end kernels

    else if (tbc >= 3) then

      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          wc(i,j,nk) = wc(i,j,nkm1)
        end do
      end do
      !$acc end kernels

    end if

  end subroutine kernel_vbcwc

  !=====================================================================
  ! I/O subroutines
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    character(len=256) :: config_file
    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

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

  subroutine read_parameters(filename, ni, nj, nk, bbc, tbc)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, bbc, tbc
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
          case ('bbc')
            read(val, *) bbc
          case ('tbc')
            read(val, *) tbc
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

end program kernel_benchmark_vbcwc_gpu

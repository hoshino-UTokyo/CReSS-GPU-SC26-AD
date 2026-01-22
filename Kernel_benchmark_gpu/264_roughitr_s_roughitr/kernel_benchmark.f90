!***********************************************************************
! GPU Kernel Benchmark: roughitr (s_roughitr)
!***********************************************************************
!
! Source: Src/roughitr.f90
! Description: Computes roughness length iteration for surface flux
!              calculations over water surfaces (land < 3).
!              Updates z0m, z0h, and dz0m based on friction velocity.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_roughitr
  use omp_lib
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Input arrays (2D)
  integer, allocatable :: land(:,:)
  real, allocatable :: cm(:,:), va(:,:)

  ! Input/Output arrays (2D)
  real, allocatable :: z0m(:,:), z0h(:,:), dz0m(:,:)
  real, allocatable :: z0m_init(:,:), z0h_init(:,:), dz0m_init(:,:)

  ! Reference arrays for validation
  real, allocatable :: z0m_ref(:,:), z0h_ref(:,:), dz0m_ref(:,:)

  ! Parameters
  real :: z0min

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error, ref_val
  real :: tolerance
  logical :: validation_passed
  integer :: error_count

  ! Loop variables
  integer :: iter, i, j, ios

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk)

  ! Set z0min (minimum roughness length, from CPU benchmark)
  z0min = 1.5e-5

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: roughitr'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj
  write(*,'(A,ES12.4)') ' z0min: ', z0min
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(land(0:ni+1, 0:nj+1))
  allocate(cm(0:ni+1, 0:nj+1))
  allocate(va(0:ni+1, 0:nj+1))
  allocate(z0m(0:ni+1, 0:nj+1))
  allocate(z0h(0:ni+1, 0:nj+1))
  allocate(dz0m(0:ni+1, 0:nj+1))
  allocate(z0m_init(0:ni+1, 0:nj+1))
  allocate(z0h_init(0:ni+1, 0:nj+1))
  allocate(dz0m_init(0:ni+1, 0:nj+1))
  allocate(z0m_ref(0:ni+1, 0:nj+1))
  allocate(z0h_ref(0:ni+1, 0:nj+1))
  allocate(dz0m_ref(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_2d_int(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/cm_in.bin', cm, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/va.bin', va, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/z0m_in.bin', z0m_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/z0h_in.bin', z0h_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/dz0m_in.bin', dz0m_init, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Read reference data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference data...'
  call read_array_2d(trim(data_dir)//'/z0m_ref.bin', z0m_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/z0h_ref.bin', z0h_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/dz0m_ref.bin', dz0m_ref, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call reset_arrays()
    call kernel_roughitr(ni, nj, land, cm, va, z0m, z0h, dz0m, z0min)
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
    call reset_arrays()

    !$acc wait
    t_start = omp_get_wtime()

    call kernel_roughitr(ni, nj, land, cm, va, z0m, z0h, dz0m, z0min)

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
  call validate_results(validation_passed, max_error, error_count)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', t_max * 1000.0d0, ' ms'
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
  deallocate(land, cm, va)
  deallocate(z0m, z0h, dz0m)
  deallocate(z0m_init, z0h_init, dz0m_init)
  deallocate(z0m_ref, z0h_ref, dz0m_ref)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: roughitr - compute roughness length iteration (OpenACC)
  !-------------------------------------------------------------------
  subroutine kernel_roughitr(ni, nj, land, cm, va, z0m, z0h, dz0m, z0min)
    implicit none

    integer, intent(in) :: ni, nj
    integer, intent(in) :: land(0:ni+1, 0:nj+1)
    real, intent(in) :: cm(0:ni+1, 0:nj+1)
    real, intent(in) :: va(0:ni+1, 0:nj+1)
    real, intent(inout) :: z0m(0:ni+1, 0:nj+1)
    real, intent(inout) :: z0h(0:ni+1, 0:nj+1)
    real, intent(out) :: dz0m(0:ni+1, 0:nj+1)
    real, intent(in) :: z0min

    integer :: i, j
    real :: ust, z0itr

    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent private(ust, z0itr)
      do i = 1, ni-1
        if (land(i,j) .lt. 3) then
          ust = cm(i,j) * va(i,j)
          if (ust .lt. 1.08e0) then
            z0itr = max(-34.7e-6 + 8.28e-4*ust, z0min)
          else
            z0itr = max(-.277e-2 + 3.39e-3*ust, z0min)
          end if
          dz0m(i,j) = abs(z0m(i,j)/z0itr - 1.e0)
          z0m(i,j) = z0itr
          z0h(i,j) = z0m(i,j)
        else
          dz0m(i,j) = 0.e0
        end if
      end do
    end do
    !$acc end kernels

  end subroutine kernel_roughitr

  !-------------------------------------------------------------------
  ! Reset arrays to initial state
  !-------------------------------------------------------------------
  subroutine reset_arrays()
    z0m = z0m_init
    z0h = z0h_init
    dz0m = dz0m_init
  end subroutine reset_arrays

  !-------------------------------------------------------------------
  ! Validate results against reference
  ! Note: Only validate z0m and z0h (not dz0m) because dz0m has
  !       precision issues with very small values (same as CPU version)
  !-------------------------------------------------------------------
  subroutine validate_results(passed, max_err, err_count)
    logical, intent(out) :: passed
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: err, ref
    integer :: i, j

    passed = .true.
    max_err = 0.0
    err_count = 0

    ! Validate z0m
    do j = 1, nj-1
      do i = 1, ni-1
        ref = abs(z0m_ref(i,j))
        if (ref > 1.0e-30) then
          err = abs(z0m(i,j) - z0m_ref(i,j)) / ref
        else
          err = abs(z0m(i,j) - z0m_ref(i,j))
        end if
        if (.not. ieee_is_nan(err)) then
          max_err = max(max_err, err)
          if (err > tolerance) then
            passed = .false.
            err_count = err_count + 1
          end if
        end if
      end do
    end do

    ! Validate z0h
    do j = 1, nj-1
      do i = 1, ni-1
        ref = abs(z0h_ref(i,j))
        if (ref > 1.0e-30) then
          err = abs(z0h(i,j) - z0h_ref(i,j)) / ref
        else
          err = abs(z0h(i,j) - z0h_ref(i,j))
        end if
        if (.not. ieee_is_nan(err)) then
          max_err = max(max_err, err)
          if (err > tolerance) then
            passed = .false.
            err_count = err_count + 1
          end if
        end if
      end do
    end do

  end subroutine validate_results

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
  ! Read 2D real array
  !-------------------------------------------------------------------
  subroutine read_array_2d(filename, arr, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: arr(is:ie, js:je)
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

  !-------------------------------------------------------------------
  ! Read 2D integer array
  !-------------------------------------------------------------------
  subroutine read_array_2d_int(filename, arr, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    integer, intent(out) :: arr(is:ie, js:je)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', &
         status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_2d_int

end program kernel_benchmark_gpu_roughitr

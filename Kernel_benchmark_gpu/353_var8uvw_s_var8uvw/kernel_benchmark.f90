!***********************************************************************
! GPU Kernel Benchmark: var8uvw (s_var8uvw)
!***********************************************************************
!
! Source: Src/var8uvw.f90
! Description: Average optional variable to u, v, and w staggered grid
!              points using simple 2-point averaging in each direction.
!
! GPU Difficulty: Easy - Simple 3D parallel loops with averaging
!
!***********************************************************************
program kernel_benchmark_var8uvw
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: wbc, ebc, exbopt

  ! Input array
  real, allocatable :: var(:,:,:)

  ! Output arrays
  real, allocatable :: var8u(:,:,:)  ! Variable at u points
  real, allocatable :: var8v(:,:,:)  ! Variable at v points
  real, allocatable :: var8w(:,:,:)  ! Variable at w points

  ! Reference outputs for validation
  real, allocatable :: var8u_ref(:,:,:)
  real, allocatable :: var8v_ref(:,:,:)
  real, allocatable :: var8w_ref(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', wbc, ebc, exbopt, ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: var8uvw'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6,A,I6)') ' wbc=', wbc, ', ebc=', ebc, ', exbopt=', exbopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(var(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8u_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8v_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8w_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/var.bin', var, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/var8u_ref.bin', var8u_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/var8v_ref.bin', var8v_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/var8w_ref.bin', var8w_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    var8u = 0.0
    var8v = 0.0
    var8w = 0.0
    call kernel_var8uvw(ni, nj, nk, var, var8u, var8v, var8w)
  end do
  !$acc wait

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    var8u = 0.0
    var8v = 0.0
    var8w = 0.0

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_var8uvw(ni, nj, nk, var, var8u, var8v, var8w)
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

  ! Validate var8u (computed for i=1:ni, j=0:nj, k=1:nk-1)
  do k = 1, nk-1
    do j = 0, nj
      do i = 1, ni
        rel_error = abs(var8u(i,j,k) - var8u_ref(i,j,k))
        if (abs(var8u_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(var8u_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate var8v (computed for i=0:ni, j=1:nj, k=1:nk-1)
  do k = 1, nk-1
    do j = 1, nj
      do i = 0, ni
        rel_error = abs(var8v(i,j,k) - var8v_ref(i,j,k))
        if (abs(var8v_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(var8v_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate var8w (computed for i=0:ni, j=0:nj, k=2:nk-1)
  do k = 2, nk-1
    do j = 0, nj
      do i = 0, ni
        rel_error = abs(var8w(i,j,k) - var8w_ref(i,j,k))
        if (abs(var8w_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(var8w_ref(i,j,k))
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
  deallocate(var, var8u, var8v, var8w, var8u_ref, var8v_ref, var8w_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: var8uvw
  ! Average variable to u, v, and w points
  !=====================================================================
  subroutine kernel_var8uvw(ni, nj, nk, var, var8u, var8v, var8w)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: var(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: var8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: var8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: var8w(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    ! Compute var8u
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = 0, nj
        !$acc loop independent
        do i = 1, ni
          var8u(i,j,k) = 0.5e0 * (var(i-1,j,k) + var(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    ! Compute var8v
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = 1, nj
        !$acc loop independent
        do i = 0, ni
          var8v(i,j,k) = 0.5e0 * (var(i,j-1,k) + var(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    ! Compute var8w
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      !$acc loop independent
      do j = 0, nj
        !$acc loop independent
        do i = 0, ni
          var8w(i,j,k) = 0.5e0 * (var(i,j,k-1) + var(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_var8uvw

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

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
  subroutine read_parameters(filename, wbc, ebc, exbopt, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: wbc, ebc, exbopt, ni, nj, nk

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    wbc = 7
    ebc = 7
    exbopt = 0
    ni = 1
    nj = 1
    nk = 1

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
          case ('wbc')
            read(val, *) wbc
          case ('ebc')
            read(val, *) ebc
          case ('exbopt')
            read(val, *) exbopt
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

end program kernel_benchmark_var8uvw

!***********************************************************************
! Kernel Benchmark: bc8u (s_bc8u)
!***********************************************************************
!
! Source: Src/bc8u.f90
! Description: Set west and east boundary conditions for optional
!              variable at u points by copying from adjacent interior
!              points based on BC type.
!
!***********************************************************************
program kernel_benchmark_bc8u
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, kmax

  ! Parameters
  integer :: wbc, ebc

  ! Input/output arrays
  real, allocatable :: var8u(:,:,:)

  ! Reference output for validation
  real, allocatable :: var8u_ref(:,:,:)

  ! Backup arrays for iteration
  real, allocatable :: var8u_in(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', wbc, ebc, ni, nj, kmax)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: bc8u'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', kmax=', kmax
  write(*,'(A,I6,A,I6)') ' Boundary conditions: wbc=', wbc, ', ebc=', ebc
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(var8u(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var8u_ref(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var8u_in(0:ni+1, 0:nj+1, 1:kmax))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/var8u_in.bin', var8u_in, 0, ni+1, 0, nj+1, 1, kmax)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/var8u_ref.bin', var8u_ref, 0, ni+1, 0, nj+1, 1, kmax)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    var8u = var8u_in
    call kernel_bc8u(wbc, ebc, ni, nj, kmax, var8u)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    var8u = var8u_in

    t_start = omp_get_wtime()
    call kernel_bc8u(wbc, ebc, ni, nj, kmax, var8u)
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

  do k = 1, kmax
    do j = 0, nj+1
      do i = 0, ni+1
        rel_error = abs(var8u(i,j,k) - var8u_ref(i,j,k))
        if (abs(var8u_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(var8u_ref(i,j,k))
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
  deallocate(var8u, var8u_ref, var8u_in, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: bc8u
  ! Set west and east boundary conditions for variable at u points
  !=====================================================================
  subroutine kernel_bc8u(wbc, ebc, ni, nj, kmax, var8u)
    implicit none

    integer, intent(in) :: wbc, ebc
    integer, intent(in) :: ni, nj, kmax
    real, intent(inout) :: var8u(0:ni+1, 0:nj+1, 1:kmax)

    integer :: j, k
    integer :: nim1, nim2

    ! MPI variables for single-process case (test_real)
    integer, parameter :: ebw = 1    ! West boundary exists
    integer, parameter :: ebe = 1    ! East boundary exists
    integer, parameter :: isub = 0   ! Subdomain index
    integer, parameter :: nisub = 1  ! Total subdomains in x

    nim1 = ni - 1
    nim2 = ni - 2

    !$omp parallel default(shared) private(k)

    ! Set the west boundary conditions
    if (ebw == 1 .and. isub == 0) then

      if (wbc == 2) then

        do k = 1, kmax
          !$omp do schedule(runtime) private(j)
          do j = 0, nj+1
            var8u(1,j,k) = var8u(3,j,k)
          end do
          !$omp end do
        end do

      else if (wbc >= 3) then

        do k = 1, kmax
          !$omp do schedule(runtime) private(j)
          do j = 0, nj+1
            var8u(1,j,k) = var8u(2,j,k)
          end do
          !$omp end do
        end do

      end if

    end if

    ! Set the east boundary conditions
    if (ebe == 1 .and. isub == nisub-1) then

      if (ebc == 2) then

        do k = 1, kmax
          !$omp do schedule(runtime) private(j)
          do j = 0, nj+1
            var8u(ni,j,k) = var8u(nim2,j,k)
          end do
          !$omp end do
        end do

      else if (ebc >= 3) then

        do k = 1, kmax
          !$omp do schedule(runtime) private(j)
          do j = 0, nj+1
            var8u(ni,j,k) = var8u(nim1,j,k)
          end do
          !$omp end do
        end do

      end if

    end if

    !$omp end parallel

  end subroutine kernel_bc8u

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
  subroutine read_parameters(filename, wbc, ebc, ni, nj, kmax)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: wbc, ebc, ni, nj, kmax

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    wbc = 3
    ebc = 3
    ni = 1
    nj = 1
    kmax = 1

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
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('kmax')
            read(val, *) kmax
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

end program kernel_benchmark_bc8u

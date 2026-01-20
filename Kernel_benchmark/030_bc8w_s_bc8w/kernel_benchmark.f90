!***********************************************************************
! Kernel Benchmark: bc8w (s_bc8w)
!***********************************************************************
!
! Source: Src/bc8w.f90
! Description: Set bottom and top boundary conditions for optional
!              variable at w points by copying from adjacent vertical
!              levels based on BC type.
!
!***********************************************************************
program kernel_benchmark_bc8w
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: bbc, tbc

  ! Input/output arrays
  real, allocatable :: var8w(:,:,:)

  ! Reference output for validation
  real, allocatable :: var8w_ref(:,:,:)

  ! Backup arrays for iteration
  real, allocatable :: var8w_in(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', bbc, tbc, ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: bc8w'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' Boundary conditions: bbc=', bbc, ', tbc=', tbc
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(var8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8w_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8w_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/var8w_in.bin', var8w_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/var8w_ref.bin', var8w_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    var8w = var8w_in
    call kernel_bc8w(bbc, tbc, ni, nj, nk, var8w)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    var8w = var8w_in

    t_start = omp_get_wtime()
    call kernel_bc8w(bbc, tbc, ni, nj, nk, var8w)
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

  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
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
  deallocate(var8w, var8w_ref, var8w_in, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: bc8w
  ! Set bottom and top boundary conditions for variable at w points
  !=====================================================================
  subroutine kernel_bc8w(bbc, tbc, ni, nj, nk, var8w)
    implicit none

    integer, intent(in) :: bbc, tbc
    integer, intent(in) :: ni, nj, nk
    real, intent(inout) :: var8w(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j
    integer :: nkm1, nkm2

    nkm1 = nk - 1
    nkm2 = nk - 2

    !$omp parallel default(shared)

    ! Set the bottom boundary conditions
    if (bbc == 2) then

      !$omp do schedule(runtime) private(i,j)
      do j = 0, nj+1
        do i = 0, ni+1
          var8w(i,j,1) = var8w(i,j,3)
        end do
      end do
      !$omp end do

    else if (bbc >= 3) then

      !$omp do schedule(runtime) private(i,j)
      do j = 0, nj+1
        do i = 0, ni+1
          var8w(i,j,1) = var8w(i,j,2)
        end do
      end do
      !$omp end do

    end if

    ! Set the top boundary conditions
    if (tbc == 2) then

      !$omp do schedule(runtime) private(i,j)
      do j = 0, nj+1
        do i = 0, ni+1
          var8w(i,j,nk) = var8w(i,j,nkm2)
        end do
      end do
      !$omp end do

    else if (tbc >= 3) then

      !$omp do schedule(runtime) private(i,j)
      do j = 0, nj+1
        do i = 0, ni+1
          var8w(i,j,nk) = var8w(i,j,nkm1)
        end do
      end do
      !$omp end do

    end if

    !$omp end parallel

  end subroutine kernel_bc8w

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
  subroutine read_parameters(filename, bbc, tbc, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: bbc, tbc, ni, nj, nk

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    bbc = 2
    tbc = 2
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
          case ('bbc')
            read(val, *) bbc
          case ('tbc')
            read(val, *) tbc
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

end program kernel_benchmark_bc8w

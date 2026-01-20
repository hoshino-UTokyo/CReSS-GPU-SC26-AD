!***********************************************************************
! Kernel Benchmark: copy3d (s_copy3d)
!***********************************************************************
!
! Source: Src/copy3d.f90
! Description: Simple 3D array copy from invar to outvar.
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_copy3d
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk
  integer :: imin, imax, jmin, jmax, kmin, kmax

  ! Input arrays
  real, allocatable :: invar(:,:,:)

  ! Output array
  real, allocatable :: outvar(:,:,:)

  ! Reference output for validation
  real, allocatable :: outvar_ref(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', &
       ni, nj, nk, imin, imax, jmin, jmax, kmin, kmax)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: copy3d'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I4,A,I4,A,I4,A,I4,A,I4,A,I4)') ' Bounds: i=[', imin, ',', imax, &
       '], j=[', jmin, ',', jmax, '], k=[', kmin, ',', kmax, ']'
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(invar(imin:imax, jmin:jmax, kmin:kmax))
  allocate(outvar(imin:imax, jmin:jmax, kmin:kmax))
  allocate(outvar_ref(imin:imax, jmin:jmax, kmin:kmax))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/invar.bin', invar, imin, imax, jmin, jmax, kmin, kmax)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/outvar_ref.bin', outvar_ref, imin, imax, jmin, jmax, kmin, kmax)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_copy3d(imin, imax, jmin, jmax, kmin, kmax, invar, outvar)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    t_start = omp_get_wtime()

    call kernel_copy3d(imin, imax, jmin, jmax, kmin, kmax, invar, outvar)

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

  do k = kmin, kmax
    do j = jmin, jmax
      do i = imin, imax
        rel_error = abs(outvar(i,j,k) - outvar_ref(i,j,k))
        if (abs(outvar_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(outvar_ref(i,j,k))
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
  deallocate(invar, outvar, outvar_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: copy3d
  ! Simple 3D array copy
  !=====================================================================
  subroutine kernel_copy3d(imin, imax, jmin, jmax, kmin, kmax, invar, outvar)
    implicit none

    integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax
    real, intent(in) :: invar(imin:imax, jmin:jmax, kmin:kmax)
    real, intent(out) :: outvar(imin:imax, jmin:jmax, kmin:kmax)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    do k = kmin, kmax
      !$omp do schedule(runtime) private(i,j)
      do j = jmin, jmax
        do i = imin, imax
          outvar(i,j,k) = invar(i,j,k)
        end do
      end do
      !$omp end do
    end do

    !$omp end parallel

  end subroutine kernel_copy3d

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    character(len=256) :: config_file
    integer :: ios
    logical :: exists

    ! Default values
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

  !=====================================================================
  ! Parameter reader
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk, imin, imax, jmin, jmax, kmin, kmax)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, imin, imax, jmin, jmax, kmin, kmax

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
          case ('imin')
            read(val, *) imin
          case ('imax')
            read(val, *) imax
          case ('jmin')
            read(val, *) jmin
          case ('jmax')
            read(val, *) jmax
          case ('kmin')
            read(val, *) kmin
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

end program kernel_benchmark_copy3d

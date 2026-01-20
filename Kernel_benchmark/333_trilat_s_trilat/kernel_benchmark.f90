!***********************************************************************
! Kernel Benchmark: trilat (s_trilat)
!***********************************************************************
!
! Source: Src/trilat.f90
! Description: Calculate Coriolis parameters (x 0.25) from latitude.
!              fc(i,j,1) is vertical component, fc(i,j,2) is horizontal.
!
!***********************************************************************
program kernel_benchmark_trilat
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj

  ! Parameters
  integer :: coropt

  ! Physical constants from m_comphy
  real, parameter :: omega = 7.292e-5  ! Angular velocity of Earth

  ! Mathematical constants from m_commath
  real, parameter :: d2r = 3.141592e0 / 180.e0  ! Degrees to radians

  ! Input array
  real, allocatable :: lat(:,:)       ! Latitude

  ! Output array
  real, allocatable :: fc(:,:,:)      ! Coriolis parameters (2 components)

  ! Reference output for validation
  real, allocatable :: fc_ref(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', coropt, ni, nj)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: trilat'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj
  write(*,'(A,I6)') ' coropt=', coropt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(lat(0:ni+1, 0:nj+1))
  allocate(fc(0:ni+1, 0:nj+1, 1:2))
  allocate(fc_ref(0:ni+1, 0:nj+1, 1:2))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_2d(trim(data_dir)//'/lat.bin', lat, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/fc_ref.bin', fc_ref, 0, ni+1, 0, nj+1, 1, 2)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    fc = 0.0
    call kernel_trilat(coropt, ni, nj, lat, fc)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    fc = 0.0

    t_start = omp_get_wtime()
    call kernel_trilat(coropt, ni, nj, lat, fc)
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output (interior points)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0

  do k = 1, 2
    do j = 1, nj-1
      do i = 1, ni-1
        rel_error = abs(fc(i,j,k) - fc_ref(i,j,k))
        if (abs(fc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(fc_ref(i,j,k))
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
  deallocate(lat, fc, fc_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: trilat
  ! Calculate Coriolis parameters from latitude
  !=====================================================================
  subroutine kernel_trilat(coropt, ni, nj, lat, fc)
    implicit none

    integer, intent(in) :: coropt
    integer, intent(in) :: ni, nj
    real, intent(in) :: lat(0:ni+1, 0:nj+1)
    real, intent(out) :: fc(0:ni+1, 0:nj+1, 1:2)

    integer :: i, j
    real :: omega5, sinlat

    omega5 = 0.5e0 * omega

    !$omp parallel default(shared)

    if (coropt == 1) then
      !$omp do schedule(runtime) private(i,j,sinlat)
      do j = 1, nj-1
        do i = 1, ni-1
          sinlat = sin(lat(i,j) * d2r)
          fc(i,j,1) = omega5 * sinlat
          fc(i,j,2) = 0.e0
        end do
      end do
      !$omp end do
    else if (coropt == 2) then
      !$omp do schedule(runtime) private(i,j,sinlat)
      do j = 1, nj-1
        do i = 1, ni-1
          sinlat = sin(lat(i,j) * d2r)
          fc(i,j,1) = omega5 * sinlat
          fc(i,j,2) = omega5 * sqrt(1.e0 - sinlat * sinlat)
        end do
      end do
      !$omp end do
    end if

    !$omp end parallel

  end subroutine kernel_trilat

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
  subroutine read_parameters(filename, coropt, ni, nj)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: coropt, ni, nj

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    coropt = 1
    ni = 1
    nj = 1

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
          case ('coropt')
            read(val, *) coropt
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
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

end program kernel_benchmark_trilat

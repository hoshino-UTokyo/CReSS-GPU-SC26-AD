!***********************************************************************
! Kernel Benchmark: coriuv (s_coriuv)
!***********************************************************************
!
! Source: Src/coriuv.f90
! Description: Calculate Coriolis force for u and v velocity equations
!
! This benchmark:
!   1. Reads parameters and input arrays from dump data
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_coriuv
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Input arrays
  real, allocatable :: fc(:,:,:)
  real, allocatable :: rst(:,:,:)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)

  ! Input/Output arrays
  real, allocatable :: ufrc(:,:,:)
  real, allocatable :: vfrc(:,:,:)
  real, allocatable :: tmp1(:,:,:)

  ! Reference output for validation
  real, allocatable :: ufrc_ref(:,:,:)
  real, allocatable :: vfrc_ref(:,:,:)
  real, allocatable :: tmp1_ref(:,:,:)

  ! Work arrays for re-initialization
  real, allocatable :: ufrc_init(:,:,:)
  real, allocatable :: vfrc_init(:,:,:)
  real, allocatable :: tmp1_init(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: coriuv'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(fc(0:ni+1, 0:nj+1, 1:2))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/fc.bin', fc, 0, ni+1, 0, nj+1, 1, 2)
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ufrc_in.bin', ufrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_in.bin', vfrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp1_in.bin', tmp1_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ufrc_ref.bin', ufrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_ref.bin', vfrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp1_ref.bin', tmp1_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ufrc = ufrc_init
    vfrc = vfrc_init
    tmp1 = tmp1_init
    call kernel_coriuv(ni, nj, nk, fc, rst, u, v, ufrc, vfrc, tmp1)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ufrc = ufrc_init
    vfrc = vfrc_init
    tmp1 = tmp1_init

    t_start = omp_get_wtime()
    call kernel_coriuv(ni, nj, nk, fc, rst, u, v, ufrc, vfrc, tmp1)
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

  ! Validate ufrc
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-1
        rel_error = abs(ufrc(i,j,k) - ufrc_ref(i,j,k))
        if (abs(ufrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(ufrc_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate vfrc
  do k = 2, nk-2
    do j = 2, nj-1
      do i = 2, ni-2
        rel_error = abs(vfrc(i,j,k) - vfrc_ref(i,j,k))
        if (abs(vfrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(vfrc_ref(i,j,k))
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
  deallocate(fc, rst, u, v)
  deallocate(ufrc, vfrc, tmp1)
  deallocate(ufrc_ref, vfrc_ref, tmp1_ref)
  deallocate(ufrc_init, vfrc_init, tmp1_init)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: coriuv
  ! Calculate Coriolis force for u and v velocity equations
  !=====================================================================
  subroutine kernel_coriuv(ni, nj, nk, fc, rst, u, v, ufrc, vfrc, tmp1)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: fc(0:ni+1, 0:nj+1, 1:2)
    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: ufrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! Calculate the Coriolis force in the x components of velocity equation.
    do k = 2, nk-2
      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-2
        do i = 1, ni-1
          tmp1(i,j,k) = fc(i,j,1) * rst(i,j,k) * (v(i,j,k) + v(i,j+1,k))
        end do
      end do
      !$omp end do
    end do

    do k = 2, nk-2
      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-2
        do i = 2, ni-1
          ufrc(i,j,k) = ufrc(i,j,k) + (tmp1(i-1,j,k) + tmp1(i,j,k))
        end do
      end do
      !$omp end do
    end do

    ! Calculate the Coriolis force in the y components of velocity equation.
    do k = 2, nk-2
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 2, ni-2
          tmp1(i,j,k) = -fc(i,j,1) * rst(i,j,k) * (u(i,j,k) + u(i+1,j,k))
        end do
      end do
      !$omp end do
    end do

    do k = 2, nk-2
      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-1
        do i = 2, ni-2
          vfrc(i,j,k) = vfrc(i,j,k) + (tmp1(i,j-1,k) + tmp1(i,j,k))
        end do
      end do
      !$omp end do
    end do

    !$omp end parallel

  end subroutine kernel_coriuv

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
  ! Parameter reader (key=value format)
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk

    character(len=256) :: line, name, val
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
        name = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))
        select case (trim(name))
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

end program kernel_benchmark_coriuv

!***********************************************************************
! Kernel Benchmark: strsten (s_strsten)
!***********************************************************************
!
! Source: Src/strsten.f90
! Description: Calculate stress tensor components
!
! NOTE: sfcopt is hardcoded based on test_real configuration (sfcopt=1)
!
!***********************************************************************
program kernel_benchmark_strsten
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Option parameter (hardcoded from test_real)
  integer, parameter :: sfcopt = 1

  ! Input arrays
  real, allocatable :: ufrc(:,:,:)
  real, allocatable :: vfrc(:,:,:)
  real, allocatable :: rkh(:,:,:)
  real, allocatable :: rkv(:,:,:)

  ! Input/Output arrays
  real, allocatable :: t11(:,:,:), t22(:,:,:), t33(:,:,:)
  real, allocatable :: t12(:,:,:), t13(:,:,:), t23(:,:,:)
  real, allocatable :: t31(:,:,:), t32(:,:,:)

  ! Initial values for re-initialization
  real, allocatable :: t11_init(:,:,:), t22_init(:,:,:), t33_init(:,:,:)
  real, allocatable :: t12_init(:,:,:), t13_init(:,:,:), t23_init(:,:,:)
  real, allocatable :: t31_init(:,:,:), t32_init(:,:,:)

  ! Reference output
  real, allocatable :: t11_ref(:,:,:), t22_ref(:,:,:), t33_ref(:,:,:)
  real, allocatable :: t12_ref(:,:,:), t13_ref(:,:,:), t23_ref(:,:,:)
  real, allocatable :: t31_ref(:,:,:), t32_ref(:,:,:)

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
  write(*,'(A)') ' Kernel Benchmark: strsten'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I3)') ' sfcopt=', sfcopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkh(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkv(0:ni+1, 0:nj+1, 1:nk))

  allocate(t11(0:ni+1, 0:nj+1, 1:nk))
  allocate(t22(0:ni+1, 0:nj+1, 1:nk))
  allocate(t33(0:ni+1, 0:nj+1, 1:nk))
  allocate(t12(0:ni+1, 0:nj+1, 1:nk))
  allocate(t13(0:ni+1, 0:nj+1, 1:nk))
  allocate(t23(0:ni+1, 0:nj+1, 1:nk))
  allocate(t31(0:ni+1, 0:nj+1, 1:nk))
  allocate(t32(0:ni+1, 0:nj+1, 1:nk))

  allocate(t11_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t22_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t33_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t12_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t13_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t23_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t31_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t32_init(0:ni+1, 0:nj+1, 1:nk))

  allocate(t11_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(t22_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(t33_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(t12_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(t13_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(t23_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(t31_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(t32_ref(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/ufrc.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rkh.bin', rkh, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rkv.bin', rkv, 0, ni+1, 0, nj+1, 1, nk)

  call read_array_3d(trim(data_dir)//'/t11_in.bin', t11_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t22_in.bin', t22_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t33_in.bin', t33_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t12_in.bin', t12_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t13_in.bin', t13_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t23_in.bin', t23_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t31_in.bin', t31_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t32_in.bin', t32_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/t11_ref.bin', t11_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t22_ref.bin', t22_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t33_ref.bin', t33_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t12_ref.bin', t12_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t13_ref.bin', t13_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t23_ref.bin', t23_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t31_ref.bin', t31_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t32_ref.bin', t32_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    t11 = t11_init; t22 = t22_init; t33 = t33_init
    t12 = t12_init; t13 = t13_init; t23 = t23_init
    t31 = t31_init; t32 = t32_init
    call kernel_strsten(ni, nj, nk, ufrc, vfrc, rkh, rkv, &
         t11, t22, t33, t12, t13, t23, t31, t32)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    t11 = t11_init; t22 = t22_init; t33 = t33_init
    t12 = t12_init; t13 = t13_init; t23 = t23_init
    t31 = t31_init; t32 = t32_init

    t_start = omp_get_wtime()
    call kernel_strsten(ni, nj, nk, ufrc, vfrc, rkh, rkv, &
         t11, t22, t33, t12, t13, t23, t31, t32)
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

  call validate_array(t11, t11_ref, ni, nj, nk, 1, ni-1, 1, nj-1, 1, nk-1, &
       tolerance, max_error, error_count)
  call validate_array(t22, t22_ref, ni, nj, nk, 1, ni-1, 1, nj-1, 1, nk-1, &
       tolerance, max_error, error_count)
  call validate_array(t33, t33_ref, ni, nj, nk, 1, ni-1, 1, nj-1, 1, nk-1, &
       tolerance, max_error, error_count)
  call validate_array(t12, t12_ref, ni, nj, nk, 2, ni-1, 2, nj-1, 1, nk-1, &
       tolerance, max_error, error_count)
  call validate_array(t13, t13_ref, ni, nj, nk, 2, ni-1, 2, nj-2, 2, nk-1, &
       tolerance, max_error, error_count)
  call validate_array(t23, t23_ref, ni, nj, nk, 2, ni-2, 2, nj-1, 2, nk-1, &
       tolerance, max_error, error_count)
  call validate_array(t31, t31_ref, ni, nj, nk, 2, ni-1, 2, nj-2, 2, nk-1, &
       tolerance, max_error, error_count)
  call validate_array(t32, t32_ref, ni, nj, nk, 2, ni-2, 2, nj-1, 2, nk-1, &
       tolerance, max_error, error_count)

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

  if (.not. validation_passed) stop 1

contains

  subroutine kernel_strsten(ni, nj, nk, ufrc, vfrc, rkh, rkv, &
       t11, t22, t33, t12, t13, t23, t31, t32)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: ufrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rkh(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rkv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t11(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t22(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t33(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t12(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t13(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t23(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t32(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! Calculate diagonal and x-y components
    do k = 1, nk-1
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          t11(i,j,k) = rkh(i,j,k) * t11(i,j,k)
          t22(i,j,k) = rkh(i,j,k) * t22(i,j,k)
          t33(i,j,k) = rkv(i,j,k) * t33(i,j,k)
        end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-1
        do i = 2, ni-1
          t12(i,j,k) = 0.25e0 * t12(i,j,k) &
               * ((rkh(i-1,j-1,k) + rkh(i,j,k)) + (rkh(i-1,j,k) + rkh(i,j-1,k)))
        end do
      end do
      !$omp end do
    end do

    ! Calculate x-z and y-z components
    do k = 2, nk-1
      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-2
        do i = 2, ni-1
          t13(i,j,k) = 0.25e0 * t13(i,j,k) &
               * ((rkv(i-1,j,k-1) + rkv(i,j,k)) + (rkv(i-1,j,k) + rkv(i,j,k-1)))
        end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-1
        do i = 2, ni-2
          t23(i,j,k) = 0.25e0 * t23(i,j,k) &
               * ((rkv(i,j-1,k-1) + rkv(i,j,k)) + (rkv(i,j-1,k) + rkv(i,j,k-1)))
        end do
      end do
      !$omp end do
    end do

    ! Surface stress (sfcopt >= 1)
    !$omp do schedule(runtime) private(i,j)
    do j = 2, nj-2
      do i = 2, ni-1
        t13(i,j,2) = ufrc(i,j,1)
      end do
    end do
    !$omp end do

    !$omp do schedule(runtime) private(i,j)
    do j = 2, nj-1
      do i = 2, ni-2
        t23(i,j,2) = vfrc(i,j,1)
      end do
    end do
    !$omp end do

    ! Calculate z-x and z-y components
    do k = 2, nk-1
      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-2
        do i = 2, ni-1
          t31(i,j,k) = 0.25e0 * t31(i,j,k) &
               * ((rkh(i-1,j,k-1) + rkh(i,j,k)) + (rkh(i-1,j,k) + rkh(i,j,k-1)))
        end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-1
        do i = 2, ni-2
          t32(i,j,k) = 0.25e0 * t32(i,j,k) &
               * ((rkh(i,j-1,k-1) + rkh(i,j,k)) + (rkh(i,j-1,k) + rkh(i,j,k-1)))
        end do
      end do
      !$omp end do
    end do

    !$omp end parallel

  end subroutine kernel_strsten

  subroutine validate_array(arr, ref, ni, nj, nk, i1, i2, j1, j2, k1, k2, &
       tol, max_err, err_cnt)
    integer, intent(in) :: ni, nj, nk, i1, i2, j1, j2, k1, k2
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(inout) :: max_err
    integer, intent(inout) :: err_cnt

    integer :: i, j, k
    real :: rel_err

    do k = k1, k2
      do j = j1, j2
        do i = i1, i2
          rel_err = abs(arr(i,j,k) - ref(i,j,k))
          if (abs(ref(i,j,k)) > 1.0e-10) then
            rel_err = rel_err / abs(ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_cnt = err_cnt + 1
        end do
      end do
    end do
  end subroutine validate_array

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

end program kernel_benchmark_strsten

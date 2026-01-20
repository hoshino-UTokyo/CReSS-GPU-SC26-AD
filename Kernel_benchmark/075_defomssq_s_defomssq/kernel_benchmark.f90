!***********************************************************************
! Kernel Benchmark: defomssq (s_defomssq)
!***********************************************************************
!
! Source: Src/defomssq.f90
! Description: Calculate magnitude of deformation tensor squared from
!              diagonal (s11, s22, s33) and off-diagonal (s12, s31, s32)
!              strain rate components using stencil averaging.
!
!***********************************************************************
program kernel_benchmark_defomssq
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Input arrays (strain rate tensor components)
  real, allocatable :: s11(:,:,:)  ! x-x component
  real, allocatable :: s22(:,:,:)  ! y-y component
  real, allocatable :: s33(:,:,:)  ! z-z component
  real, allocatable :: s12(:,:,:)  ! x-y component
  real, allocatable :: s31(:,:,:)  ! z-x component
  real, allocatable :: s32(:,:,:)  ! z-y component

  ! Output array
  real, allocatable :: ssq(:,:,:)  ! Deformation magnitude squared

  ! Reference output for validation
  real, allocatable :: ssq_ref(:,:,:)

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
  write(*,'(A)') ' Kernel Benchmark: defomssq'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(s11(0:ni+1, 0:nj+1, 1:nk))
  allocate(s22(0:ni+1, 0:nj+1, 1:nk))
  allocate(s33(0:ni+1, 0:nj+1, 1:nk))
  allocate(s12(0:ni+1, 0:nj+1, 1:nk))
  allocate(s31(0:ni+1, 0:nj+1, 1:nk))
  allocate(s32(0:ni+1, 0:nj+1, 1:nk))
  allocate(ssq(0:ni+1, 0:nj+1, 1:nk))
  allocate(ssq_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/s11.bin', s11, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s22.bin', s22, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s33.bin', s33, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s12.bin', s12, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s31.bin', s31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s32.bin', s32, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ssq_ref.bin', ssq_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ssq = 0.0
    call kernel_defomssq(ni, nj, nk, s11, s22, s33, s12, s31, s32, ssq)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ssq = 0.0

    t_start = omp_get_wtime()
    call kernel_defomssq(ni, nj, nk, s11, s22, s33, s12, s31, s32, ssq)
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

  do k = 1, nk-1
    do j = 1, nj-1
      do i = 1, ni-1
        rel_error = abs(ssq(i,j,k) - ssq_ref(i,j,k))
        if (abs(ssq_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(ssq_ref(i,j,k))
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
  deallocate(s11, s22, s33, s12, s31, s32, ssq, ssq_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: defomssq
  ! Calculate deformation magnitude squared
  !=====================================================================
  subroutine kernel_defomssq(ni, nj, nk, s11, s22, s33, s12, s31, s32, ssq)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: s11(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: s22(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: s33(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: s12(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: s31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: s32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: ssq(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    real :: s128s, s318s, s328s

    !$omp parallel default(shared) private(k)

    do k = 1, nk-1
      !$omp do schedule(runtime) private(i,j,s128s,s318s,s328s)
      do j = 1, nj-1
        do i = 1, ni-1
          s128s = (s12(i,j,k) + s12(i+1,j+1,k)) + (s12(i+1,j,k) + s12(i,j+1,k))
          s318s = (s31(i,j,k) + s31(i+1,j,k+1)) + (s31(i,j,k+1) + s31(i+1,j,k))
          s328s = (s32(i,j,k) + s32(i,j+1,k+1)) + (s32(i,j+1,k) + s32(i,j,k+1))

          ssq(i,j,k) = 0.5e0 * (s33(i,j,k)*s33(i,j,k) &
                       + (s11(i,j,k)*s11(i,j,k) + s22(i,j,k)*s22(i,j,k))) &
                       + 0.0625e0 * (s128s*s128s + (s318s*s318s + s328s*s328s))
        end do
      end do
      !$omp end do
    end do

    !$omp end parallel

  end subroutine kernel_defomssq

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
  subroutine read_parameters(filename, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
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

end program kernel_benchmark_defomssq

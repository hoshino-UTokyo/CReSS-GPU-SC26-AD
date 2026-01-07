!***********************************************************************
! Kernel Benchmark: advbspi (s_advbspi)
!***********************************************************************
!
! Source: Src/advbspi.f90
! Description: Calculate base state pressure advection for horizontally
!              explicit and vertically implicit method
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_advbspi
  use omp_lib
  implicit none

  ! Physical constants (read from params.txt)
  real :: g

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  real :: weicoe
  character(len=4) :: fproc

  ! Input arrays
  real, allocatable :: rst(:,:,:)
  real, allocatable :: w(:,:,:)
  real, allocatable :: pfrc(:,:,:)
  real, allocatable :: phdiv(:,:,:)
  real, allocatable :: pvdiv(:,:,:)

  ! Output array
  real, allocatable :: fp(:,:,:)

  ! Reference output for validation
  real, allocatable :: fp_ref(:,:,:)

  ! Work array for re-initialization
  real, allocatable :: fp_init(:,:,:)

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
       ni, nj, nk, weicoe, g, fproc)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: advbspi'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,F10.6)') ' weicoe=', weicoe
  write(*,'(A,F10.6)') ' g=', g
  write(*,'(A,A)') ' fproc=', fproc
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(w(0:ni+1, 0:nj+1, 1:nk))
  allocate(pfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(phdiv(0:ni+1, 0:nj+1, 1:nk))
  allocate(pvdiv(0:ni+1, 0:nj+1, 1:nk))
  allocate(fp(0:ni+1, 0:nj+1, 1:nk))
  allocate(fp_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(fp_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pfrc.bin', pfrc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/phdiv.bin', phdiv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pvdiv.bin', pvdiv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/fp_in.bin', fp_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/fp_ref.bin', fp_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Debug output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Debug: checking input values...'
  write(*,'(A,A,A)') ' fproc = [', fproc, ']'
  write(*,'(A,3ES15.7)') ' fp_init(100,100,10-12) = ', fp_init(100,100,10), fp_init(100,100,11), fp_init(100,100,12)
  write(*,'(A,3ES15.7)') ' fp_ref(100,100,10-12)  = ', fp_ref(100,100,10), fp_ref(100,100,11), fp_ref(100,100,12)
  write(*,'(A,3ES15.7)') ' rst(100,100,10-12)     = ', rst(100,100,10), rst(100,100,11), rst(100,100,12)
  write(*,'(A,3ES15.7)') ' w(100,100,10-12)       = ', w(100,100,10), w(100,100,11), w(100,100,12)
  write(*,'(A,3ES15.7)') ' pvdiv(100,100,10-12)   = ', pvdiv(100,100,10), pvdiv(100,100,11), pvdiv(100,100,12)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    fp = fp_init
    call kernel_advbspi(fproc, weicoe, ni, nj, nk, &
         rst, w, pfrc, phdiv, pvdiv, fp)
  end do

  write(*,'(A,3ES15.7)') ' fp after warmup(100,100,10-12) = ', fp(100,100,10), fp(100,100,11), fp(100,100,12)

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    fp = fp_init

    t_start = omp_get_wtime()

    call kernel_advbspi(fproc, weicoe, ni, nj, nk, &
         rst, w, pfrc, phdiv, pvdiv, fp)

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

  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(fp(i,j,k) - fp_ref(i,j,k))
        if (abs(fp_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(fp_ref(i,j,k))
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
  deallocate(rst, w, pfrc, phdiv, pvdiv)
  deallocate(fp, fp_ref, fp_init)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: advbspi
  ! Calculate base state pressure advection
  !=====================================================================
  subroutine kernel_advbspi(fproc, weicoe, ni, nj, nk, &
       rst, w, pfrc, phdiv, pvdiv, fp)
    implicit none

    character(len=4), intent(in) :: fproc
    real, intent(in) :: weicoe
    integer, intent(in) :: ni, nj, nk

    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: pfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: phdiv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: pvdiv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: fp(0:ni+1, 0:nj+1, 1:nk)

    ! Local variables
    real :: weic1m, g05
    integer :: i, j, k

    weic1m = 1.0e0 - weicoe
    g05 = 0.5e0 * g

    !$omp parallel default(shared) private(k)

    if (fproc(1:4) == 'back') then

      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-2
            fp(i,j,k) = pfrc(i,j,k) + phdiv(i,j,k) &
                 + weic1m * (g05 * (w(i,j,k) + w(i,j,k+1)) * rst(i,j,k) + pvdiv(i,j,k))
          end do
        end do
        !$omp end do
      end do

    else if (fproc(1:4) == 'fore') then

      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-2
            fp(i,j,k) = fp(i,j,k) &
                 + weicoe * (g05 * (w(i,j,k) + w(i,j,k+1)) * rst(i,j,k) + pvdiv(i,j,k))
          end do
        end do
        !$omp end do
      end do

    end if

    !$omp end parallel

  end subroutine kernel_advbspi

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
  ! Parameter reader (key=value format)
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk, weicoe, g, fproc)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    real, intent(out) :: weicoe, g
    character(len=4), intent(out) :: fproc

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
          case ('weicoe')
            read(val, *) weicoe
          case ('g')
            read(val, *) g
          case ('fproc')
            fproc = val(1:4)
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

end program kernel_benchmark_advbspi

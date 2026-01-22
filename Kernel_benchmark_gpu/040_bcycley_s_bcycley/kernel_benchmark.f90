!***********************************************************************
! GPU Kernel Benchmark: bcycley (s_bcycley)
!***********************************************************************
!
! Source: Src/bcycley.f90
! Description: Set periodic boundary conditions in y direction by copying
!              values between south and north boundaries for cyclic domains.
!              Only active when sbc==-1 and nbc==-1 (cyclic BC).
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_bcycley
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, kmax

  ! Parameters
  integer :: sbc, nbc, jssnd, jsrcv, jnsnd, jnrcv

  ! Input/output array
  real, allocatable :: var(:,:,:)

  ! Reference output for validation
  real, allocatable :: var_ref(:,:,:)

  ! Backup array for iteration
  real, allocatable :: var_in(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', sbc, nbc, &
                       jssnd, jsrcv, jnsnd, jnrcv, ni, nj, kmax)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: bcycley'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', kmax=', kmax
  write(*,'(A,I6,A,I6)') ' Boundary conditions: sbc=', sbc, ', nbc=', nbc
  write(*,'(A,I6,A,I6,A,I6,A,I6)') ' Indices: jssnd=', jssnd, ', jsrcv=', jsrcv, &
                                   ', jnsnd=', jnsnd, ', jnrcv=', jnrcv
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(var(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var_ref(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var_in(0:ni+1, 0:nj+1, 1:kmax))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/var_in.bin', var_in, 0, ni+1, 0, nj+1, 1, kmax)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/var_ref.bin', var_ref, 0, ni+1, 0, nj+1, 1, kmax)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    var = var_in
    call kernel_bcycley(sbc, nbc, jssnd, jsrcv, jnsnd, jnrcv, ni, nj, kmax, var)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    var = var_in

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_bcycley(sbc, nbc, jssnd, jsrcv, jnsnd, jnrcv, ni, nj, kmax, var)
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

  do k = 1, kmax
    do j = 0, nj+1
      do i = 0, ni+1
        rel_error = abs(var(i,j,k) - var_ref(i,j,k))
        if (abs(var_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(var_ref(i,j,k))
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
  deallocate(var, var_ref, var_in, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: bcycley (OpenACC version)
  ! Set periodic boundary conditions in y direction
  !=====================================================================
  subroutine kernel_bcycley(sbc, nbc, jssnd, jsrcv, jnsnd, jnrcv, ni, nj, kmax, var)
    implicit none

    integer, intent(in) :: sbc, nbc
    integer, intent(in) :: jssnd, jsrcv, jnsnd, jnrcv
    integer, intent(in) :: ni, nj, kmax
    real, intent(inout) :: var(0:ni+1, 0:nj+1, 1:kmax)

    integer :: i, k

    ! MPI variable for single-process case
    integer, parameter :: njsub = 1  ! Number of subdomains in y

    if (njsub == 1) then
      if (sbc == -1 .and. nbc == -1) then
        ! South boundary receives from north interior
        !$acc kernels
        !$acc loop independent
        do k = 1, kmax
          !$acc loop independent
          do i = 0, ni+1
            var(i,jsrcv,k) = var(i,jnsnd,k)
          end do
        end do
        !$acc end kernels

        ! North boundary receives from south interior
        !$acc kernels
        !$acc loop independent
        do k = 1, kmax
          !$acc loop independent
          do i = 0, ni+1
            var(i,jnrcv,k) = var(i,jssnd,k)
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_bcycley

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
  subroutine read_parameters(filename, sbc, nbc, jssnd, jsrcv, jnsnd, jnrcv, ni, nj, kmax)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: sbc, nbc, jssnd, jsrcv, jnsnd, jnrcv, ni, nj, kmax

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    sbc = 7
    nbc = 7
    jssnd = 4
    jsrcv = 0
    jnsnd = 896
    jnrcv = 900
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
          case ('sbc')
            read(val, *) sbc
          case ('nbc')
            read(val, *) nbc
          case ('jssnd')
            read(val, *) jssnd
          case ('jsrcv')
            read(val, *) jsrcv
          case ('jnsnd')
            read(val, *) jnsnd
          case ('jnrcv')
            read(val, *) jnrcv
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

end program kernel_benchmark_gpu_bcycley

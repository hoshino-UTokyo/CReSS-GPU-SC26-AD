!***********************************************************************
! GPU Kernel Benchmark: bcycle (s_bcycle)
!***********************************************************************
!
! Source: Src/bcycle.f90
! Description: Sets periodic boundary conditions by copying values between
!              west/east and south/north boundaries for cyclic domains.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_bcycle
  use omp_lib
  implicit none

  integer :: ni, nj, kmax
  integer :: wbc, ebc, sbc, nbc
  integer :: nisub, njsub
  integer :: iwsnd, iwrcv, iesnd, iercv
  integer :: jssnd, jsrcv, jnsnd, jnrcv

  real, allocatable :: var(:,:,:)
  real, allocatable :: var_ref(:,:,:)
  real, allocatable :: var_input(:,:,:)

  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  real :: max_error, rel_error, tolerance
  integer :: error_count
  logical :: validation_passed
  integer :: iter, i, j, k

  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, kmax, &
       wbc, ebc, sbc, nbc, nisub, njsub, &
       iwsnd, iwrcv, iesnd, iercv, jssnd, jsrcv, jnsnd, jnrcv)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: bcycle'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', kmax=', kmax
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  allocate(var(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var_ref(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var_input(0:ni+1, 0:nj+1, 1:kmax))
  allocate(times(num_iterations))

  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/var_in.bin', var_input, 0, ni+1, 0, nj+1, 1, kmax)
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/var_ref.bin', var_ref, 0, ni+1, 0, nj+1, 1, kmax)

  ! Warmup
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    var = var_input
    call kernel_bcycle(wbc, ebc, sbc, nbc, nisub, njsub, &
         iwsnd, iwrcv, iesnd, iercv, jssnd, jsrcv, jnsnd, jnrcv, &
         ni, nj, kmax, var)
    !$acc wait
  end do

  ! Benchmark
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  do iter = 1, num_iterations
    var = var_input
    !$acc wait
    t_start = omp_get_wtime()
    call kernel_bcycle(wbc, ebc, sbc, nbc, nisub, njsub, &
         iwsnd, iwrcv, iesnd, iercv, jssnd, jsrcv, jnsnd, jnrcv, &
         ni, nj, kmax, var)
    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do
  t_avg = t_total / dble(num_iterations)

  ! Validate
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0
  do k = 1, kmax
    do j = 0, nj+1
      do i = 0, ni+1
        rel_error = abs(var(i,j,k) - var_ref(i,j,k))
        if (abs(var_ref(i,j,k)) > 1.0e-20) rel_error = rel_error / abs(var_ref(i,j,k))
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do
  validation_passed = (error_count == 0)

  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,ES12.4)') ' Max error:    ', max_error
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  deallocate(var, var_ref, var_input, times)
  if (.not. validation_passed) stop 1

contains

  subroutine kernel_bcycle(wbc, ebc, sbc, nbc, nisub, njsub, &
       iwsnd, iwrcv, iesnd, iercv, jssnd, jsrcv, jnsnd, jnrcv, &
       ni, nj, kmax, var)
    integer, intent(in) :: wbc, ebc, sbc, nbc, nisub, njsub
    integer, intent(in) :: iwsnd, iwrcv, iesnd, iercv
    integer, intent(in) :: jssnd, jsrcv, jnsnd, jnrcv
    integer, intent(in) :: ni, nj, kmax
    real, intent(inout) :: var(0:ni+1, 0:nj+1, 1:kmax)

    integer :: i, j, k

    if (nisub == 1) then
      if (wbc == -1 .and. ebc == -1) then
        !$acc kernels
        !$acc loop independent
        do k = 1, kmax
          !$acc loop independent
          do j = 0, nj+1
            var(iwrcv,j,k) = var(iesnd,j,k)
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k = 1, kmax
          !$acc loop independent
          do j = 0, nj+1
            var(iercv,j,k) = var(iwsnd,j,k)
          end do
        end do
        !$acc end kernels
      end if
    end if

    if (njsub == 1) then
      if (sbc == -1 .and. nbc == -1) then
        !$acc kernels
        !$acc loop independent
        do k = 1, kmax
          !$acc loop independent
          do i = 0, ni+1
            var(i,jsrcv,k) = var(i,jnsnd,k)
          end do
        end do
        !$acc end kernels

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

  end subroutine kernel_bcycle

  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    character(len=256) :: arg
    integer :: ios, nargs
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    nargs = command_argument_count()
    if (nargs >= 1) then
      call get_command_argument(1, arg)
      data_dir = trim(arg)
      return
    end if

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

  subroutine read_parameters(filename, ni, nj, kmax, &
       wbc, ebc, sbc, nbc, nisub, njsub, &
       iwsnd, iwrcv, iesnd, iercv, jssnd, jsrcv, jnsnd, jnrcv)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, kmax
    integer, intent(out) :: wbc, ebc, sbc, nbc, nisub, njsub
    integer, intent(out) :: iwsnd, iwrcv, iesnd, iercv
    integer, intent(out) :: jssnd, jsrcv, jnsnd, jnrcv
    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) stop 'Cannot open params file'
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
          case ('kmax')
            read(val, *) kmax
          case ('wbc')
            read(val, *) wbc
          case ('ebc')
            read(val, *) ebc
          case ('sbc')
            read(val, *) sbc
          case ('nbc')
            read(val, *) nbc
          case ('nisub')
            read(val, *) nisub
          case ('njsub')
            read(val, *) njsub
          case ('iwsnd')
            read(val, *) iwsnd
          case ('iwrcv')
            read(val, *) iwrcv
          case ('iesnd')
            read(val, *) iesnd
          case ('iercv')
            read(val, *) iercv
          case ('jssnd')
            read(val, *) jssnd
          case ('jsrcv')
            read(val, *) jsrcv
          case ('jnsnd')
            read(val, *) jnsnd
          case ('jnrcv')
            read(val, *) jnrcv
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

    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) stop 'Cannot open data file'
    read(10) arr
    close(10)
  end subroutine read_array_3d

end program kernel_benchmark_bcycle

!***********************************************************************
! Kernel Benchmark: advs (s_advs)
!***********************************************************************
!
! Source: Src/advs.f90
! Description: Calculate scalar advection using 2nd/4th order schemes
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_advs
  use omp_lib
  implicit none

  ! Math constants
  real, parameter :: oned24 = 1.0e0 / 24.0e0
  real, parameter :: fourd3 = 4.0e0 / 3.0e0

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: advopt
  integer :: iwest, ieast, jsouth, jnorth
  real :: dxiv, dyiv, dziv

  ! Input arrays
  real, allocatable :: rstxu(:,:,:)
  real, allocatable :: rstxv(:,:,:)
  real, allocatable :: rstxwc(:,:,:)
  real, allocatable :: s(:,:,:)

  ! Output array
  real, allocatable :: sfrc(:,:,:)

  ! Work arrays
  real, allocatable :: vadv(:,:,:)
  real, allocatable :: tmp1(:,:,:)
  real, allocatable :: tmp2(:,:,:)
  real, allocatable :: tmp3(:,:,:)

  ! Reference output for validation
  real, allocatable :: sfrc_ref(:,:,:)

  ! Initial value for re-initialization
  real, allocatable :: sfrc_init(:,:,:)

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
       ni, nj, nk, advopt, iwest, ieast, jsouth, jnorth, dxiv, dyiv, dziv)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: advs'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' advopt=', advopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rstxu(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxv(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxwc(0:ni+1, 0:nj+1, 1:nk))
  allocate(s(0:ni+1, 0:nj+1, 1:nk))
  allocate(sfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(sfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(sfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vadv(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rstxu.bin', rstxu, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rstxv.bin', rstxv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rstxwc.bin', rstxwc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s.bin', s, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/sfrc_in.bin', sfrc_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/sfrc_ref.bin', sfrc_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    sfrc = sfrc_init
    call kernel_advs(advopt, iwest, ieast, jsouth, jnorth, &
         dxiv, dyiv, dziv, ni, nj, nk, &
         rstxu, rstxv, rstxwc, s, sfrc, vadv, tmp1, tmp2, tmp3)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    sfrc = sfrc_init

    t_start = omp_get_wtime()

    call kernel_advs(advopt, iwest, ieast, jsouth, jnorth, &
         dxiv, dyiv, dziv, ni, nj, nk, &
         rstxu, rstxv, rstxwc, s, sfrc, vadv, tmp1, tmp2, tmp3)

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
        rel_error = abs(sfrc(i,j,k) - sfrc_ref(i,j,k))
        if (abs(sfrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(sfrc_ref(i,j,k))
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
  deallocate(rstxu, rstxv, rstxwc, s)
  deallocate(sfrc, sfrc_init, sfrc_ref)
  deallocate(vadv, tmp1, tmp2, tmp3)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: advs
  ! Calculate scalar advection
  !=====================================================================
  subroutine kernel_advs(advopt, iwest, ieast, jsouth, jnorth, &
       dxiv, dyiv, dziv, ni, nj, nk, &
       rstxu, rstxv, rstxwc, s, sfrc, vadv, tmp1, tmp2, tmp3)
    implicit none

    integer, intent(in) :: advopt
    integer, intent(in) :: iwest, ieast, jsouth, jnorth
    real, intent(in) :: dxiv, dyiv, dziv
    integer, intent(in) :: ni, nj, nk

    real, intent(in) :: rstxu(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rstxv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rstxwc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: s(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: sfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vadv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp3(0:ni+1, 0:nj+1, 1:nk)

    ! Local variables
    real :: dxv05n, dyv05n, dzv05n
    real :: dxv24, dyv24, dzv24
    integer :: i, j, k

    dxv05n = -0.5e0 * dxiv
    dyv05n = -0.5e0 * dyiv
    dzv05n = -0.5e0 * dziv

    dxv24 = oned24 * dxiv
    dyv24 = oned24 * dyiv
    dzv24 = oned24 * dziv

    !$omp parallel default(shared) private(k)

    if (advopt <= 3) then

      ! Calculate the 2nd order scalar advection
      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-1
            tmp1(i,j,k) = rstxu(i,j,k) * (s(i,j,k) - s(i-1,j,k)) * dxv05n
          end do
        end do
        !$omp end do

        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-1
          do i = 2, ni-2
            tmp2(i,j,k) = rstxv(i,j,k) * (s(i,j,k) - s(i,j-1,k)) * dyv05n
          end do
        end do
        !$omp end do
      end do

      do k = 2, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-2
            tmp3(i,j,k) = rstxwc(i,j,k) * (s(i,j,k) - s(i,j,k-1)) * dzv05n
          end do
        end do
        !$omp end do
      end do

      if (advopt == 1) then

        do k = 2, nk-2
          !$omp do schedule(runtime) private(i,j)
          do j = 2, nj-2
            do i = 2, ni-2
              sfrc(i,j,k) = (tmp3(i,j,k) + tmp3(i,j,k+1)) &
                   + ((tmp1(i,j,k) + tmp1(i+1,j,k)) &
                   + (tmp2(i,j,k) + tmp2(i,j+1,k)))
            end do
          end do
          !$omp end do
        end do

      else

        if (advopt == 2) then

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                sfrc(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) &
                     + (tmp2(i,j,k) + tmp2(i,j+1,k))
              end do
            end do
            !$omp end do
          end do

        else if (advopt == 3) then

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                sfrc(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) &
                     + (tmp2(i,j,k) + tmp2(i,j+1,k))
                vadv(i,j,k) = tmp3(i,j,k) + tmp3(i,j,k+1)
              end do
            end do
            !$omp end do
          end do

        end if

      end if

      ! Calculate the 4th order scalar advection
      if (advopt == 2 .or. advopt == 3) then

        do k = 2, nk-2
          !$omp do schedule(runtime) private(i,j)
          do j = 2+jsouth, nj-2-jnorth
            do i = 1+iwest, ni-1-ieast
              tmp1(i,j,k) = (rstxu(i,j,k) + rstxu(i+1,j,k)) &
                   * (s(i+1,j,k) - s(i-1,j,k)) * dxv24
            end do
          end do
          !$omp end do

          !$omp do schedule(runtime) private(i,j)
          do j = 1+jsouth, nj-1-jnorth
            do i = 2+iwest, ni-2-ieast
              tmp2(i,j,k) = (rstxv(i,j,k) + rstxv(i,j+1,k)) &
                   * (s(i,j+1,k) - s(i,j-1,k)) * dyv24
            end do
          end do
          !$omp end do
        end do

        do k = 2, nk-2
          !$omp do schedule(runtime) private(i,j)
          do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-2-ieast
              sfrc(i,j,k) = fourd3 * sfrc(i,j,k) &
                   + ((tmp1(i-1,j,k) + tmp1(i+1,j,k)) &
                   + (tmp2(i,j-1,k) + tmp2(i,j+1,k)))
            end do
          end do
          !$omp end do
        end do

        if (advopt == 2) then

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                sfrc(i,j,k) = sfrc(i,j,k) + (tmp3(i,j,k) + tmp3(i,j,k+1))
              end do
            end do
            !$omp end do
          end do

        else if (advopt == 3) then

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2+jsouth, nj-2-jnorth
              do i = 2+iwest, ni-2-ieast
                tmp3(i,j,k) = (rstxwc(i,j,k) + rstxwc(i,j,k+1)) &
                     * (s(i,j,k+1) - s(i,j,k-1)) * dzv24
              end do
            end do
            !$omp end do
          end do

          do k = 3, nk-3
            !$omp do schedule(runtime) private(i,j)
            do j = 2+jsouth, nj-2-jnorth
              do i = 2+iwest, ni-2-ieast
                vadv(i,j,k) = fourd3 * vadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k+1))
              end do
            end do
            !$omp end do
          end do

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                sfrc(i,j,k) = sfrc(i,j,k) + vadv(i,j,k)
              end do
            end do
            !$omp end do
          end do

        end if

      end if

    else

      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-2
            sfrc(i,j,k) = 0.0e0
          end do
        end do
        !$omp end do
      end do

    end if

    !$omp end parallel

  end subroutine kernel_advs

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
  ! Parameter reader
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk, advopt, &
       iwest, ieast, jsouth, jnorth, dxiv, dyiv, dziv)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, advopt
    integer, intent(out) :: iwest, ieast, jsouth, jnorth
    real, intent(out) :: dxiv, dyiv, dziv

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
          case ('advopt')
            read(val, *) advopt
          case ('iwest')
            read(val, *) iwest
          case ('ieast')
            read(val, *) ieast
          case ('jsouth')
            read(val, *) jsouth
          case ('jnorth')
            read(val, *) jnorth
          case ('dxiv')
            read(val, *) dxiv
          case ('dyiv')
            read(val, *) dyiv
          case ('dziv')
            read(val, *) dziv
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

end program kernel_benchmark_advs

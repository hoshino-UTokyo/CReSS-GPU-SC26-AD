!***********************************************************************
! GPU Kernel Benchmark: diver2d (s_diver2d)
!***********************************************************************
!
! Description: Calculate 2D horizontal negative divergence using velocity
!              components u,v weighted by Jacobian and map scale factors
!              GPU version using OpenACC
!
!***********************************************************************
program kernel_benchmark_diver2d
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: mpopt, mfcopt
  real :: dxiv, dyiv

  ! Input arrays
  real, allocatable :: mf(:,:)
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: rmf8u(:,:,:)
  real, allocatable :: rmf8v(:,:,:)
  real, allocatable :: var8u(:,:,:)
  real, allocatable :: var8v(:,:,:)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)

  ! Output array
  real, allocatable :: div2d(:,:,:)

  ! Work arrays
  real, allocatable :: tmp1(:,:,:)
  real, allocatable :: tmp2(:,:,:)

  ! Reference output for validation
  real, allocatable :: div2d_ref(:,:,:)

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
       ni, nj, nk, mpopt, mfcopt, dxiv, dyiv)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: diver2d'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' Options: mpopt=', mpopt, ', mfcopt=', mfcopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rmf8u(0:ni+1, 0:nj+1, 1:3))
  allocate(rmf8v(0:ni+1, 0:nj+1, 1:3))
  allocate(var8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(div2d(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(div2d_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_2d(trim(data_dir)//'/mf.bin', mf, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call read_array_3d(trim(data_dir)//'/rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call read_array_3d(trim(data_dir)//'/rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  call read_array_3d(trim(data_dir)//'/var8u.bin', var8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/var8v.bin', var8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/div2d_ref.bin', div2d_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_diver2d(mpopt, mfcopt, dxiv, dyiv, ni, nj, nk, &
         mf, rmf, rmf8u, rmf8v, var8u, var8v, u, v, div2d, tmp1, tmp2)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    call cpu_time(t_start)

    call kernel_diver2d(mpopt, mfcopt, dxiv, dyiv, ni, nj, nk, &
         mf, rmf, rmf8u, rmf8v, var8u, var8v, u, v, div2d, tmp1, tmp2)
    !$acc wait

    call cpu_time(t_end)
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

  do k = 1, nk-1
    do j = 1, nj-1
      do i = 1, ni-1
        rel_error = abs(div2d(i,j,k) - div2d_ref(i,j,k))
        if (abs(div2d_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(div2d_ref(i,j,k))
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
  deallocate(mf, rmf, rmf8u, rmf8v)
  deallocate(var8u, var8v, u, v)
  deallocate(div2d, tmp1, tmp2)
  deallocate(div2d_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: diver2d - GPU version using OpenACC
  !=====================================================================
  subroutine kernel_diver2d(mpopt, mfcopt, dxiv, dyiv, ni, nj, nk, &
       mf, rmf, rmf8u, rmf8v, var8u, var8v, u, v, div2d, tmp1, tmp2)
    implicit none

    integer, intent(in) :: mpopt, mfcopt
    real, intent(in) :: dxiv, dyiv
    integer, intent(in) :: ni, nj, nk

    real, intent(in) :: mf(0:ni+1, 0:nj+1)
    real, intent(in) :: rmf(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: rmf8u(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: rmf8v(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: var8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: var8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)

    real, intent(out) :: div2d(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    ! Phase 1: Multiply velocities by optional variables and map factors

    if (mfcopt == 0) then

      !$acc kernels
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni
            tmp1(i,j,k) = var8u(i,j,k) * u(i,j,k)
          end do
        end do
        !$acc loop independent
        do j = 1, nj
          !$acc loop independent
          do i = 1, ni-1
            tmp2(i,j,k) = var8v(i,j,k) * v(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

    else

      if (mpopt == 0 .or. mpopt == 10) then

        !$acc kernels
        do k = 1, nk-1
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni
              tmp1(i,j,k) = var8u(i,j,k) * u(i,j,k)
            end do
          end do
          !$acc loop independent
          do j = 1, nj
            !$acc loop independent
            do i = 1, ni-1
              tmp2(i,j,k) = rmf8v(i,j,2) * var8v(i,j,k) * v(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else if (mpopt == 5) then

        !$acc kernels
        do k = 1, nk-1
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni
              tmp1(i,j,k) = rmf8u(i,j,2) * var8u(i,j,k) * u(i,j,k)
            end do
          end do
          !$acc loop independent
          do j = 1, nj
            !$acc loop independent
            do i = 1, ni-1
              tmp2(i,j,k) = var8v(i,j,k) * v(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        do k = 1, nk-1
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni
              tmp1(i,j,k) = rmf8u(i,j,2) * var8u(i,j,k) * u(i,j,k)
            end do
          end do
          !$acc loop independent
          do j = 1, nj
            !$acc loop independent
            do i = 1, ni-1
              tmp2(i,j,k) = rmf8v(i,j,2) * var8v(i,j,k) * v(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      end if

    end if

    ! Phase 2: Calculate divergence

    if (mfcopt == 0) then

      !$acc kernels
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            div2d(i,j,k) = (tmp1(i,j,k) - tmp1(i+1,j,k)) * dxiv &
                 + (tmp2(i,j,k) - tmp2(i,j+1,k)) * dyiv
          end do
        end do
      end do
      !$acc end kernels

    else

      if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then

        !$acc kernels
        do k = 1, nk-1
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni-1
              div2d(i,j,k) = mf(i,j) * ((tmp1(i,j,k) - tmp1(i+1,j,k)) * dxiv &
                   + (tmp2(i,j,k) - tmp2(i,j+1,k)) * dyiv)
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        do k = 1, nk-1
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni-1
              div2d(i,j,k) = rmf(i,j,1) * ((tmp1(i,j,k) - tmp1(i+1,j,k)) * dxiv &
                   + (tmp2(i,j,k) - tmp2(i,j+1,k)) * dyiv)
            end do
          end do
        end do
        !$acc end kernels

      end if

    end if

  end subroutine kernel_diver2d

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    integer :: ios
    logical :: exists

    ! Default values
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
  subroutine read_parameters(filename, ni, nj, nk, mpopt, mfcopt, dxiv, dyiv)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, mpopt, mfcopt
    real, intent(out) :: dxiv, dyiv

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
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
          case ('dxiv')
            read(val, *) dxiv
          case ('dyiv')
            read(val, *) dyiv
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

end program kernel_benchmark_diver2d

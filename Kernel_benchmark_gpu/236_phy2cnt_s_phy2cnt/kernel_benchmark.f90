!***********************************************************************
! GPU Kernel Benchmark: phy2cnt (s_phy2cnt)
!***********************************************************************
!
! Source: Src/phy2cnt.f90
! Description: Calculate zeta components of contravariant velocity from
!              physical velocity components with terrain-following coords.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_phy2cnt
  use omp_lib
  implicit none

  ! Grid dimensions
  integer :: ni, nj, nk

  ! Physics options
  integer :: sthopt   ! Vertical grid stretching option
  integer :: trnopt   ! Terrain height option
  integer :: mpopt    ! Map projection option
  integer :: mfcopt   ! Map scale factor option

  ! Arrays
  real, allocatable :: j31(:,:,:)     ! z-x Jacobian
  real, allocatable :: j32(:,:,:)     ! z-y Jacobian
  real, allocatable :: jcb8w(:,:,:)   ! Jacobian at w points
  real, allocatable :: mf(:,:)        ! Map scale factors
  real, allocatable :: u(:,:,:)       ! x velocity
  real, allocatable :: v(:,:,:)       ! y velocity
  real, allocatable :: w(:,:,:)       ! z velocity
  real, allocatable :: wc(:,:,:)      ! contravariant velocity (output)
  real, allocatable :: wc_ref(:,:,:)  ! Reference output

  ! Work arrays
  real, allocatable :: mf25(:,:)      ! 0.25 x mf
  real, allocatable :: j31u2(:,:,:)   ! 2.0 x j31 x u
  real, allocatable :: j32v2(:,:,:)   ! 2.0 x j32 x v

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
       ni, nj, nk, sthopt, trnopt, mpopt, mfcopt)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: phy2cnt'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' sthopt: ', sthopt
  write(*,'(A,I6)') ' trnopt: ', trnopt
  write(*,'(A,I6)') ' mpopt:  ', mpopt
  write(*,'(A,I6)') ' mfcopt: ', mfcopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(w(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf25(0:ni+1, 0:nj+1))
  allocate(j31u2(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32v2(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8w.bin', jcb8w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/mf.bin', mf, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/w.bin', w, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/wc_ref.bin', wc_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    wc = 0.0
    j31u2 = 0.0
    j32v2 = 0.0
    mf25 = 0.0
    call kernel_phy2cnt(sthopt, trnopt, mpopt, mfcopt, ni, nj, nk, &
         j31, j32, jcb8w, mf, u, v, w, wc, mf25, j31u2, j32v2)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    wc = 0.0
    j31u2 = 0.0
    j32v2 = 0.0
    mf25 = 0.0

    !$acc wait
    t_start = omp_get_wtime()

    call kernel_phy2cnt(sthopt, trnopt, mpopt, mfcopt, ni, nj, nk, &
         j31, j32, jcb8w, mf, u, v, w, wc, mf25, j31u2, j32v2)

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

  do k = 2, nk-1
    do j = 1, nj-1
      do i = 1, ni-1
        rel_error = abs(wc(i,j,k) - wc_ref(i,j,k))
        if (abs(wc_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(wc_ref(i,j,k))
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
  deallocate(j31, j32, jcb8w, mf, u, v, w, wc, wc_ref)
  deallocate(mf25, j31u2, j32v2, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: phy2cnt (OpenACC version)
  !=====================================================================
  subroutine kernel_phy2cnt(sthopt, trnopt, mpopt, mfcopt, ni, nj, nk, &
       j31, j32, jcb8w, mf, u, v, w, wc, mf25, j31u2, j32v2)
    implicit none

    integer, intent(in) :: sthopt, trnopt, mpopt, mfcopt
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: j31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: mf(0:ni+1, 0:nj+1)
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: wc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: mf25(0:ni+1, 0:nj+1)
    real, intent(inout) :: j31u2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: j32v2(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    if(trnopt.eq.0) then

      if(sthopt.eq.0) then

        !$acc kernels
        !$acc loop independent collapse(3)
        do k=2,nk-1
          do j=1,nj-1
            do i=1,ni-1
              wc(i,j,k)=w(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else if(sthopt.ge.1) then

        !$acc kernels
        !$acc loop independent collapse(3)
        do k=2,nk-1
          do j=1,nj-1
            do i=1,ni-1
              wc(i,j,k)=w(i,j,k)/jcb8w(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      end if

    else

      !$acc kernels
      !$acc loop independent collapse(3)
      do k=2,nk-1
        do j=1,nj-1
          do i=1,ni
            j31u2(i,j,k)=(u(i,j,k-1)+u(i,j,k))*j31(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k=2,nk-1
        do j=1,nj
          do i=1,ni-1
            j32v2(i,j,k)=(v(i,j,k-1)+v(i,j,k))*j32(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      if(mfcopt.eq.0) then

        !$acc kernels
        !$acc loop independent collapse(3)
        do k=2,nk-1
          do j=1,nj-1
            do i=1,ni-1
              wc(i,j,k)=(.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k)) &
                   +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else

        if(mpopt.eq.0.or.mpopt.eq.10) then

          !$acc kernels
          !$acc loop independent collapse(3)
          do k=2,nk-1
            do j=1,nj-1
              do i=1,ni-1
                wc(i,j,k)=(.25e0*(mf(i,j)*(j31u2(i,j,k)+j31u2(i+1,j,k)) &
                     +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

        else if(mpopt.eq.5) then

          !$acc kernels
          !$acc loop independent collapse(3)
          do k=2,nk-1
            do j=1,nj-1
              do i=1,ni-1
                wc(i,j,k)=(.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k)) &
                     +mf(i,j)*(j32v2(i,j,k)+j32v2(i,j+1,k))) &
                     +w(i,j,k))/jcb8w(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

        else

          !$acc kernels
          !$acc loop independent collapse(2)
          do j=1,nj-1
            do i=1,ni-1
              mf25(i,j)=.25e0*mf(i,j)
            end do
          end do
          !$acc end kernels

          !$acc kernels
          !$acc loop independent collapse(3)
          do k=2,nk-1
            do j=1,nj-1
              do i=1,ni-1
                wc(i,j,k)=(mf25(i,j)*((j31u2(i,j,k)+j31u2(i+1,j,k)) &
                     +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

        end if

      end if

    end if

  end subroutine kernel_phy2cnt

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
  subroutine read_parameters(filename, ni, nj, nk, sthopt, trnopt, mpopt, mfcopt)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    integer, intent(out) :: sthopt, trnopt, mpopt, mfcopt

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
          case ('sthopt')
            read(val, *) sthopt
          case ('trnopt')
            read(val, *) trnopt
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary 2D array reader
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

  !=====================================================================
  ! Binary 3D array reader
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

end program kernel_benchmark_gpu_phy2cnt

!***********************************************************************
! GPU Kernel Benchmark: stepwi (s_stepwi) - Matrix setup part
!***********************************************************************
!
! Source: Src/stepwi.f90
! Description: Prepares coefficient matrices for vertical implicit solver
!              of w-equation, computing tridiagonal matrix elements.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_stepwi
  use omp_lib
  implicit none

  ! Physical constant (read from params.txt)
  real :: g  ! Gravitational acceleration

  ! Grid dimensions
  integer :: ni, nj, nk

  ! Parameters
  real :: dts      ! Small time step
  real :: dziv     ! Inverse of dz
  real :: weicoe   ! Weighting coefficient
  integer :: buyopt ! Buoyancy option

  ! Arrays
  real, allocatable :: rst8w(:,:,:), rbr(:,:,:), rcsq(:,:,:)
  real, allocatable :: jcb(:,:,:), rst(:,:,:)
  real, allocatable :: fw(:,:,:), wf(:,:,:), wc(:,:,:)
  real, allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:)
  real, allocatable :: fw_input(:,:,:), wf_input(:,:,:)
  real, allocatable :: tmp1_ref(:,:,:), tmp2_ref(:,:,:), tmp3_ref(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, dts, dziv, weicoe, buyopt, g)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: stepwi (matrix setup)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dts:    ', dts
  write(*,'(A,ES12.4)') ' dziv:   ', dziv
  write(*,'(A,ES12.4)') ' weicoe: ', weicoe
  write(*,'(A,ES12.4)') ' g:      ', g
  write(*,'(A,I6)') ' buyopt: ', buyopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rst8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rcsq(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(fw(0:ni+1, 0:nj+1, 1:nk))
  allocate(wf(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
  allocate(fw_input(0:ni+1, 0:nj+1, 1:nk))
  allocate(wf_input(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rst8w.bin', rst8w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rcsq.bin', rcsq, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/fw_in.bin', fw_input, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wf_in.bin', wf_input, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/tmp1_ref.bin', tmp1_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp2_ref.bin', tmp2_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp3_ref.bin', tmp3_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    fw = fw_input
    wf = wf_input
    tmp1 = 0.0; tmp2 = 0.0; tmp3 = 0.0; wc = 0.0
    call kernel_stepwi(buyopt, dts, dziv, weicoe, ni, nj, nk, &
         rst8w, rbr, rcsq, jcb, rst, fw, wf, wc, tmp1, tmp2, tmp3)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    fw = fw_input
    wf = wf_input
    tmp1 = 0.0; tmp2 = 0.0; tmp3 = 0.0; wc = 0.0

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_stepwi(buyopt, dts, dziv, weicoe, ni, nj, nk, &
         rst8w, rbr, rcsq, jcb, rst, fw, wf, wc, tmp1, tmp2, tmp3)
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

  ! Validate tmp1, tmp2, tmp3 (tridiagonal matrix coefficients)
  do k = 3, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        ! tmp1
        rel_error = abs(tmp1(i,j,k) - tmp1_ref(i,j,k))
        if (abs(tmp1_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(tmp1_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1

        ! tmp2
        rel_error = abs(tmp2(i,j,k) - tmp2_ref(i,j,k))
        if (abs(tmp2_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(tmp2_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1

        ! tmp3
        rel_error = abs(tmp3(i,j,k) - tmp3_ref(i,j,k))
        if (abs(tmp3_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(tmp3_ref(i,j,k))
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
  deallocate(rst8w, rbr, rcsq, jcb, rst, fw, wf, wc)
  deallocate(tmp1, tmp2, tmp3, fw_input, wf_input)
  deallocate(tmp1_ref, tmp2_ref, tmp3_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! GPU Kernel: stepwi (matrix coefficient setup) - OpenACC
  !=====================================================================
  subroutine kernel_stepwi(buyopt, dts, dziv, weicoe, ni, nj, nk, &
       rst8w, rbr, rcsq, jcb, rst, fw, wf, wc, tmp1, tmp2, tmp3)
    implicit none

    integer, intent(in) :: buyopt
    real, intent(in) :: dts, dziv, weicoe
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: rst8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rcsq(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: fw(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: wf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: wc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: tmp3(0:ni+1, 0:nj+1, 1:nk)

    real :: g05, sbsqzi, sbsqg5
    real :: rstiv, mm, nn, a, b
    integer :: i, j, k

    g05 = 0.5 * g
    sbsqzi = dts * dts * weicoe * weicoe * dziv
    sbsqg5 = 0.5 * g * dts * dts * weicoe * weicoe

    ! Step 1: Update wf
    !$acc kernels
    !$acc loop independent
    do k = 3, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          wf(i,j,k) = wf(i,j,k) + dts * fw(i,j,k) / rst8w(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Step 2: Compute intermediate values
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          tmp1(i,j,k) = g05 * rbr(i,j,k) / rcsq(i,j,k)
          tmp2(i,j,k) = dziv / jcb(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          fw(i,j,k) = tmp1(i,j,k) + tmp2(i,j,k)
          wc(i,j,k) = tmp1(i,j,k) - tmp2(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Step 3: Compute tridiagonal matrix coefficients
    if (buyopt == 0) then

      !$acc kernels
      !$acc loop independent
      do k = 3, nk-2
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            a = sbsqzi / rst8w(i,j,k)
            mm = a * rcsq(i,j,k-1)
            nn = a * rcsq(i,j,k)
            tmp1(i,j,k) = -mm * fw(i,j,k-1)
            tmp2(i,j,k) = 1.0 + (nn * fw(i,j,k) - mm * wc(i,j,k-1))
            tmp3(i,j,k) = nn * wc(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

    else if (buyopt == 1) then

      !$acc kernels
      !$acc loop independent
      do k = 3, nk-2
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            rstiv = 1.0 / rst8w(i,j,k)
            a = rstiv * sbsqzi
            b = rstiv * sbsqg5
            mm = b * rst(i,j,k-1) - a * rcsq(i,j,k-1)
            nn = b * rst(i,j,k) + a * rcsq(i,j,k)
            tmp1(i,j,k) = mm * fw(i,j,k-1)
            tmp2(i,j,k) = 1.0 + (mm * wc(i,j,k-1) + nn * fw(i,j,k))
            tmp3(i,j,k) = nn * wc(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

    end if

  end subroutine kernel_stepwi

  !=====================================================================
  ! I/O subroutines
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

  subroutine read_parameters(filename, ni, nj, nk, dts, dziv, weicoe, buyopt, g)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, buyopt
    real, intent(out) :: dts, dziv, weicoe, g
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
          case ('dts')
            read(val, *) dts
          case ('dziv')
            read(val, *) dziv
          case ('weicoe')
            read(val, *) weicoe
          case ('buyopt')
            read(val, *) buyopt
          case ('g')
            read(val, *) g
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
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

end program kernel_benchmark_gpu_stepwi

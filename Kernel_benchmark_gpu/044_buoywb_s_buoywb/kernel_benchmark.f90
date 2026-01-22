!***********************************************************************
! GPU Kernel Benchmark: buoywb (s_buoywb)
!***********************************************************************
!
! Description: Calculate buoyancy forcing for vertical velocity in large
!              time step integration. GPU version using OpenACC.
!
!***********************************************************************
program kernel_benchmark_buoywb
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: gwmopt, cphopt
  real, parameter :: g = 9.8        ! Gravity acceleration
  real, parameter :: epsva = 0.622  ! rd/rv ratio

  ! Input arrays
  real, allocatable :: ptbr(:,:,:), qvbr(:,:,:), rst(:,:,:)
  real, allocatable :: ptp(:,:,:), qv(:,:,:), qall(:,:,:)

  ! Input/output arrays
  real, allocatable :: wfrc(:,:,:), qvd(:,:,:), wb8s(:,:,:)

  ! Reference output for validation
  real, allocatable :: wfrc_ref(:,:,:), qvd_ref(:,:,:), wb8s_ref(:,:,:)

  ! Backup arrays for iteration
  real, allocatable :: wfrc_in(:,:,:), qvd_in(:,:,:), wb8s_in(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, arr_error
  real :: tolerance
  integer :: error_count
  logical :: validation_passed

  ! Loop variables
  integer :: iter

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', gwmopt, cphopt, ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: buoywb'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' Options: gwmopt=', gwmopt, ', cphopt=', cphopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk), qvbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk), ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk), qall(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc(0:ni+1, 0:nj+1, 1:nk), qvd(0:ni+1, 0:nj+1, 1:nk))
  allocate(wb8s(0:ni+1, 0:nj+1, 1:nk))

  allocate(wfrc_ref(0:ni+1, 0:nj+1, 1:nk), qvd_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(wb8s_ref(0:ni+1, 0:nj+1, 1:nk))

  allocate(wfrc_in(0:ni+1, 0:nj+1, 1:nk), qvd_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(wb8s_in(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvbr.bin', qvbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qall.bin', qall, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wfrc_in.bin', wfrc_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvd_in.bin', qvd_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wb8s_in.bin', wb8s_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/wfrc_ref.bin', wfrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvd_ref.bin', qvd_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wb8s_ref.bin', wb8s_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    wfrc = wfrc_in; qvd = qvd_in; wb8s = wb8s_in
    call kernel_buoywb(gwmopt, cphopt, ni, nj, nk, g, epsva, &
                       ptbr, qvbr, rst, ptp, qv, qall, wfrc, qvd, wb8s)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    wfrc = wfrc_in; qvd = qvd_in; wb8s = wb8s_in

    call cpu_time(t_start)
    call kernel_buoywb(gwmopt, cphopt, ni, nj, nk, g, epsva, &
                       ptbr, qvbr, rst, ptp, qv, qall, wfrc, qvd, wb8s)
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

  call validate_array(wfrc, wfrc_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error
  call validate_array(qvd, qvd_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error
  call validate_array(wb8s, wb8s_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error

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
  deallocate(ptbr, qvbr, rst, ptp, qv, qall, wfrc, qvd, wb8s)
  deallocate(wfrc_ref, qvd_ref, wb8s_ref)
  deallocate(wfrc_in, qvd_in, wb8s_in)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: buoywb - GPU version using OpenACC
  !=====================================================================
  subroutine kernel_buoywb(gwmopt, cphopt, ni, nj, nk, g, epsva, &
                           ptbr, qvbr, rst, ptp, qv, qall, wfrc, qvd, wb8s)
    implicit none

    integer, intent(in) :: gwmopt, cphopt
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: g, epsva
    real, intent(in) :: ptbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qvbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qall(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: wfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: qvd(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: wb8s(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    real :: g05

    g05 = 0.5 * g

    ! For moist air case - compute qvd first
    !$acc kernels
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          qvd(i,j,k) = qv(i,j,k) - qvbr(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Compute wb8s based on cphopt and gwmopt
    if (abs(cphopt) == 0) then
      ! No cloud microphysics
      if (gwmopt == 0) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-2
              wb8s(i,j,k) = g05 * rst(i,j,k) &
                * (qvd(i,j,k)/(epsva + qvbr(i,j,k)) &
                - qvd(i,j,k)/(1.0 + qvbr(i,j,k)) + ptp(i,j,k)/ptbr(i,j,k))
            end do
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-2
              wb8s(i,j,k) = g05 * rst(i,j,k) &
                * (qvd(i,j,k)/(epsva + qvbr(i,j,k)) &
                - qvd(i,j,k)/(1.0 + qvbr(i,j,k)))
            end do
          end do
        end do
        !$acc end kernels
      end if
    else
      ! With cloud microphysics
      if (gwmopt == 0) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-2
              wb8s(i,j,k) = g05 * rst(i,j,k) &
                * (qvd(i,j,k)/(epsva + qvbr(i,j,k)) &
                - (qvd(i,j,k) + qall(i,j,k))/(1.0 + qvbr(i,j,k)) &
                + ptp(i,j,k)/ptbr(i,j,k))
            end do
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-2
              wb8s(i,j,k) = g05 * rst(i,j,k) &
                * (qvd(i,j,k)/(epsva + qvbr(i,j,k)) &
                - (qvd(i,j,k) + qall(i,j,k))/(1.0 + qvbr(i,j,k)))
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! Finally be averaged vertically
    !$acc kernels
    do k = 3, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          wfrc(i,j,k) = wfrc(i,j,k) + (wb8s(i,j,k-1) + wb8s(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_buoywb

  !=====================================================================
  ! Validate array
  !=====================================================================
  subroutine validate_array(arr, arr_ref, ni, nj, nk, tol, max_err, err_count)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: arr_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(inout) :: err_count

    integer :: i, j, k
    real :: rel_error

    max_err = 0.0
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          rel_error = abs(arr(i,j,k) - arr_ref(i,j,k))
          if (abs(arr_ref(i,j,k)) > 1.0e-10) then
            rel_error = rel_error / abs(arr_ref(i,j,k))
          end if
          if (rel_error > max_err) max_err = rel_error
          if (rel_error > tol) err_count = err_count + 1
        end do
      end do
    end do
  end subroutine validate_array

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
  subroutine read_parameters(filename, gwmopt, cphopt, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: gwmopt, cphopt, ni, nj, nk

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    gwmopt = 0
    cphopt = 0
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
          case ('gwmopt')
            read(val, *) gwmopt
          case ('cphopt')
            read(val, *) cphopt
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

end program kernel_benchmark_buoywb

!***********************************************************************
! GPU Kernel Benchmark: buoywsi (s_buoywsi)
!***********************************************************************
!
! Source: Src/buoywsi.f90
! Description: Calculate buoyancy for small time steps integration
!              using horizontally explicit and vertically implicit method
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_buoywsi
  use omp_lib
  implicit none

  ! Physical constants (read from params.txt)
  real :: g

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: gwmopt
  real :: weicoe, dts

  ! Input arrays
  real, allocatable :: rbr(:,:,:)
  real, allocatable :: ptbr(:,:,:)
  real, allocatable :: rst(:,:,:)
  real, allocatable :: rcsq(:,:,:)
  real, allocatable :: pp(:,:,:)
  real, allocatable :: ptp(:,:,:)
  real, allocatable :: fp(:,:,:)

  ! Output array
  real, allocatable :: fw(:,:,:)

  ! Work array
  real, allocatable :: wb8s(:,:,:)

  ! Reference output for validation
  real, allocatable :: fw_ref(:,:,:)

  ! Initial value for re-initialization
  real, allocatable :: fw_init(:,:,:)

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
       ni, nj, nk, gwmopt, weicoe, dts, g)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: buoywsi'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' gwmopt=', gwmopt
  write(*,'(A,F10.6)') ' weicoe=', weicoe
  write(*,'(A,F10.6)') ' g=', g
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(rcsq(0:ni+1, 0:nj+1, 1:nk))
  allocate(pp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(fp(0:ni+1, 0:nj+1, 1:nk))
  allocate(fw(0:ni+1, 0:nj+1, 1:nk))
  allocate(fw_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(fw_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(wb8s(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rcsq.bin', rcsq, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pp.bin', pp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/fp.bin', fp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/fw_in.bin', fw_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/fw_ref.bin', fw_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    fw = fw_init
    call kernel_buoywsi(gwmopt, weicoe, dts, ni, nj, nk, &
         rbr, ptbr, rst, rcsq, pp, ptp, fp, fw, wb8s)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    fw = fw_init

    !$acc wait
    t_start = omp_get_wtime()

    call kernel_buoywsi(gwmopt, weicoe, dts, ni, nj, nk, &
         rbr, ptbr, rst, rcsq, pp, ptp, fp, fw, wb8s)

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

  do k = 3, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(fw(i,j,k) - fw_ref(i,j,k))
        if (abs(fw_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(fw_ref(i,j,k))
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
  deallocate(rbr, ptbr, rst, rcsq, pp, ptp, fp)
  deallocate(fw, fw_init, fw_ref, wb8s)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: buoywsi (OpenACC version)
  ! Calculate buoyancy for small time steps
  !=====================================================================
  subroutine kernel_buoywsi(gwmopt, weicoe, dts, ni, nj, nk, &
       rbr, ptbr, rst, rcsq, pp, ptp, fp, fw, wb8s)
    implicit none

    integer, intent(in) :: gwmopt
    real, intent(in) :: weicoe, dts
    integer, intent(in) :: ni, nj, nk

    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rcsq(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: pp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: fp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: fw(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: wb8s(0:ni+1, 0:nj+1, 1:nk)

    ! Local variables
    real :: dtw, g05n
    integer :: i, j, k

    dtw = dts * weicoe
    g05n = -0.5e0 * g

    if (gwmopt == 0) then

      !$acc kernels
      !$acc loop independent
      do k = 2, nk-2
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            wb8s(i,j,k) = (pp(i,j,k) * rst(i,j,k) + fp(i,j,k) * rbr(i,j,k) * dtw) &
                 / rcsq(i,j,k) * g05n
          end do
        end do
      end do
      !$acc end kernels

    else

      !$acc kernels
      !$acc loop independent
      do k = 2, nk-2
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            wb8s(i,j,k) = ((pp(i,j,k) * rst(i,j,k) + fp(i,j,k) * rbr(i,j,k) * dtw) &
                 / rcsq(i,j,k) - ptp(i,j,k) * rst(i,j,k) / ptbr(i,j,k)) * g05n
          end do
        end do
      end do
      !$acc end kernels

    end if

    !$acc kernels
    !$acc loop independent
    do k = 3, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          fw(i,j,k) = fw(i,j,k) + (wb8s(i,j,k-1) + wb8s(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_buoywsi

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
  ! Parameter reader (key=value format)
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk, gwmopt, weicoe, dts, g)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, gwmopt
    real, intent(out) :: weicoe, dts, g

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
          case ('gwmopt')
            read(val, *) gwmopt
          case ('weicoe')
            read(val, *) weicoe
          case ('dts')
            read(val, *) dts
          case ('g')
            read(val, *) g
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

end program kernel_benchmark_gpu_buoywsi

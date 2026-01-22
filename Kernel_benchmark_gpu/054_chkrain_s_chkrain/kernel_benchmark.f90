!***********************************************************************
! GPU Kernel Benchmark: chkrain (s_chkrain)
!***********************************************************************
!
! Source: Src/chkrain.f90
! Description: Check precipitation on surface by comparing precipitation
!              rates against threshold (prmin). Sets fall flag to 1.0 if
!              precipitation detected, -1.0 otherwise.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_chkrain
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nqw, nqi

  ! Parameters
  integer :: cphopt, haiopt

  ! Physical constants from m_comphy
  real, parameter :: prmin = 1.389e-7  ! Minimum precipitation rate

  ! Input arrays
  real, allocatable :: prwtr(:,:,:,:)  ! Precipitation for water
  real, allocatable :: price(:,:,:,:)  ! Precipitation for ice

  ! Output array
  real, allocatable :: fall(:,:)

  ! Reference output for validation
  real, allocatable :: fall_ref(:,:)

  ! Backup arrays for iteration
  real, allocatable :: prwtr_in(:,:,:,:)
  real, allocatable :: price_in(:,:,:,:)

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
  integer :: iter, i, j

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', cphopt, haiopt, &
                       ni, nj, nqw, nqi)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: chkrain'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj
  write(*,'(A,I6,A,I6)') ' Categories: nqw=', nqw, ', nqi=', nqi
  write(*,'(A,I6,A,I6)') ' Options: cphopt=', cphopt, ', haiopt=', haiopt
  write(*,'(A,ES12.4)') ' prmin=', prmin
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(prwtr(0:ni+1, 0:nj+1, 1:2, 1:nqw))
  allocate(price(0:ni+1, 0:nj+1, 1:2, 1:nqi))
  allocate(fall(0:ni+1, 0:nj+1))
  allocate(fall_ref(0:ni+1, 0:nj+1))
  allocate(prwtr_in(0:ni+1, 0:nj+1, 1:2, 1:nqw))
  allocate(price_in(0:ni+1, 0:nj+1, 1:2, 1:nqi))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_4d(trim(data_dir)//'/prwtr.bin', prwtr_in, 0, ni+1, 0, nj+1, 1, 2, 1, nqw)
  call read_array_4d(trim(data_dir)//'/price.bin', price_in, 0, ni+1, 0, nj+1, 1, 2, 1, nqi)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_2d(trim(data_dir)//'/fall_ref.bin', fall_ref, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    prwtr = prwtr_in
    price = price_in
    fall = 0.0
    call kernel_chkrain(cphopt, haiopt, ni, nj, nqw, nqi, prwtr, price, fall)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    prwtr = prwtr_in
    price = price_in
    fall = 0.0

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_chkrain(cphopt, haiopt, ni, nj, nqw, nqi, prwtr, price, fall)
    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output (interior points only - kernel doesn't set boundaries)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0

  do j = 1, nj-1
    do i = 1, ni-1
      rel_error = abs(fall(i,j) - fall_ref(i,j))
      if (abs(fall_ref(i,j)) > 1.0e-10) then
        rel_error = rel_error / abs(fall_ref(i,j))
      end if
      if (rel_error > max_error) max_error = rel_error
      if (rel_error > tolerance) error_count = error_count + 1
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
  deallocate(prwtr, price, fall, fall_ref, prwtr_in, price_in, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: chkrain (OpenACC version)
  ! Check precipitation on surface - moist run with bulk method
  !=====================================================================
  subroutine kernel_chkrain(cphopt, haiopt, ni, nj, nqw, nqi, prwtr, price, fall)
    implicit none

    integer, intent(in) :: cphopt, haiopt
    integer, intent(in) :: ni, nj, nqw, nqi
    real, intent(in) :: prwtr(0:ni+1, 0:nj+1, 1:2, 1:nqw)
    real, intent(in) :: price(0:ni+1, 0:nj+1, 1:2, 1:nqi)
    real, intent(out) :: fall(0:ni+1, 0:nj+1)

    integer :: i, j

    if (abs(cphopt) >= 2 .and. abs(cphopt) < 10) then
      if (haiopt == 0) then
        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            if (prwtr(i,j,1,1) > prmin .or. &
                prwtr(i,j,1,2) > prmin .or. &
                price(i,j,1,1) > prmin .or. &
                price(i,j,1,2) > prmin .or. &
                price(i,j,1,3) > prmin) then
              fall(i,j) = 1.0e0
            else
              fall(i,j) = -1.0e0
            end if
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            if (prwtr(i,j,1,1) > prmin .or. &
                prwtr(i,j,1,2) > prmin .or. &
                price(i,j,1,1) > prmin .or. &
                price(i,j,1,2) > prmin .or. &
                price(i,j,1,3) > prmin .or. &
                price(i,j,1,4) > prmin) then
              fall(i,j) = 1.0e0
            else
              fall(i,j) = -1.0e0
            end if
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_chkrain

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
  subroutine read_parameters(filename, cphopt, haiopt, ni, nj, nqw, nqi)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: cphopt, haiopt, ni, nj, nqw, nqi

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    cphopt = 0
    haiopt = 0
    ni = 1
    nj = 1
    nqw = 2
    nqi = 3

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
          case ('cphopt')
            read(val, *) cphopt
          case ('haiopt')
            read(val, *) haiopt
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nqw')
            read(val, *) nqw
          case ('nqi')
            read(val, *) nqi
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

  subroutine read_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_4d

end program kernel_benchmark_gpu_chkrain

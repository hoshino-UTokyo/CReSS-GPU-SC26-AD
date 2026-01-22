!***********************************************************************
! GPU Kernel Benchmark: chksat (s_chksat)
!***********************************************************************
!
! Source: Src/chksat.f90
! Description: Check and avoid super saturation mixing ratio by
!              limiting qv to saturation value computed from pressure
!              and temperature fields.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_chksat
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: imin, imax, jmin, jmax, kmin, kmax

  ! Physical constants from m_comphy
  real :: t0, epsva, es0
  real, parameter :: rd = 287.e0
  real, parameter :: cp = 1004.e0
  real, parameter :: p0 = 1.e5
  real, parameter :: tlow = 233.16e0

  ! Input arrays
  real, allocatable :: pbr(:,:,:)
  real, allocatable :: ptbr(:,:,:)
  real, allocatable :: pp(:,:,:)
  real, allocatable :: ptp(:,:,:)

  ! Input/output array
  real, allocatable :: qv(:,:,:)

  ! Reference output for validation
  real, allocatable :: qv_ref(:,:,:)

  ! Backup array for iteration
  real, allocatable :: qv_in(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', imin, imax, jmin, jmax, &
                       kmin, kmax, t0, epsva, es0)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: chksat'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6,A,I6,A,I6,A,I6)') ' Array bounds: i=', imin, ':', imax, &
       ', j=', jmin, ':', jmax, ', k=', kmin, ':', kmax
  write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') ' t0=', t0, ', epsva=', epsva, ', es0=', es0
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(pbr(imin:imax, jmin:jmax, kmin:kmax))
  allocate(ptbr(imin:imax, jmin:jmax, kmin:kmax))
  allocate(pp(imin:imax, jmin:jmax, kmin:kmax))
  allocate(ptp(imin:imax, jmin:jmax, kmin:kmax))
  allocate(qv(imin:imax, jmin:jmax, kmin:kmax))
  allocate(qv_ref(imin:imax, jmin:jmax, kmin:kmax))
  allocate(qv_in(imin:imax, jmin:jmax, kmin:kmax))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/pbr.bin', pbr, imin, imax, jmin, jmax, kmin, kmax)
  call read_array_3d(trim(data_dir)//'/ptbr.bin', ptbr, imin, imax, jmin, jmax, kmin, kmax)
  call read_array_3d(trim(data_dir)//'/pp.bin', pp, imin, imax, jmin, jmax, kmin, kmax)
  call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, imin, imax, jmin, jmax, kmin, kmax)
  call read_array_3d(trim(data_dir)//'/qv_in.bin', qv_in, imin, imax, jmin, jmax, kmin, kmax)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/qv_ref.bin', qv_ref, imin, imax, jmin, jmax, kmin, kmax)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    qv = qv_in
    call kernel_chksat(imin, imax, jmin, jmax, kmin, kmax, &
                       t0, epsva, es0, pbr, ptbr, pp, ptp, qv)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    qv = qv_in

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_chksat(imin, imax, jmin, jmax, kmin, kmax, &
                       t0, epsva, es0, pbr, ptbr, pp, ptp, qv)
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

  do k = kmin, kmax
    do j = jmin, jmax-1
      do i = imin, imax-1
        rel_error = abs(qv(i,j,k) - qv_ref(i,j,k))
        if (abs(qv_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(qv_ref(i,j,k))
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
  deallocate(pbr, ptbr, pp, ptp, qv, qv_ref, qv_in, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: chksat (OpenACC version)
  ! Check and avoid super saturation
  !=====================================================================
  subroutine kernel_chksat(imin, imax, jmin, jmax, kmin, kmax, &
                           t0, epsva, es0, pbr, ptbr, pp, ptp, qv)
    implicit none

    integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax
    real, intent(in) :: t0, epsva, es0
    real, intent(in) :: pbr(imin:imax, jmin:jmax, kmin:kmax)
    real, intent(in) :: ptbr(imin:imax, jmin:jmax, kmin:kmax)
    real, intent(in) :: pp(imin:imax, jmin:jmax, kmin:kmax)
    real, intent(in) :: ptp(imin:imax, jmin:jmax, kmin:kmax)
    real, intent(inout) :: qv(imin:imax, jmin:jmax, kmin:kmax)

    ! Loop bounds
    integer :: istr, iend, jstr, jend, kstr, kend

    ! Derived constants
    real :: rddvcp, p0iv

    ! Local variables
    integer :: i, j, k
    real :: t, es, qvs, pres

    ! Set loop bounds
    istr = imin
    iend = imax - 1
    jstr = jmin
    jend = jmax - 1
    kstr = kmin
    kend = kmax

    ! Set derived constants
    rddvcp = rd / cp
    p0iv = 1.e0 / p0

    !$acc kernels
    !$acc loop independent
    do k = kstr, kend
      !$acc loop independent
      do j = jstr, jend
        !$acc loop independent
        do i = istr, iend
          pres = pbr(i,j,k) + pp(i,j,k)
          t = (ptbr(i,j,k) + ptp(i,j,k)) * exp(rddvcp * log(p0iv * pres))

          if (t > tlow) then
            es = es0 * exp(17.269e0 * (t - t0) / (t - 35.86e0))
          else
            es = es0 * exp(21.875e0 * (t - t0) / (t - 7.66e0))
          end if

          qvs = epsva * es / (pres - es)
          qv(i,j,k) = min(qv(i,j,k), qvs)
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_chksat

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
  subroutine read_parameters(filename, imin, imax, jmin, jmax, kmin, kmax, &
                             t0, epsva, es0)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: imin, imax, jmin, jmax, kmin, kmax
    real, intent(out) :: t0, epsva, es0

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    imin = 0
    imax = 1
    jmin = 0
    jmax = 1
    kmin = 1
    kmax = 1
    t0 = 273.16
    epsva = 0.622
    es0 = 610.78

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
          case ('imin')
            read(val, *) imin
          case ('imax')
            read(val, *) imax
          case ('jmin')
            read(val, *) jmin
          case ('jmax')
            read(val, *) jmax
          case ('kmin')
            read(val, *) kmin
          case ('kmax')
            read(val, *) kmax
          case ('t0')
            read(val, *) t0
          case ('epsva')
            read(val, *) epsva
          case ('es0')
            read(val, *) es0
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

end program kernel_benchmark_gpu_chksat

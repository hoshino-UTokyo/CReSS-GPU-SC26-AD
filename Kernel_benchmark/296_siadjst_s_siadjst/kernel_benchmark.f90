!***********************************************************************
! Kernel Benchmark: siadjst (s_siadjst)
!***********************************************************************
!
! Source: Src/siadjst.f90
! Description: Saturation adjustment for ice - converts between water vapor
!              and cloud ice based on saturation conditions at low temperatures
!
!***********************************************************************
program kernel_benchmark_siadjst
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Physical constants
  real :: thresq, cp, t0, epsva, es0, lv0, lf0, cwmci, mi0iv
  real, parameter :: tlow = 233.16e0

  ! Input arrays
  real, allocatable :: ptbr(:,:,:)  ! Base state potential temperature
  real, allocatable :: pi(:,:,:)    ! Exner function
  real, allocatable :: p(:,:,:)     ! Pressure

  ! Input/output arrays
  real, allocatable :: ptp(:,:,:)   ! Potential temperature perturbation
  real, allocatable :: qv(:,:,:)    ! Water vapor mixing ratio
  real, allocatable :: qi(:,:,:)    ! Cloud ice mixing ratio
  real, allocatable :: nci(:,:,:)   ! Concentration of cloud ice

  ! Reference outputs
  real, allocatable :: ptp_ref(:,:,:)
  real, allocatable :: qv_ref(:,:,:)
  real, allocatable :: qi_ref(:,:,:)
  real, allocatable :: nci_ref(:,:,:)

  ! Input copies for repeated runs
  real, allocatable :: ptp_in(:,:,:)
  real, allocatable :: qv_in(:,:,:)
  real, allocatable :: qi_in(:,:,:)
  real, allocatable :: nci_in(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error, tolerance
  integer :: error_count, total_errors
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
       ni, nj, nk, thresq, cp, t0, epsva, es0, lv0, lf0, cwmci, mi0iv)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: siadjst'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(pi(0:ni+1, 0:nj+1, 1:nk))
  allocate(p(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pi.bin', pi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/p.bin', p, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptp_in.bin', ptp_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv_in.bin', qv_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qi_in.bin', qi_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/nci_in.bin', nci_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ptp_ref.bin', ptp_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv_ref.bin', qv_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qi_ref.bin', qi_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/nci_ref.bin', nci_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ptp = ptp_in
    qv = qv_in
    qi = qi_in
    nci = nci_in
    call kernel_siadjst(ni, nj, nk, thresq, cp, t0, epsva, es0, lv0, lf0, &
         cwmci, mi0iv, ptbr, pi, p, ptp, qv, qi, nci)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ! Reset arrays to initial state
    ptp = ptp_in
    qv = qv_in
    qi = qi_in
    nci = nci_in

    t_start = omp_get_wtime()

    call kernel_siadjst(ni, nj, nk, thresq, cp, t0, epsva, es0, lv0, lf0, &
         cwmci, mi0iv, ptbr, pi, p, ptp, qv, qi, nci)

    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  total_errors = 0

  ! Validate ptp
  call validate_array(ptp, ptp_ref, tolerance, 'ptp', max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate qv
  call validate_array(qv, qv_ref, tolerance, 'qv', max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate qi
  call validate_array(qi, qi_ref, tolerance, 'qi', max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate nci
  call validate_array(nci, nci_ref, tolerance, 'nci', max_error, error_count)
  total_errors = total_errors + error_count

  validation_passed = (total_errors == 0)

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
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Total error count:  ', total_errors
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(ptbr, pi, p, ptp, qv, qi, nci)
  deallocate(ptp_ref, qv_ref, qi_ref, nci_ref)
  deallocate(ptp_in, qv_in, qi_in, nci_in, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: siadjst
  !=====================================================================
  subroutine kernel_siadjst(ni, nj, nk, thresq, cp, t0, epsva, es0, lv0, lf0, &
       cwmci, mi0iv, ptbr, pi, p, ptp, qv, qi, nci)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: thresq, cp, t0, epsva, es0, lv0, lf0, cwmci, mi0iv
    real, intent(in) :: ptbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: pi(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: p(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: ptp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: qv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: qi(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: nci(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    real :: t, tcel, esi, qvsi, lscpi, dqi, a, b

    !$omp parallel default(shared) private(k)

    do k = 1, nk-1

      !$omp do schedule(runtime) private(i,j,t,tcel,esi,qvsi,lscpi,dqi,a,b)

      do j = 1, nj-1
        do i = 1, ni-1
          t = (ptbr(i,j,k) + ptp(i,j,k)) * pi(i,j,k)

          if (t <= tlow) then

            tcel = t - t0
            a = 1.0e0 / (t - 7.66e0)
            b = a * tcel
            esi = es0 * exp(21.875e0 * b)
            qvsi = epsva * esi / (p(i,j,k) - esi)

            if (qi(i,j,k) > thresq .or. qv(i,j,k) > qvsi) then

              lscpi = (lv0 * exp((0.167e0 + 3.67e-4*t) * log(t0/t)) &
                   + (lf0 + cwmci*tcel)) / (cp * pi(i,j,k))

              dqi = (qvsi - qv(i,j,k)) &
                   / (1.0e0 + 21.875e0*a*(1.0e0-b)*qvsi*lscpi*pi(i,j,k))

              if (qi(i,j,k) > dqi) then
                if (qi(i,j,k) > thresq) then
                  nci(i,j,k) = nci(i,j,k) - dqi * nci(i,j,k) / qi(i,j,k)
                else
                  nci(i,j,k) = nci(i,j,k) - dqi * mi0iv
                end if
                ptp(i,j,k) = ptp(i,j,k) - dqi * lscpi
                qv(i,j,k) = qv(i,j,k) + dqi
                qi(i,j,k) = qi(i,j,k) - dqi
              else
                nci(i,j,k) = 0.0e0
                ptp(i,j,k) = ptp(i,j,k) - qi(i,j,k) * lscpi
                qv(i,j,k) = qv(i,j,k) + qi(i,j,k)
                qi(i,j,k) = 0.0e0
              end if

            end if

            ! Second iteration
            t = (ptbr(i,j,k) + ptp(i,j,k)) * pi(i,j,k)

            if (t <= tlow) then
              tcel = t - t0
              a = 1.0e0 / (t - 7.66e0)
              b = a * tcel
              esi = es0 * exp(21.875e0 * b)
              qvsi = epsva * esi / (p(i,j,k) - esi)

              if (qi(i,j,k) > thresq .or. qv(i,j,k) > qvsi) then

                lscpi = (lv0 * exp((0.167e0 + 3.67e-4*t) * log(t0/t)) &
                     + (lf0 + cwmci*tcel)) / (cp * pi(i,j,k))

                dqi = (qvsi - qv(i,j,k)) &
                     / (1.0e0 + 21.875e0*a*(1.0e0-b)*qvsi*lscpi*pi(i,j,k))

                if (qi(i,j,k) > dqi) then
                  if (qi(i,j,k) > thresq) then
                    nci(i,j,k) = nci(i,j,k) - dqi * nci(i,j,k) / qi(i,j,k)
                  else
                    nci(i,j,k) = nci(i,j,k) - dqi * mi0iv
                  end if
                  ptp(i,j,k) = ptp(i,j,k) - dqi * lscpi
                  qv(i,j,k) = qv(i,j,k) + dqi
                  qi(i,j,k) = qi(i,j,k) - dqi
                else
                  nci(i,j,k) = 0.0e0
                  ptp(i,j,k) = ptp(i,j,k) - qi(i,j,k) * lscpi
                  qv(i,j,k) = qv(i,j,k) + qi(i,j,k)
                  qi(i,j,k) = 0.0e0
                end if

              end if
            end if

          end if

        end do
      end do

      !$omp end do

    end do

    !$omp end parallel

  end subroutine kernel_siadjst

  !=====================================================================
  ! Validation helper
  !=====================================================================
  subroutine validate_array(arr, ref, tol, name, max_err, err_count)
    real, intent(in) :: arr(0:,0:,1:), ref(0:,0:,1:)
    real, intent(in) :: tol
    character(len=*), intent(in) :: name
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j, k, ni_l, nj_l, nk_l

    ni_l = ubound(arr,1) - 1
    nj_l = ubound(arr,2) - 1
    nk_l = ubound(arr,3)

    max_err = 0.0
    err_count = 0

    do k = 1, nk_l-1
      do j = 1, nj_l-1
        do i = 1, ni_l-1
          rel_err = abs(arr(i,j,k) - ref(i,j,k))
          if (abs(ref(i,j,k)) > 1.0e-20) then
            rel_err = rel_err / abs(ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

    write(*,'(A,A,A,I8,A,ES12.4)') '   ', name, ': errors=', err_count, ', max_rel_err=', max_err

  end subroutine validate_array

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
  subroutine read_parameters(filename, ni, nj, nk, thresq, cp, t0, epsva, &
       es0, lv0, lf0, cwmci, mi0iv)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    real, intent(out) :: thresq, cp, t0, epsva, es0, lv0, lf0, cwmci, mi0iv

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
          case ('thresq')
            read(val, *) thresq
          case ('cp')
            read(val, *) cp
          case ('t0')
            read(val, *) t0
          case ('epsva')
            read(val, *) epsva
          case ('es0')
            read(val, *) es0
          case ('lv0')
            read(val, *) lv0
          case ('lf0')
            read(val, *) lf0
          case ('cwmci')
            read(val, *) cwmci
          case ('mi0iv')
            read(val, *) mi0iv
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

end program kernel_benchmark_siadjst

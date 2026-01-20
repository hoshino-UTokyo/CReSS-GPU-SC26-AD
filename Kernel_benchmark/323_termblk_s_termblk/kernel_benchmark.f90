!***********************************************************************
! Kernel Benchmark: termblk
!***********************************************************************
! Source: Src/termblk.f90
! Subroutine: s_termblk
! Profile ID: 323
! Description: Calculate terminal velocities for cloud water, rain, ice,
!              snow, graupel, and hail based on microphysics options
!***********************************************************************
program kernel_benchmark_termblk
  use omp_lib
  implicit none

  ! === Grid dimensions ===
  integer :: cphopt, haiopt, ni, nj, nk

  ! === Scalar parameters ===
  real :: thresq, r0
  real :: buc3, bur3, bui3, bus3, bug3, buh3
  real :: cucq, curq, cusq, cugq, cuhq
  real :: cucn, curn, cusn, cugn, cuhn
  real :: ccrw6, ccri6
  real :: cdiaqc, cdiaqr, cdiaqs, cdiaqg, cdiaqh

  ! === Compile-time constants (from temparam.f90) ===
  integer, parameter :: flqcqi_opt = 1    ! Terminal velocity option
  real, parameter :: ucqcst = 0.1e0       ! Constant cloud water velocity
  real, parameter :: ucncst = 0.1e0       ! Constant cloud water conc velocity
  real, parameter :: uiqcst = 0.0001e0    ! Constant cloud ice velocity
  real, parameter :: uincst = 0.0001e0    ! Constant cloud ice conc velocity

  ! === Module constants (from m_comphy) ===
  real, parameter :: guc = 1.e0
  real, parameter :: gur = 0.5e0
  real, parameter :: gui = 0.33e0
  real, parameter :: gus = 0.5e0
  real, parameter :: gug = 0.5e0
  real, parameter :: guh = 0.5e0
  real, parameter :: auc = 2.98e7
  real, parameter :: aui = 770.e0

  ! === Input arrays ===
  real, allocatable :: rbv(:,:,:)
  real, allocatable :: qc(:,:,:), qr(:,:,:), qi(:,:,:)
  real, allocatable :: qs(:,:,:), qg(:,:,:), qh(:,:,:)
  real, allocatable :: ncc(:,:,:), ncr(:,:,:), nci(:,:,:)
  real, allocatable :: ncs(:,:,:), ncg(:,:,:), nch(:,:,:)

  ! === Output arrays ===
  real, allocatable :: ucq(:,:,:), urq(:,:,:), uiq(:,:,:)
  real, allocatable :: usq(:,:,:), ugq(:,:,:), uhq(:,:,:)
  real, allocatable :: ucn(:,:,:), urn(:,:,:), uin(:,:,:)
  real, allocatable :: usn(:,:,:), ugn(:,:,:), uhn(:,:,:)

  ! === Reference arrays for validation ===
  real, allocatable :: ucq_ref(:,:,:), urq_ref(:,:,:), uiq_ref(:,:,:)
  real, allocatable :: usq_ref(:,:,:), ugq_ref(:,:,:), uhq_ref(:,:,:)
  real, allocatable :: ucn_ref(:,:,:), urn_ref(:,:,:), uin_ref(:,:,:)
  real, allocatable :: usn_ref(:,:,:), ugn_ref(:,:,:), uhn_ref(:,:,:)

  ! === Benchmark control ===
  integer :: num_iterations, warmup_iterations, iter
  real :: tolerance
  character(len=512) :: data_dir

  ! === Timing ===
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max

  ! === Validation ===
  logical :: passed
  real :: max_err, max_err_total
  integer :: err_count, err_count_total

  ! === Main program flow ===

  ! 1. Read configuration
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! 2. Read parameters
  call read_parameters(trim(data_dir)//'/params.txt')

  ! 3. Allocate arrays
  allocate(rbv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qc(0:ni+1, 0:nj+1, 1:nk))
  allocate(qr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(qg(0:ni+1, 0:nj+1, 1:nk))
  allocate(qh(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncc(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncr(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncs(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncg(0:ni+1, 0:nj+1, 1:nk))
  allocate(nch(0:ni+1, 0:nj+1, 1:nk))

  allocate(ucq(0:ni+1, 0:nj+1, 1:nk))
  allocate(urq(0:ni+1, 0:nj+1, 1:nk))
  allocate(uiq(0:ni+1, 0:nj+1, 1:nk))
  allocate(usq(0:ni+1, 0:nj+1, 1:nk))
  allocate(ugq(0:ni+1, 0:nj+1, 1:nk))
  allocate(uhq(0:ni+1, 0:nj+1, 1:nk))
  allocate(ucn(0:ni+1, 0:nj+1, 1:nk))
  allocate(urn(0:ni+1, 0:nj+1, 1:nk))
  allocate(uin(0:ni+1, 0:nj+1, 1:nk))
  allocate(usn(0:ni+1, 0:nj+1, 1:nk))
  allocate(ugn(0:ni+1, 0:nj+1, 1:nk))
  allocate(uhn(0:ni+1, 0:nj+1, 1:nk))

  allocate(ucq_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(urq_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(uiq_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(usq_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ugq_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(uhq_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ucn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(urn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(uin_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(usn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ugn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(uhn_ref(0:ni+1, 0:nj+1, 1:nk))

  ! 4. Read input data
  call read_array_3d(trim(data_dir)//'/rbv.bin', rbv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qc.bin', qc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qr.bin', qr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qi.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qg.bin', qg, 0, ni+1, 0, nj+1, 1, nk)
  if (haiopt /= 0) then
    call read_array_3d(trim(data_dir)//'/qh.bin', qh, 0, ni+1, 0, nj+1, 1, nk)
  else
    qh = 0.0
  end if
  call read_array_3d(trim(data_dir)//'/ncc.bin', ncc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncr.bin', ncr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/nci.bin', nci, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncs.bin', ncs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncg.bin', ncg, 0, ni+1, 0, nj+1, 1, nk)
  if (haiopt /= 0) then
    call read_array_3d(trim(data_dir)//'/nch.bin', nch, 0, ni+1, 0, nj+1, 1, nk)
  else
    nch = 0.0
  end if

  ! 5. Read reference output
  call read_array_3d(trim(data_dir)//'/ucq_ref.bin', ucq_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/urq_ref.bin', urq_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/uiq_ref.bin', uiq_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/usq_ref.bin', usq_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ugq_ref.bin', ugq_ref, 0, ni+1, 0, nj+1, 1, nk)
  if (haiopt /= 0) then
    call read_array_3d(trim(data_dir)//'/uhq_ref.bin', uhq_ref, 0, ni+1, 0, nj+1, 1, nk)
  else
    uhq_ref = 0.0
  end if
  call read_array_3d(trim(data_dir)//'/ucn_ref.bin', ucn_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/urn_ref.bin', urn_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/uin_ref.bin', uin_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/usn_ref.bin', usn_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ugn_ref.bin', ugn_ref, 0, ni+1, 0, nj+1, 1, nk)
  if (haiopt /= 0) then
    call read_array_3d(trim(data_dir)//'/uhn_ref.bin', uhn_ref, 0, ni+1, 0, nj+1, 1, nk)
  else
    uhn_ref = 0.0
  end if

  ! Print header
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: termblk'
  write(*,'(A)') '=================================================='
  write(*,'(A,I8,A,I8,A,I8)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I8)') ' cphopt: ', cphopt
  write(*,'(A,I8)') ' haiopt: ', haiopt
  write(*,'(A,I8)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  ! 6. Warmup iterations
  do iter = 1, warmup_iterations
    call kernel_termblk(cphopt, haiopt, thresq, ni, nj, nk, r0, &
         buc3, bur3, bui3, bus3, bug3, buh3, &
         cucq, curq, cusq, cugq, cuhq, &
         cucn, curn, cusn, cugn, cuhn, &
         ccrw6, ccri6, cdiaqc, cdiaqr, cdiaqs, cdiaqg, cdiaqh, &
         rbv, qc, qr, qi, qs, qg, qh, ncc, ncr, nci, ncs, ncg, nch, &
         ucq, urq, uiq, usq, ugq, uhq, ucn, urn, uin, usn, ugn, uhn)
  end do

  ! 7. Benchmark iterations (timed)
  t_total = 0.0d0
  t_min = huge(1.0d0)
  t_max = 0.0d0

  do iter = 1, num_iterations
    t_start = omp_get_wtime()
    call kernel_termblk(cphopt, haiopt, thresq, ni, nj, nk, r0, &
         buc3, bur3, bui3, bus3, bug3, buh3, &
         cucq, curq, cusq, cugq, cuhq, &
         cucn, curn, cusn, cugn, cuhn, &
         ccrw6, ccri6, cdiaqc, cdiaqr, cdiaqs, cdiaqg, cdiaqh, &
         rbv, qc, qr, qi, qs, qg, qh, ncc, ncr, nci, ncs, ncg, nch, &
         ucq, urq, uiq, usq, ugq, uhq, ucn, urn, uin, usn, ugn, uhn)
    t_end = omp_get_wtime()

    t_total = t_total + (t_end - t_start)
    if (t_end - t_start < t_min) t_min = t_end - t_start
    if (t_end - t_start > t_max) t_max = t_end - t_start
  end do
  t_avg = t_total / dble(num_iterations)

  ! 8. Validate output (only interior points: 1:ni-1, 1:nj-1, 1:nk-1)
  max_err_total = 0.0
  err_count_total = 0

  call validate_output_3d(ucq, ucq_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(urq, urq_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(uiq, uiq_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(usq, usq_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(ugq, ugq_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(uhq, uhq_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(ucn, ucn_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(urn, urn_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(uin, uin_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(usn, usn_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(ugn, ugn_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(uhn, uhn_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, nk-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  passed = (err_count_total == 0)

  ! 9. Report results
  write(*,'(A,F12.3,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:     ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:     ', t_max * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_err_total
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Total error count:  ', err_count_total
  if (passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

contains

  !-----------------------------------------------------------------
  ! The kernel (extracted from original source)
  !-----------------------------------------------------------------
  subroutine kernel_termblk(cphopt, haiopt, thresq, ni, nj, nk, r0, &
       buc3, bur3, bui3, bus3, bug3, buh3, &
       cucq, curq, cusq, cugq, cuhq, &
       cucn, curn, cusn, cugn, cuhn, &
       ccrw6, ccri6, cdiaqc, cdiaqr, cdiaqs, cdiaqg, cdiaqh, &
       rbv, qc, qr, qi, qs, qg, qh, ncc, ncr, nci, ncs, ncg, nch, &
       ucq, urq, uiq, usq, ugq, uhq, ucn, urn, uin, usn, ugn, uhn)
    implicit none

    integer, intent(in) :: cphopt, haiopt, ni, nj, nk
    real, intent(in) :: thresq, r0
    real, intent(in) :: buc3, bur3, bui3, bus3, bug3, buh3
    real, intent(in) :: cucq, curq, cusq, cugq, cuhq
    real, intent(in) :: cucn, curn, cusn, cugn, cuhn
    real, intent(in) :: ccrw6, ccri6, cdiaqc, cdiaqr, cdiaqs, cdiaqg, cdiaqh
    real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qh(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncc(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: nci(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: nch(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: ucq(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: urq(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: uiq(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: usq(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: ugq(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: uhq(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: ucn(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: urn(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: uin(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: usn(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: ugn(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: uhn(0:ni+1,0:nj+1,1:nk)

    integer :: i, j, k
    real :: lnr0r

    !$omp parallel default(shared) private(k)

    ! Set the terminal velocity of cloud water and cloud ice
    ! Using flqcqi_opt = 1 (constant velocities)
    if (flqcqi_opt == 1) then
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qc(i,j,k) > thresq) then
              ucq(i,j,k) = ucqcst
              ucn(i,j,k) = ucncst
            else
              ucq(i,j,k) = 0.e0
              ucn(i,j,k) = 0.e0
            end if

            if (qi(i,j,k) > thresq) then
              uiq(i,j,k) = uiqcst
              uin(i,j,k) = uincst
            else
              uiq(i,j,k) = 0.e0
              uin(i,j,k) = 0.e0
            end if
          end do
        end do
        !$omp end do
      end do

    else if (flqcqi_opt == 2) then
      if (abs(cphopt) <= 3) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j,lnr0r)
          do j = 1, nj-1
            do i = 1, ni-1
              lnr0r = log(r0 * rbv(i,j,k))

              if (qc(i,j,k) > thresq) then
                ucq(i,j,k) = auc * exp(guc * lnr0r) &
                     * exp(buc3 * log(ccrw6 * qc(i,j,k) / ncc(i,j,k)))
                ucn(i,j,k) = ucq(i,j,k)
              else
                ucq(i,j,k) = 0.e0
                ucn(i,j,k) = 0.e0
              end if

              if (qi(i,j,k) > thresq) then
                uiq(i,j,k) = aui * exp(gui * lnr0r) &
                     * exp(bui3 * log(ccri6 * qi(i,j,k) / nci(i,j,k)))
                uin(i,j,k) = uiq(i,j,k)
              else
                uiq(i,j,k) = 0.e0
                uin(i,j,k) = 0.e0
              end if
            end do
          end do
          !$omp end do
        end do
      else
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j,lnr0r)
          do j = 1, nj-1
            do i = 1, ni-1
              lnr0r = log(r0 * rbv(i,j,k))

              if (qc(i,j,k) > thresq) then
                ucq(i,j,k) = cucq * exp(guc * lnr0r) &
                     * exp(buc3 * log(cdiaqc * qc(i,j,k) / ncc(i,j,k)))
                ucn(i,j,k) = cucn * ucq(i,j,k)
              else
                ucq(i,j,k) = 0.e0
                ucn(i,j,k) = 0.e0
              end if

              if (qi(i,j,k) > thresq) then
                uiq(i,j,k) = aui * exp(gui * lnr0r) &
                     * exp(bui3 * log(ccri6 * qi(i,j,k) / nci(i,j,k)))
                uin(i,j,k) = uiq(i,j,k)
              else
                uiq(i,j,k) = 0.e0
                uin(i,j,k) = 0.e0
              end if
            end do
          end do
          !$omp end do
        end do
      end if
    end if

    ! Calculate terminal velocity of rain, snow, graupel (and hail if haiopt /= 0)
    if (haiopt == 0) then
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j,lnr0r)
        do j = 1, nj-1
          do i = 1, ni-1
            lnr0r = log(r0 * rbv(i,j,k))

            if (qr(i,j,k) > thresq) then
              urq(i,j,k) = curq * exp(gur * lnr0r) &
                   * exp(bur3 * log(cdiaqr * qr(i,j,k) / ncr(i,j,k)))
              urn(i,j,k) = curn * urq(i,j,k)
            else
              urq(i,j,k) = 0.e0
              urn(i,j,k) = 0.e0
            end if

            if (qs(i,j,k) > thresq) then
              usq(i,j,k) = cusq * exp(gus * lnr0r) &
                   * exp(bus3 * log(cdiaqs * qs(i,j,k) / ncs(i,j,k)))
              usn(i,j,k) = cusn * usq(i,j,k)
            else
              usq(i,j,k) = 0.e0
              usn(i,j,k) = 0.e0
            end if

            if (qg(i,j,k) > thresq) then
              ugq(i,j,k) = cugq * exp(gug * lnr0r) &
                   * exp(bug3 * log(cdiaqg * qg(i,j,k) / ncg(i,j,k)))
              ugn(i,j,k) = cugn * ugq(i,j,k)
            else
              ugq(i,j,k) = 0.e0
              ugn(i,j,k) = 0.e0
            end if

            ! No hail when haiopt == 0
            uhq(i,j,k) = 0.e0
            uhn(i,j,k) = 0.e0
          end do
        end do
        !$omp end do
      end do
    else
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j,lnr0r)
        do j = 1, nj-1
          do i = 1, ni-1
            lnr0r = log(r0 * rbv(i,j,k))

            if (qr(i,j,k) > thresq) then
              urq(i,j,k) = curq * exp(gur * lnr0r) &
                   * exp(bur3 * log(cdiaqr * qr(i,j,k) / ncr(i,j,k)))
              urn(i,j,k) = curn * urq(i,j,k)
            else
              urq(i,j,k) = 0.e0
              urn(i,j,k) = 0.e0
            end if

            if (qs(i,j,k) > thresq) then
              usq(i,j,k) = cusq * exp(gus * lnr0r) &
                   * exp(bus3 * log(cdiaqs * qs(i,j,k) / ncs(i,j,k)))
              usn(i,j,k) = cusn * usq(i,j,k)
            else
              usq(i,j,k) = 0.e0
              usn(i,j,k) = 0.e0
            end if

            if (qg(i,j,k) > thresq) then
              ugq(i,j,k) = cugq * exp(gug * lnr0r) &
                   * exp(bug3 * log(cdiaqg * qg(i,j,k) / ncg(i,j,k)))
              ugn(i,j,k) = cugn * ugq(i,j,k)
            else
              ugq(i,j,k) = 0.e0
              ugn(i,j,k) = 0.e0
            end if

            if (qh(i,j,k) > thresq) then
              uhq(i,j,k) = cuhq * exp(guh * lnr0r) &
                   * exp(buh3 * log(cdiaqh * qh(i,j,k) / nch(i,j,k)))
              uhn(i,j,k) = cuhn * uhq(i,j,k)
            else
              uhq(i,j,k) = 0.e0
              uhn(i,j,k) = 0.e0
            end if
          end do
        end do
        !$omp end do
      end do
    end if

    !$omp end parallel

  end subroutine kernel_termblk

  !-----------------------------------------------------------------
  ! I/O utilities
  !-----------------------------------------------------------------
  subroutine read_config(data_dir, num_iterations, warmup_iterations, tolerance)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iterations, warmup_iterations
    real, intent(out) :: tolerance

    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iterations = 10
    warmup_iterations = 2
    tolerance = 1.0e-5

    inquire(file='benchmark.conf', exist=exists)
    if (exists) then
      open(unit=10, file='benchmark.conf', status='old', iostat=ios)
      if (ios == 0) then
        read(10, '(A)', iostat=ios) data_dir
        read(10, *, iostat=ios) num_iterations
        read(10, *, iostat=ios) warmup_iterations
        read(10, *, iostat=ios) tolerance
        close(10)
      end if
    end if

    data_dir = trim(adjustl(data_dir))
  end subroutine read_config

  subroutine read_parameters(filename)
    implicit none
    character(len=*), intent(in) :: filename

    character(len=256) :: line
    character(len=64) :: varname
    integer :: ios, eq_pos

    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error opening parameter file: ', trim(filename)
      stop 1
    end if

    do while (.true.)
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        varname = trim(adjustl(line(1:eq_pos-1)))
        select case (trim(varname))
        case ('cphopt'); read(line(eq_pos+1:), *) cphopt
        case ('haiopt'); read(line(eq_pos+1:), *) haiopt
        case ('ni'); read(line(eq_pos+1:), *) ni
        case ('nj'); read(line(eq_pos+1:), *) nj
        case ('nk'); read(line(eq_pos+1:), *) nk
        case ('thresq'); read(line(eq_pos+1:), *) thresq
        case ('r0'); read(line(eq_pos+1:), *) r0
        case ('buc3'); read(line(eq_pos+1:), *) buc3
        case ('bur3'); read(line(eq_pos+1:), *) bur3
        case ('bui3'); read(line(eq_pos+1:), *) bui3
        case ('bus3'); read(line(eq_pos+1:), *) bus3
        case ('bug3'); read(line(eq_pos+1:), *) bug3
        case ('buh3'); read(line(eq_pos+1:), *) buh3
        case ('ccri6'); read(line(eq_pos+1:), *) ccri6
        case ('ccrw6'); read(line(eq_pos+1:), *) ccrw6
        case ('cdiaqc'); read(line(eq_pos+1:), *) cdiaqc
        case ('cdiaqg'); read(line(eq_pos+1:), *) cdiaqg
        case ('cdiaqh'); read(line(eq_pos+1:), *) cdiaqh
        case ('cdiaqr'); read(line(eq_pos+1:), *) cdiaqr
        case ('cdiaqs'); read(line(eq_pos+1:), *) cdiaqs
        case ('cucn'); read(line(eq_pos+1:), *) cucn
        case ('cucq'); read(line(eq_pos+1:), *) cucq
        case ('cugn'); read(line(eq_pos+1:), *) cugn
        case ('cugq'); read(line(eq_pos+1:), *) cugq
        case ('cuhn'); read(line(eq_pos+1:), *) cuhn
        case ('cuhq'); read(line(eq_pos+1:), *) cuhq
        case ('curn'); read(line(eq_pos+1:), *) curn
        case ('curq'); read(line(eq_pos+1:), *) curq
        case ('cusn'); read(line(eq_pos+1:), *) cusn
        case ('cusq'); read(line(eq_pos+1:), *) cusq
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: ios

    open(unit=20, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(20, iostat=ios) arr
    if (ios /= 0) then
      write(*,*) 'Error reading file: ', trim(filename)
      stop 1
    end if
    close(20)
  end subroutine read_array_3d

  subroutine validate_output_3d(output, reference, tol, max_err, err_count, is, ie, js, je, ks, ke)
    implicit none
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(in) :: output(0:,0:,1:)
    real, intent(in) :: reference(0:,0:,1:)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err, abs_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = ks, ke
      do j = js, je
        do i = is, ie
          abs_err = abs(output(i,j,k) - reference(i,j,k))
          if (abs(reference(i,j,k)) > 1.0e-30) then
            rel_err = abs_err / abs(reference(i,j,k))
          else
            rel_err = abs_err
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do
  end subroutine validate_output_3d

end program kernel_benchmark_termblk

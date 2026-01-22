program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup, num_iterations
  real(8) :: tolerance

  ! Dimensions and options
  integer :: cphopt, ni, nj, nk

  ! Physical constants from params.txt
  real :: dtb, thresq, t0, r0, rhow, rhos
  real :: aur, aus, aug, bur, bus, bug
  real :: eir, eis, eig, ers, erg, esr, esg
  real :: cc, oned9, gf3bur, gf3bus, gf3bug

  ! Additional physical constant (hardcoded from m_comphy)
  real, parameter :: tlow = 233.16

  ! Derived values (computed before parallel region)
  real :: rwdv9, bur2, bus2, bug2, esgiv
  real :: cclcr, cclcs, cclcg, cclri, cclrs, cclrg
  real :: cclis, cclig, cclsr, cclsg, cclsg3
  real :: cclrsn, cclsgn

  ! Input arrays
  real, allocatable :: t(:,:,:), rbr(:,:,:)
  real, allocatable :: qc(:,:,:), qr(:,:,:), qi(:,:,:), qs(:,:,:), qg(:,:,:)
  real, allocatable :: ncr(:,:,:), ncs(:,:,:), ncg(:,:,:)
  real, allocatable :: urq(:,:,:), usq(:,:,:), ugq(:,:,:)
  real, allocatable :: urn(:,:,:), usn(:,:,:), ugn(:,:,:)
  real, allocatable :: mu(:,:,:), mi(:,:,:)
  real, allocatable :: diaqc(:,:,:), diaqr(:,:,:), diaqs(:,:,:), diaqg(:,:,:)

  ! Output arrays
  real, allocatable :: clcr(:,:,:), clcs(:,:,:), clcg(:,:,:)
  real, allocatable :: clri(:,:,:), clrs(:,:,:), clrg(:,:,:)
  real, allocatable :: clir(:,:,:), clis(:,:,:), clig(:,:,:)
  real, allocatable :: clsr(:,:,:), clsg(:,:,:)
  real, allocatable :: clrin(:,:,:), clrsn(:,:,:), clsrn(:,:,:), clsgn(:,:,:)
  real, allocatable :: ecs(:,:,:)

  ! Reference arrays
  real, allocatable :: clcr_ref(:,:,:), clcs_ref(:,:,:), clcg_ref(:,:,:)
  real, allocatable :: clri_ref(:,:,:), clrs_ref(:,:,:), clrg_ref(:,:,:)
  real, allocatable :: clir_ref(:,:,:), clis_ref(:,:,:), clig_ref(:,:,:)
  real, allocatable :: clsr_ref(:,:,:), clsg_ref(:,:,:)
  real, allocatable :: clrin_ref(:,:,:), clrsn_ref(:,:,:), clsrn_ref(:,:,:), clsgn_ref(:,:,:)
  real, allocatable :: ecs_ref(:,:,:)

  ! Timing variables
  real(8) :: start_time, end_time, elapsed_time
  real(8) :: total_time, min_time, max_time, avg_time
  real(8), allocatable :: times(:)

  ! Validation variables
  integer :: error_count, total_errors
  real(8) :: max_rel_error

  ! Loop variables
  integer :: iter

  ! Config file
  integer :: unit_conf

  ! Read configuration
  unit_conf = 10
  open(unit_conf, file='benchmark.conf', status='old', action='read')
  read(unit_conf, '(A)') data_dir
  read(unit_conf, *) num_warmup
  read(unit_conf, *) num_iterations
  read(unit_conf, *) tolerance
  close(unit_conf)

  data_dir = trim(adjustl(data_dir))
  print '(A)', '=== Collect Kernel Benchmark (GPU) ==='
  print '(A,A)', 'Data directory: ', trim(data_dir)
  print '(A,I0)', 'Warmup iterations: ', num_warmup
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,ES10.2)', 'Tolerance: ', tolerance

  ! Read parameters
  call read_params(trim(data_dir) // '/params.txt')

  print '(A)', ''
  print '(A)', '--- Parameters ---'
  print '(A,I0)', 'cphopt = ', cphopt
  print '(A,I0,A,I0,A,I0)', 'ni x nj x nk = ', ni, ' x ', nj, ' x ', nk

  ! Compute derived values
  rwdv9 = oned9 * rhow
  bur2 = 2.0 + bur
  bus2 = 2.0 + bus
  bug2 = 2.0 + bug
  esgiv = 1.0 / esg
  cclcr = 0.25 * aur * gf3bur * cc * dtb
  cclcs = 0.25 * aus * gf3bus * cc * dtb
  cclcg = 0.25 * aug * gf3bug * cc * dtb
  cclri = 0.25 * aur * gf3bur * eir * cc * dtb
  cclrs = ers * cc * cc * rhow * dtb
  cclrg = erg * cc * cc * rhow * dtb
  cclis = 0.25 * aus * gf3bus * eis * cc * dtb
  cclig = 0.25 * aug * gf3bug * eig * cc * dtb
  cclsr = esr * cc * cc * rhos * dtb
  cclsg = esg * cc * cc * rhos * dtb
  cclsg3 = cc * cc * rhos * dtb
  cclrsn = 0.5 * ers * cc * dtb
  cclsgn = 0.5 * esg * cc * dtb

  ! Allocate input arrays
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qc(0:ni+1, 0:nj+1, 1:nk))
  allocate(qr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(qg(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncs(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncg(0:ni+1, 0:nj+1, 1:nk))
  allocate(urq(0:ni+1, 0:nj+1, 1:nk))
  allocate(usq(0:ni+1, 0:nj+1, 1:nk))
  allocate(ugq(0:ni+1, 0:nj+1, 1:nk))
  allocate(urn(0:ni+1, 0:nj+1, 1:nk))
  allocate(usn(0:ni+1, 0:nj+1, 1:nk))
  allocate(ugn(0:ni+1, 0:nj+1, 1:nk))
  allocate(mu(0:ni+1, 0:nj+1, 1:nk))
  allocate(mi(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqc(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqr(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqs(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqg(0:ni+1, 0:nj+1, 1:nk))

  ! Allocate output arrays
  allocate(clcr(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clri(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clir(0:ni+1, 0:nj+1, 1:nk))
  allocate(clis(0:ni+1, 0:nj+1, 1:nk))
  allocate(clig(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsr(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrin(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrsn(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsrn(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsgn(0:ni+1, 0:nj+1, 1:nk))
  allocate(ecs(0:ni+1, 0:nj+1, 1:nk))

  ! Allocate reference arrays
  allocate(clcr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcs_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clri_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrs_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clir_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clis_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clig_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrin_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrsn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsrn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsgn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ecs_ref(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  ! Read input arrays
  print '(A)', ''
  print '(A)', '--- Loading input data ---'
  call read_array_3d(trim(data_dir) // '/t.bin', t, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/rbr.bin', rbr, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qc.bin', qc, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qr.bin', qr, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qi.bin', qi, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qs.bin', qs, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qg.bin', qg, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ncr.bin', ncr, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ncs.bin', ncs, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ncg.bin', ncg, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/urq.bin', urq, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/usq.bin', usq, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ugq.bin', ugq, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/urn.bin', urn, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/usn.bin', usn, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ugn.bin', ugn, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/mu.bin', mu, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/mi.bin', mi, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/diaqc.bin', diaqc, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/diaqr.bin', diaqr, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/diaqs.bin', diaqs, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/diaqg.bin', diaqg, ni, nj, nk)

  ! Read reference arrays
  call read_array_3d(trim(data_dir) // '/clcr_ref.bin', clcr_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clcs_ref.bin', clcs_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clcg_ref.bin', clcg_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clri_ref.bin', clri_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clrs_ref.bin', clrs_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clrg_ref.bin', clrg_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clir_ref.bin', clir_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clis_ref.bin', clis_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clig_ref.bin', clig_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clsr_ref.bin', clsr_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clsg_ref.bin', clsg_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clrin_ref.bin', clrin_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clrsn_ref.bin', clrsn_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clsrn_ref.bin', clsrn_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clsgn_ref.bin', clsgn_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ecs_ref.bin', ecs_ref, ni, nj, nk)
  print '(A)', 'Data loaded successfully'

  ! Warmup iterations
  print '(A)', ''
  print '(A)', '--- Warmup ---'
  do iter = 1, num_warmup
    call run_kernel()
    !$acc wait
  end do
  print '(A,I0,A)', 'Completed ', num_warmup, ' warmup iterations'

  ! Benchmark iterations
  print '(A)', ''
  print '(A)', '--- Benchmark ---'
  total_time = 0.0d0
  min_time = huge(1.0d0)
  max_time = 0.0d0

  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()
    call run_kernel()
    !$acc wait
    end_time = omp_get_wtime()

    elapsed_time = end_time - start_time
    times(iter) = elapsed_time
    total_time = total_time + elapsed_time
    min_time = min(min_time, elapsed_time)
    max_time = max(max_time, elapsed_time)
  end do

  avg_time = total_time / num_iterations

  print '(A)', ''
  print '(A)', '=== Timing Results ==='
  print '(A,ES12.4,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,ES12.4,A)', 'Min time:     ', min_time * 1000.0d0, ' ms'
  print '(A,ES12.4,A)', 'Max time:     ', max_time * 1000.0d0, ' ms'
  print '(A,ES12.4,A)', 'Total time:   ', total_time * 1000.0d0, ' ms'

  ! Validation
  print '(A)', ''
  print '(A)', '=== Validation ==='
  total_errors = 0

  call validate_array('clcr', clcr, clcr_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clcs', clcs, clcs_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clcg', clcg, clcg_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clri', clri, clri_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clrs', clrs, clrs_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clrg', clrg, clrg_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clir', clir, clir_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clis', clis, clis_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clig', clig, clig_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clsr', clsr, clsr_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clsg', clsg, clsg_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clrin', clrin, clrin_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clrsn', clrsn, clrsn_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clsrn', clsrn, clsrn_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('clsgn', clsgn, clsgn_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('ecs', ecs, ecs_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  print '(A)', ''
  if (total_errors == 0) then
    print '(A)', '*** VALIDATION PASSED ***'
  else
    print '(A,I0,A)', '*** VALIDATION FAILED: ', total_errors, ' errors ***'
  end if

contains

  subroutine run_kernel()
    integer :: i, j, k
    real :: cstk, r0rsq, qr2b, qs2b, qg2b
    real :: diaqr2, diaqs2, diaqg2, diaqr3, diaqs3
    real :: sink, a, b, c

    ! In the case nk = 1
    if (nk == 1) then
      ! Case abs(cphopt) >= 3
      if (abs(cphopt) >= 3) then
        !$acc kernels
        !$acc loop independent collapse(2) &
        !$acc& private(i,j,cstk,r0rsq,qr2b,qs2b,qg2b) &
        !$acc& private(diaqr2,diaqs2,diaqg2,diaqr3,diaqs3,sink,a,b,c)
        do j = 1, nj-1
          do i = 1, ni-1
            ! Set common variables
            r0rsq = sqrt(r0 * rbr(i,j,1))

            if (qr(i,j,1) > thresq) then
              diaqr2 = diaqr(i,j,1) * diaqr(i,j,1)
              diaqr3 = diaqr(i,j,1) * diaqr2
              qr2b = r0rsq * exp(bur2 * log(diaqr(i,j,1)))
            else
              diaqr2 = 0.0
              diaqr3 = 0.0
              qr2b = 0.0
            end if

            if (qs(i,j,1) > thresq) then
              diaqs2 = diaqs(i,j,1) * diaqs(i,j,1)
              diaqs3 = diaqs(i,j,1) * diaqs2
              qs2b = r0rsq * exp(bus2 * log(diaqs(i,j,1)))
            else
              diaqs2 = 0.0
              diaqs3 = 0.0
              qs2b = 0.0
            end if

            if (qg(i,j,1) > thresq) then
              diaqg2 = diaqg(i,j,1) * diaqg(i,j,1)
              qg2b = r0rsq * exp(bug2 * log(diaqg(i,j,1)))
            else
              diaqg2 = 0.0
              qg2b = 0.0
            end if

            ! Cloud water and rain/snow/graupel collection
            if (t(i,j,1) > tlow) then
              if (qc(i,j,1) > thresq) then
                cstk = rwdv9 * diaqc(i,j,1) * diaqc(i,j,1)

                if (qr(i,j,1) > thresq) then
                  a = cstk * urq(i,j,1) / (diaqr(i,j,1) * mu(i,j,1))
                  b = a + 0.5
                  a = a * a / (b * b)
                  clcr(i,j,1) = cclcr * a * qr2b * ncr(i,j,1) * qc(i,j,1)
                else
                  clcr(i,j,1) = 0.0
                end if

                if (qs(i,j,1) > thresq) then
                  a = cstk * usq(i,j,1) / (diaqs(i,j,1) * mu(i,j,1))
                  b = a + 0.5
                  ecs(i,j,1) = a * a / (b * b)
                  clcs(i,j,1) = cclcs * qs2b * ecs(i,j,1) * ncs(i,j,1) * qc(i,j,1)
                else
                  clcs(i,j,1) = 0.0
                  ecs(i,j,1) = 0.0
                end if

                if (qg(i,j,1) > thresq) then
                  a = cstk * ugq(i,j,1) / (diaqg(i,j,1) * mu(i,j,1))
                  b = a + 0.5
                  a = a * a / (b * b)
                  clcg(i,j,1) = cclcg * a * qg2b * ncg(i,j,1) * qc(i,j,1)
                else
                  clcg(i,j,1) = 0.0
                end if

                sink = clcr(i,j,1) + clcs(i,j,1) + clcg(i,j,1)
                if (qc(i,j,1) < sink) then
                  a = qc(i,j,1) / sink
                  clcr(i,j,1) = clcr(i,j,1) * a
                  clcs(i,j,1) = clcs(i,j,1) * a
                  clcg(i,j,1) = clcg(i,j,1) * a
                end if
              else
                clcr(i,j,1) = 0.0
                clcs(i,j,1) = 0.0
                clcg(i,j,1) = 0.0
                ecs(i,j,1) = 0.0
              end if

              ! Rain water and cloud ice/snow/graupel
              if (qr(i,j,1) > thresq) then
                if (qi(i,j,1) > thresq) then
                  if (t(i,j,1) < t0) then
                    clri(i,j,1) = cclri * qr2b * ncr(i,j,1) * qi(i,j,1)
                    clir(i,j,1) = clri(i,j,1)
                    clrin(i,j,1) = clri(i,j,1) / mi(i,j,1)
                  else
                    clri(i,j,1) = 0.0
                    clir(i,j,1) = 0.0
                    clrin(i,j,1) = 0.0
                  end if
                else
                  clri(i,j,1) = 0.0
                  clir(i,j,1) = 0.0
                  clrin(i,j,1) = 0.0
                end if

                if (qs(i,j,1) > thresq) then
                  a = 2.0 * diaqr(i,j,1) * diaqs(i,j,1)
                  b = ncr(i,j,1) * ncs(i,j,1) * rbr(i,j,1)
                  c = urq(i,j,1) - usq(i,j,1)
                  c = sqrt(c * c + 0.04 * urq(i,j,1) * usq(i,j,1))
                  clrs(i,j,1) = b * c * cclrs * (5.0 * diaqr2 + a + 0.5 * diaqs2) * diaqr3
                  clsr(i,j,1) = b * c * cclsr * (5.0 * diaqs2 + a + 0.5 * diaqr2) * diaqs3
                  c = urn(i,j,1) - usn(i,j,1)
                  clrsn(i,j,1) = b * cclrsn * (diaqr2 + 0.5 * a + diaqs2) &
                               * sqrt(c * c + 0.04 * urn(i,j,1) * usn(i,j,1))
                  clsrn(i,j,1) = clrsn(i,j,1)
                else
                  clrs(i,j,1) = 0.0
                  clsr(i,j,1) = 0.0
                  clrsn(i,j,1) = 0.0
                  clsrn(i,j,1) = 0.0
                end if

                if (qg(i,j,1) > thresq) then
                  a = 5.0 * diaqr2 + 2.0 * diaqr(i,j,1) * diaqg(i,j,1) + 0.5 * diaqg2
                  b = urq(i,j,1) - ugq(i,j,1)
                  clrg(i,j,1) = a * cclrg * diaqr3 * ncr(i,j,1) * ncg(i,j,1) &
                              * rbr(i,j,1) * sqrt(b * b + 0.04 * urq(i,j,1) * ugq(i,j,1))
                else
                  clrg(i,j,1) = 0.0
                end if

                sink = clri(i,j,1) + clrs(i,j,1) + clrg(i,j,1)
                if (qr(i,j,1) < sink) then
                  a = qr(i,j,1) / sink
                  clri(i,j,1) = clri(i,j,1) * a
                  clrs(i,j,1) = clrs(i,j,1) * a
                  clrg(i,j,1) = clrg(i,j,1) * a
                end if

                sink = clrin(i,j,1) + clrsn(i,j,1)
                if (ncr(i,j,1) < sink) then
                  a = ncr(i,j,1) / sink
                  clrin(i,j,1) = clrin(i,j,1) * a
                  clrsn(i,j,1) = clrsn(i,j,1) * a
                end if
              else
                clri(i,j,1) = 0.0
                clrs(i,j,1) = 0.0
                clrg(i,j,1) = 0.0
                clir(i,j,1) = 0.0
                clsr(i,j,1) = 0.0
                clrin(i,j,1) = 0.0
                clrsn(i,j,1) = 0.0
                clsrn(i,j,1) = 0.0
              end if
            else
              ! Temperature below tlow
              clcr(i,j,1) = 0.0
              clcs(i,j,1) = 0.0
              clcg(i,j,1) = 0.0
              clri(i,j,1) = 0.0
              clrs(i,j,1) = 0.0
              clrg(i,j,1) = 0.0
              clir(i,j,1) = 0.0
              clsr(i,j,1) = 0.0
              clrin(i,j,1) = 0.0
              clrsn(i,j,1) = 0.0
              clsrn(i,j,1) = 0.0
              ecs(i,j,1) = 0.0
            end if

            ! Cloud ice and snow/graupel
            if (qi(i,j,1) > thresq) then
              if (t(i,j,1) < t0) then
                if (qs(i,j,1) > thresq) then
                  clis(i,j,1) = cclis * qs2b * ncs(i,j,1) * qi(i,j,1)
                else
                  clis(i,j,1) = 0.0
                end if
                if (qg(i,j,1) > thresq) then
                  clig(i,j,1) = cclig * qg2b * ncg(i,j,1) * qi(i,j,1)
                else
                  clig(i,j,1) = 0.0
                end if
                sink = clir(i,j,1) + clis(i,j,1) + clig(i,j,1)
                if (qi(i,j,1) < sink) then
                  a = qi(i,j,1) / sink
                  clir(i,j,1) = clir(i,j,1) * a
                  clis(i,j,1) = clis(i,j,1) * a
                  clig(i,j,1) = clig(i,j,1) * a
                end if
              else
                clis(i,j,1) = 0.0
                clig(i,j,1) = 0.0
              end if
            else
              clis(i,j,1) = 0.0
              clig(i,j,1) = 0.0
            end if

            ! Snow and graupel
            if (qs(i,j,1) > thresq) then
              if (qg(i,j,1) > thresq) then
                if (t(i,j,1) < t0) then
                  a = diaqs(i,j,1) * diaqg(i,j,1)
                  b = ncs(i,j,1) * ncg(i,j,1) * rbr(i,j,1)
                  c = usq(i,j,1) - ugq(i,j,1)
                  clsg(i,j,1) = b * cclsg * (5.0 * diaqs2 + 2.0 * a + 0.5 * diaqg2) &
                              * diaqs3 * sqrt(c * c + 0.04 * usq(i,j,1) * ugq(i,j,1))
                  c = usn(i,j,1) - ugn(i,j,1)
                  clsgn(i,j,1) = b * cclsgn * (diaqs2 + a + diaqg2) &
                               * sqrt(c * c + 0.04 * usn(i,j,1) * ugn(i,j,1))
                else
                  a = diaqs(i,j,1) * diaqg(i,j,1)
                  b = esgiv * ncs(i,j,1) * ncg(i,j,1) * rbr(i,j,1)
                  c = usq(i,j,1) - ugq(i,j,1)
                  clsg(i,j,1) = b * cclsg * (5.0 * diaqs2 + 2.0 * a + 0.5 * diaqg2) &
                              * diaqs3 * sqrt(c * c + 0.04 * usq(i,j,1) * ugq(i,j,1))
                  c = usn(i,j,1) - ugn(i,j,1)
                  clsgn(i,j,1) = b * cclsgn * (diaqs2 + a + diaqg2) &
                               * sqrt(c * c + 0.04 * usn(i,j,1) * ugn(i,j,1))
                end if
              else
                clsg(i,j,1) = 0.0
                clsgn(i,j,1) = 0.0
              end if

              sink = clsr(i,j,1) + clsg(i,j,1)
              if (qs(i,j,1) < sink) then
                a = qs(i,j,1) / sink
                clsr(i,j,1) = clsr(i,j,1) * a
                clsg(i,j,1) = clsg(i,j,1) * a
              end if

              sink = clsrn(i,j,1) + clsgn(i,j,1)
              if (ncs(i,j,1) < sink) then
                a = ncs(i,j,1) / sink
                clsrn(i,j,1) = clsrn(i,j,1) * a
                clsgn(i,j,1) = clsgn(i,j,1) * a
              end if
            else
              clsg(i,j,1) = 0.0
              clsgn(i,j,1) = 0.0
            end if
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine run_kernel

  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    integer :: unit_num, ios
    character(len=256) :: line
    character(len=64) :: key
    character(len=128) :: value_str
    integer :: eq_pos

    unit_num = 20
    open(unit_num, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening params file: ', trim(filename)
      stop 1
    end if

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      eq_pos = index(line, '=')
      if (eq_pos == 0) cycle
      key = adjustl(line(1:eq_pos-1))
      value_str = adjustl(line(eq_pos+1:))

      select case (trim(key))
        case ('cphopt'); read(value_str, *) cphopt
        case ('ni'); read(value_str, *) ni
        case ('nj'); read(value_str, *) nj
        case ('nk'); read(value_str, *) nk
        case ('dtb'); read(value_str, *) dtb
        case ('thresq'); read(value_str, *) thresq
        case ('t0'); read(value_str, *) t0
        case ('r0'); read(value_str, *) r0
        case ('rhow'); read(value_str, *) rhow
        case ('rhos'); read(value_str, *) rhos
        case ('aur'); read(value_str, *) aur
        case ('aus'); read(value_str, *) aus
        case ('aug'); read(value_str, *) aug
        case ('bur'); read(value_str, *) bur
        case ('bus'); read(value_str, *) bus
        case ('bug'); read(value_str, *) bug
        case ('eir'); read(value_str, *) eir
        case ('eis'); read(value_str, *) eis
        case ('eig'); read(value_str, *) eig
        case ('ers'); read(value_str, *) ers
        case ('erg'); read(value_str, *) erg
        case ('esr'); read(value_str, *) esr
        case ('esg'); read(value_str, *) esg
        case ('cc'); read(value_str, *) cc
        case ('oned9'); read(value_str, *) oned9
        case ('gf3bur'); read(value_str, *) gf3bur
        case ('gf3bus'); read(value_str, *) gf3bus
        case ('gf3bug'); read(value_str, *) gf3bug
      end select
    end do
    close(unit_num)
  end subroutine read_params

  subroutine read_array_3d(filename, arr, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni, nj, nk
    real, intent(out) :: arr(0:ni+1, 0:nj+1, 1:nk)
    integer :: unit_num, ios

    unit_num = 30
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', convert='big_endian', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_3d

  subroutine validate_array(name, arr, ref, tol, err_count, max_err)
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:nk)
    real(8), intent(in) :: tol
    integer, intent(out) :: err_count
    real(8), intent(out) :: max_err

    integer :: i, j, k
    real(8) :: rel_err, abs_val

    err_count = 0
    max_err = 0.0d0

    do k = 1, nk
      do j = 1, nj-1
        do i = 1, ni-1
          abs_val = abs(dble(ref(i,j,k)))
          if (abs_val > 1.0d-30) then
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k))) / abs_val
          else
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k)))
          end if
          max_err = max(max_err, rel_err)
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

    print '(A,A,A,I0,A,ES12.4)', 'Array ', trim(name), ': errors=', err_count, &
          ', max_rel_err=', max_err
  end subroutine validate_array

end program kernel_benchmark

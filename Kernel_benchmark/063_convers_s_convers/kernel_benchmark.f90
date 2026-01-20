program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup, num_iterations
  real(8) :: tolerance

  ! Dimensions and options
  integer :: cphopt, ni, nj, nk

  ! Physical constants (from params.txt)
  real :: dtb, thresq, t0, rhow, rhoi, rhos, rhog, ms0, diaqcm, cc, oned3, oned6

  ! Additional physical constants (hardcoded from m_comphy and m_commath)
  real, parameter :: tlow = 233.16
  real, parameter :: r0 = 1.225
  real, parameter :: xl = 0.25
  real, parameter :: aui = 770.0
  real, parameter :: aus = 27.0
  real, parameter :: bus = 0.5
  real, parameter :: ecc = 0.55
  real, parameter :: eii = 0.1
  real, parameter :: g = 9.8
  real, parameter :: gfbus = 1.77245
  real, parameter :: sevnd3 = 7.0 / 3.0

  ! Derived values (computed before parallel region)
  real :: busm1, ms05, qccm, diaqs0, cagin, ccncr, ccnsg, ccnsgn

  ! Input arrays
  real, allocatable :: t(:,:,:), rbr(:,:,:), rbv(:,:,:)
  real, allocatable :: qc(:,:,:), qi(:,:,:), qs(:,:,:)
  real, allocatable :: ncc(:,:,:), ncs(:,:,:)
  real, allocatable :: mu(:,:,:), mi(:,:,:)
  real, allocatable :: diaqi(:,:,:), diaqs(:,:,:)
  real, allocatable :: clcs(:,:,:), vdvi(:,:,:), vdvs(:,:,:), ecs(:,:,:)

  ! Output arrays
  real, allocatable :: cncr(:,:,:), cnis(:,:,:), cnsg(:,:,:), cnsgn(:,:,:)

  ! Reference arrays
  real, allocatable :: cncr_ref(:,:,:), cnis_ref(:,:,:)
  real, allocatable :: cnsg_ref(:,:,:), cnsgn_ref(:,:,:)

  ! Timing variables
  real(8) :: start_time, end_time, elapsed_time
  real(8) :: total_time, min_time, max_time, avg_time
  real(8), allocatable :: times(:)

  ! Validation variables
  integer :: error_count, total_errors
  real(8) :: max_rel_error

  ! Loop variables
  integer :: iter, i, j, k
  real :: qcm

  ! Config file
  integer :: unit_conf
  character(len=256) :: line

  ! Read configuration
  unit_conf = 10
  open(unit_conf, file='benchmark.conf', status='old', action='read')
  read(unit_conf, '(A)') data_dir
  read(unit_conf, *) num_warmup
  read(unit_conf, *) num_iterations
  read(unit_conf, *) tolerance
  close(unit_conf)

  data_dir = trim(adjustl(data_dir))
  print '(A)', '=== Convers Kernel Benchmark ==='
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
  print '(A,ES12.4)', 'dtb = ', dtb
  print '(A,ES12.4)', 'thresq = ', thresq

  ! Compute derived values (matching the original code)
  busm1 = bus - 1.0
  ms05 = 0.5 * ms0
  qccm = oned6 * rhow * cc * diaqcm * diaqcm * diaqcm
  diaqs0 = exp(oned3 * log(6.0 * ms0 / (cc * rhos)))
  cagin = r0 * exp(3.0 * log(0.5 * oned3 * aui * eii * xl / rhoi * dtb))
  ccncr = 0.104 * g * ecc * exp(-oned3 * log(rhow)) * dtb
  ccnsgn = 1.0 / (rhog - rhos)
  ccnsg = rhog * ccnsgn
  ccnsgn = 1.5 * aus * gfbus * sqrt(r0) * ccnsgn * dtb

  ! Allocate arrays
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qc(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncc(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncs(0:ni+1, 0:nj+1, 1:nk))
  allocate(mu(0:ni+1, 0:nj+1, 1:nk))
  allocate(mi(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqi(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcs(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvi(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvs(0:ni+1, 0:nj+1, 1:nk))
  allocate(ecs(0:ni+1, 0:nj+1, 1:nk))

  allocate(cncr(0:ni+1, 0:nj+1, 1:nk))
  allocate(cnis(0:ni+1, 0:nj+1, 1:nk))
  allocate(cnsg(0:ni+1, 0:nj+1, 1:nk))
  allocate(cnsgn(0:ni+1, 0:nj+1, 1:nk))

  allocate(cncr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(cnis_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(cnsg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(cnsgn_ref(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  ! Read input arrays
  print '(A)', ''
  print '(A)', '--- Loading input data ---'
  call read_array_3d(trim(data_dir) // '/t.bin', t, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/rbr.bin', rbr, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/rbv.bin', rbv, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qc.bin', qc, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qi.bin', qi, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qs.bin', qs, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ncc.bin', ncc, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ncs.bin', ncs, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/mu.bin', mu, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/mi.bin', mi, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/diaqi.bin', diaqi, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/diaqs.bin', diaqs, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/clcs.bin', clcs, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/vdvi.bin', vdvi, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/vdvs.bin', vdvs, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ecs.bin', ecs, ni, nj, nk)

  ! Read reference arrays
  call read_array_3d(trim(data_dir) // '/cncr_ref.bin', cncr_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/cnis_ref.bin', cnis_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/cnsg_ref.bin', cnsg_ref, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/cnsgn_ref.bin', cnsgn_ref, ni, nj, nk)
  print '(A)', 'Data loaded successfully'

  ! Warmup iterations
  print '(A)', ''
  print '(A)', '--- Warmup ---'
  do iter = 1, num_warmup
    call run_kernel()
  end do
  print '(A,I0,A)', 'Completed ', num_warmup, ' warmup iterations'

  ! Benchmark iterations
  print '(A)', ''
  print '(A)', '--- Benchmark ---'
  total_time = 0.0d0
  min_time = huge(1.0d0)
  max_time = 0.0d0

  do iter = 1, num_iterations
    start_time = omp_get_wtime()
    call run_kernel()
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

  call validate_array('cncr', cncr, cncr_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  call validate_array('cnis', cnis, cnis_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  call validate_array('cnsg', cnsg, cnsg_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  call validate_array('cnsgn', cnsgn, cnsgn_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  print '(A)', ''
  if (total_errors == 0) then
    print '(A)', '*** VALIDATION PASSED ***'
  else
    print '(A,I0,A)', '*** VALIDATION FAILED: ', total_errors, ' errors ***'
  end if

  ! Cleanup
  deallocate(t, rbr, rbv, qc, qi, qs, ncc, ncs, mu, mi)
  deallocate(diaqi, diaqs, clcs, vdvi, vdvs, ecs)
  deallocate(cncr, cnis, cnsg, cnsgn)
  deallocate(cncr_ref, cnis_ref, cnsg_ref, cnsgn_ref)
  deallocate(times)

contains

  subroutine run_kernel()
    integer :: i, j, k
    real :: qcm

    !$omp parallel default(shared) private(k)

    ! In the case nk = 1
    if (nk == 1) then
      ! Case abs(cphopt) >= 3
      if (abs(cphopt) >= 3) then
        !$omp do schedule(runtime) private(i,j,qcm)
        do j = 1, nj-1
          do i = 1, ni-1
            ! Calculate conversion rate from cloud water to rain water
            if (qc(i,j,1) > thresq) then
              qcm = qccm * ncc(i,j,1)
              if (t(i,j,1) > tlow .and. qc(i,j,1) >= qcm) then
                cncr(i,j,1) = ccncr * rbr(i,j,1) * exp(-oned3 * log(ncc(i,j,1))) &
                            * exp(sevnd3 * log(qc(i,j,1))) / mu(i,j,1)
              else
                cncr(i,j,1) = 0.0
              end if
            else
              cncr(i,j,1) = 0.0
            end if

            ! Calculate conversion rates from cloud ice to snow and snow to graupel
            if (t(i,j,1) < t0) then
              ! Cloud ice to snow
              if (qi(i,j,1) > thresq) then
                cnis(i,j,1) = rbr(i,j,1) * exp(oned3 * log(cagin * rbv(i,j,1))) &
                            * qi(i,j,1) * qi(i,j,1) / log(diaqs0 / diaqi(i,j,1))
                if (vdvi(i,j,1) >= 0.0) then
                  if (mi(i,j,1) < ms05) then
                    cnis(i,j,1) = cnis(i,j,1) + mi(i,j,1) / (ms0 - mi(i,j,1)) * vdvi(i,j,1)
                  else
                    cnis(i,j,1) = cnis(i,j,1) + (vdvi(i,j,1) + (1.0 - ms05 / mi(i,j,1)) * qi(i,j,1))
                  end if
                end if
                cnis(i,j,1) = min(cnis(i,j,1), qi(i,j,1))
              else
                cnis(i,j,1) = 0.0
              end if

              ! Snow to graupel
              if (qs(i,j,1) > thresq) then
                if (vdvs(i,j,1) > 0.0 .and. vdvs(i,j,1) < clcs(i,j,1)) then
                  cnsg(i,j,1) = min(ccnsg * clcs(i,j,1), qs(i,j,1))
                  cnsgn(i,j,1) = min(ccnsgn * ecs(i,j,1) * qc(i,j,1) * ncs(i,j,1) &
                               * sqrt(rbr(i,j,1)) * exp(busm1 * log(diaqs(i,j,1))), ncs(i,j,1))
                else
                  cnsg(i,j,1) = 0.0
                  cnsgn(i,j,1) = 0.0
                end if
              else
                cnsg(i,j,1) = 0.0
                cnsgn(i,j,1) = 0.0
              end if
            else
              cnis(i,j,1) = 0.0
              cnsg(i,j,1) = 0.0
              cnsgn(i,j,1) = 0.0
            end if
          end do
        end do
        !$omp end do

      ! Case abs(cphopt) == 2
      else if (abs(cphopt) == 2) then
        !$omp do schedule(runtime) private(i,j,qcm)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qc(i,j,1) > thresq) then
              qcm = qccm * ncc(i,j,1)
              if (t(i,j,1) > tlow .and. qc(i,j,1) >= qcm) then
                cncr(i,j,1) = ccncr * rbr(i,j,1) * exp(-oned3 * log(ncc(i,j,1))) &
                            * exp(sevnd3 * log(qc(i,j,1))) / mu(i,j,1)
              else
                cncr(i,j,1) = 0.0
              end if
            else
              cncr(i,j,1) = 0.0
            end if

            if (t(i,j,1) < t0) then
              if (qi(i,j,1) > thresq) then
                cnis(i,j,1) = rbr(i,j,1) * exp(oned3 * log(cagin * rbv(i,j,1))) &
                            * qi(i,j,1) * qi(i,j,1) / log(diaqs0 / diaqi(i,j,1))
                if (vdvi(i,j,1) >= 0.0) then
                  if (mi(i,j,1) < ms05) then
                    cnis(i,j,1) = cnis(i,j,1) + mi(i,j,1) / (ms0 - mi(i,j,1)) * vdvi(i,j,1)
                  else
                    cnis(i,j,1) = cnis(i,j,1) + (vdvi(i,j,1) + (1.0 - ms05 / mi(i,j,1)) * qi(i,j,1))
                  end if
                end if
                cnis(i,j,1) = min(cnis(i,j,1), qi(i,j,1))
              else
                cnis(i,j,1) = 0.0
              end if

              if (qs(i,j,1) > thresq) then
                if (vdvs(i,j,1) > 0.0 .and. vdvs(i,j,1) < clcs(i,j,1)) then
                  cnsg(i,j,1) = min(ccnsg * clcs(i,j,1), qs(i,j,1))
                else
                  cnsg(i,j,1) = 0.0
                end if
              else
                cnsg(i,j,1) = 0.0
              end if
            else
              cnis(i,j,1) = 0.0
              cnsg(i,j,1) = 0.0
            end if
          end do
        end do
        !$omp end do
      end if

    ! In the case nk > 1
    else
      if (abs(cphopt) >= 3) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j,qcm)
          do j = 1, nj-1
            do i = 1, ni-1
              if (qc(i,j,k) > thresq) then
                qcm = qccm * ncc(i,j,k)
                if (t(i,j,k) > tlow .and. qc(i,j,k) >= qcm) then
                  cncr(i,j,k) = ccncr * rbr(i,j,k) * exp(-oned3 * log(ncc(i,j,k))) &
                              * exp(sevnd3 * log(qc(i,j,k))) / mu(i,j,k)
                else
                  cncr(i,j,k) = 0.0
                end if
              else
                cncr(i,j,k) = 0.0
              end if

              if (t(i,j,k) < t0) then
                if (qi(i,j,k) > thresq) then
                  cnis(i,j,k) = rbr(i,j,k) * exp(oned3 * log(cagin * rbv(i,j,k))) &
                              * qi(i,j,k) * qi(i,j,k) / log(diaqs0 / diaqi(i,j,k))
                  if (vdvi(i,j,k) >= 0.0) then
                    if (mi(i,j,k) < ms05) then
                      cnis(i,j,k) = cnis(i,j,k) + mi(i,j,k) / (ms0 - mi(i,j,k)) * vdvi(i,j,k)
                    else
                      cnis(i,j,k) = cnis(i,j,k) + (vdvi(i,j,k) + (1.0 - ms05 / mi(i,j,k)) * qi(i,j,k))
                    end if
                  end if
                  cnis(i,j,k) = min(cnis(i,j,k), qi(i,j,k))
                else
                  cnis(i,j,k) = 0.0
                end if

                if (qs(i,j,k) > thresq) then
                  if (vdvs(i,j,k) > 0.0 .and. vdvs(i,j,k) < clcs(i,j,k)) then
                    cnsg(i,j,k) = min(ccnsg * clcs(i,j,k), qs(i,j,k))
                    cnsgn(i,j,k) = min(ccnsgn * ecs(i,j,k) * qc(i,j,k) * ncs(i,j,k) &
                                 * sqrt(rbr(i,j,k)) * exp(busm1 * log(diaqs(i,j,k))), ncs(i,j,k))
                  else
                    cnsg(i,j,k) = 0.0
                    cnsgn(i,j,k) = 0.0
                  end if
                else
                  cnsg(i,j,k) = 0.0
                  cnsgn(i,j,k) = 0.0
                end if
              else
                cnis(i,j,k) = 0.0
                cnsg(i,j,k) = 0.0
                cnsgn(i,j,k) = 0.0
              end if
            end do
          end do
          !$omp end do
        end do

      else if (abs(cphopt) == 2) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j,qcm)
          do j = 1, nj-1
            do i = 1, ni-1
              if (qc(i,j,k) > thresq) then
                qcm = qccm * ncc(i,j,k)
                if (t(i,j,k) > tlow .and. qc(i,j,k) >= qcm) then
                  cncr(i,j,k) = ccncr * rbr(i,j,k) * exp(-oned3 * log(ncc(i,j,k))) &
                              * exp(sevnd3 * log(qc(i,j,k))) / mu(i,j,k)
                else
                  cncr(i,j,k) = 0.0
                end if
              else
                cncr(i,j,k) = 0.0
              end if

              if (t(i,j,k) < t0) then
                if (qi(i,j,k) > thresq) then
                  cnis(i,j,k) = rbr(i,j,k) * exp(oned3 * log(cagin * rbv(i,j,k))) &
                              * qi(i,j,k) * qi(i,j,k) / log(diaqs0 / diaqi(i,j,k))
                  if (vdvi(i,j,k) >= 0.0) then
                    if (mi(i,j,k) < ms05) then
                      cnis(i,j,k) = cnis(i,j,k) + mi(i,j,k) / (ms0 - mi(i,j,k)) * vdvi(i,j,k)
                    else
                      cnis(i,j,k) = cnis(i,j,k) + (vdvi(i,j,k) + (1.0 - ms05 / mi(i,j,k)) * qi(i,j,k))
                    end if
                  end if
                  cnis(i,j,k) = min(cnis(i,j,k), qi(i,j,k))
                else
                  cnis(i,j,k) = 0.0
                end if

                if (qs(i,j,k) > thresq) then
                  if (vdvs(i,j,k) > 0.0 .and. vdvs(i,j,k) < clcs(i,j,k)) then
                    cnsg(i,j,k) = min(ccnsg * clcs(i,j,k), qs(i,j,k))
                  else
                    cnsg(i,j,k) = 0.0
                  end if
                else
                  cnsg(i,j,k) = 0.0
                end if
              else
                cnis(i,j,k) = 0.0
                cnsg(i,j,k) = 0.0
              end if
            end do
          end do
          !$omp end do
        end do
      end if
    end if

    !$omp end parallel
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
        case ('cphopt')
          read(value_str, *) cphopt
        case ('ni')
          read(value_str, *) ni
        case ('nj')
          read(value_str, *) nj
        case ('nk')
          read(value_str, *) nk
        case ('dtb')
          read(value_str, *) dtb
        case ('thresq')
          read(value_str, *) thresq
        case ('t0')
          read(value_str, *) t0
        case ('rhow')
          read(value_str, *) rhow
        case ('rhoi')
          read(value_str, *) rhoi
        case ('rhos')
          read(value_str, *) rhos
        case ('rhog')
          read(value_str, *) rhog
        case ('ms0')
          read(value_str, *) ms0
        case ('diaqcm')
          read(value_str, *) diaqcm
        case ('cc')
          read(value_str, *) cc
        case ('oned3')
          read(value_str, *) oned3
        case ('oned6')
          read(value_str, *) oned6
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
          if (rel_err > tol) then
            err_count = err_count + 1
            if (err_count <= 5) then
              print '(A,A,A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4,A,ES12.4)', &
                'Error in ', trim(name), ' at (', i, ',', j, ',', k, &
                '): computed=', arr(i,j,k), ' ref=', ref(i,j,k), ' rel_err=', rel_err
            end if
          end if
        end do
      end do
    end do

    print '(A,A,A,I0,A,ES12.4)', 'Array ', trim(name), ': errors=', err_count, &
          ', max_rel_err=', max_err
  end subroutine validate_array

end program kernel_benchmark

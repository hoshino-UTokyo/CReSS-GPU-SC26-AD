!> GPU Kernel Benchmark: nuc1stc (s_nuc1stc)
!> Compute ice nucleation rates (freezing, contact, homogeneous)
!> GPU Difficulty: Hard - Complex conditionals with many private variables
program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup, num_iterations
  real(8) :: tolerance

  ! Dimensions
  integer :: ni, nj, nk

  ! Physical constants from params.txt
  real :: dtb, thresq, t0, rv, rhow
  real :: lmb20, p20, t20, diasl, boltz, kpa, cc

  ! Additional physical constant (hardcoded from m_comphy)
  real, parameter :: tlow = 233.16

  ! Derived values (computed before parallel region)
  real :: rwdt2, cc45, cknd, cdar, kpa25, kpa50

  ! Input arrays
  real, allocatable :: t(:,:,:), p(:,:,:), qc(:,:,:), ncc(:,:,:)
  real, allocatable :: tcel(:,:,:), lv(:,:,:), kp(:,:,:), mu(:,:,:)
  real, allocatable :: diaqc(:,:,:)

  ! Output array
  real, allocatable :: nuci(:,:,:)

  ! Reference array
  real, allocatable :: nuci_ref(:,:,:)

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
  print '(A)', '=== GPU Kernel Benchmark: nuc1stc ==='
  print '(A,A)', 'Data directory: ', trim(data_dir)
  print '(A,I0)', 'Warmup iterations: ', num_warmup
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,ES10.2)', 'Tolerance: ', tolerance

  ! Read parameters
  call read_params(trim(data_dir) // '/params.txt')

  print '(A)', ''
  print '(A)', '--- Parameters ---'
  print '(A,I0,A,I0,A,I0)', 'ni x nj x nk = ', ni, ' x ', nj, ' x ', nk
  print '(A,ES12.4)', 'dtb = ', dtb

  ! Compute derived values
  rwdt2 = (100.0 / rhow) * dtb
  cc45 = 4.0e5 * cc
  cknd = 2.0 * lmb20 * p20 / (t20 * diasl)
  cdar = boltz / (3.0 * cc * diasl)
  kpa25 = 2.5 * kpa
  kpa50 = 5.0 * kpa

  ! Allocate arrays
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(p(0:ni+1, 0:nj+1, 1:nk))
  allocate(qc(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tcel(0:ni+1, 0:nj+1, 1:nk))
  allocate(lv(0:ni+1, 0:nj+1, 1:nk))
  allocate(kp(0:ni+1, 0:nj+1, 1:nk))
  allocate(mu(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqc(0:ni+1, 0:nj+1, 1:nk))
  allocate(nuci(0:ni+1, 0:nj+1, 1:nk))
  allocate(nuci_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  ! Read input arrays
  print '(A)', ''
  print '(A)', '--- Loading input data ---'
  call read_array_3d(trim(data_dir) // '/t.bin', t, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/p.bin', p, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qc.bin', qc, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/ncc.bin', ncc, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/tcel.bin', tcel, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/lv.bin', lv, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/kp.bin', kp, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/mu.bin', mu, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/diaqc.bin', diaqc, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/nuci_ref.bin', nuci_ref, ni, nj, nk)
  print '(A)', 'Data loaded successfully'

  ! Warmup iterations
  print '(A)', ''
  print '(A)', '--- Warmup ---'
  do iter = 1, num_warmup
    call run_kernel()
  end do
  !$acc wait
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

  call validate_array('nuci', nuci, nuci_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  print '(A)', ''
  if (total_errors == 0) then
    print '(A)', '*** VALIDATION PASSED ***'
  else
    print '(A,I0,A)', '*** VALIDATION FAILED: ', total_errors, ' errors ***'
  end if

  ! Cleanup
  deallocate(t, p, qc, ncc, tcel, lv, kp, mu, diaqc)
  deallocate(nuci, nuci_ref, times)

contains

  subroutine run_kernel()
    integer :: i, j, k
    real :: tc, piv, knd, dar, f1, f2, ft
    real :: nufci, nucci, nuhci

    if (nk == 1) then
      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent private(tc, piv, knd, dar, f1, f2, ft, nufci, nucci, nuhci)
        do i = 1, ni-1
          if (qc(i,j,1) > thresq) then
            ! Condensation nucleation
            if (t(i,j,1) > tlow .and. t(i,j,1) < t0) then
              nufci = min(rwdt2 * (exp(-0.66 * tcel(i,j,1)) - 1.0) &
                        * qc(i,j,1) * qc(i,j,1) / ncc(i,j,1), qc(i,j,1))
            else
              nufci = 0.0
            end if

            ! Contact nucleation
            tc = t(i,j,1)
            if (tc > tlow .and. tc < 270.16) then
              piv = 1.0 / p(i,j,1)
              knd = cknd * piv * t(i,j,1)
              dar = cdar * tc * (1.0 + knd) / mu(i,j,1)
              f1 = cc45 * exp(1.3 * log(270.16 - tc)) * diaqc(i,j,1)

              if (tc < t(i,j,1)) then
                f2 = kpa * piv * (t(i,j,1) - tc)
                ft = (kp(i,j,1) + kpa25 * knd) &
                   * (0.4 + 0.58 * knd + 0.16 * exp(-1.0 / knd)) &
                   / ((1.0 + 3.0 * knd) * (2.0 * kp(i,j,1) + kpa50 * knd + kpa))
                nucci = f1 * f2 * (rv * t(i,j,1) / lv(i,j,1) + ft)
              else
                nucci = 0.0
              end if
              nucci = (nucci + f1 * dar) * qc(i,j,1) * dtb
            else
              nucci = 0.0
            end if

            ! Homogeneous nucleation
            if (t(i,j,1) <= tlow) then
              nuhci = qc(i,j,1)
            else
              nuhci = 0.0
            end if

            ! Total nucleation rate
            nuci(i,j,1) = nufci + nucci + nuhci
          else
            nuci(i,j,1) = 0.0
          end if
        end do
      end do
      !$acc end kernels

    else
      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(tc, piv, knd, dar, f1, f2, ft, nufci, nucci, nuhci)
          do i = 1, ni-1
            if (qc(i,j,k) > thresq) then
              if (t(i,j,k) > tlow .and. t(i,j,k) < t0) then
                nufci = min(rwdt2 * (exp(-0.66 * tcel(i,j,k)) - 1.0) &
                          * qc(i,j,k) * qc(i,j,k) / ncc(i,j,k), qc(i,j,k))
              else
                nufci = 0.0
              end if

              tc = t(i,j,k)
              if (tc > tlow .and. tc < 270.16) then
                piv = 1.0 / p(i,j,k)
                knd = cknd * piv * t(i,j,k)
                dar = cdar * tc * (1.0 + knd) / mu(i,j,k)
                f1 = cc45 * exp(1.3 * log(270.16 - tc)) * diaqc(i,j,k)

                if (tc < t(i,j,k)) then
                  f2 = kpa * piv * (t(i,j,k) - tc)
                  ft = (kp(i,j,k) + kpa25 * knd) &
                     * (0.4 + 0.58 * knd + 0.16 * exp(-1.0 / knd)) &
                     / ((1.0 + 3.0 * knd) * (2.0 * kp(i,j,k) + kpa50 * knd + kpa))
                  nucci = f1 * f2 * (rv * t(i,j,k) / lv(i,j,k) + ft)
                else
                  nucci = 0.0
                end if
                nucci = (nucci + f1 * dar) * qc(i,j,k) * dtb
              else
                nucci = 0.0
              end if

              if (t(i,j,k) <= tlow) then
                nuhci = qc(i,j,k)
              else
                nuhci = 0.0
              end if

              nuci(i,j,k) = nufci + nucci + nuhci
            else
              nuci(i,j,k) = 0.0
            end if
          end do
        end do
      end do
      !$acc end kernels
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
        case ('ni'); read(value_str, *) ni
        case ('nj'); read(value_str, *) nj
        case ('nk'); read(value_str, *) nk
        case ('dtb'); read(value_str, *) dtb
        case ('thresq'); read(value_str, *) thresq
        case ('t0'); read(value_str, *) t0
        case ('rv'); read(value_str, *) rv
        case ('rhow'); read(value_str, *) rhow
        case ('lmb20'); read(value_str, *) lmb20
        case ('p20'); read(value_str, *) p20
        case ('t20'); read(value_str, *) t20
        case ('diasl'); read(value_str, *) diasl
        case ('boltz'); read(value_str, *) boltz
        case ('kpa'); read(value_str, *) kpa
        case ('cc'); read(value_str, *) cc
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

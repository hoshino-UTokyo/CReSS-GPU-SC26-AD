!***********************************************************************
! GPU Kernel Benchmark: setblk (s_setblk)
!***********************************************************************
!
! Source: Src/setblk.f90
! Description: Calculate thermodynamic properties (T, saturation, latent heat),
!              air properties (viscosity, conductivity), and hydrometeor parameters
! GPU Port: OpenACC with Unified Memory
!
! Note: The nu(i,j) array is updated in the first loop section and used in the
!       second loop section. We use separate acc kernel blocks to handle this
!       dependency properly.
!
!***********************************************************************
program kernel_benchmark_gpu_setblk
  use omp_lib
  implicit none

  ! === Grid dimensions ===
  integer :: cphopt, ni, nj, nk

  ! === Scalar parameters ===
  real :: thresq, epsva, es0, lv0, lf0, r0
  real :: ccri6, ccrw6, cdiaqc, cdiaqg, cdiaqr, cdiaqs
  real :: cdv, cnu, cvntg, cvntr, cvnts, cwmci, epses0
  real :: pdiaqg, pdiaqr, pdiaqs, t23iv

  ! === Physical constants (from m_comphy and m_commath) ===
  real, parameter :: t0 = 273.16e0    ! Melting point temperature (K)
  real, parameter :: p20 = 101325.e0  ! Standard pressure (Pa)
  real, parameter :: kp0 = 10.79e0    ! Thermal conductivity coefficient
  real, parameter :: oned3 = 1.e0/3.e0  ! 1/3

  ! === Input arrays ===
  real, allocatable :: ptbr(:,:,:), rbr(:,:,:), rbv(:,:,:)
  real, allocatable :: pi(:,:,:), p(:,:,:), ptp(:,:,:)
  real, allocatable :: qv(:,:,:), qc(:,:,:), qr(:,:,:)
  real, allocatable :: qi(:,:,:), qs(:,:,:), qg(:,:,:)
  real, allocatable :: ncc(:,:,:), ncr(:,:,:), nci(:,:,:)
  real, allocatable :: ncs(:,:,:), ncg(:,:,:)

  ! === Output arrays ===
  real, allocatable :: t(:,:,:), tcel(:,:,:)
  real, allocatable :: qvsst0(:,:,:), qvsw(:,:,:), qvsi(:,:,:)
  real, allocatable :: lv(:,:,:), ls(:,:,:), lf(:,:,:)
  real, allocatable :: kp(:,:,:), mu(:,:,:), dv(:,:,:)
  real, allocatable :: mi(:,:,:)
  real, allocatable :: diaqc(:,:,:), diaqr(:,:,:), diaqi(:,:,:)
  real, allocatable :: diaqs(:,:,:), diaqg(:,:,:)
  real, allocatable :: vntr(:,:,:), vnts(:,:,:), vntg(:,:,:)

  ! === 2D working array (inout) ===
  real, allocatable :: nu(:,:), nu_in(:,:)

  ! === Reference arrays for validation ===
  real, allocatable :: tcel_ref(:,:,:)
  real, allocatable :: qvsst0_ref(:,:,:), qvsw_ref(:,:,:), qvsi_ref(:,:,:)
  real, allocatable :: lv_ref(:,:,:), ls_ref(:,:,:), lf_ref(:,:,:)
  real, allocatable :: kp_ref(:,:,:), mu_ref(:,:,:), dv_ref(:,:,:)
  real, allocatable :: mi_ref(:,:,:)
  real, allocatable :: diaqc_ref(:,:,:), diaqr_ref(:,:,:), diaqi_ref(:,:,:)
  real, allocatable :: diaqs_ref(:,:,:), diaqg_ref(:,:,:)
  real, allocatable :: vntr_ref(:,:,:), vnts_ref(:,:,:), vntg_ref(:,:,:)
  real, allocatable :: nu_ref(:,:)

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
  ! Input 3D arrays
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbv(0:ni+1, 0:nj+1, 1:nk))
  allocate(pi(0:ni+1, 0:nj+1, 1:nk))
  allocate(p(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qc(0:ni+1, 0:nj+1, 1:nk))
  allocate(qr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(qg(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncc(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncr(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncs(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncg(0:ni+1, 0:nj+1, 1:nk))

  ! Output 3D arrays
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(tcel(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsst0(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsw(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsi(0:ni+1, 0:nj+1, 1:nk))
  allocate(lv(0:ni+1, 0:nj+1, 1:nk))
  allocate(ls(0:ni+1, 0:nj+1, 1:nk))
  allocate(lf(0:ni+1, 0:nj+1, 1:nk))
  allocate(kp(0:ni+1, 0:nj+1, 1:nk))
  allocate(mu(0:ni+1, 0:nj+1, 1:nk))
  allocate(dv(0:ni+1, 0:nj+1, 1:nk))
  allocate(mi(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqc(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqr(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqi(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqs(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqg(0:ni+1, 0:nj+1, 1:nk))
  allocate(vntr(0:ni+1, 0:nj+1, 1:nk))
  allocate(vnts(0:ni+1, 0:nj+1, 1:nk))
  allocate(vntg(0:ni+1, 0:nj+1, 1:nk))

  ! 2D working array
  allocate(nu(0:ni+1, 0:nj+1))
  allocate(nu_in(0:ni+1, 0:nj+1))

  ! Reference arrays
  allocate(tcel_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsst0_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsw_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsi_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(lv_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ls_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(lf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(kp_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(mu_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(dv_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(mi_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqi_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqs_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vntr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vnts_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vntg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(nu_ref(0:ni+1, 0:nj+1))

  ! 4. Read input data
  call read_array_3d(trim(data_dir)//'/ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rbv.bin', rbv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pi.bin', pi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/p.bin', p, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qc.bin', qc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qr.bin', qr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qi.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qg.bin', qg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncc.bin', ncc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncr.bin', ncr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/nci.bin', nci, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncs.bin', ncs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncg.bin', ncg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/nu_in.bin', nu_in, 0, ni+1, 0, nj+1)

  ! 5. Read reference output
  call read_array_3d(trim(data_dir)//'/tcel_ref.bin', tcel_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvsst0_ref.bin', qvsst0_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvsw_ref.bin', qvsw_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvsi_ref.bin', qvsi_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/lv_ref.bin', lv_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ls_ref.bin', ls_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/lf_ref.bin', lf_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/kp_ref.bin', kp_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/mu_ref.bin', mu_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/dv_ref.bin', dv_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/mi_ref.bin', mi_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/diaqc_ref.bin', diaqc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/diaqr_ref.bin', diaqr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/diaqi_ref.bin', diaqi_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/diaqs_ref.bin', diaqs_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/diaqg_ref.bin', diaqg_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vntr_ref.bin', vntr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vnts_ref.bin', vnts_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vntg_ref.bin', vntg_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/nu_ref.bin', nu_ref, 0, ni+1, 0, nj+1)

  ! Print header
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: setblk'
  write(*,'(A)') '=================================================='
  write(*,'(A,I8,A,I8,A,I8)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I8)') ' cphopt: ', cphopt
  write(*,'(A,I8)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I8)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  ! 6. Warmup iterations
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    nu = nu_in  ! Reset 2D working array
    call kernel_setblk(cphopt, thresq, ni, nj, nk, &
         epsva, es0, lv0, lf0, r0, cwmci, epses0, t23iv, cnu, cdv, &
         pdiaqr, pdiaqs, pdiaqg, ccrw6, ccri6, cdiaqc, cdiaqr, cdiaqs, cdiaqg, &
         cvntr, cvnts, cvntg, &
         ptbr, rbr, rbv, pi, p, ptp, qv, qc, qr, qi, qs, qg, &
         ncc, ncr, nci, ncs, ncg, &
         t, tcel, qvsst0, qvsw, qvsi, lv, ls, lf, kp, mu, dv, mi, &
         diaqc, diaqr, diaqi, diaqs, diaqg, vntr, vnts, vntg, nu)
    !$acc wait
  end do

  ! 7. Benchmark iterations (timed)
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = huge(1.0d0)
  t_max = 0.0d0

  do iter = 1, num_iterations
    nu = nu_in  ! Reset 2D working array

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_setblk(cphopt, thresq, ni, nj, nk, &
         epsva, es0, lv0, lf0, r0, cwmci, epses0, t23iv, cnu, cdv, &
         pdiaqr, pdiaqs, pdiaqg, ccrw6, ccri6, cdiaqc, cdiaqr, cdiaqs, cdiaqg, &
         cvntr, cvnts, cvntg, &
         ptbr, rbr, rbv, pi, p, ptp, qv, qc, qr, qi, qs, qg, &
         ncc, ncr, nci, ncs, ncg, &
         t, tcel, qvsst0, qvsw, qvsi, lv, ls, lf, kp, mu, dv, mi, &
         diaqc, diaqr, diaqi, diaqs, diaqg, vntr, vnts, vntg, nu)
    !$acc wait
    t_end = omp_get_wtime()

    t_total = t_total + (t_end - t_start)
    if (t_end - t_start < t_min) t_min = t_end - t_start
    if (t_end - t_start > t_max) t_max = t_end - t_start
  end do
  t_avg = t_total / dble(num_iterations)

  ! 8. Validate output (only interior points: 1:ni-1, 1:nj-1)
  max_err_total = 0.0
  err_count_total = 0

  call validate_output_3d(tcel, tcel_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(qvsst0, qvsst0_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(qvsw, qvsw_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(qvsi, qvsi_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(lv, lv_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(ls, ls_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(lf, lf_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(kp, kp_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(mu, mu_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(dv, dv_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(mi, mi_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(diaqc, diaqc_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(diaqr, diaqr_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(diaqi, diaqi_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(diaqs, diaqs_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(diaqg, diaqg_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(vntr, vntr_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(vnts, vnts_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(vntg, vntg_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_2d(nu, nu_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1)
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  passed = (err_count_total == 0)

  ! 9. Report results
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', t_max * 1000.0d0, ' ms'
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

  if (.not. passed) stop 1

contains

  !-----------------------------------------------------------------
  ! The kernel (GPU version with OpenACC)
  !-----------------------------------------------------------------
  subroutine kernel_setblk(cphopt, thresq, ni, nj, nk, &
       epsva, es0, lv0, lf0, r0, cwmci, epses0, t23iv, cnu, cdv, &
       pdiaqr, pdiaqs, pdiaqg, ccrw6, ccri6, cdiaqc, cdiaqr, cdiaqs, cdiaqg, &
       cvntr, cvnts, cvntg, &
       ptbr, rbr, rbv, pi, p, ptp, qv, qc, qr, qi, qs, qg, &
       ncc, ncr, nci, ncs, ncg, &
       t, tcel, qvsst0, qvsw, qvsi, lv, ls, lf, kp, mu, dv, mi, &
       diaqc, diaqr, diaqi, diaqs, diaqg, vntr, vnts, vntg, nu)
    implicit none

    integer, intent(in) :: cphopt, ni, nj, nk
    real, intent(in) :: thresq, epsva, es0, lv0, lf0, r0
    real, intent(in) :: cwmci, epses0, t23iv, cnu, cdv
    real, intent(in) :: pdiaqr, pdiaqs, pdiaqg, ccrw6, ccri6
    real, intent(in) :: cdiaqc, cdiaqr, cdiaqs, cdiaqg
    real, intent(in) :: cvntr, cvnts, cvntg
    real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncc(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: nci(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncg(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: t(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: tcel(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: qvsst0(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: qvsw(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: qvsi(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: lv(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: ls(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: lf(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: kp(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: mu(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: dv(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: mi(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: diaqc(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: diaqr(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: diaqi(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: diaqs(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: diaqg(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: vntr(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: vnts(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: vntg(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: nu(0:ni+1,0:nj+1)

    integer :: i, j, k
    real :: esw, esi, p20dvp, cvnt

    if (nk == 1) then
      ! nk = 1 case

      ! First loop section: calculate thermodynamic properties and update nu(i,j)
      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent private(esw, esi, p20dvp)
        do i = 1, ni-1
          ! Calculate air temperature
          t(i,j,1) = (ptbr(i,j,1) + ptp(i,j,1)) * pi(i,j,1)
          tcel(i,j,1) = t(i,j,1) - t0

          ! Calculate saturation mixing ratio
          esw = es0 * exp(17.269e0 * tcel(i,j,1) / (t(i,j,1) - 35.86e0))
          esi = es0 * exp(21.875e0 * tcel(i,j,1) / (t(i,j,1) - 7.66e0))

          qvsst0(i,j,1) = qv(i,j,1) - epses0 / (p(i,j,1) - es0)
          qvsw(i,j,1) = epsva * esw / (p(i,j,1) - esw)
          qvsi(i,j,1) = epsva * esi / (p(i,j,1) - esi)

          ! Calculate latent heat
          lv(i,j,1) = lv0 * exp((.167e0 + 3.67e-4 * t(i,j,1)) * log(t0 / t(i,j,1)))
          lf(i,j,1) = lf0 + cwmci * tcel(i,j,1)
          ls(i,j,1) = lv(i,j,1) + lf(i,j,1)

          ! Calculate viscosity and diffusivity
          p20dvp = p20 / p(i,j,1)
          nu(i,j) = p20dvp * exp(1.754e0 * log(cnu * t(i,j,1)))
          kp(i,j,1) = (kp0 / (120.e0 + t(i,j,1))) * exp(1.5e0 * log(t23iv * t(i,j,1)))
          mu(i,j,1) = rbr(i,j,1) * nu(i,j)
          dv(i,j,1) = p20dvp * exp(1.81e0 * log(cdv * t(i,j,1)))
        end do
      end do
      !$acc end kernels

      ! Second loop section: calculate hydrometeor parameters using nu(i,j)
      ! Separate kernel block to ensure nu(i,j) is updated before use
      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent private(cvnt)
        do i = 1, ni-1
          cvnt = sqrt(sqrt(r0 * rbv(i,j,1)) / nu(i,j))

          if (qc(i,j,1) > thresq) then
            diaqc(i,j,1) = exp(oned3 * log(ccrw6 * qc(i,j,1) / ncc(i,j,1)))
          else
            diaqc(i,j,1) = 0.e0
          end if

          if (qr(i,j,1) > thresq) then
            diaqr(i,j,1) = exp(oned3 * log(cdiaqr * qr(i,j,1) / ncr(i,j,1)))
            vntr(i,j,1) = rbr(i,j,1) * ncr(i,j,1) * (.78e0 * diaqr(i,j,1) &
                 + cvnt * cvntr * exp(pdiaqr * log(diaqr(i,j,1))))
          else
            diaqr(i,j,1) = 0.e0
            vntr(i,j,1) = 0.e0
          end if

          if (qi(i,j,1) > thresq) then
            mi(i,j,1) = qi(i,j,1) / nci(i,j,1)
            diaqi(i,j,1) = exp(oned3 * log(ccri6 * mi(i,j,1)))
          else
            mi(i,j,1) = 0.e0
            diaqi(i,j,1) = 0.e0
          end if

          if (qs(i,j,1) > thresq) then
            diaqs(i,j,1) = exp(oned3 * log(cdiaqs * qs(i,j,1) / ncs(i,j,1)))
            vnts(i,j,1) = rbr(i,j,1) * ncs(i,j,1) * (.78e0 * diaqs(i,j,1) &
                 + cvnt * cvnts * exp(pdiaqs * log(diaqs(i,j,1))))
          else
            diaqs(i,j,1) = 0.e0
            vnts(i,j,1) = 0.e0
          end if

          if (qg(i,j,1) > thresq) then
            diaqg(i,j,1) = exp(oned3 * log(cdiaqg * qg(i,j,1) / ncg(i,j,1)))
            vntg(i,j,1) = rbr(i,j,1) * ncg(i,j,1) * (.78e0 * diaqg(i,j,1) &
                 + cvnt * cvntg * exp(pdiaqg * log(diaqg(i,j,1))))
          else
            diaqg(i,j,1) = 0.e0
            vntg(i,j,1) = 0.e0
          end if
        end do
      end do
      !$acc end kernels

    else
      ! nk > 1 case
      ! For GPU: parallelize all three dimensions (k, j, i)
      ! Since nu(i,j) is updated at each k and used in the same k iteration,
      ! we need separate kernel blocks for the two loop sections

      ! First loop section: calculate thermodynamic properties and update nu(i,j)
      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(esw, esi, p20dvp)
          do i = 1, ni-1
            t(i,j,k) = (ptbr(i,j,k) + ptp(i,j,k)) * pi(i,j,k)
            tcel(i,j,k) = t(i,j,k) - t0

            esw = es0 * exp(17.269e0 * tcel(i,j,k) / (t(i,j,k) - 35.86e0))
            esi = es0 * exp(21.875e0 * tcel(i,j,k) / (t(i,j,k) - 7.66e0))

            qvsst0(i,j,k) = qv(i,j,k) - epses0 / (p(i,j,k) - es0)
            qvsw(i,j,k) = epsva * esw / (p(i,j,k) - esw)
            qvsi(i,j,k) = epsva * esi / (p(i,j,k) - esi)

            lv(i,j,k) = lv0 * exp((.167e0 + 3.67e-4 * t(i,j,k)) * log(t0 / t(i,j,k)))
            lf(i,j,k) = lf0 + cwmci * tcel(i,j,k)
            ls(i,j,k) = lv(i,j,k) + lf(i,j,k)

            p20dvp = p20 / p(i,j,k)
            kp(i,j,k) = (kp0 / (120.e0 + t(i,j,k))) * exp(1.5e0 * log(t23iv * t(i,j,k)))
            dv(i,j,k) = p20dvp * exp(1.81e0 * log(cdv * t(i,j,k)))
          end do
        end do
      end do
      !$acc end kernels

      ! Compute nu(i,j) and mu(i,j,k) for the last k level (nk-1) only
      ! This is because nu is 2D and gets overwritten at each k in the original code
      ! The final nu value is from k=nk-1
      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent private(p20dvp)
        do i = 1, ni-1
          p20dvp = p20 / p(i,j,nk-1)
          nu(i,j) = p20dvp * exp(1.754e0 * log(cnu * t(i,j,nk-1)))
        end do
      end do
      !$acc end kernels

      ! Compute mu for all k levels using the final nu value
      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(p20dvp)
          do i = 1, ni-1
            ! Recompute nu for this k level to get correct mu
            p20dvp = p20 / p(i,j,k)
            mu(i,j,k) = rbr(i,j,k) * (p20dvp * exp(1.754e0 * log(cnu * t(i,j,k))))
          end do
        end do
      end do
      !$acc end kernels

      ! Second loop section: calculate hydrometeor parameters
      ! For cvnt calculation, we need nu at each (i,j) for the corresponding k
      ! In the original code, nu(i,j) is computed at each k and used immediately
      ! For GPU, we compute nu inline for each k level
      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(cvnt, p20dvp)
          do i = 1, ni-1
            ! Recompute nu for this k level (same as in original sequential code)
            p20dvp = p20 / p(i,j,k)
            cvnt = sqrt(sqrt(r0 * rbv(i,j,k)) / (p20dvp * exp(1.754e0 * log(cnu * t(i,j,k)))))

            if (qc(i,j,k) > thresq) then
              diaqc(i,j,k) = exp(oned3 * log(ccrw6 * qc(i,j,k) / ncc(i,j,k)))
            else
              diaqc(i,j,k) = 0.e0
            end if

            if (qr(i,j,k) > thresq) then
              diaqr(i,j,k) = exp(oned3 * log(cdiaqr * qr(i,j,k) / ncr(i,j,k)))
              vntr(i,j,k) = rbr(i,j,k) * ncr(i,j,k) * (.78e0 * diaqr(i,j,k) &
                   + cvnt * cvntr * exp(pdiaqr * log(diaqr(i,j,k))))
            else
              diaqr(i,j,k) = 0.e0
              vntr(i,j,k) = 0.e0
            end if

            if (qi(i,j,k) > thresq) then
              mi(i,j,k) = qi(i,j,k) / nci(i,j,k)
              diaqi(i,j,k) = exp(oned3 * log(ccri6 * mi(i,j,k)))
            else
              mi(i,j,k) = 0.e0
              diaqi(i,j,k) = 0.e0
            end if

            if (qs(i,j,k) > thresq) then
              diaqs(i,j,k) = exp(oned3 * log(cdiaqs * qs(i,j,k) / ncs(i,j,k)))
              vnts(i,j,k) = rbr(i,j,k) * ncs(i,j,k) * (.78e0 * diaqs(i,j,k) &
                   + cvnt * cvnts * exp(pdiaqs * log(diaqs(i,j,k))))
            else
              diaqs(i,j,k) = 0.e0
              vnts(i,j,k) = 0.e0
            end if

            if (qg(i,j,k) > thresq) then
              diaqg(i,j,k) = exp(oned3 * log(cdiaqg * qg(i,j,k) / ncg(i,j,k)))
              vntg(i,j,k) = rbr(i,j,k) * ncg(i,j,k) * (.78e0 * diaqg(i,j,k) &
                   + cvnt * cvntg * exp(pdiaqg * log(diaqg(i,j,k))))
            else
              diaqg(i,j,k) = 0.e0
              vntg(i,j,k) = 0.e0
            end if
          end do
        end do
      end do
      !$acc end kernels
    end if

  end subroutine kernel_setblk

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
        case ('ni'); read(line(eq_pos+1:), *) ni
        case ('nj'); read(line(eq_pos+1:), *) nj
        case ('nk'); read(line(eq_pos+1:), *) nk
        case ('thresq'); read(line(eq_pos+1:), *) thresq
        case ('epsva'); read(line(eq_pos+1:), *) epsva
        case ('es0'); read(line(eq_pos+1:), *) es0
        case ('lv0'); read(line(eq_pos+1:), *) lv0
        case ('lf0'); read(line(eq_pos+1:), *) lf0
        case ('r0'); read(line(eq_pos+1:), *) r0
        case ('ccri6'); read(line(eq_pos+1:), *) ccri6
        case ('ccrw6'); read(line(eq_pos+1:), *) ccrw6
        case ('cdiaqc'); read(line(eq_pos+1:), *) cdiaqc
        case ('cdiaqg'); read(line(eq_pos+1:), *) cdiaqg
        case ('cdiaqr'); read(line(eq_pos+1:), *) cdiaqr
        case ('cdiaqs'); read(line(eq_pos+1:), *) cdiaqs
        case ('cdv'); read(line(eq_pos+1:), *) cdv
        case ('cnu'); read(line(eq_pos+1:), *) cnu
        case ('cvntg'); read(line(eq_pos+1:), *) cvntg
        case ('cvntr'); read(line(eq_pos+1:), *) cvntr
        case ('cvnts'); read(line(eq_pos+1:), *) cvnts
        case ('cwmci'); read(line(eq_pos+1:), *) cwmci
        case ('epses0'); read(line(eq_pos+1:), *) epses0
        case ('pdiaqg'); read(line(eq_pos+1:), *) pdiaqg
        case ('pdiaqr'); read(line(eq_pos+1:), *) pdiaqr
        case ('pdiaqs'); read(line(eq_pos+1:), *) pdiaqs
        case ('t23iv'); read(line(eq_pos+1:), *) t23iv
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

  subroutine read_array_2d(filename, arr, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: arr(is:ie, js:je)
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
  end subroutine read_array_2d

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

  subroutine validate_output_2d(output, reference, tol, max_err, err_count, is, ie, js, je)
    implicit none
    integer, intent(in) :: is, ie, js, je
    real, intent(in) :: output(0:,0:)
    real, intent(in) :: reference(0:,0:)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err, abs_err
    integer :: i, j

    max_err = 0.0
    err_count = 0

    do j = js, je
      do i = is, ie
        abs_err = abs(output(i,j) - reference(i,j))
        if (abs(reference(i,j)) > 1.0e-30) then
          rel_err = abs_err / abs(reference(i,j))
        else
          rel_err = abs_err
        end if
        if (rel_err > max_err) max_err = rel_err
        if (rel_err > tol) err_count = err_count + 1
      end do
    end do
  end subroutine validate_output_2d

end program kernel_benchmark_gpu_setblk

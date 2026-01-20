!***********************************************************************
! Kernel Benchmark: depsit
!***********************************************************************
! Source: Src/depsit.f90
! Subroutine: s_depsit
! Profile ID: 078
! Description: Calculate evaporation rate from rain to vapor and
!              deposition rates from vapor to ice hydrometeors
!***********************************************************************
program kernel_benchmark_depsit
  use omp_lib
  implicit none

  ! === Grid dimensions ===
  integer :: ni, nj, nk

  ! === Scalar parameters ===
  real :: dtb, thresq, rv
  real :: ccdtb4, t27311

  ! === Physical constants ===
  real, parameter :: t0 = 273.16e0      ! Melting point temperature [K]
  real, parameter :: tlow = 233.16e0    ! Lowest super-cooled temperature [K]
  real, parameter :: cc = 3.141592e0    ! Pi
  real, parameter :: ln10 = 2.302585e0  ! log(10)

  ! === Koenig lookup tables ===
  real :: ckoe(0:40), pkoe(0:40)

  ! === Arrays (allocatable) ===
  real, allocatable :: rbr(:,:,:), rbv(:,:,:)
  real, allocatable :: qv(:,:,:), qr(:,:,:), qi(:,:,:), qs(:,:,:), qg(:,:,:)
  real, allocatable :: nci(:,:,:), t(:,:,:), tcel(:,:,:)
  real, allocatable :: qvsst0(:,:,:), qvsw(:,:,:), qvsi(:,:,:)
  real, allocatable :: lv(:,:,:), ls(:,:,:), lf(:,:,:)
  real, allocatable :: kp(:,:,:), dv(:,:,:), mi(:,:,:)
  real, allocatable :: vntr(:,:,:), vnts(:,:,:), vntg(:,:,:)
  real, allocatable :: clcs(:,:,:), clcg(:,:,:)
  real, allocatable :: mlsr(:,:,:), mlgr(:,:,:)

  ! === Output arrays ===
  real, allocatable :: vdvr(:,:,:), vdvi(:,:,:), vdvs(:,:,:), vdvg(:,:,:)
  real, allocatable :: vdvr_ref(:,:,:), vdvi_ref(:,:,:), vdvs_ref(:,:,:), vdvg_ref(:,:,:)

  ! === Benchmark control ===
  integer :: num_iterations, warmup_iterations, iter
  real :: tolerance
  character(len=512) :: data_dir

  ! === Timing ===
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max

  ! === Validation ===
  logical :: passed
  real :: max_err_vdvr, max_err_vdvi, max_err_vdvs, max_err_vdvg
  integer :: err_count_vdvr, err_count_vdvi, err_count_vdvs, err_count_vdvg

  ! === Main program flow ===

  ! 1. Read configuration
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! 2. Initialize Koenig tables
  call init_koenig_tables(ckoe, pkoe)

  ! 3. Read parameters
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, dtb, thresq, rv, ccdtb4, t27311)

  ! 4. Allocate arrays
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(qg(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci(0:ni+1, 0:nj+1, 1:nk))
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(tcel(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsst0(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsw(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsi(0:ni+1, 0:nj+1, 1:nk))
  allocate(lv(0:ni+1, 0:nj+1, 1:nk))
  allocate(ls(0:ni+1, 0:nj+1, 1:nk))
  allocate(lf(0:ni+1, 0:nj+1, 1:nk))
  allocate(kp(0:ni+1, 0:nj+1, 1:nk))
  allocate(dv(0:ni+1, 0:nj+1, 1:nk))
  allocate(mi(0:ni+1, 0:nj+1, 1:nk))
  allocate(vntr(0:ni+1, 0:nj+1, 1:nk))
  allocate(vnts(0:ni+1, 0:nj+1, 1:nk))
  allocate(vntg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcg(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlsr(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlgr(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvr(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvi(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvs(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvg(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvi_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvs_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvg_ref(0:ni+1, 0:nj+1, 1:nk))

  ! 5. Read input data
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rbv.bin', rbv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qr.bin', qr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qi.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qg.bin', qg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/nci.bin', nci, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tcel.bin', tcel, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvsst0.bin', qvsst0, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvsw.bin', qvsw, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvsi.bin', qvsi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/lv.bin', lv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ls.bin', ls, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/lf.bin', lf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/kp.bin', kp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/dv.bin', dv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/mi.bin', mi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vntr.bin', vntr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vnts.bin', vnts, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vntg.bin', vntg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clcs.bin', clcs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clcg.bin', clcg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/mlsr.bin', mlsr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/mlgr.bin', mlgr, 0, ni+1, 0, nj+1, 1, nk)

  ! 6. Read reference output
  call read_array_3d(trim(data_dir)//'/vdvr_ref.bin', vdvr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vdvi_ref.bin', vdvi_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vdvs_ref.bin', vdvs_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vdvg_ref.bin', vdvg_ref, 0, ni+1, 0, nj+1, 1, nk)

  ! Print header
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: depsit'
  write(*,'(A)') '=================================================='
  write(*,'(A,I8,A,I8,A,I8)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I8)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  ! 7. Warmup iterations
  do iter = 1, warmup_iterations
    call kernel_depsit(ni, nj, nk, dtb, thresq, rv, ccdtb4, t27311, &
         rbr, rbv, qv, qr, qi, qs, qg, nci, t, tcel, &
         qvsst0, qvsw, qvsi, lv, ls, lf, kp, dv, mi, &
         vntr, vnts, vntg, clcs, clcg, mlsr, mlgr, &
         vdvr, vdvi, vdvs, vdvg, ckoe, pkoe)
  end do

  ! 8. Benchmark iterations (timed)
  t_total = 0.0d0
  t_min = huge(1.0d0)
  t_max = 0.0d0

  do iter = 1, num_iterations
    t_start = omp_get_wtime()
    call kernel_depsit(ni, nj, nk, dtb, thresq, rv, ccdtb4, t27311, &
         rbr, rbv, qv, qr, qi, qs, qg, nci, t, tcel, &
         qvsst0, qvsw, qvsi, lv, ls, lf, kp, dv, mi, &
         vntr, vnts, vntg, clcs, clcg, mlsr, mlgr, &
         vdvr, vdvi, vdvs, vdvg, ckoe, pkoe)
    t_end = omp_get_wtime()
    t_total = t_total + (t_end - t_start)
    if (t_end - t_start < t_min) t_min = t_end - t_start
    if (t_end - t_start > t_max) t_max = t_end - t_start
  end do
  t_avg = t_total / dble(num_iterations)

  ! 9. Validate output
  call validate_output_3d(vdvr, vdvr_ref, tolerance, max_err_vdvr, err_count_vdvr, &
       0, ni+1, 0, nj+1, 1, nk)
  call validate_output_3d(vdvi, vdvi_ref, tolerance, max_err_vdvi, err_count_vdvi, &
       0, ni+1, 0, nj+1, 1, nk)
  call validate_output_3d(vdvs, vdvs_ref, tolerance, max_err_vdvs, err_count_vdvs, &
       0, ni+1, 0, nj+1, 1, nk)
  call validate_output_3d(vdvg, vdvg_ref, tolerance, max_err_vdvg, err_count_vdvg, &
       0, ni+1, 0, nj+1, 1, nk)

  passed = (err_count_vdvr == 0 .and. err_count_vdvi == 0 .and. &
            err_count_vdvs == 0 .and. err_count_vdvg == 0)

  ! 10. Report results
  write(*,'(A,F12.3,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:     ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:     ', t_max * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max error (vdvr): ', max_err_vdvr
  write(*,'(A,ES12.4)') ' Max error (vdvi): ', max_err_vdvi
  write(*,'(A,ES12.4)') ' Max error (vdvs): ', max_err_vdvs
  write(*,'(A,ES12.4)') ' Max error (vdvg): ', max_err_vdvg
  write(*,'(A,ES12.4)') ' Tolerance:        ', tolerance
  write(*,'(A,I12)') ' Error count (vdvr): ', err_count_vdvr
  write(*,'(A,I12)') ' Error count (vdvi): ', err_count_vdvi
  write(*,'(A,I12)') ' Error count (vdvs): ', err_count_vdvs
  write(*,'(A,I12)') ' Error count (vdvg): ', err_count_vdvg
  if (passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! 11. Cleanup
  deallocate(rbr, rbv, qv, qr, qi, qs, qg, nci, t, tcel)
  deallocate(qvsst0, qvsw, qvsi, lv, ls, lf, kp, dv, mi)
  deallocate(vntr, vnts, vntg, clcs, clcg, mlsr, mlgr)
  deallocate(vdvr, vdvi, vdvs, vdvg)
  deallocate(vdvr_ref, vdvi_ref, vdvs_ref, vdvg_ref)

contains

  !-----------------------------------------------------------------
  ! The kernel (extracted from original source)
  !-----------------------------------------------------------------
  subroutine kernel_depsit(ni, nj, nk, dtb, thresq, rv, ccdtb4, t27311, &
       rbr, rbv, qv, qr, qi, qs, qg, nci, t, tcel, &
       qvsst0, qvsw, qvsi, lv, ls, lf, kp, dv, mi, &
       vntr, vnts, vntg, clcs, clcg, mlsr, mlgr, &
       vdvr, vdvi, vdvs, vdvg, ckoe, pkoe)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: dtb, thresq, rv, ccdtb4, t27311
    real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: nci(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: tcel(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qvsst0(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qvsw(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qvsi(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: lv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ls(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: lf(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: kp(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: dv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: mi(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: vntr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: vnts(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: vntg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: clcs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: clcg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: mlsr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: mlgr(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: vdvr(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: vdvi(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: vdvs(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: vdvg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ckoe(0:40), pkoe(0:40)

    ! Local variables
    integer :: i, j, k, itsc
    real :: qvssw, qvssi, gi, cvdvx1, cvdvx2, cvdvx3, cvdvx4, sink, a, b, c

    !$omp parallel default(shared) private(k)

    ! Case nk = 1
    if (nk == 1) then

      !$omp do schedule(runtime) private(i,j,itsc) &
      !$omp&   private(qvssw,qvssi,gi,cvdvx1,cvdvx2,cvdvx3,cvdvx4,sink,a,b,c)
      do j = 1, nj-1
        do i = 1, ni-1

          if (t(i,j,1) > tlow) then

            qvssw = qv(i,j,1) - qvsw(i,j,1)
            qvssi = qv(i,j,1) - qvsi(i,j,1)

            a = 1.e0 / (rv * kp(i,j,1) * t(i,j,1) * t(i,j,1))
            b = dv(i,j,1) * rbr(i,j,1)
            c = ccdtb4 * rbv(i,j,1)

            gi = 1.e0 / (ls(i,j,1)*ls(i,j,1)*a + 1.e0/(b*qvsi(i,j,1)))

            cvdvx1 = c * (qv(i,j,1)/qvsw(i,j,1) - 1.e0) &
                 / (lv(i,j,1)*lv(i,j,1)*a + 1.e0/(b*qvsw(i,j,1)))

            cvdvx2 = ccdtb4 * dv(i,j,1) * qvsst0(i,j,1)
            cvdvx3 = c * (qv(i,j,1)/qvsi(i,j,1) - 1.e0)
            cvdvx4 = a * ls(i,j,1) * lf(i,j,1)

            ! Evaporation from rain
            if (qr(i,j,1) > thresq) then
              if (qv(i,j,1) < qvsw(i,j,1)) then
                vdvr(i,j,1) = max(cvdvx1*vntr(i,j,1), qvssw)
              else
                vdvr(i,j,1) = 0.e0
              end if
            else
              vdvr(i,j,1) = 0.e0
            end if

            ! Deposition to cloud ice
            if (qi(i,j,1) > thresq) then
              if (t(i,j,1) < t27311) then
                itsc = int(-tcel(i,j,1))
                vdvi(i,j,1) = ckoe(itsc) * (qv(i,j,1)-qvsi(i,j,1)) &
                     / (qvsw(i,j,1)-qvsi(i,j,1)) * nci(i,j,1) &
                     * exp(pkoe(itsc)*log(mi(i,j,1))) * dtb
              else
                vdvi(i,j,1) = 0.e0
              end if
            else
              vdvi(i,j,1) = 0.e0
            end if

            ! Deposition to snow
            if (qs(i,j,1) > thresq) then
              if (t(i,j,1) > t0) then
                if (mlsr(i,j,1) > 0.e0) then
                  vdvs(i,j,1) = cvdvx2 * vnts(i,j,1)
                else
                  vdvs(i,j,1) = cvdvx1 * vnts(i,j,1)
                end if
              else
                vdvs(i,j,1) = gi * (cvdvx3*vnts(i,j,1) - cvdvx4*clcs(i,j,1))
              end if
            else
              vdvs(i,j,1) = 0.e0
            end if

            ! Deposition to graupel
            if (qg(i,j,1) > thresq) then
              if (t(i,j,1) > t0) then
                if (mlgr(i,j,1) > 0.e0) then
                  vdvg(i,j,1) = cvdvx2 * vntg(i,j,1)
                else
                  vdvg(i,j,1) = cvdvx1 * vntg(i,j,1)
                end if
              else
                vdvg(i,j,1) = gi * (cvdvx3*vntg(i,j,1) - cvdvx4*clcg(i,j,1))
              end if
            else
              vdvg(i,j,1) = 0.e0
            end if

            ! Adjust deposition rate
            sink = vdvi(i,j,1) + vdvs(i,j,1) + vdvg(i,j,1)
            if ((qvssi < sink .and. qvssi > 0.e0) .or. &
                (qvssi > sink .and. qvssi < 0.e0)) then
              a = qvssi / sink
              vdvi(i,j,1) = vdvi(i,j,1) * a
              vdvs(i,j,1) = vdvs(i,j,1) * a
              vdvg(i,j,1) = vdvg(i,j,1) * a
            end if

          else
            vdvr(i,j,1) = 0.e0
            vdvi(i,j,1) = 0.e0
            vdvs(i,j,1) = 0.e0
            vdvg(i,j,1) = 0.e0
          end if

        end do
      end do
      !$omp end do

    else
      ! Case nk > 1
      do k = 1, nk-1

        !$omp do schedule(runtime) private(i,j,itsc) &
        !$omp&   private(qvssw,qvssi,gi,cvdvx1,cvdvx2,cvdvx3,cvdvx4,sink,a,b,c)
        do j = 1, nj-1
          do i = 1, ni-1

            if (t(i,j,k) > tlow) then

              qvssw = qv(i,j,k) - qvsw(i,j,k)
              qvssi = qv(i,j,k) - qvsi(i,j,k)

              a = 1.e0 / (rv * kp(i,j,k) * t(i,j,k) * t(i,j,k))
              b = dv(i,j,k) * rbr(i,j,k)
              c = ccdtb4 * rbv(i,j,k)

              gi = 1.e0 / (ls(i,j,k)*ls(i,j,k)*a + 1.e0/(b*qvsi(i,j,k)))

              cvdvx1 = c * (qv(i,j,k)/qvsw(i,j,k) - 1.e0) &
                   / (lv(i,j,k)*lv(i,j,k)*a + 1.e0/(b*qvsw(i,j,k)))

              cvdvx2 = ccdtb4 * dv(i,j,k) * qvsst0(i,j,k)
              cvdvx3 = c * (qv(i,j,k)/qvsi(i,j,k) - 1.e0)
              cvdvx4 = a * ls(i,j,k) * lf(i,j,k)

              ! Evaporation from rain
              if (qr(i,j,k) > thresq) then
                if (qv(i,j,k) < qvsw(i,j,k)) then
                  vdvr(i,j,k) = max(cvdvx1*vntr(i,j,k), qvssw)
                else
                  vdvr(i,j,k) = 0.e0
                end if
              else
                vdvr(i,j,k) = 0.e0
              end if

              ! Deposition to cloud ice
              if (qi(i,j,k) > thresq) then
                if (t(i,j,k) < t27311) then
                  itsc = int(-tcel(i,j,k))
                  vdvi(i,j,k) = ckoe(itsc) * (qv(i,j,k)-qvsi(i,j,k)) &
                       / (qvsw(i,j,k)-qvsi(i,j,k)) * nci(i,j,k) &
                       * exp(pkoe(itsc)*log(mi(i,j,k))) * dtb
                else
                  vdvi(i,j,k) = 0.e0
                end if
              else
                vdvi(i,j,k) = 0.e0
              end if

              ! Deposition to snow
              if (qs(i,j,k) > thresq) then
                if (t(i,j,k) > t0) then
                  if (mlsr(i,j,k) > 0.e0) then
                    vdvs(i,j,k) = cvdvx2 * vnts(i,j,k)
                  else
                    vdvs(i,j,k) = cvdvx1 * vnts(i,j,k)
                  end if
                else
                  vdvs(i,j,k) = gi * (cvdvx3*vnts(i,j,k) - cvdvx4*clcs(i,j,k))
                end if
              else
                vdvs(i,j,k) = 0.e0
              end if

              ! Deposition to graupel
              if (qg(i,j,k) > thresq) then
                if (t(i,j,k) > t0) then
                  if (mlgr(i,j,k) > 0.e0) then
                    vdvg(i,j,k) = cvdvx2 * vntg(i,j,k)
                  else
                    vdvg(i,j,k) = cvdvx1 * vntg(i,j,k)
                  end if
                else
                  vdvg(i,j,k) = gi * (cvdvx3*vntg(i,j,k) - cvdvx4*clcg(i,j,k))
                end if
              else
                vdvg(i,j,k) = 0.e0
              end if

              ! Adjust deposition rate
              sink = vdvi(i,j,k) + vdvs(i,j,k) + vdvg(i,j,k)
              if ((qvssi < sink .and. qvssi > 0.e0) .or. &
                  (qvssi > sink .and. qvssi < 0.e0)) then
                a = qvssi / sink
                vdvi(i,j,k) = vdvi(i,j,k) * a
                vdvs(i,j,k) = vdvs(i,j,k) * a
                vdvg(i,j,k) = vdvg(i,j,k) * a
              end if

            else
              vdvr(i,j,k) = 0.e0
              vdvi(i,j,k) = 0.e0
              vdvs(i,j,k) = 0.e0
              vdvg(i,j,k) = 0.e0
            end if

          end do
        end do
        !$omp end do

      end do
    end if

    !$omp end parallel

  end subroutine kernel_depsit

  !-----------------------------------------------------------------
  ! Initialize Koenig tables
  !-----------------------------------------------------------------
  subroutine init_koenig_tables(ckoe, pkoe)
    implicit none
    real, intent(out) :: ckoe(0:40), pkoe(0:40)

    real, parameter :: ln10 = 2.302585e0
    integer :: i

    ! Initial values from setref.f90
    pkoe(0)=0.e0
    pkoe(1)=.4006e0
    pkoe(2)=.4831e0
    pkoe(3)=.5320e0
    pkoe(4)=.5307e0
    pkoe(5)=.5319e0
    pkoe(6)=.5249e0
    pkoe(7)=.4888e0
    pkoe(8)=.3894e0
    pkoe(9)=.4047e0
    pkoe(10)=.4318e0
    pkoe(11)=.4771e0
    pkoe(12)=.5183e0
    pkoe(13)=.5463e0
    pkoe(14)=.5651e0
    pkoe(15)=.5813e0
    pkoe(16)=.5655e0
    pkoe(17)=.5478e0
    pkoe(18)=.5203e0
    pkoe(19)=.4906e0
    pkoe(20)=.4447e0
    pkoe(21)=.4126e0
    pkoe(22)=.3960e0
    pkoe(23)=.4149e0
    pkoe(24)=.4320e0
    pkoe(25)=.4506e0
    pkoe(26)=.4483e0
    pkoe(27)=.4460e0
    pkoe(28)=.4433e0
    pkoe(29)=.4413e0
    pkoe(30)=.4382e0
    pkoe(31)=.4361e0
    pkoe(32)=.4340e0
    pkoe(33)=.4319e0
    pkoe(34)=.4298e0
    pkoe(35)=.4277e0
    pkoe(36)=.4256e0
    pkoe(37)=.4235e0
    pkoe(38)=.4214e0
    pkoe(39)=.4193e0
    pkoe(40)=.4172e0

    ckoe(0)=0.e0
    ckoe(1)=.7939e-10
    ckoe(2)=.7841e-9
    ckoe(3)=.3369e-8
    ckoe(4)=.4336e-8
    ckoe(5)=.5285e-8
    ckoe(6)=.3728e-8
    ckoe(7)=.1852e-8
    ckoe(8)=.2991e-9
    ckoe(9)=.4248e-9
    ckoe(10)=.7434e-9
    ckoe(11)=.1812e-8
    ckoe(12)=.4394e-8
    ckoe(13)=.9145e-8
    ckoe(14)=.1725e-7
    ckoe(15)=.3348e-7
    ckoe(16)=.1725e-7
    ckoe(17)=.9175e-8
    ckoe(18)=.4412e-8
    ckoe(19)=.2252e-8
    ckoe(20)=.9115e-9
    ckoe(21)=.4876e-9
    ckoe(22)=.3473e-9
    ckoe(23)=.4758e-9
    ckoe(24)=.6306e-9
    ckoe(25)=.8573e-9
    ckoe(26)=.7868e-9
    ckoe(27)=.7192e-9
    ckoe(28)=.6513e-9
    ckoe(29)=.5956e-9
    ckoe(30)=.5333e-9
    ckoe(31)=.4834e-9
    ckoe(32)=.4335e-9
    ckoe(33)=.3836e-9
    ckoe(34)=.3337e-9
    ckoe(35)=.2838e-9
    ckoe(36)=.2339e-9
    ckoe(37)=.1840e-9
    ckoe(38)=.1341e-9
    ckoe(39)=.8420e-10
    ckoe(40)=.3430e-10

    ! Apply transformation: ckoe(i) = ckoe(i) * exp(3*pkoe(i)*ln10)
    do i = 0, 40
      ckoe(i) = ckoe(i) * exp(3.e0 * pkoe(i) * ln10)
    end do

  end subroutine init_koenig_tables

  !-----------------------------------------------------------------
  ! I/O utilities
  !-----------------------------------------------------------------
  subroutine read_config(data_dir, num_iterations, warmup_iterations, tolerance)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iterations, warmup_iterations
    real, intent(out) :: tolerance

    character(len=512) :: line
    integer :: ios
    logical :: exists

    ! Default values
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

  subroutine read_parameters(filename, ni, nj, nk, dtb, thresq, rv, ccdtb4, t27311)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    real, intent(out) :: dtb, thresq, rv, ccdtb4, t27311

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
        case ('ni')
          read(line(eq_pos+1:), *) ni
        case ('nj')
          read(line(eq_pos+1:), *) nj
        case ('nk')
          read(line(eq_pos+1:), *) nk
        case ('dtb')
          read(line(eq_pos+1:), *) dtb
        case ('thresq')
          read(line(eq_pos+1:), *) thresq
        case ('rv')
          read(line(eq_pos+1:), *) rv
        case ('ccdtb4')
          read(line(eq_pos+1:), *) ccdtb4
        case ('t27311')
          read(line(eq_pos+1:), *) t27311
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

    open(unit=20, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
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

  subroutine validate_output_3d(output, reference, tol, max_err, err_count, &
       is, ie, js, je, ks, ke)
    implicit none
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(in) :: output(is:ie, js:je, ks:ke)
    real, intent(in) :: reference(is:ie, js:je, ks:ke)
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

end program kernel_benchmark_depsit

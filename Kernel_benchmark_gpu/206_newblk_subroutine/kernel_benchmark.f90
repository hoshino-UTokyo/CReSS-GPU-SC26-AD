!***********************************************************************
! GPU Kernel Benchmark: newblk (s_newblk)
!***********************************************************************
!
! Source: Src/newblk.f90
! Description: Solves microphysics budget equations for potential temperature,
!              mixing ratios (qv, qc, qr, qi, qs, qg), and concentrations.
!              This benchmark is for cphopt=3, nk=1 case.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_newblk
  use omp_lib
  implicit none

  ! Parameters
  integer :: cphopt, ni, nj, nk
  real :: thresq, cp_const, mi0iv, mr0iv, ms0iv

  ! Input arrays (read-only)
  real, allocatable :: pi(:,:,:)
  real, allocatable :: qcp(:,:,:), qrp(:,:,:), qip(:,:,:), qsp(:,:,:), qgp(:,:,:)
  real, allocatable :: nccp(:,:,:), ncsp(:,:,:), ncgp(:,:,:)
  real, allocatable :: lv(:,:,:), ls(:,:,:), lf(:,:,:), mi(:,:,:)
  real, allocatable :: nuvi(:,:,:), nuci(:,:,:)
  real, allocatable :: clcr(:,:,:), clcs(:,:,:), clcg(:,:,:)
  real, allocatable :: clri(:,:,:), clrs(:,:,:), clrg(:,:,:)
  real, allocatable :: clir(:,:,:), clis(:,:,:), clig(:,:,:)
  real, allocatable :: clsr(:,:,:), clsg(:,:,:), clrsg(:,:,:)
  real, allocatable :: clrin(:,:,:), clrsn(:,:,:), clsrn(:,:,:), clsgn(:,:,:)
  real, allocatable :: agcn(:,:,:), agrn(:,:,:), agin(:,:,:), agsn(:,:,:)
  real, allocatable :: vdvr(:,:,:), vdvi(:,:,:), vdvs(:,:,:), vdvg(:,:,:)
  real, allocatable :: cncr(:,:,:), cnis(:,:,:), cnsg(:,:,:), cnsgn(:,:,:)
  real, allocatable :: spsi(:,:,:), spgi(:,:,:)
  real, allocatable :: mlic(:,:,:), mlsr(:,:,:), mlgr(:,:,:)
  real, allocatable :: frrg(:,:,:), frrgn(:,:,:)
  real, allocatable :: shsr(:,:,:), shgr(:,:,:)

  ! Input/output arrays
  real, allocatable :: ptpf(:,:,:), qvf(:,:,:), qcf(:,:,:), qrf(:,:,:)
  real, allocatable :: qif(:,:,:), qsf(:,:,:), qgf(:,:,:)
  real, allocatable :: ncif(:,:,:), ncsf(:,:,:), ncgf(:,:,:)

  ! Backup and reference arrays
  real, allocatable :: ptpf_in(:,:,:), qvf_in(:,:,:), qcf_in(:,:,:), qrf_in(:,:,:)
  real, allocatable :: qif_in(:,:,:), qsf_in(:,:,:), qgf_in(:,:,:)
  real, allocatable :: ncif_in(:,:,:), ncsf_in(:,:,:), ncgf_in(:,:,:)
  real, allocatable :: ptpf_ref(:,:,:), qvf_ref(:,:,:), qcf_ref(:,:,:), qrf_ref(:,:,:)
  real, allocatable :: qif_ref(:,:,:), qsf_ref(:,:,:), qgf_ref(:,:,:)
  real, allocatable :: ncif_ref(:,:,:), ncsf_ref(:,:,:), ncgf_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir
  real :: tolerance
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  integer :: error_count, total_errors
  logical :: validation_passed

  integer :: iter

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt')

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: newblk'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' cphopt=', cphopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(pi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qcp(0:ni+1, 0:nj+1, 1:nk), qrp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qip(0:ni+1, 0:nj+1, 1:nk), qsp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qgp(0:ni+1, 0:nj+1, 1:nk))
  allocate(nccp(0:ni+1, 0:nj+1, 1:nk), ncsp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncgp(0:ni+1, 0:nj+1, 1:nk))
  allocate(lv(0:ni+1, 0:nj+1, 1:nk), ls(0:ni+1, 0:nj+1, 1:nk))
  allocate(lf(0:ni+1, 0:nj+1, 1:nk), mi(0:ni+1, 0:nj+1, 1:nk))
  allocate(nuvi(0:ni+1, 0:nj+1, 1:nk), nuci(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcr(0:ni+1, 0:nj+1, 1:nk), clcs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcg(0:ni+1, 0:nj+1, 1:nk), clri(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrs(0:ni+1, 0:nj+1, 1:nk), clrg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clir(0:ni+1, 0:nj+1, 1:nk), clis(0:ni+1, 0:nj+1, 1:nk))
  allocate(clig(0:ni+1, 0:nj+1, 1:nk), clsr(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsg(0:ni+1, 0:nj+1, 1:nk), clrsg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrin(0:ni+1, 0:nj+1, 1:nk), clrsn(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsrn(0:ni+1, 0:nj+1, 1:nk), clsgn(0:ni+1, 0:nj+1, 1:nk))
  allocate(agcn(0:ni+1, 0:nj+1, 1:nk), agrn(0:ni+1, 0:nj+1, 1:nk))
  allocate(agin(0:ni+1, 0:nj+1, 1:nk), agsn(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvr(0:ni+1, 0:nj+1, 1:nk), vdvi(0:ni+1, 0:nj+1, 1:nk))
  allocate(vdvs(0:ni+1, 0:nj+1, 1:nk), vdvg(0:ni+1, 0:nj+1, 1:nk))
  allocate(cncr(0:ni+1, 0:nj+1, 1:nk), cnis(0:ni+1, 0:nj+1, 1:nk))
  allocate(cnsg(0:ni+1, 0:nj+1, 1:nk), cnsgn(0:ni+1, 0:nj+1, 1:nk))
  allocate(spsi(0:ni+1, 0:nj+1, 1:nk), spgi(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlic(0:ni+1, 0:nj+1, 1:nk), mlsr(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlgr(0:ni+1, 0:nj+1, 1:nk), frrg(0:ni+1, 0:nj+1, 1:nk))
  allocate(frrgn(0:ni+1, 0:nj+1, 1:nk), shsr(0:ni+1, 0:nj+1, 1:nk))
  allocate(shgr(0:ni+1, 0:nj+1, 1:nk))

  allocate(ptpf(0:ni+1, 0:nj+1, 1:nk), qvf(0:ni+1, 0:nj+1, 1:nk))
  allocate(qcf(0:ni+1, 0:nj+1, 1:nk), qrf(0:ni+1, 0:nj+1, 1:nk))
  allocate(qif(0:ni+1, 0:nj+1, 1:nk), qsf(0:ni+1, 0:nj+1, 1:nk))
  allocate(qgf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncif(0:ni+1, 0:nj+1, 1:nk), ncsf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncgf(0:ni+1, 0:nj+1, 1:nk))

  allocate(ptpf_in(0:ni+1, 0:nj+1, 1:nk), qvf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(qcf_in(0:ni+1, 0:nj+1, 1:nk), qrf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(qif_in(0:ni+1, 0:nj+1, 1:nk), qsf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(qgf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncif_in(0:ni+1, 0:nj+1, 1:nk), ncsf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncgf_in(0:ni+1, 0:nj+1, 1:nk))

  allocate(ptpf_ref(0:ni+1, 0:nj+1, 1:nk), qvf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qcf_ref(0:ni+1, 0:nj+1, 1:nk), qrf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qif_ref(0:ni+1, 0:nj+1, 1:nk), qsf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qgf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncif_ref(0:ni+1, 0:nj+1, 1:nk), ncsf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncgf_ref(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/pi.bin', pi)
  call read_array_3d(trim(data_dir)//'/qcp.bin', qcp)
  call read_array_3d(trim(data_dir)//'/qrp.bin', qrp)
  call read_array_3d(trim(data_dir)//'/qip.bin', qip)
  call read_array_3d(trim(data_dir)//'/qsp.bin', qsp)
  call read_array_3d(trim(data_dir)//'/qgp.bin', qgp)
  call read_array_3d(trim(data_dir)//'/nccp.bin', nccp)
  call read_array_3d(trim(data_dir)//'/ncsp.bin', ncsp)
  call read_array_3d(trim(data_dir)//'/ncgp.bin', ncgp)
  call read_array_3d(trim(data_dir)//'/lv.bin', lv)
  call read_array_3d(trim(data_dir)//'/ls.bin', ls)
  call read_array_3d(trim(data_dir)//'/lf.bin', lf)
  call read_array_3d(trim(data_dir)//'/mi.bin', mi)
  call read_array_3d(trim(data_dir)//'/nuvi.bin', nuvi)
  call read_array_3d(trim(data_dir)//'/nuci.bin', nuci)
  call read_array_3d(trim(data_dir)//'/clcr.bin', clcr)
  call read_array_3d(trim(data_dir)//'/clcs.bin', clcs)
  call read_array_3d(trim(data_dir)//'/clcg.bin', clcg)
  call read_array_3d(trim(data_dir)//'/clri.bin', clri)
  call read_array_3d(trim(data_dir)//'/clrs.bin', clrs)
  call read_array_3d(trim(data_dir)//'/clrg.bin', clrg)
  call read_array_3d(trim(data_dir)//'/clir.bin', clir)
  call read_array_3d(trim(data_dir)//'/clis.bin', clis)
  call read_array_3d(trim(data_dir)//'/clig.bin', clig)
  call read_array_3d(trim(data_dir)//'/clsr.bin', clsr)
  call read_array_3d(trim(data_dir)//'/clsg.bin', clsg)
  call read_array_3d(trim(data_dir)//'/clrsg.bin', clrsg)
  call read_array_3d(trim(data_dir)//'/clrin.bin', clrin)
  call read_array_3d(trim(data_dir)//'/clrsn.bin', clrsn)
  call read_array_3d(trim(data_dir)//'/clsrn.bin', clsrn)
  call read_array_3d(trim(data_dir)//'/clsgn.bin', clsgn)
  call read_array_3d(trim(data_dir)//'/agcn.bin', agcn)
  call read_array_3d(trim(data_dir)//'/agrn.bin', agrn)
  call read_array_3d(trim(data_dir)//'/agin.bin', agin)
  call read_array_3d(trim(data_dir)//'/agsn.bin', agsn)
  call read_array_3d(trim(data_dir)//'/vdvr.bin', vdvr)
  call read_array_3d(trim(data_dir)//'/vdvi.bin', vdvi)
  call read_array_3d(trim(data_dir)//'/vdvs.bin', vdvs)
  call read_array_3d(trim(data_dir)//'/vdvg.bin', vdvg)
  call read_array_3d(trim(data_dir)//'/cncr.bin', cncr)
  call read_array_3d(trim(data_dir)//'/cnis.bin', cnis)
  call read_array_3d(trim(data_dir)//'/cnsg.bin', cnsg)
  call read_array_3d(trim(data_dir)//'/cnsgn.bin', cnsgn)
  call read_array_3d(trim(data_dir)//'/spsi.bin', spsi)
  call read_array_3d(trim(data_dir)//'/spgi.bin', spgi)
  call read_array_3d(trim(data_dir)//'/mlic.bin', mlic)
  call read_array_3d(trim(data_dir)//'/mlsr.bin', mlsr)
  call read_array_3d(trim(data_dir)//'/mlgr.bin', mlgr)
  call read_array_3d(trim(data_dir)//'/frrg.bin', frrg)
  call read_array_3d(trim(data_dir)//'/frrgn.bin', frrgn)
  call read_array_3d(trim(data_dir)//'/shsr.bin', shsr)
  call read_array_3d(trim(data_dir)//'/shgr.bin', shgr)

  call read_array_3d(trim(data_dir)//'/ptpf_in.bin', ptpf_in)
  call read_array_3d(trim(data_dir)//'/qvf_in.bin', qvf_in)
  call read_array_3d(trim(data_dir)//'/qcf_in.bin', qcf_in)
  call read_array_3d(trim(data_dir)//'/qrf_in.bin', qrf_in)
  call read_array_3d(trim(data_dir)//'/qif_in.bin', qif_in)
  call read_array_3d(trim(data_dir)//'/qsf_in.bin', qsf_in)
  call read_array_3d(trim(data_dir)//'/qgf_in.bin', qgf_in)
  call read_array_3d(trim(data_dir)//'/ncif_in.bin', ncif_in)
  call read_array_3d(trim(data_dir)//'/ncsf_in.bin', ncsf_in)
  call read_array_3d(trim(data_dir)//'/ncgf_in.bin', ncgf_in)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ptpf_ref.bin', ptpf_ref)
  call read_array_3d(trim(data_dir)//'/qvf_ref.bin', qvf_ref)
  call read_array_3d(trim(data_dir)//'/qcf_ref.bin', qcf_ref)
  call read_array_3d(trim(data_dir)//'/qrf_ref.bin', qrf_ref)
  call read_array_3d(trim(data_dir)//'/qif_ref.bin', qif_ref)
  call read_array_3d(trim(data_dir)//'/qsf_ref.bin', qsf_ref)
  call read_array_3d(trim(data_dir)//'/qgf_ref.bin', qgf_ref)
  call read_array_3d(trim(data_dir)//'/ncif_ref.bin', ncif_ref)
  call read_array_3d(trim(data_dir)//'/ncsf_ref.bin', ncsf_ref)
  call read_array_3d(trim(data_dir)//'/ncgf_ref.bin', ncgf_ref)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ptpf = ptpf_in; qvf = qvf_in; qcf = qcf_in; qrf = qrf_in
    qif = qif_in; qsf = qsf_in; qgf = qgf_in
    ncif = ncif_in; ncsf = ncsf_in; ncgf = ncgf_in
    call kernel_newblk()
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ptpf = ptpf_in; qvf = qvf_in; qcf = qcf_in; qrf = qrf_in
    qif = qif_in; qsf = qsf_in; qgf = qgf_in
    ncif = ncif_in; ncsf = ncsf_in; ncgf = ncgf_in

    !$acc wait
    t_start = omp_get_wtime()

    call kernel_newblk()

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
  total_errors = 0

  call validate_3d(ptpf, ptpf_ref, 'ptpf', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(qvf, qvf_ref, 'qvf', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(qcf, qcf_ref, 'qcf', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(qrf, qrf_ref, 'qrf', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(qif, qif_ref, 'qif', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(qsf, qsf_ref, 'qsf', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(qgf, qgf_ref, 'qgf', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(ncif, ncif_ref, 'ncif', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(ncsf, ncsf_ref, 'ncsf', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(ncgf, ncgf_ref, 'ncgf', tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
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
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', total_errors
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: newblk for cphopt=3, nk=1 (OpenACC version)
  !=====================================================================
  subroutine kernel_newblk()
    implicit none
    integer :: i, j
    real :: cppiv, clix, spxi, cxcr, cxsg
    real :: nmci, mssr, msgr, nuvdvx, clnmci, cfmsxr
    real :: qvsink, qxsink, dqv

    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent private(cppiv,clix,spxi,cxcr,cxsg,nmci,mssr,msgr,nuvdvx,clnmci,cfmsxr,qvsink,qxsink,dqv)
      do i = 1, ni-1
        cppiv = 1.0e0 / (cp_const * pi(i,j,1))
        clix = clir(i,j,1) + clis(i,j,1) + clig(i,j,1)
        spxi = spsi(i,j,1) + spgi(i,j,1)
        cxcr = clcr(i,j,1) + cncr(i,j,1)
        cxsg = clsr(i,j,1) + clsg(i,j,1) + cnsg(i,j,1)
        nmci = nuci(i,j,1) - mlic(i,j,1)
        mssr = mlsr(i,j,1) + shsr(i,j,1)
        msgr = mlgr(i,j,1) + shgr(i,j,1)
        nuvdvx = vdvi(i,j,1) + vdvs(i,j,1) + vdvg(i,j,1) + nuvi(i,j,1)
        clnmci = clcs(i,j,1) + clcg(i,j,1) + nmci
        cfmsxr = clri(i,j,1) + clrs(i,j,1) + clrg(i,j,1) + clrsg(i,j,1) &
               + frrg(i,j,1) - mssr - msgr

        ! New potential temperature
        ptpf(i,j,1) = ptpf(i,j,1) + (vdvr(i,j,1)*lv(i,j,1) &
                    + (clnmci + cfmsxr)*lf(i,j,1) + nuvdvx*ls(i,j,1)) * cppiv

        ! Mixing ratios
        qvsink = nuvdvx + vdvr(i,j,1)
        qxsink = cfmsxr - cxcr - vdvr(i,j,1)

        if (qrf(i,j,1) < qxsink) then
          qxsink = cxcr + clnmci + (qxsink - qrf(i,j,1))
          qrf(i,j,1) = 0.0e0
        else
          qrf(i,j,1) = qrf(i,j,1) - qxsink
          qxsink = cxcr + clnmci
        end if

        if (qcf(i,j,1) < qxsink) then
          dqv = qxsink - qcf(i,j,1)
          ptpf(i,j,1) = ptpf(i,j,1) + cppiv * dqv * lv(i,j,1)
          qvsink = qvsink + dqv
          qcf(i,j,1) = 0.0e0
        else
          qcf(i,j,1) = qcf(i,j,1) - qxsink
        end if

        qxsink = clix - spxi - nmci + cnis(i,j,1) - vdvi(i,j,1) - nuvi(i,j,1)

        if (qif(i,j,1) < qxsink) then
          dqv = qxsink - qif(i,j,1)
          ptpf(i,j,1) = ptpf(i,j,1) + cppiv * dqv * ls(i,j,1)
          qvsink = qvsink + dqv
          qif(i,j,1) = 0.0e0
        else
          qif(i,j,1) = qif(i,j,1) - qxsink
        end if

        qvf(i,j,1) = max(qvf(i,j,1) - qvsink, 0.0e0)

        ! Snow and graupel
        qsf(i,j,1) = max(qsf(i,j,1) - (cxsg + mssr + spsi(i,j,1) - vdvs(i,j,1) &
                   - clcs(i,j,1) - clrs(i,j,1) - clis(i,j,1) - cnis(i,j,1)), 0.0e0)

        qgf(i,j,1) = max(qgf(i,j,1) + (cxsg - msgr - spgi(i,j,1) + vdvg(i,j,1) &
                   + clri(i,j,1) + clir(i,j,1) + clcg(i,j,1) &
                   + clrg(i,j,1) + clig(i,j,1) + clrsg(i,j,1) + frrg(i,j,1)), 0.0e0)

        ! Cloud ice concentrations
        if (qcp(i,j,1) > thresq) then
          ncif(i,j,1) = ncif(i,j,1) + nuci(i,j,1) * nccp(i,j,1) / qcp(i,j,1)
        end if

        if (qip(i,j,1) > thresq) then
          if (vdvi(i,j,1) < 0.0e0) then
            ncif(i,j,1) = ncif(i,j,1) - (agin(i,j,1) &
                        + (clix + mlic(i,j,1) - vdvi(i,j,1)) / mi(i,j,1) &
                        - (spxi + nuvi(i,j,1)) * mi0iv + cnis(i,j,1) * ms0iv)
          else
            ncif(i,j,1) = ncif(i,j,1) - (agin(i,j,1) &
                        + (clix + mlic(i,j,1)) / mi(i,j,1) &
                        - (spxi + nuvi(i,j,1)) * mi0iv + cnis(i,j,1) * ms0iv)
          end if
        else
          ncif(i,j,1) = ncif(i,j,1) + (spxi + nuvi(i,j,1)) * mi0iv
        end if

        ! Snow concentrations
        if (qsp(i,j,1) > thresq) then
          if (vdvs(i,j,1) < 0.0e0) then
            ncsf(i,j,1) = ncsf(i,j,1) - (agsn(i,j,1) - cnis(i,j,1) * ms0iv &
                        + clsrn(i,j,1) + clsgn(i,j,1) + cnsgn(i,j,1) &
                        + (mlsr(i,j,1) - vdvs(i,j,1)) * ncsp(i,j,1) / qsp(i,j,1))
          else
            ncsf(i,j,1) = ncsf(i,j,1) - (agsn(i,j,1) - cnis(i,j,1) * ms0iv &
                        + clsrn(i,j,1) + clsgn(i,j,1) + cnsgn(i,j,1) &
                        + mlsr(i,j,1) * ncsp(i,j,1) / qsp(i,j,1))
          end if
        else
          ncsf(i,j,1) = ncsf(i,j,1) + cnis(i,j,1) * ms0iv
        end if

        ! Graupel concentrations
        if (qgp(i,j,1) > thresq) then
          if (vdvg(i,j,1) < 0.0e0) then
            ncgf(i,j,1) = ncgf(i,j,1) + ((vdvg(i,j,1) - mlgr(i,j,1)) * ncgp(i,j,1) / qgp(i,j,1) &
                        + clrsn(i,j,1) + clrin(i,j,1) + cnsgn(i,j,1) + frrgn(i,j,1))
          else
            ncgf(i,j,1) = ncgf(i,j,1) - (mlgr(i,j,1) * ncgp(i,j,1) / qgp(i,j,1) &
                        - clrsn(i,j,1) - clrin(i,j,1) - cnsgn(i,j,1) - frrgn(i,j,1))
          end if
        else
          ncgf(i,j,1) = ncgf(i,j,1) + (clrsn(i,j,1) + clrin(i,j,1) + cnsgn(i,j,1) + frrgn(i,j,1))
        end if
      end do
    end do
    !$acc end kernels

  end subroutine kernel_newblk

  !=====================================================================
  ! Validation
  !=====================================================================
  subroutine validate_3d(arr, ref, name, tol, max_err, err_count)
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:nk)
    character(len=*), intent(in) :: name
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = 1, nk
      do j = 1, nj-1
        do i = 1, ni-1
          rel_err = abs(arr(i,j,k) - ref(i,j,k))
          if (abs(ref(i,j,k)) > 1.0e-10) then
            rel_err = rel_err / abs(ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

  end subroutine validate_3d

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
  subroutine read_parameters(filename)
    character(len=*), intent(in) :: filename

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    cphopt = 3
    cp_const = 1004.0

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
          case ('cphopt'); read(val, *) cphopt
          case ('ni');     read(val, *) ni
          case ('nj');     read(val, *) nj
          case ('nk');     read(val, *) nk
          case ('thresq'); read(val, *) thresq
          case ('cp');     read(val, *) cp_const
          case ('mi0iv');  read(val, *) mi0iv
          case ('mr0iv');  read(val, *) mr0iv
          case ('ms0iv');  read(val, *) ms0iv
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary array reader
  !=====================================================================
  subroutine read_array_3d(filename, arr)
    character(len=*), intent(in) :: filename
    real, intent(out) :: arr(0:ni+1, 0:nj+1, 1:nk)

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

end program kernel_benchmark_gpu_newblk

!***********************************************************************
      module m_newblk
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2001/05/29, 2001/10/18,
!                   2001/11/20, 2002/01/15, 2002/03/29, 2002/04/02,
!                   2002/12/02, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/06/10, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2005/01/07, 2005/04/04, 2005/10/05,
!                   2006/02/13, 2006/04/03, 2006/09/30, 2007/10/19,
!                   2007/11/26, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2011/03/11, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the new potential temperature perturbation, the mixing ratio
!     and concentrations.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: newblk, s_newblk

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface newblk

        module procedure s_newblk

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_newblk(cphopt,thresq,ni,nj,nk,pi,qcp,qrp,qip,qsp,qgp,&
     &                    nccp,ncsp,ncgp,lv,ls,lf,mi,nuvi,nuci,         &
     &                    clcr,clcs,clcg,clri,clrs,clrg,clir,clis,clig, &
     &                    clsr,clsg,clrsg,clrin,clrsn,clsrn,clsgn,      &
     &                    agcn,agrn,agin,agsn,vdvr,vdvi,vdvs,vdvg,      &
     &                    cncr,cnis,cnsg,cnsgn,spsi,spgi,mlic,mlsr,mlgr,&
     &                    frrg,frrgn,shsr,shgr,ptpf,qvf,qcf,qrf,        &
     &                    qif,qsf,qgf,nccf,ncrf,ncif,ncsf,ncgf)
!***********************************************************************

! Input variables

      integer, intent(in) :: cphopt
                       ! Option for cloud micro physics

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exner function

      real, intent(in) :: qcp(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at past

      real, intent(in) :: qrp(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at past

      real, intent(in) :: qip(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio at past

      real, intent(in) :: qsp(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio at past

      real, intent(in) :: qgp(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio at past

      real, intent(in) :: nccp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water at past

      real, intent(in) :: ncsp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow at past

      real, intent(in) :: ncgp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel at past

      real, intent(in) :: lv(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of evapolation

      real, intent(in) :: ls(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of sublimation

      real, intent(in) :: lf(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of fusion

      real, intent(in) :: mi(0:ni+1,0:nj+1,1:nk)
                       ! Mean mass of cloud ice

      real, intent(in) :: nuvi(0:ni+1,0:nj+1,1:nk)
                       ! Nucleation rate of deposition or sorption

      real, intent(in) :: nuci(0:ni+1,0:nj+1,1:nk)
                       ! Nucleation rate
                       ! of condensation, contact and homogeneous

      real, intent(in) :: clcr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and rain water

      real, intent(in) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and snow

      real, intent(in) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and graupel

      real, intent(in) :: clri(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, intent(in) :: clrs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! from rain water to snow

      real, intent(in) :: clrg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and graupel

      real, intent(in) :: clir(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, intent(in) :: clis(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud ice and snow

      real, intent(in) :: clig(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud ice and graupel

      real, intent(in) :: clsr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! from snow to rain water

      real, intent(in) :: clsg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between snow and graupel

      real, intent(in) :: clrsg(0:ni+1,0:nj+1,1:nk)
                       ! Production rate of graupel
                       ! from collection rate form rain to snow

      real, intent(in) :: clrin(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and cloud ice

      real, intent(in) :: clrsn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, intent(in) :: clsrn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, intent(in) :: clsgn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between snow and graupel

      real, intent(in) :: agcn(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for cloud water

      real, intent(in) :: agrn(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for rain water

      real, intent(in) :: agin(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for cloud ice

      real, intent(in) :: agsn(0:ni+1,0:nj+1,1:nk)
                       ! Aggregation rate for snow

      real, intent(in) :: vdvr(0:ni+1,0:nj+1,1:nk)
                       ! Evaporation rate from rain water to water vapor

      real, intent(in) :: vdvi(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to cloud ice

      real, intent(in) :: vdvs(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to snow

      real, intent(in) :: vdvg(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to graupel

      real, intent(in) :: cncr(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from cloud water to rain water

      real, intent(in) :: cnis(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from cloud ice to snow

      real, intent(in) :: cnsg(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from snow to graupel

      real, intent(in) :: cnsgn(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate for concentrations
                       ! from snow to graupel

      real, intent(in) :: spsi(0:ni+1,0:nj+1,1:nk)
                       ! Secondary nucleation rate from snow

      real, intent(in) :: spgi(0:ni+1,0:nj+1,1:nk)
                       ! Secondary nucleation rate from graupel

      real, intent(in) :: mlic(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from cloud ice to cloud water

      real, intent(in) :: mlsr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from snow to rain water

      real, intent(in) :: mlgr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from graupel to rain water

      real, intent(in) :: frrg(0:ni+1,0:nj+1,1:nk)
                       ! Freezing rate from rain water to graupel

      real, intent(in) :: frrgn(0:ni+1,0:nj+1,1:nk)
                       ! Freezing rate for concentrations
                       ! from rain water to graupel

      real, intent(in) :: shsr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from snow

      real, intent(in) :: shgr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from graupel

! Input and output variables

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qcf(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at future

      real, intent(inout) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at future

      real, intent(inout) :: qif(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio at future

      real, intent(inout) :: qsf(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio at future

      real, intent(inout) :: qgf(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio at future

      real, intent(inout) :: nccf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water at future

      real, intent(inout) :: ncrf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water at future

      real, intent(inout) :: ncif(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice at future

      real, intent(inout) :: ncsf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow at future

      real, intent(inout) :: ncgf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel at future

! Internal shared variables

      real mr0iv       ! Inverse of mr0
      real mi0iv       ! Inverse of mi0
      real ms0iv       ! Inverse of ms0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real mciv        ! Inverse of mean mass of cloud water
      real miiv        ! Inverse of mean mass of cloud ice
      real msiv        ! Inverse of mean mass of snow
      real mgiv        ! Inverse of mean mass of graupel

      real cppiv       ! 1.0 / (cp x pi)

      real clix        ! clir + clis + clig

      real spxi        ! spsi + spgi

      real cxcr        ! clcr + cncr

      real cxsg        ! clsr + clsg + cnsg

      real nmci        ! nuci - mlic

      real mssr        ! mlsr + shsr
      real msgr        ! mlgr + shgr

      real nuvdvx      ! vdvi + vdvs + vdvg + nuvi

      real clnmci      ! clcs + clcg + nmci

      real cfmsxr      ! clri + clrs + clrg + clrsg + frrg - mssr - msgr

      real qvsink      ! Sink amount of water vapor mixing ratio
      real qxsink      ! Sink amount of optional mixing ratio

      real dqv         ! Over estimated sink amount
                       ! of cloud water or cloud ice


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_newblk = 0
      integer, parameter :: DUMP_TARGET_newblk = 45720
      logical, save :: dump_done_newblk = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variables.

      mr0iv=1.e0/mr0
      mi0iv=1.e0/mi0
      ms0iv=1.e0/ms0

! -----

!!!! Solve the new potential temperature perturbation, the mixing ratio
!!!! and concentrations.

!@llm start meta_info ----------------------------------------------------
! Location: newblk.f90 :: subroutine s_newblk
! Summary : Solves microphysics budget equations for potential temperature,
!           mixing ratios (qv, qc, qr, qi, qs, qg), and concentrations.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constants (cp, mr0, mi0, ms0) from comphy.
!   - No synchronization constructs.
!   - Complex conditionals on cphopt (cloud physics option: 2, 3, or 4).
!   - Separate code paths for nk=1 (2D) vs nk>1 (3D).
!   - Many private variables for microphysics rate calculations.
!   - Contains threshold checks (qxp > thresq) with conditional updates.
!   - All grid points are independent (embarrassingly parallel).
!   - Uses intrinsic max() and abs() - GPU compatible.
! Next:
!   - Select code path based on cphopt outside kernel.
!   - OpenACC kernels with collapse(2) or collapse(3) for 3D case.
!   - Large number of input arrays - ensure efficient data movement.
!   - Consider kernel fusion for related calculations.
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 17.986s (0.60%)
!   - AvgTime: 0.393ms
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('newblk.f90', 's_newblk', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_newblk = dump_call_count_newblk + 1
if (dump_call_count_newblk == DUMP_TARGET_newblk .and. .not. dump_done_newblk) then
  call dump_init('newblk')
  ! FIXME: cphopt is array - call dump_scalar_i('cphopt', cphopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('thresq', thresq)
  ! FIXME: cp is array - call dump_scalar_r('cp', cp)
  call dump_array_3d('pi.bin', pi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qcp.bin', qcp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qrp.bin', qrp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qip.bin', qip, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qsp.bin', qsp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qgp.bin', qgp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nccp.bin', nccp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncsp.bin', ncsp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncgp.bin', ncgp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('lv.bin', lv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ls.bin', ls, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('lf.bin', lf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mi.bin', mi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nuvi.bin', nuvi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nuci.bin', nuci, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcr.bin', clcr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcs.bin', clcs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcg.bin', clcg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clri.bin', clri, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrs.bin', clrs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrg.bin', clrg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clir.bin', clir, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clis.bin', clis, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clig.bin', clig, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsr.bin', clsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsg.bin', clsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrsg.bin', clrsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrin.bin', clrin, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrsn.bin', clrsn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsrn.bin', clsrn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsgn.bin', clsgn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('agcn.bin', agcn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('agrn.bin', agrn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('agin.bin', agin, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('agsn.bin', agsn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvr.bin', vdvr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvi.bin', vdvi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvs.bin', vdvs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvg.bin', vdvg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cncr.bin', cncr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cnis.bin', cnis, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cnsg.bin', cnsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cnsgn.bin', cnsgn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('spsi.bin', spsi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('spgi.bin', spgi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mlic.bin', mlic, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mlsr.bin', mlsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mlgr.bin', mlgr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('frrg.bin', frrg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('frrgn.bin', frrgn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('shsr.bin', shsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('shgr.bin', shgr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptpf_in.bin', ptpf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvf_in.bin', qvf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qcf_in.bin', qcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qrf_in.bin', qrf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qif_in.bin', qif, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qsf_in.bin', qsf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qgf_in.bin', qgf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nccf_in.bin', nccf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncrf_in.bin', ncrf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncif_in.bin', ncif, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncsf_in.bin', ncsf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncgf_in.bin', ncgf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('mi0iv', mi0iv)
  call dump_scalar_r('mr0iv', mr0iv)
  call dump_scalar_r('ms0iv', ms0iv)
end if

!$omp parallel default(shared) private(k)

!!! In the case nk = 1.

      if(nk.eq.1) then

!! Solve the new potential temperature perturbation, the mixing ratio
!! and cloud ice concentrations.

        if(abs(cphopt).eq.2) then

!$omp do schedule(runtime) private(i,j,cppiv,clix,spxi,cxcr,cxsg)       &
!$omp&   private(nmci,mssr,msgr,nuvdvx,clnmci,cfmsxr,qvsink,qxsink,dqv)

          do j=1,nj-1
          do i=1,ni-1

! Set the common used variables.

            cppiv=1.e0/(cp*pi(i,j,1))

            clix=clir(i,j,1)+clis(i,j,1)+clig(i,j,1)

            spxi=spsi(i,j,1)+spgi(i,j,1)

            cxcr=clcr(i,j,1)+cncr(i,j,1)

            cxsg=clsr(i,j,1)+clsg(i,j,1)+cnsg(i,j,1)

            nmci=nuci(i,j,1)-mlic(i,j,1)

            mssr=mlsr(i,j,1)+shsr(i,j,1)
            msgr=mlgr(i,j,1)+shgr(i,j,1)

            nuvdvx=vdvi(i,j,1)+vdvs(i,j,1)+vdvg(i,j,1)+nuvi(i,j,1)

            clnmci=clcs(i,j,1)+clcg(i,j,1)+nmci

            cfmsxr=clri(i,j,1)+clrs(i,j,1)+clrg(i,j,1)+clrsg(i,j,1)     &
     &        +frrg(i,j,1)-mssr-msgr

! -----

! Get the new potential temperature perturbation.

            ptpf(i,j,1)=ptpf(i,j,1)+(vdvr(i,j,1)*lv(i,j,1)              &
     &        +(clnmci+cfmsxr)*lf(i,j,1)+nuvdvx*ls(i,j,1))*cppiv

! -----

! Get the new water vapor, the cloud water, the rain water and the cloud
! ice mixing ratio.

            qvsink=nuvdvx+vdvr(i,j,1)

            qxsink=cfmsxr-cxcr-vdvr(i,j,1)

            if(qrf(i,j,1).lt.qxsink) then

              qxsink=cxcr+clnmci+(qxsink-qrf(i,j,1))

              qrf(i,j,1)=0.e0

            else

              qrf(i,j,1)=qrf(i,j,1)-qxsink

              qxsink=cxcr+clnmci

            end if

            if(qcf(i,j,1).lt.qxsink) then

              dqv=qxsink-qcf(i,j,1)

              ptpf(i,j,1)=ptpf(i,j,1)+cppiv*dqv*lv(i,j,1)

              qvsink=qvsink+dqv

              qcf(i,j,1)=0.e0

            else

              qcf(i,j,1)=qcf(i,j,1)-qxsink

            end if

            qxsink=clix-spxi-nmci+cnis(i,j,1)-vdvi(i,j,1)-nuvi(i,j,1)

            if(qif(i,j,1).lt.qxsink) then

              dqv=qxsink-qif(i,j,1)

              ptpf(i,j,1)=ptpf(i,j,1)+cppiv*dqv*ls(i,j,1)

              qvsink=qvsink+dqv

              qif(i,j,1)=0.e0

            else

              qif(i,j,1)=qif(i,j,1)-qxsink

            end if

            qvf(i,j,1)=max(qvf(i,j,1)-qvsink,0.e0)

! -----

! Get the new snow and the graupel mixing ratio.

            qsf(i,j,1)=max(qsf(i,j,1)                                   &
     &        -(cxsg+mssr+spsi(i,j,1)-vdvs(i,j,1)                       &
     &        -clcs(i,j,1)-clrs(i,j,1)-clis(i,j,1)-cnis(i,j,1)),0.e0)

            qgf(i,j,1)=max(qgf(i,j,1)                                   &
     &        +(cxsg-msgr-spgi(i,j,1)+vdvg(i,j,1)                       &
     &        +clri(i,j,1)+clir(i,j,1)+clcg(i,j,1)                      &
     &        +clrg(i,j,1)+clig(i,j,1)+clrsg(i,j,1)+frrg(i,j,1)),0.e0)

! -----

! Get the new cloud ice concentrations.

            if(qcp(i,j,1).gt.thresq) then

              ncif(i,j,1)=ncif(i,j,1)+nuci(i,j,1)*nccp(i,j,1)/qcp(i,j,1)

            end if

            if(qip(i,j,1).gt.thresq) then

              if(vdvi(i,j,1).lt.0.e0) then

                ncif(i,j,1)=ncif(i,j,1)-(agin(i,j,1)                    &
     &            +(clix+mlic(i,j,1)-vdvi(i,j,1))/mi(i,j,1)             &
     &            -(spxi+nuvi(i,j,1))*mi0iv+cnis(i,j,1)*ms0iv)

              else

                ncif(i,j,1)=ncif(i,j,1)                                 &
     &            -(agin(i,j,1)+(clix+mlic(i,j,1))/mi(i,j,1)            &
     &            -(spxi+nuvi(i,j,1))*mi0iv+cnis(i,j,1)*ms0iv)

              end if

            else

              ncif(i,j,1)=ncif(i,j,1)+(spxi+nuvi(i,j,1))*mi0iv

            end if

! -----

          end do
          end do

!$omp end do

!! -----

!! Solve the new potential temperature perturbation, the mixing ratio
!! and ice concentrations.

        else if(abs(cphopt).eq.3) then

!$omp do schedule(runtime) private(i,j,cppiv,clix,spxi,cxcr,cxsg)       &
!$omp&   private(nmci,mssr,msgr,nuvdvx,clnmci,cfmsxr,qvsink,qxsink,dqv)

          do j=1,nj-1
          do i=1,ni-1

! Set the common used variables.

            cppiv=1.e0/(cp*pi(i,j,1))

            clix=clir(i,j,1)+clis(i,j,1)+clig(i,j,1)

            spxi=spsi(i,j,1)+spgi(i,j,1)

            cxcr=clcr(i,j,1)+cncr(i,j,1)

            cxsg=clsr(i,j,1)+clsg(i,j,1)+cnsg(i,j,1)

            nmci=nuci(i,j,1)-mlic(i,j,1)

            mssr=mlsr(i,j,1)+shsr(i,j,1)
            msgr=mlgr(i,j,1)+shgr(i,j,1)

            nuvdvx=vdvi(i,j,1)+vdvs(i,j,1)+vdvg(i,j,1)+nuvi(i,j,1)

            clnmci=clcs(i,j,1)+clcg(i,j,1)+nmci

            cfmsxr=clri(i,j,1)+clrs(i,j,1)+clrg(i,j,1)+clrsg(i,j,1)     &
     &        +frrg(i,j,1)-mssr-msgr

! -----

! Get the new potential temperature perturbation.

            ptpf(i,j,1)=ptpf(i,j,1)+(vdvr(i,j,1)*lv(i,j,1)              &
     &        +(clnmci+cfmsxr)*lf(i,j,1)+nuvdvx*ls(i,j,1))*cppiv

! -----

! Get the new water vapor, the cloud water, the rain water and the cloud
! ice mixing ratio.

            qvsink=nuvdvx+vdvr(i,j,1)

            qxsink=cfmsxr-cxcr-vdvr(i,j,1)

            if(qrf(i,j,1).lt.qxsink) then

              qxsink=cxcr+clnmci+(qxsink-qrf(i,j,1))

              qrf(i,j,1)=0.e0

            else

              qrf(i,j,1)=qrf(i,j,1)-qxsink

              qxsink=cxcr+clnmci

            end if

            if(qcf(i,j,1).lt.qxsink) then

              dqv=qxsink-qcf(i,j,1)

              ptpf(i,j,1)=ptpf(i,j,1)+cppiv*dqv*lv(i,j,1)

              qvsink=qvsink+dqv

              qcf(i,j,1)=0.e0

            else

              qcf(i,j,1)=qcf(i,j,1)-qxsink

            end if

            qxsink=clix-spxi-nmci+cnis(i,j,1)-vdvi(i,j,1)-nuvi(i,j,1)

            if(qif(i,j,1).lt.qxsink) then

              dqv=qxsink-qif(i,j,1)

              ptpf(i,j,1)=ptpf(i,j,1)+cppiv*dqv*ls(i,j,1)

              qvsink=qvsink+dqv

              qif(i,j,1)=0.e0

            else

              qif(i,j,1)=qif(i,j,1)-qxsink

            end if

            qvf(i,j,1)=max(qvf(i,j,1)-qvsink,0.e0)

! -----

! Get the new snow and the graupel mixing ratio.

            qsf(i,j,1)=max(qsf(i,j,1)                                   &
     &        -(cxsg+mssr+spsi(i,j,1)-vdvs(i,j,1)                       &
     &        -clcs(i,j,1)-clrs(i,j,1)-clis(i,j,1)-cnis(i,j,1)),0.e0)

            qgf(i,j,1)=max(qgf(i,j,1)                                   &
     &        +(cxsg-msgr-spgi(i,j,1)+vdvg(i,j,1)                       &
     &        +clri(i,j,1)+clir(i,j,1)+clcg(i,j,1)                      &
     &        +clrg(i,j,1)+clig(i,j,1)+clrsg(i,j,1)+frrg(i,j,1)),0.e0)

! -----

! Get the new cloud ice concentrations.

            if(qcp(i,j,1).gt.thresq) then

              ncif(i,j,1)=ncif(i,j,1)+nuci(i,j,1)*nccp(i,j,1)/qcp(i,j,1)

            end if

            if(qip(i,j,1).gt.thresq) then

              if(vdvi(i,j,1).lt.0.e0) then

                ncif(i,j,1)=ncif(i,j,1)-(agin(i,j,1)                    &
     &            +(clix+mlic(i,j,1)-vdvi(i,j,1))/mi(i,j,1)             &
     &            -(spxi+nuvi(i,j,1))*mi0iv+cnis(i,j,1)*ms0iv)

              else

                ncif(i,j,1)=ncif(i,j,1)                                 &
     &            -(agin(i,j,1)+(clix+mlic(i,j,1))/mi(i,j,1)            &
     &            -(spxi+nuvi(i,j,1))*mi0iv+cnis(i,j,1)*ms0iv)

              end if

            else

              ncif(i,j,1)=ncif(i,j,1)+(spxi+nuvi(i,j,1))*mi0iv

            end if

! -----

! Get the new snow concentrations.

            if(qsp(i,j,1).gt.thresq) then

              if(vdvs(i,j,1).lt.0.e0) then

                ncsf(i,j,1)=ncsf(i,j,1)                                 &
     &            -(agsn(i,j,1)-cnis(i,j,1)*ms0iv                       &
     &            +clsrn(i,j,1)+clsgn(i,j,1)+cnsgn(i,j,1)               &
     &            +(mlsr(i,j,1)-vdvs(i,j,1))*ncsp(i,j,1)/qsp(i,j,1))

              else

                ncsf(i,j,1)=ncsf(i,j,1)-(agsn(i,j,1)-cnis(i,j,1)*ms0iv  &
     &            +clsrn(i,j,1)+clsgn(i,j,1)+cnsgn(i,j,1)               &
     &            +mlsr(i,j,1)*ncsp(i,j,1)/qsp(i,j,1))

              end if

            else

              ncsf(i,j,1)=ncsf(i,j,1)+cnis(i,j,1)*ms0iv

            end if

! -----

! Get the new graupel concentrations.

            if(qgp(i,j,1).gt.thresq) then

              if(vdvg(i,j,1).lt.0.e0) then

                ncgf(i,j,1)=ncgf(i,j,1)                                 &
     &            +((vdvg(i,j,1)-mlgr(i,j,1))*ncgp(i,j,1)/qgp(i,j,1)    &
     &            +clrsn(i,j,1)+clrin(i,j,1)+cnsgn(i,j,1)+frrgn(i,j,1))

              else

                ncgf(i,j,1)=ncgf(i,j,1)                                 &
     &            -(mlgr(i,j,1)*ncgp(i,j,1)/qgp(i,j,1)                  &
     &            -clrsn(i,j,1)-clrin(i,j,1)-cnsgn(i,j,1)-frrgn(i,j,1))

              end if

            else

              ncgf(i,j,1)=ncgf(i,j,1)                                   &
     &          +(clrsn(i,j,1)+clrin(i,j,1)+cnsgn(i,j,1)+frrgn(i,j,1))

            end if

! -----

          end do
          end do

!$omp end do

!! -----

!! Solve the new potential temperature perturbation, the mixing ratio
!! and water and ice concentrations.

        else if(abs(cphopt).eq.4) then

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,mciv,miiv,msiv,mgiv,cppiv,clix,spxi,cxcr,cxsg)     &
!$omp&   private(nmci,mssr,msgr,nuvdvx,clnmci,cfmsxr,qvsink,qxsink,dqv)

          do j=1,nj-1
          do i=1,ni-1

! Set the common used variables.

            cppiv=1.e0/(cp*pi(i,j,1))

            clix=clir(i,j,1)+clis(i,j,1)+clig(i,j,1)

            spxi=spsi(i,j,1)+spgi(i,j,1)

            cxcr=clcr(i,j,1)+cncr(i,j,1)

            cxsg=clsr(i,j,1)+clsg(i,j,1)+cnsg(i,j,1)

            nmci=nuci(i,j,1)-mlic(i,j,1)

            mssr=mlsr(i,j,1)+shsr(i,j,1)
            msgr=mlgr(i,j,1)+shgr(i,j,1)

            nuvdvx=vdvi(i,j,1)+vdvs(i,j,1)+vdvg(i,j,1)+nuvi(i,j,1)

            clnmci=clcs(i,j,1)+clcg(i,j,1)+nmci

            cfmsxr=clri(i,j,1)+clrs(i,j,1)+clrg(i,j,1)+clrsg(i,j,1)     &
     &        +frrg(i,j,1)-mssr-msgr

! -----

! Get the new potential temperature perturbation.

            ptpf(i,j,1)=ptpf(i,j,1)+(vdvr(i,j,1)*lv(i,j,1)              &
     &        +(clnmci+cfmsxr)*lf(i,j,1)+nuvdvx*ls(i,j,1))*cppiv

! -----

! Get the new water vapor, the cloud water, the rain water and the cloud
! ice mixing ratio.

            qvsink=nuvdvx+vdvr(i,j,1)

            qxsink=cfmsxr-cxcr-vdvr(i,j,1)

            if(qrf(i,j,1).lt.qxsink) then

              qxsink=cxcr+clnmci+(qxsink-qrf(i,j,1))

              qrf(i,j,1)=0.e0

            else

              qrf(i,j,1)=qrf(i,j,1)-qxsink

              qxsink=cxcr+clnmci

            end if

            if(qcf(i,j,1).lt.qxsink) then

              dqv=qxsink-qcf(i,j,1)

              ptpf(i,j,1)=ptpf(i,j,1)+cppiv*dqv*lv(i,j,1)

              qvsink=qvsink+dqv

              qcf(i,j,1)=0.e0

            else

              qcf(i,j,1)=qcf(i,j,1)-qxsink

            end if

            qxsink=clix-spxi-nmci+cnis(i,j,1)-vdvi(i,j,1)-nuvi(i,j,1)

            if(qif(i,j,1).lt.qxsink) then

              dqv=qxsink-qif(i,j,1)

              ptpf(i,j,1)=ptpf(i,j,1)+cppiv*dqv*ls(i,j,1)

              qvsink=qvsink+dqv

              qif(i,j,1)=0.e0

            else

              qif(i,j,1)=qif(i,j,1)-qxsink

            end if

            qvf(i,j,1)=max(qvf(i,j,1)-qvsink,0.e0)

! -----

! Get the new snow and the graupel mixing ratio.

            qsf(i,j,1)=max(qsf(i,j,1)                                   &
     &        -(cxsg+mssr+spsi(i,j,1)-vdvs(i,j,1)                       &
     &        -clcs(i,j,1)-clrs(i,j,1)-clis(i,j,1)-cnis(i,j,1)),0.e0)

            qgf(i,j,1)=max(qgf(i,j,1)                                   &
     &        +(cxsg-msgr-spgi(i,j,1)+vdvg(i,j,1)                       &
     &        +clri(i,j,1)+clir(i,j,1)+clcg(i,j,1)                      &
     &        +clrg(i,j,1)+clig(i,j,1)+clrsg(i,j,1)+frrg(i,j,1)),0.e0)

! -----

! Get the new cloud water and cloud ice concentrations.

            if(qcp(i,j,1).gt.thresq) then

              mciv=nccp(i,j,1)/qcp(i,j,1)

              nccf(i,j,1)=nccf(i,j,1)-(agcn(i,j,1)+cncr(i,j,1)*mr0iv    &
     &          +(nuci(i,j,1)+clcr(i,j,1)+clcs(i,j,1)+clcg(i,j,1))*mciv)

              ncif(i,j,1)=ncif(i,j,1)+nuci(i,j,1)*mciv

            end if

            if(qip(i,j,1).gt.thresq) then

              miiv=1.e0/mi(i,j,1)

              nccf(i,j,1)=nccf(i,j,1)+mlic(i,j,1)*miiv

              if(vdvi(i,j,1).lt.0.e0) then

                ncif(i,j,1)=ncif(i,j,1)-(agin(i,j,1)                    &
     &            +(clix+mlic(i,j,1)-vdvi(i,j,1))*miiv                  &
     &            -(spxi+nuvi(i,j,1))*mi0iv+cnis(i,j,1)*ms0iv)

              else

                ncif(i,j,1)=ncif(i,j,1)                                 &
     &            -(agin(i,j,1)+(clix+mlic(i,j,1))*miiv                 &
     &            -(spxi+nuvi(i,j,1))*mi0iv+cnis(i,j,1)*ms0iv)

              end if

            else

              ncif(i,j,1)=ncif(i,j,1)+(spxi+nuvi(i,j,1))*mi0iv

            end if

! -----

! Get the new rain water, snow and graupel concentrations.

            if(qrp(i,j,1).gt.thresq) then

              ncrf(i,j,1)=ncrf(i,j,1)-(clrin(i,j,1)+clrsn(i,j,1)        &
     &          +frrgn(i,j,1)+agrn(i,j,1)-cncr(i,j,1)*mr0iv)

            else

              ncrf(i,j,1)=ncrf(i,j,1)+cncr(i,j,1)*mr0iv

            end if

            if(qsp(i,j,1).gt.thresq) then

              msiv=ncsp(i,j,1)/qsp(i,j,1)

              ncrf(i,j,1)=ncrf(i,j,1)+mlsr(i,j,1)*msiv

              if(vdvs(i,j,1).lt.0.e0) then

                ncsf(i,j,1)=ncsf(i,j,1)                                 &
     &            -(agsn(i,j,1)-cnis(i,j,1)*ms0iv                       &
     &            +(mlsr(i,j,1)-vdvs(i,j,1))*msiv                       &
     &            +clsrn(i,j,1)+clsgn(i,j,1)+cnsgn(i,j,1))

              else

                ncsf(i,j,1)=ncsf(i,j,1)                                 &
     &            -(agsn(i,j,1)-cnis(i,j,1)*ms0iv+mlsr(i,j,1)*msiv      &
     &            +clsrn(i,j,1)+clsgn(i,j,1)+cnsgn(i,j,1))

              end if

            else

              ncsf(i,j,1)=ncsf(i,j,1)+cnis(i,j,1)*ms0iv

            end if

            if(qgp(i,j,1).gt.thresq) then

              mgiv=ncgp(i,j,1)/qgp(i,j,1)

              ncrf(i,j,1)=ncrf(i,j,1)+mlgr(i,j,1)*mgiv

              if(vdvg(i,j,1).lt.0.e0) then

                ncgf(i,j,1)=ncgf(i,j,1)                                 &
     &            +((vdvg(i,j,1)-mlgr(i,j,1))*mgiv                      &
     &            +clrsn(i,j,1)+clrin(i,j,1)+cnsgn(i,j,1)+frrgn(i,j,1))

              else

                ncgf(i,j,1)=ncgf(i,j,1)-(mlgr(i,j,1)*mgiv               &
     &            -clrsn(i,j,1)-clrin(i,j,1)-cnsgn(i,j,1)-frrgn(i,j,1))

              end if

            else

              ncgf(i,j,1)=ncgf(i,j,1)                                   &
     &          +(clrsn(i,j,1)+clrin(i,j,1)+cnsgn(i,j,1)+frrgn(i,j,1))

            end if

! -----

          end do
          end do

!$omp end do

        end if

!! -----

!!! -----

!!! In the case nk > 1.

      else

!! Solve the new potential temperature perturbation, the mixing ratio
!! and cloud ice concentrations.

        if(abs(cphopt).eq.2) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,cppiv,clix,spxi,cxcr,cxsg)       &
!$omp&   private(nmci,mssr,msgr,nuvdvx,clnmci,cfmsxr,qvsink,qxsink,dqv)

            do j=1,nj-1
            do i=1,ni-1

! Set the common used variables.

              cppiv=1.e0/(cp*pi(i,j,k))

              clix=clir(i,j,k)+clis(i,j,k)+clig(i,j,k)

              spxi=spsi(i,j,k)+spgi(i,j,k)

              cxcr=clcr(i,j,k)+cncr(i,j,k)

              cxsg=clsr(i,j,k)+clsg(i,j,k)+cnsg(i,j,k)

              nmci=nuci(i,j,k)-mlic(i,j,k)

              mssr=mlsr(i,j,k)+shsr(i,j,k)
              msgr=mlgr(i,j,k)+shgr(i,j,k)

              nuvdvx=vdvi(i,j,k)+vdvs(i,j,k)+vdvg(i,j,k)+nuvi(i,j,k)

              clnmci=clcs(i,j,k)+clcg(i,j,k)+nmci

              cfmsxr=clri(i,j,k)+clrs(i,j,k)+clrg(i,j,k)+clrsg(i,j,k)   &
     &          +frrg(i,j,k)-mssr-msgr

! -----

! Get the new potential temperature perturbation.

              ptpf(i,j,k)=ptpf(i,j,k)+(vdvr(i,j,k)*lv(i,j,k)            &
     &          +(clnmci+cfmsxr)*lf(i,j,k)+nuvdvx*ls(i,j,k))*cppiv

! -----

! Get the new water vapor, the cloud water, the rain water and the cloud
! ice mixing ratio.

              qvsink=nuvdvx+vdvr(i,j,k)

              qxsink=cfmsxr-cxcr-vdvr(i,j,k)

              if(qrf(i,j,k).lt.qxsink) then

                qxsink=cxcr+clnmci+(qxsink-qrf(i,j,k))

                qrf(i,j,k)=0.e0

              else

                qrf(i,j,k)=qrf(i,j,k)-qxsink

                qxsink=cxcr+clnmci

              end if

              if(qcf(i,j,k).lt.qxsink) then

                dqv=qxsink-qcf(i,j,k)

                ptpf(i,j,k)=ptpf(i,j,k)+cppiv*dqv*lv(i,j,k)

                qvsink=qvsink+dqv

                qcf(i,j,k)=0.e0

              else

                qcf(i,j,k)=qcf(i,j,k)-qxsink

              end if

              qxsink=clix-spxi-nmci+cnis(i,j,k)-vdvi(i,j,k)-nuvi(i,j,k)

              if(qif(i,j,k).lt.qxsink) then

                dqv=qxsink-qif(i,j,k)

                ptpf(i,j,k)=ptpf(i,j,k)+cppiv*dqv*ls(i,j,k)

                qvsink=qvsink+dqv

                qif(i,j,k)=0.e0

              else

                qif(i,j,k)=qif(i,j,k)-qxsink

              end if

              qvf(i,j,k)=max(qvf(i,j,k)-qvsink,0.e0)

! -----

! Get the new snow and the graupel mixing ratio.

              qsf(i,j,k)=max(qsf(i,j,k)                                 &
     &          -(cxsg+mssr+spsi(i,j,k)-vdvs(i,j,k)                     &
     &          -clcs(i,j,k)-clrs(i,j,k)-clis(i,j,k)-cnis(i,j,k)),0.e0)

              qgf(i,j,k)=max(qgf(i,j,k)                                 &
     &          +(cxsg-msgr-spgi(i,j,k)+vdvg(i,j,k)                     &
     &          +clri(i,j,k)+clir(i,j,k)+clcg(i,j,k)                    &
     &          +clrg(i,j,k)+clig(i,j,k)+clrsg(i,j,k)+frrg(i,j,k)),0.e0)

! -----

! Get the new cloud ice concentrations.

              if(qcp(i,j,k).gt.thresq) then

                ncif(i,j,k)=ncif(i,j,k)                                 &
     &            +nuci(i,j,k)*nccp(i,j,k)/qcp(i,j,k)

              end if

              if(qip(i,j,k).gt.thresq) then

                if(vdvi(i,j,k).lt.0.e0) then

                  ncif(i,j,k)=ncif(i,j,k)-(agin(i,j,k)                  &
     &              +(clix+mlic(i,j,k)-vdvi(i,j,k))/mi(i,j,k)           &
     &              -(spxi+nuvi(i,j,k))*mi0iv+cnis(i,j,k)*ms0iv)

                else

                  ncif(i,j,k)=ncif(i,j,k)                               &
     &              -(agin(i,j,k)+(clix+mlic(i,j,k))/mi(i,j,k)          &
     &              -(spxi+nuvi(i,j,k))*mi0iv+cnis(i,j,k)*ms0iv)

                end if

              else

                ncif(i,j,k)=ncif(i,j,k)+(spxi+nuvi(i,j,k))*mi0iv

              end if

! -----

            end do
            end do

!$omp end do

          end do

!! -----

!! Solve the new potential temperature perturbation, the mixing ratio
!! and ice concentrations.

        else if(abs(cphopt).eq.3) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,cppiv,clix,spxi,cxcr,cxsg)       &
!$omp&   private(nmci,mssr,msgr,nuvdvx,clnmci,cfmsxr,qvsink,qxsink,dqv)

            do j=1,nj-1
            do i=1,ni-1

! Set the common used variables.

              cppiv=1.e0/(cp*pi(i,j,k))

              clix=clir(i,j,k)+clis(i,j,k)+clig(i,j,k)

              spxi=spsi(i,j,k)+spgi(i,j,k)

              cxcr=clcr(i,j,k)+cncr(i,j,k)

              cxsg=clsr(i,j,k)+clsg(i,j,k)+cnsg(i,j,k)

              nmci=nuci(i,j,k)-mlic(i,j,k)

              mssr=mlsr(i,j,k)+shsr(i,j,k)
              msgr=mlgr(i,j,k)+shgr(i,j,k)

              nuvdvx=vdvi(i,j,k)+vdvs(i,j,k)+vdvg(i,j,k)+nuvi(i,j,k)

              clnmci=clcs(i,j,k)+clcg(i,j,k)+nmci

              cfmsxr=clri(i,j,k)+clrs(i,j,k)+clrg(i,j,k)+clrsg(i,j,k)   &
     &          +frrg(i,j,k)-mssr-msgr

! -----

! Get the new potential temperature perturbation.

              ptpf(i,j,k)=ptpf(i,j,k)+(vdvr(i,j,k)*lv(i,j,k)            &
     &          +(clnmci+cfmsxr)*lf(i,j,k)+nuvdvx*ls(i,j,k))*cppiv

! -----

! Get the new water vapor, the cloud water, the rain water and the cloud
! ice mixing ratio.

              qvsink=nuvdvx+vdvr(i,j,k)

              qxsink=cfmsxr-cxcr-vdvr(i,j,k)

              if(qrf(i,j,k).lt.qxsink) then

                qxsink=cxcr+clnmci+(qxsink-qrf(i,j,k))

                qrf(i,j,k)=0.e0

              else

                qrf(i,j,k)=qrf(i,j,k)-qxsink

                qxsink=cxcr+clnmci

              end if

              if(qcf(i,j,k).lt.qxsink) then

                dqv=qxsink-qcf(i,j,k)

                ptpf(i,j,k)=ptpf(i,j,k)+cppiv*dqv*lv(i,j,k)

                qvsink=qvsink+dqv

                qcf(i,j,k)=0.e0

              else

                qcf(i,j,k)=qcf(i,j,k)-qxsink

              end if

              qxsink=clix-spxi-nmci+cnis(i,j,k)-vdvi(i,j,k)-nuvi(i,j,k)

              if(qif(i,j,k).lt.qxsink) then

                dqv=qxsink-qif(i,j,k)

                ptpf(i,j,k)=ptpf(i,j,k)+cppiv*dqv*ls(i,j,k)

                qvsink=qvsink+dqv

                qif(i,j,k)=0.e0

              else

                qif(i,j,k)=qif(i,j,k)-qxsink

              end if

              qvf(i,j,k)=max(qvf(i,j,k)-qvsink,0.e0)

! -----

! Get the new snow and the graupel mixing ratio.

              qsf(i,j,k)=max(qsf(i,j,k)                                 &
     &          -(cxsg+mssr+spsi(i,j,k)-vdvs(i,j,k)                     &
     &          -clcs(i,j,k)-clrs(i,j,k)-clis(i,j,k)-cnis(i,j,k)),0.e0)

              qgf(i,j,k)=max(qgf(i,j,k)                                 &
     &          +(cxsg-msgr-spgi(i,j,k)+vdvg(i,j,k)                     &
     &          +clri(i,j,k)+clir(i,j,k)+clcg(i,j,k)                    &
     &          +clrg(i,j,k)+clig(i,j,k)+clrsg(i,j,k)+frrg(i,j,k)),0.e0)

! -----

! Get the new cloud ice concentrations.

              if(qcp(i,j,k).gt.thresq) then

                ncif(i,j,k)=ncif(i,j,k)                                 &
     &            +nuci(i,j,k)*nccp(i,j,k)/qcp(i,j,k)

              end if

              if(qip(i,j,k).gt.thresq) then

                if(vdvi(i,j,k).lt.0.e0) then

                  ncif(i,j,k)=ncif(i,j,k)-(agin(i,j,k)                  &
     &              +(clix+mlic(i,j,k)-vdvi(i,j,k))/mi(i,j,k)           &
     &              -(spxi+nuvi(i,j,k))*mi0iv+cnis(i,j,k)*ms0iv)

                else

                  ncif(i,j,k)=ncif(i,j,k)                               &
     &              -(agin(i,j,k)+(clix+mlic(i,j,k))/mi(i,j,k)          &
     &              -(spxi+nuvi(i,j,k))*mi0iv+cnis(i,j,k)*ms0iv)

                end if

              else

                ncif(i,j,k)=ncif(i,j,k)+(spxi+nuvi(i,j,k))*mi0iv

              end if

! -----

! Get the new snow concentrations.

              if(qsp(i,j,k).gt.thresq) then

                if(vdvs(i,j,k).lt.0.e0) then

                  ncsf(i,j,k)=ncsf(i,j,k)                               &
     &              -(agsn(i,j,k)-cnis(i,j,k)*ms0iv                     &
     &              +clsrn(i,j,k)+clsgn(i,j,k)+cnsgn(i,j,k)             &
     &              +(mlsr(i,j,k)-vdvs(i,j,k))*ncsp(i,j,k)/qsp(i,j,k))

                else

                  ncsf(i,j,k)=ncsf(i,j,k)-(agsn(i,j,k)-cnis(i,j,k)*ms0iv&
     &              +clsrn(i,j,k)+clsgn(i,j,k)+cnsgn(i,j,k)             &
     &              +mlsr(i,j,k)*ncsp(i,j,k)/qsp(i,j,k))

                end if

              else

                ncsf(i,j,k)=ncsf(i,j,k)+cnis(i,j,k)*ms0iv

              end if

! -----

! Get the new graupel concentrations.

              if(qgp(i,j,k).gt.thresq) then

                if(vdvg(i,j,k).lt.0.e0) then

                 ncgf(i,j,k)=ncgf(i,j,k)                                &
     &             +((vdvg(i,j,k)-mlgr(i,j,k))*ncgp(i,j,k)/qgp(i,j,k)   &
     &             +clrsn(i,j,k)+clrin(i,j,k)+cnsgn(i,j,k)+frrgn(i,j,k))

                else

                 ncgf(i,j,k)=ncgf(i,j,k)                                &
     &             -(mlgr(i,j,k)*ncgp(i,j,k)/qgp(i,j,k)                 &
     &             -clrsn(i,j,k)-clrin(i,j,k)-cnsgn(i,j,k)-frrgn(i,j,k))

                end if

              else

                ncgf(i,j,k)=ncgf(i,j,k)                                 &
     &            +(clrsn(i,j,k)+clrin(i,j,k)+cnsgn(i,j,k)+frrgn(i,j,k))

              end if

! -----

            end do
            end do

!$omp end do

          end do

!! -----

!! Solve the new potential temperature perturbation, the mixing ratio
!! and water and ice concentrations.

        else if(abs(cphopt).eq.4) then

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,mciv,miiv,msiv,mgiv,cppiv,clix,spxi,cxcr,cxsg)     &
!$omp&   private(nmci,mssr,msgr,nuvdvx,clnmci,cfmsxr,qvsink,qxsink,dqv)

            do j=1,nj-1
            do i=1,ni-1

! Set the common used variables.

              cppiv=1.e0/(cp*pi(i,j,k))

              clix=clir(i,j,k)+clis(i,j,k)+clig(i,j,k)

              spxi=spsi(i,j,k)+spgi(i,j,k)

              cxcr=clcr(i,j,k)+cncr(i,j,k)

              cxsg=clsr(i,j,k)+clsg(i,j,k)+cnsg(i,j,k)

              nmci=nuci(i,j,k)-mlic(i,j,k)

              mssr=mlsr(i,j,k)+shsr(i,j,k)
              msgr=mlgr(i,j,k)+shgr(i,j,k)

              nuvdvx=vdvi(i,j,k)+vdvs(i,j,k)+vdvg(i,j,k)+nuvi(i,j,k)

              clnmci=clcs(i,j,k)+clcg(i,j,k)+nmci

              cfmsxr=clri(i,j,k)+clrs(i,j,k)+clrg(i,j,k)+clrsg(i,j,k)   &
     &          +frrg(i,j,k)-mssr-msgr

! -----

! Get the new potential temperature perturbation.

              ptpf(i,j,k)=ptpf(i,j,k)+(vdvr(i,j,k)*lv(i,j,k)            &
     &          +(clnmci+cfmsxr)*lf(i,j,k)+nuvdvx*ls(i,j,k))*cppiv

! -----

! Get the new water vapor, the cloud water, the rain water and the cloud
! ice mixing ratio.

              qvsink=nuvdvx+vdvr(i,j,k)

              qxsink=cfmsxr-cxcr-vdvr(i,j,k)

              if(qrf(i,j,k).lt.qxsink) then

                qxsink=cxcr+clnmci+(qxsink-qrf(i,j,k))

                qrf(i,j,k)=0.e0

              else

                qrf(i,j,k)=qrf(i,j,k)-qxsink

                qxsink=cxcr+clnmci

              end if

              if(qcf(i,j,k).lt.qxsink) then

                dqv=qxsink-qcf(i,j,k)

                ptpf(i,j,k)=ptpf(i,j,k)+cppiv*dqv*lv(i,j,k)

                qvsink=qvsink+dqv

                qcf(i,j,k)=0.e0

              else

                qcf(i,j,k)=qcf(i,j,k)-qxsink

              end if

              qxsink=clix-spxi-nmci+cnis(i,j,k)-vdvi(i,j,k)-nuvi(i,j,k)

              if(qif(i,j,k).lt.qxsink) then

                dqv=qxsink-qif(i,j,k)

                ptpf(i,j,k)=ptpf(i,j,k)+cppiv*dqv*ls(i,j,k)

                qvsink=qvsink+dqv

                qif(i,j,k)=0.e0

              else

                qif(i,j,k)=qif(i,j,k)-qxsink

              end if

              qvf(i,j,k)=max(qvf(i,j,k)-qvsink,0.e0)

! -----

! Get the new snow and the graupel mixing ratio.

              qsf(i,j,k)=max(qsf(i,j,k)                                 &
     &          -(cxsg+mssr+spsi(i,j,k)-vdvs(i,j,k)                     &
     &          -clcs(i,j,k)-clrs(i,j,k)-clis(i,j,k)-cnis(i,j,k)),0.e0)

              qgf(i,j,k)=max(qgf(i,j,k)                                 &
     &          +(cxsg-msgr-spgi(i,j,k)+vdvg(i,j,k)                     &
     &          +clri(i,j,k)+clir(i,j,k)+clcg(i,j,k)                    &
     &          +clrg(i,j,k)+clig(i,j,k)+clrsg(i,j,k)+frrg(i,j,k)),0.e0)

! -----

! Get the new cloud water and cloud ice concentrations.

              if(qcp(i,j,k).gt.thresq) then

                mciv=nccp(i,j,k)/qcp(i,j,k)

                nccf(i,j,k)=nccf(i,j,k)-(agcn(i,j,k)+cncr(i,j,k)*mr0iv  &
     &            +(nuci(i,j,k)+clcr(i,j,k)+clcs(i,j,k)                 &
     &            +clcg(i,j,k))*mciv)

                ncif(i,j,k)=ncif(i,j,k)+nuci(i,j,k)*mciv

              end if

              if(qip(i,j,k).gt.thresq) then

                miiv=1.e0/mi(i,j,k)

                nccf(i,j,k)=nccf(i,j,k)+mlic(i,j,k)*miiv

                if(vdvi(i,j,k).lt.0.e0) then

                  ncif(i,j,k)=ncif(i,j,k)-(agin(i,j,k)                  &
     &              +(clix+mlic(i,j,k)-vdvi(i,j,k))*miiv                &
     &              -(spxi+nuvi(i,j,k))*mi0iv+cnis(i,j,k)*ms0iv)

                else

                  ncif(i,j,k)=ncif(i,j,k)                               &
     &              -(agin(i,j,k)+(clix+mlic(i,j,k))*miiv               &
     &              -(spxi+nuvi(i,j,k))*mi0iv+cnis(i,j,k)*ms0iv)

                end if

              else

                ncif(i,j,k)=ncif(i,j,k)+(spxi+nuvi(i,j,k))*mi0iv

              end if

! -----

! Get the new rain water, snow and graupel concentrations.

              if(qrp(i,j,k).gt.thresq) then

                ncrf(i,j,k)=ncrf(i,j,k)-(clrin(i,j,k)+clrsn(i,j,k)      &
     &            +frrgn(i,j,k)+agrn(i,j,k)-cncr(i,j,k)*mr0iv)

              else

                ncrf(i,j,k)=ncrf(i,j,k)+cncr(i,j,k)*mr0iv

              end if

              if(qsp(i,j,k).gt.thresq) then

                msiv=ncsp(i,j,k)/qsp(i,j,k)

                ncrf(i,j,k)=ncrf(i,j,k)+mlsr(i,j,k)*msiv

                if(vdvs(i,j,k).lt.0.e0) then

                  ncsf(i,j,k)=ncsf(i,j,k)                               &
     &              -(agsn(i,j,k)-cnis(i,j,k)*ms0iv                     &
     &              +(mlsr(i,j,k)-vdvs(i,j,k))*msiv                     &
     &              +clsrn(i,j,k)+clsgn(i,j,k)+cnsgn(i,j,k))

                else

                  ncsf(i,j,k)=ncsf(i,j,k)                               &
     &              -(agsn(i,j,k)-cnis(i,j,k)*ms0iv+mlsr(i,j,k)*msiv    &
     &              +clsrn(i,j,k)+clsgn(i,j,k)+cnsgn(i,j,k))

                end if

              else

                ncsf(i,j,k)=ncsf(i,j,k)+cnis(i,j,k)*ms0iv

              end if

              if(qgp(i,j,k).gt.thresq) then

                mgiv=ncgp(i,j,k)/qgp(i,j,k)

                ncrf(i,j,k)=ncrf(i,j,k)+mlgr(i,j,k)*mgiv

                if(vdvg(i,j,k).lt.0.e0) then

                 ncgf(i,j,k)=ncgf(i,j,k)                                &
     &             +((vdvg(i,j,k)-mlgr(i,j,k))*mgiv                     &
     &             +clrsn(i,j,k)+clrin(i,j,k)+cnsgn(i,j,k)+frrgn(i,j,k))

                else

                 ncgf(i,j,k)=ncgf(i,j,k)-(mlgr(i,j,k)*mgiv              &
     &             -clrsn(i,j,k)-clrin(i,j,k)-cnsgn(i,j,k)-frrgn(i,j,k))

                end if

              else

                ncgf(i,j,k)=ncgf(i,j,k)                                 &
     &            +(clrsn(i,j,k)+clrin(i,j,k)+cnsgn(i,j,k)+frrgn(i,j,k))

              end if

! -----

            end do
            end do

!$omp end do

          end do

        end if

!! -----

      end if

!!! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_newblk == DUMP_TARGET_newblk .and. .not. dump_done_newblk) then
  call dump_array_3d('ptpf_ref.bin', ptpf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvf_ref.bin', qvf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qcf_ref.bin', qcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qrf_ref.bin', qrf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qif_ref.bin', qif, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qsf_ref.bin', qsf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qgf_ref.bin', qgf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nccf_ref.bin', nccf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncrf_ref.bin', ncrf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncif_ref.bin', ncif, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncsf_ref.bin', ncsf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncgf_ref.bin', ncgf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_newblk = .true.
end if


call profile_stop(prof_id1, loop_len)

!!!! -----

      end subroutine s_newblk

!-----7--------------------------------------------------------------7--

      end module m_newblk

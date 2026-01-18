!***********************************************************************
      module m_more0q
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2001/10/18, 2001/11/20,
!                   2002/01/15, 2002/04/02, 2002/04/09, 2002/12/06,
!                   2003/04/30, 2003/05/19, 2003/12/12, 2004/08/01,
!                   2004/09/01, 2004/09/25, 2004/10/12, 2004/12/17,
!                   2005/04/04, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13, 2013/10/08

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the mixing ratio to not being less than 0.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_temparam
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: more0q, s_more0q

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface more0q

        module procedure s_more0q

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_more0q(thresq,ni,nj,nk,qcp,qrp,qip,qsp,qgp,          &
     &                    qcf,qrf,qif,qsf,qgf,nuvi,nuci,clcr,clcs,clcg, &
     &                    clri,clrs,clrg,clir,clis,clig,clsr,clsg,clrsg,&
     &                    vdvr,vdvi,vdvs,vdvg,cncr,cnis,cnsg,spsi,spgi, &
     &                    mlic,mlsr,mlgr,frrg,shsr,shgr)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

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

      real, intent(in) :: qcf(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at future

      real, intent(in) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at future

      real, intent(in) :: qif(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio at future

      real, intent(in) :: qsf(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio at future

      real, intent(in) :: qgf(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio at future

! Input and output variables

      real, intent(inout) :: nuvi(0:ni+1,0:nj+1,1:nk)
                       ! Nucleation rate of deposition or sorption

      real, intent(inout) :: nuci(0:ni+1,0:nj+1,1:nk)
                       ! Nucleation rate
                       ! of condensation, contact and homogeneous

      real, intent(inout) :: clcr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and rain water

      real, intent(inout) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and snow

      real, intent(inout) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and graupel

      real, intent(inout) :: clri(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, intent(inout) :: clrs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! from rain water to snow

      real, intent(inout) :: clrg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and graupel

      real, intent(inout) :: clir(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and cloud ice

      real, intent(inout) :: clis(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud ice and snow

      real, intent(inout) :: clig(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud ice and graupel

      real, intent(inout) :: clsr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! from snow to rain water

      real, intent(inout) :: clsg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between snow and graupel

      real, intent(inout) :: clrsg(0:ni+1,0:nj+1,1:nk)
                       ! Production rate of graupel
                       ! from collection rate form rain to snow

      real, intent(inout) :: vdvr(0:ni+1,0:nj+1,1:nk)
                       ! Evaporation rate from rain water to water vapor

      real, intent(inout) :: vdvi(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to cloud ice

      real, intent(inout) :: vdvs(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to snow

      real, intent(inout) :: vdvg(0:ni+1,0:nj+1,1:nk)
                       ! Deposiotn rate from water vapor to graupel

      real, intent(inout) :: cncr(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from cloud water to rain water

      real, intent(inout) :: cnis(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from cloud ice to snow

      real, intent(inout) :: cnsg(0:ni+1,0:nj+1,1:nk)
                       ! Conversion rate from snow to graupel

      real, intent(inout) :: spsi(0:ni+1,0:nj+1,1:nk)
                       ! Secondary nucleation rate from snow

      real, intent(inout) :: spgi(0:ni+1,0:nj+1,1:nk)
                       ! Secondary nucleation rate from graupel

      real, intent(inout) :: mlic(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from cloud ice to cloud water

      real, intent(inout) :: mlsr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from snow to rain water

      real, intent(inout) :: mlgr(0:ni+1,0:nj+1,1:nk)
                       ! Melting rate from graupel to rain water

      real, intent(inout) :: frrg(0:ni+1,0:nj+1,1:nk)
                       ! Freezing rate from rain water to graupel

      real, intent(inout) :: shsr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from snow

      real, intent(inout) :: shgr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from graupel

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real sink        ! Total sink amount

      real handle      ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_more0q = 0
      integer, parameter :: DUMP_TARGET_more0q = 45720
      logical, save :: dump_done_more0q = .false.


!-----7--------------------------------------------------------------7--

!!!! Force the mixing ratio to not being less than 0.

!@llm start meta_info ----------------------------------------------------
! Location: more0q.f90 :: s_more0q
! Summary : Force mixing ratios to be non-negative by scaling microphysical process
!           rates when sink exceeds available mass for cloud water, ice, rain, snow, graupel
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic max() function - GPU compatible
!   - Uses module variable evapor_opt from m_temparam
!   - Complex conditional logic with multiple branching per hydrometeor type
!   - Writes to many microphysical rate arrays (nuvi, nuci, clcr, clcs, ... etc)
!   - Branching on nk==1 for 2D vs 3D handling
!   - No synchronization constructs besides implicit barriers at !$omp end do
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - May need to restructure conditionals for better GPU branch divergence
!   - Collapse nested i,j loops for better GPU occupancy
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 2.650s (0.09%)
!   - AvgTime: 0.058ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('more0q.f90', 's_more0q', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_more0q = dump_call_count_more0q + 1
if (dump_call_count_more0q == DUMP_TARGET_more0q .and. .not. dump_done_more0q) then
  call dump_init('more0q')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('thresq', thresq)
  call dump_array_3d('qcp.bin', qcp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qrp.bin', qrp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qip.bin', qip, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qsp.bin', qsp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qgp.bin', qgp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qcf.bin', qcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qrf.bin', qrf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qif.bin', qif, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qsf.bin', qsf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qgf.bin', qgf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nuvi_in.bin', nuvi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nuci_in.bin', nuci, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcr_in.bin', clcr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcs_in.bin', clcs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcg_in.bin', clcg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clri_in.bin', clri, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrs_in.bin', clrs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrg_in.bin', clrg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clir_in.bin', clir, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clis_in.bin', clis, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clig_in.bin', clig, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsr_in.bin', clsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsg_in.bin', clsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrsg_in.bin', clrsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvr_in.bin', vdvr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvi_in.bin', vdvi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvs_in.bin', vdvs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvg_in.bin', vdvg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cncr_in.bin', cncr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cnis_in.bin', cnis, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cnsg_in.bin', cnsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('spsi_in.bin', spsi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('spgi_in.bin', spgi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mlic_in.bin', mlic, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mlsr_in.bin', mlsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mlgr_in.bin', mlgr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('frrg_in.bin', frrg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('shsr_in.bin', shsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('shgr_in.bin', shgr, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

!!! In the case nk = 1.

      if(nk.eq.1) then

! In the case option, evapor_opt = 0, evaporation rate set to zero.

        if(evapor_opt.eq.0) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(qrp(i,j,1).gt.thresq) then
              vdvr(i,j,1)=max(vdvr(i,j,1),0.e0)
            end if

            if(qip(i,j,1).gt.thresq) then
              vdvi(i,j,1)=max(vdvi(i,j,1),0.e0)
            end if

            if(qsp(i,j,1).gt.thresq) then
              vdvs(i,j,1)=max(vdvs(i,j,1),0.e0)
            end if

            if(qgp(i,j,1).gt.thresq) then
              vdvg(i,j,1)=max(vdvg(i,j,1),0.e0)
            end if

          end do
          end do

!$omp end do

        end if

! -----

!! Perform forcing for each grid.

!$omp do schedule(runtime) private(i,j,sink,handle)

        do j=1,nj-1
        do i=1,ni-1

! For the cloud water mixing ratio.

          if(qcp(i,j,1).gt.thresq) then

            sink=nuci(i,j,1)+clcr(i,j,1)+clcs(i,j,1)+clcg(i,j,1)        &
     &        +cncr(i,j,1)-mlic(i,j,1)

            if(qcf(i,j,1).lt.sink) then

              handle=qcf(i,j,1)/sink

              nuci(i,j,1)=nuci(i,j,1)*handle
              clcr(i,j,1)=clcr(i,j,1)*handle
              clcs(i,j,1)=clcs(i,j,1)*handle
              clcg(i,j,1)=clcg(i,j,1)*handle
              cncr(i,j,1)=cncr(i,j,1)*handle

              mlic(i,j,1)=mlic(i,j,1)*handle

            end if

          end if

! -----

! For the cloud ice mixing ratio.

          if(qip(i,j,1).gt.thresq) then

            sink=clir(i,j,1)+clis(i,j,1)+clig(i,j,1)+cnis(i,j,1)        &
     &        +mlic(i,j,1)-vdvi(i,j,1)-nuvi(i,j,1)-nuci(i,j,1)          &
     &        -spsi(i,j,1)-spgi(i,j,1)

            if(qif(i,j,1).lt.sink) then

              handle=qif(i,j,1)/sink

              clir(i,j,1)=clir(i,j,1)*handle
              clis(i,j,1)=clis(i,j,1)*handle
              clig(i,j,1)=clig(i,j,1)*handle
              cnis(i,j,1)=cnis(i,j,1)*handle
              mlic(i,j,1)=mlic(i,j,1)*handle

              vdvi(i,j,1)=vdvi(i,j,1)*handle
              nuvi(i,j,1)=nuvi(i,j,1)*handle
              nuci(i,j,1)=nuci(i,j,1)*handle
              spsi(i,j,1)=spsi(i,j,1)*handle
              spgi(i,j,1)=spgi(i,j,1)*handle

            end if

          end if

! -----

! For the rain water mixing ratio.

          if(qrp(i,j,1).gt.thresq) then

            sink=clri(i,j,1)+clrs(i,j,1)+clrg(i,j,1)+clrsg(i,j,1)       &
     &        +frrg(i,j,1)-vdvr(i,j,1)-clcr(i,j,1)-cncr(i,j,1)          &
     &        -mlsr(i,j,1)-mlgr(i,j,1)-shsr(i,j,1)-shgr(i,j,1)

            if(qrf(i,j,1).lt.sink) then

              handle=qrf(i,j,1)/sink

              clri(i,j,1)=clri(i,j,1)*handle
              clrs(i,j,1)=clrs(i,j,1)*handle
              clrg(i,j,1)=clrg(i,j,1)*handle
              clrsg(i,j,1)=clrsg(i,j,1)*handle
              frrg(i,j,1)=frrg(i,j,1)*handle

              vdvr(i,j,1)=vdvr(i,j,1)*handle
              clcr(i,j,1)=clcr(i,j,1)*handle
              cncr(i,j,1)=cncr(i,j,1)*handle
              mlsr(i,j,1)=mlsr(i,j,1)*handle
              mlgr(i,j,1)=mlgr(i,j,1)*handle
              shsr(i,j,1)=shsr(i,j,1)*handle
              shgr(i,j,1)=shgr(i,j,1)*handle

            end if

          end if

! -----

! For the snow mixing ratio.

          if(qsp(i,j,1).gt.thresq) then

            sink=clsr(i,j,1)+clsg(i,j,1)+cnsg(i,j,1)+spsi(i,j,1)        &
     &        +mlsr(i,j,1)+shsr(i,j,1)-vdvs(i,j,1)-clcs(i,j,1)          &
     &        -clrs(i,j,1)-clis(i,j,1)-cnis(i,j,1)

            if(qsf(i,j,1).lt.sink) then

              handle=qsf(i,j,1)/sink

              clsr(i,j,1)=clsr(i,j,1)*handle
              clsg(i,j,1)=clsg(i,j,1)*handle
              cnsg(i,j,1)=cnsg(i,j,1)*handle
              spsi(i,j,1)=spsi(i,j,1)*handle
              mlsr(i,j,1)=mlsr(i,j,1)*handle
              shsr(i,j,1)=shsr(i,j,1)*handle

              vdvs(i,j,1)=vdvs(i,j,1)*handle
              clcs(i,j,1)=clcs(i,j,1)*handle
              clrs(i,j,1)=clrs(i,j,1)*handle
              clis(i,j,1)=clis(i,j,1)*handle
              cnis(i,j,1)=cnis(i,j,1)*handle

            end if

          end if

! -----

! For the graupel mixing ratio.

          if(qgp(i,j,1).gt.thresq) then

            sink=spgi(i,j,1)+mlgr(i,j,1)+shgr(i,j,1)-vdvg(i,j,1)        &
     &        -clri(i,j,1)-clir(i,j,1)-clsr(i,j,1)-clcg(i,j,1)          &
     &        -clrg(i,j,1)-clig(i,j,1)-clsg(i,j,1)-clrsg(i,j,1)         &
     &        -cnsg(i,j,1)-frrg(i,j,1)

            if(qgf(i,j,1).lt.sink) then

              handle=qgf(i,j,1)/sink

              spgi(i,j,1)=spgi(i,j,1)*handle
              mlgr(i,j,1)=mlgr(i,j,1)*handle
              shgr(i,j,1)=shgr(i,j,1)*handle

              vdvg(i,j,1)=vdvg(i,j,1)*handle
              clri(i,j,1)=clri(i,j,1)*handle
              clir(i,j,1)=clir(i,j,1)*handle
              clsr(i,j,1)=clsr(i,j,1)*handle
              clcg(i,j,1)=clcg(i,j,1)*handle
              clrg(i,j,1)=clrg(i,j,1)*handle
              clig(i,j,1)=clig(i,j,1)*handle
              clsg(i,j,1)=clsg(i,j,1)*handle
              clrsg(i,j,1)=clrsg(i,j,1)*handle
              cnsg(i,j,1)=cnsg(i,j,1)*handle
              frrg(i,j,1)=frrg(i,j,1)*handle

            end if

          end if

! -----

        end do
        end do

!$omp end do

!! -----

!!! -----

!!! In the case nk > 1.

      else

! In the case option, evapor_opt = 0, evaporation rate set to zero.

        if(evapor_opt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(qrp(i,j,k).gt.thresq) then
                vdvr(i,j,k)=max(vdvr(i,j,k),0.e0)
              end if

              if(qip(i,j,k).gt.thresq) then
                vdvi(i,j,k)=max(vdvi(i,j,k),0.e0)
              end if

              if(qsp(i,j,k).gt.thresq) then
                vdvs(i,j,k)=max(vdvs(i,j,k),0.e0)
              end if

              if(qgp(i,j,k).gt.thresq) then
                vdvg(i,j,k)=max(vdvg(i,j,k),0.e0)
              end if

            end do
            end do

!$omp end do

          end do

        end if

! -----

!! Perform forcing for each grid.

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,sink,handle)

          do j=1,nj-1
          do i=1,ni-1

! For the cloud water mixing ratio.

            if(qcp(i,j,k).gt.thresq) then

              sink=nuci(i,j,k)+clcr(i,j,k)+clcs(i,j,k)+clcg(i,j,k)      &
     &          +cncr(i,j,k)-mlic(i,j,k)

              if(qcf(i,j,k).lt.sink) then

                handle=qcf(i,j,k)/sink

                nuci(i,j,k)=nuci(i,j,k)*handle
                clcr(i,j,k)=clcr(i,j,k)*handle
                clcs(i,j,k)=clcs(i,j,k)*handle
                clcg(i,j,k)=clcg(i,j,k)*handle
                cncr(i,j,k)=cncr(i,j,k)*handle

                mlic(i,j,k)=mlic(i,j,k)*handle

              end if

            end if

! -----

! For the cloud ice mixing ratio.

            if(qip(i,j,k).gt.thresq) then

              sink=clir(i,j,k)+clis(i,j,k)+clig(i,j,k)+cnis(i,j,k)      &
     &          +mlic(i,j,k)-vdvi(i,j,k)-nuvi(i,j,k)-nuci(i,j,k)        &
     &          -spsi(i,j,k)-spgi(i,j,k)

              if(qif(i,j,k).lt.sink) then

                handle=qif(i,j,k)/sink

                clir(i,j,k)=clir(i,j,k)*handle
                clis(i,j,k)=clis(i,j,k)*handle
                clig(i,j,k)=clig(i,j,k)*handle
                cnis(i,j,k)=cnis(i,j,k)*handle
                mlic(i,j,k)=mlic(i,j,k)*handle

                vdvi(i,j,k)=vdvi(i,j,k)*handle
                nuvi(i,j,k)=nuvi(i,j,k)*handle
                nuci(i,j,k)=nuci(i,j,k)*handle
                spsi(i,j,k)=spsi(i,j,k)*handle
                spgi(i,j,k)=spgi(i,j,k)*handle

              end if

            end if

! -----

! For the rain water mixing ratio.

            if(qrp(i,j,k).gt.thresq) then

              sink=clri(i,j,k)+clrs(i,j,k)+clrg(i,j,k)+clrsg(i,j,k)     &
     &          +frrg(i,j,k)-vdvr(i,j,k)-clcr(i,j,k)-cncr(i,j,k)        &
     &          -mlsr(i,j,k)-mlgr(i,j,k)-shsr(i,j,k)-shgr(i,j,k)

              if(qrf(i,j,k).lt.sink) then

                handle=qrf(i,j,k)/sink

                clri(i,j,k)=clri(i,j,k)*handle
                clrs(i,j,k)=clrs(i,j,k)*handle
                clrg(i,j,k)=clrg(i,j,k)*handle
                clrsg(i,j,k)=clrsg(i,j,k)*handle
                frrg(i,j,k)=frrg(i,j,k)*handle

                vdvr(i,j,k)=vdvr(i,j,k)*handle
                clcr(i,j,k)=clcr(i,j,k)*handle
                cncr(i,j,k)=cncr(i,j,k)*handle
                mlsr(i,j,k)=mlsr(i,j,k)*handle
                mlgr(i,j,k)=mlgr(i,j,k)*handle
                shsr(i,j,k)=shsr(i,j,k)*handle
                shgr(i,j,k)=shgr(i,j,k)*handle

              end if

            end if

! -----

! For the snow mixing ratio.

            if(qsp(i,j,k).gt.thresq) then

              sink=clsr(i,j,k)+clsg(i,j,k)+cnsg(i,j,k)+spsi(i,j,k)      &
     &          +mlsr(i,j,k)+shsr(i,j,k)-vdvs(i,j,k)-clcs(i,j,k)        &
     &          -clrs(i,j,k)-clis(i,j,k)-cnis(i,j,k)

              if(qsf(i,j,k).lt.sink) then

                handle=qsf(i,j,k)/sink

                clsr(i,j,k)=clsr(i,j,k)*handle
                clsg(i,j,k)=clsg(i,j,k)*handle
                cnsg(i,j,k)=cnsg(i,j,k)*handle
                spsi(i,j,k)=spsi(i,j,k)*handle
                mlsr(i,j,k)=mlsr(i,j,k)*handle
                shsr(i,j,k)=shsr(i,j,k)*handle

                vdvs(i,j,k)=vdvs(i,j,k)*handle
                clcs(i,j,k)=clcs(i,j,k)*handle
                clrs(i,j,k)=clrs(i,j,k)*handle
                clis(i,j,k)=clis(i,j,k)*handle
                cnis(i,j,k)=cnis(i,j,k)*handle

              end if

            end if

! -----

! For the graupel mixing ratio.

            if(qgp(i,j,k).gt.thresq) then

              sink=spgi(i,j,k)+mlgr(i,j,k)+shgr(i,j,k)-vdvg(i,j,k)      &
     &          -clri(i,j,k)-clir(i,j,k)-clsr(i,j,k)-clcg(i,j,k)        &
     &          -clrg(i,j,k)-clig(i,j,k)-clsg(i,j,k)-clrsg(i,j,k)       &
     &          -cnsg(i,j,k)-frrg(i,j,k)

              if(qgf(i,j,k).lt.sink) then

                handle=qgf(i,j,k)/sink

                spgi(i,j,k)=spgi(i,j,k)*handle
                mlgr(i,j,k)=mlgr(i,j,k)*handle
                shgr(i,j,k)=shgr(i,j,k)*handle

                vdvg(i,j,k)=vdvg(i,j,k)*handle
                clri(i,j,k)=clri(i,j,k)*handle
                clir(i,j,k)=clir(i,j,k)*handle
                clsr(i,j,k)=clsr(i,j,k)*handle
                clcg(i,j,k)=clcg(i,j,k)*handle
                clrg(i,j,k)=clrg(i,j,k)*handle
                clig(i,j,k)=clig(i,j,k)*handle
                clsg(i,j,k)=clsg(i,j,k)*handle
                clrsg(i,j,k)=clrsg(i,j,k)*handle
                cnsg(i,j,k)=cnsg(i,j,k)*handle
                frrg(i,j,k)=frrg(i,j,k)*handle

              end if

            end if

! -----

          end do
          end do

!$omp end do

        end do

!! -----

      end if

!!! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_more0q == DUMP_TARGET_more0q .and. .not. dump_done_more0q) then
  call dump_array_3d('nuvi_ref.bin', nuvi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nuci_ref.bin', nuci, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcr_ref.bin', clcr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcs_ref.bin', clcs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcg_ref.bin', clcg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clri_ref.bin', clri, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrs_ref.bin', clrs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrg_ref.bin', clrg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clir_ref.bin', clir, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clis_ref.bin', clis, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clig_ref.bin', clig, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsr_ref.bin', clsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsg_ref.bin', clsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrsg_ref.bin', clrsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvr_ref.bin', vdvr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvi_ref.bin', vdvi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvs_ref.bin', vdvs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vdvg_ref.bin', vdvg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cncr_ref.bin', cncr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cnis_ref.bin', cnis, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('cnsg_ref.bin', cnsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('spsi_ref.bin', spsi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('spgi_ref.bin', spgi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mlic_ref.bin', mlic, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mlsr_ref.bin', mlsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('mlgr_ref.bin', mlgr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('frrg_ref.bin', frrg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('shsr_ref.bin', shsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('shgr_ref.bin', shgr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_more0q = .true.
end if


call profile_stop(prof_id1, loop_len)

!!!! -----

      end subroutine s_more0q

!-----7--------------------------------------------------------------7--

      end module m_more0q

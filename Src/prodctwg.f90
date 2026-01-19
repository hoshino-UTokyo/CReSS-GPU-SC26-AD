!***********************************************************************
      module m_prodctwg
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2001/06/29, 2001/10/18,
!                   2001/11/20, 2002/01/15, 2002/04/02, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/12/12, 2004/03/22,
!                   2004/04/15, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2005/01/07, 2005/04/04, 2006/03/06,
!                   2006/04/03, 2007/06/27, 2007/10/19, 2007/11/26,
!                   2008/01/11, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the graupel production rate.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_dump_kernel
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: prodctwg, s_prodctwg

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface prodctwg

        module procedure s_prodctwg

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_prodctwg(cphopt,dtb,thresq,ni,nj,nk,rbr,             &
     &                      qi,qs,qg,ncs,tcel,qvsst0,lv,lf,kp,dv,       &
     &                      vntg,clcg,clrg,clir,clis,clig,clsr,clsg,    &
     &                      clsrn,clsgn,pgwet)
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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow

      real, intent(in) :: tcel(0:ni+1,0:nj+1,1:nk)
                       ! Ambient air temperature

      real, intent(in) :: qvsst0(0:ni+1,0:nj+1,1:nk)
                       ! Super saturation mixing ratio at melting point

      real, intent(in) :: lv(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of evapolation

      real, intent(in) :: lf(0:ni+1,0:nj+1,1:nk)
                       ! Latent heat of fusion

      real, intent(in) :: kp(0:ni+1,0:nj+1,1:nk)
                       ! Thermal conductivity of air

      real, intent(in) :: dv(0:ni+1,0:nj+1,1:nk)
                       ! Molecular diffusivity of water

      real, intent(in) :: vntg(0:ni+1,0:nj+1,1:nk)
                       ! Ventilation factor for graupel

      real, intent(in) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between cloud water and graupel

      real, intent(in) :: clrg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate
                       ! between rain water and graupel

! Input and output variables

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

      real, intent(inout) :: clsrn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, intent(inout) :: clsgn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between snow and graupel

! Output variable

      real, intent(out) :: pgwet(0:ni+1,0:nj+1,1:nk)
                       ! Graupel production rate for moist process

! Internal shared variables

      real cc2dtn      ! - 2.0 x cc x dtb

      real eigiv       ! Inverse of eig
      real esgiv       ! Inverse of esg

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real pgdry       ! Graupel production rate for dry process

      real lfice       ! lf + cw x tcel

      real cligw       ! eigiv x clig
      real clsgw       ! esgiv x clsg

      real sink        ! Sink amount in collection

      real a           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_prodctwg = 0
      integer, parameter :: DUMP_TARGET_prodctwg = 45720
      logical, save :: dump_done_prodctwg = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variable.

      cc2dtn=-2.e0*cc*dtb

      eigiv=1.e0/eig
      esgiv=1.e0/esg

! -----

!!! Calculate the graupel production rate.

!@llm start meta_info ----------------------------------------------------
! Location: prodctwg.f90 :: s_prodctwg
! Summary : Calculate graupel production rate for wet/dry growth processes
!           in cloud microphysics scheme.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Complex conditional logic with many if/else branches
!   - Reads/writes to multiple 3D arrays (clir, clis, clig, clsr, clsg, etc.)
!   - Many private local variables (pgdry, lfice, cligw, clsgw, sink, a)
!   - Data-dependent branching may cause thread divergence on GPU
!   - No explicit barriers but implicit at !$omp end do
! Next:
!   - Use OpenACC kernels or parallel loop with private clause for locals
!   - Thread divergence from conditionals may reduce GPU efficiency
!   - Consider restructuring conditionals to minimize divergence
!   - Ensure all read/write arrays are in data region
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 0.986s (0.03%)
!   - AvgTime: 0.022ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('prodctwg.f90', 's_prodctwg', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_prodctwg = dump_call_count_prodctwg + 1
if (dump_call_count_prodctwg == DUMP_TARGET_prodctwg .and. .not. dump_done_prodctwg) then
  call dump_init('prodctwg')
  call dump_scalar_i('cphopt', cphopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('dtb', dtb)
  call dump_scalar_r('thresq', thresq)
  call dump_scalar_r('ci', ci)
  call dump_scalar_r('cw', cw)
  call dump_array_3d('rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qi.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qg.bin', qg, 0, ni+1, 0, nj+1, 1, nk)
  if (abs(cphopt) >= 3) then
    call dump_array_3d('ncs.bin', ncs, 0, ni+1, 0, nj+1, 1, nk)
  end if
  call dump_array_3d('tcel.bin', tcel, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvsst0.bin', qvsst0, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('lv.bin', lv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('lf.bin', lf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('kp.bin', kp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('dv.bin', dv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vntg.bin', vntg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcg.bin', clcg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrg.bin', clrg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clir_in.bin', clir, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clis_in.bin', clis, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clig_in.bin', clig, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsr_in.bin', clsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsg_in.bin', clsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsrn_in.bin', clsrn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsgn_in.bin', clsgn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('cc2dtn', cc2dtn)
  call dump_scalar_r('eigiv', eigiv)
  call dump_scalar_r('esgiv', esgiv)
end if

!$omp parallel default(shared) private(k)

!! In the case nk = 1.

      if(nk.eq.1) then

! Perform calculating in the case the option abs(cphopt) is equal to 2.

        if(abs(cphopt).eq.2) then

!$omp do schedule(runtime) private(i,j,pgdry,lfice,cligw,clsgw,sink,a)

          do j=1,nj-1
          do i=1,ni-1

            if(qg(i,j,1).gt.thresq) then

              if(tcel(i,j,1).lt.t0cel) then

                pgdry=clcg(i,j,1)+clrg(i,j,1)+clig(i,j,1)+clsg(i,j,1)

                lfice=lf(i,j,1)+cw*tcel(i,j,1)

                cligw=eigiv*clig(i,j,1)
                clsgw=esgiv*clsg(i,j,1)

                pgwet(i,j,1)=cc2dtn*vntg(i,j,1)                         &
     &            *(lv(i,j,1)*dv(i,j,1)*rbr(i,j,1)*qvsst0(i,j,1)        &
     &            +kp(i,j,1)*tcel(i,j,1))/(lfice*rbr(i,j,1))            &
     &            +(cligw+clsgw)*(1.e0-ci*tcel(i,j,1)/lfice)

                if(pgwet(i,j,1).gt.0.e0.and.pgwet(i,j,1).lt.pgdry) then

                  clig(i,j,1)=cligw
                  clsg(i,j,1)=clsgw

                  sink=clir(i,j,1)+clis(i,j,1)+clig(i,j,1)

                  if(qi(i,j,1).lt.sink) then

                    a=qi(i,j,1)/sink

                    clir(i,j,1)=clir(i,j,1)*a
                    clis(i,j,1)=clis(i,j,1)*a
                    clig(i,j,1)=clig(i,j,1)*a

                  end if

                  sink=clsr(i,j,1)+clsg(i,j,1)

                  if(qs(i,j,1).lt.sink) then

                    a=qs(i,j,1)/sink

                    clsr(i,j,1)=clsr(i,j,1)*a
                    clsg(i,j,1)=clsg(i,j,1)*a

                  end if

                else

                  pgwet(i,j,1)=-1.e0

                end if

              else

                pgwet(i,j,1)=-1.e0

              end if

            else

              pgwet(i,j,1)=-1.e0

            end if

          end do
          end do

!$omp end do

! -----

! Perform calculating in the case the option abs(cphopt) is greater
! than 2.

        else if(abs(cphopt).ge.3) then

!$omp do schedule(runtime) private(i,j,pgdry,lfice,cligw,clsgw,sink,a)

          do j=1,nj-1
          do i=1,ni-1

            if(qg(i,j,1).gt.thresq) then

              if(tcel(i,j,1).lt.t0cel) then

                pgdry=clcg(i,j,1)+clrg(i,j,1)+clig(i,j,1)+clsg(i,j,1)

                lfice=lf(i,j,1)+cw*tcel(i,j,1)

                cligw=eigiv*clig(i,j,1)
                clsgw=esgiv*clsg(i,j,1)

                pgwet(i,j,1)=cc2dtn*vntg(i,j,1)                         &
     &            *(lv(i,j,1)*dv(i,j,1)*rbr(i,j,1)*qvsst0(i,j,1)        &
     &            +kp(i,j,1)*tcel(i,j,1))/(lfice*rbr(i,j,1))            &
     &            +(cligw+clsgw)*(1.e0-ci*tcel(i,j,1)/lfice)

                if(pgwet(i,j,1).gt.0.e0.and.pgwet(i,j,1).lt.pgdry) then

                  clig(i,j,1)=cligw
                  clsg(i,j,1)=clsgw

                  clsgn(i,j,1)=esgiv*clsgn(i,j,1)

                  sink=clir(i,j,1)+clis(i,j,1)+clig(i,j,1)

                  if(qi(i,j,1).lt.sink) then

                    a=qi(i,j,1)/sink

                    clir(i,j,1)=clir(i,j,1)*a
                    clis(i,j,1)=clis(i,j,1)*a
                    clig(i,j,1)=clig(i,j,1)*a

                  end if

                  sink=clsr(i,j,1)+clsg(i,j,1)

                  if(qs(i,j,1).lt.sink) then

                    a=qs(i,j,1)/sink

                    clsr(i,j,1)=clsr(i,j,1)*a
                    clsg(i,j,1)=clsg(i,j,1)*a

                  end if

                  sink=clsrn(i,j,1)+clsgn(i,j,1)

                  if(ncs(i,j,1).lt.sink) then

                    a=ncs(i,j,1)/sink

                    clsrn(i,j,1)=clsrn(i,j,1)*a
                    clsgn(i,j,1)=clsgn(i,j,1)*a

                  end if

                else

                  pgwet(i,j,1)=-1.e0

                end if

              else

                pgwet(i,j,1)=-1.e0

              end if

            else

              pgwet(i,j,1)=-1.e0

            end if

          end do
          end do

!$omp end do

        end if

! -----

!! -----

!! In the case nk > 1.

      else

! Perform calculating in the case the option abs(cphopt) is equal to 2.

        if(abs(cphopt).eq.2) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,pgdry,lfice,cligw,clsgw,sink,a)

            do j=1,nj-1
            do i=1,ni-1

              if(qg(i,j,k).gt.thresq) then

                if(tcel(i,j,k).lt.t0cel) then

                  pgdry=clcg(i,j,k)+clrg(i,j,k)+clig(i,j,k)+clsg(i,j,k)

                  lfice=lf(i,j,k)+cw*tcel(i,j,k)

                  cligw=eigiv*clig(i,j,k)
                  clsgw=esgiv*clsg(i,j,k)

                  pgwet(i,j,k)=cc2dtn*vntg(i,j,k)                       &
     &              *(lv(i,j,k)*dv(i,j,k)*rbr(i,j,k)*qvsst0(i,j,k)      &
     &              +kp(i,j,k)*tcel(i,j,k))/(lfice*rbr(i,j,k))          &
     &              +(cligw+clsgw)*(1.e0-ci*tcel(i,j,k)/lfice)

                  if(pgwet(i,j,k).gt.0.e0                               &
     &              .and.pgwet(i,j,k).lt.pgdry) then

                    clig(i,j,k)=cligw
                    clsg(i,j,k)=clsgw

                    sink=clir(i,j,k)+clis(i,j,k)+clig(i,j,k)

                    if(qi(i,j,k).lt.sink) then

                      a=qi(i,j,k)/sink

                      clir(i,j,k)=clir(i,j,k)*a
                      clis(i,j,k)=clis(i,j,k)*a
                      clig(i,j,k)=clig(i,j,k)*a

                    end if

                    sink=clsr(i,j,k)+clsg(i,j,k)

                    if(qs(i,j,k).lt.sink) then

                      a=qs(i,j,k)/sink

                      clsr(i,j,k)=clsr(i,j,k)*a
                      clsg(i,j,k)=clsg(i,j,k)*a

                    end if

                  else

                    pgwet(i,j,k)=-1.e0

                  end if

                else

                  pgwet(i,j,k)=-1.e0

                end if

              else

                pgwet(i,j,k)=-1.e0

              end if

            end do
            end do

!$omp end do

          end do

! -----

! Perform calculating in the case the option abs(cphopt) is greater
! than 2.

        else if(abs(cphopt).ge.3) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,pgdry,lfice,cligw,clsgw,sink,a)

            do j=1,nj-1
            do i=1,ni-1

              if(qg(i,j,k).gt.thresq) then

                if(tcel(i,j,k).lt.t0cel) then

                  pgdry=clcg(i,j,k)+clrg(i,j,k)+clig(i,j,k)+clsg(i,j,k)

                  lfice=lf(i,j,k)+cw*tcel(i,j,k)

                  cligw=eigiv*clig(i,j,k)
                  clsgw=esgiv*clsg(i,j,k)

                  pgwet(i,j,k)=cc2dtn*vntg(i,j,k)                       &
     &              *(lv(i,j,k)*dv(i,j,k)*rbr(i,j,k)*qvsst0(i,j,k)      &
     &              +kp(i,j,k)*tcel(i,j,k))/(lfice*rbr(i,j,k))          &
     &              +(cligw+clsgw)*(1.e0-ci*tcel(i,j,k)/lfice)

                  if(pgwet(i,j,k).gt.0.e0                               &
     &              .and.pgwet(i,j,k).lt.pgdry) then

                    clig(i,j,k)=cligw
                    clsg(i,j,k)=clsgw

                    clsgn(i,j,k)=esgiv*clsgn(i,j,k)

                    sink=clir(i,j,k)+clis(i,j,k)+clig(i,j,k)

                    if(qi(i,j,k).lt.sink) then

                      a=qi(i,j,k)/sink

                      clir(i,j,k)=clir(i,j,k)*a
                      clis(i,j,k)=clis(i,j,k)*a
                      clig(i,j,k)=clig(i,j,k)*a

                    end if

                    sink=clsr(i,j,k)+clsg(i,j,k)

                    if(qs(i,j,k).lt.sink) then

                      a=qs(i,j,k)/sink

                      clsr(i,j,k)=clsr(i,j,k)*a
                      clsg(i,j,k)=clsg(i,j,k)*a

                    end if

                    sink=clsrn(i,j,k)+clsgn(i,j,k)

                    if(ncs(i,j,k).lt.sink) then

                      a=ncs(i,j,k)/sink

                      clsrn(i,j,k)=clsrn(i,j,k)*a
                      clsgn(i,j,k)=clsgn(i,j,k)*a

                    end if

                  else

                    pgwet(i,j,k)=-1.e0

                  end if

                else

                  pgwet(i,j,k)=-1.e0

                end if

              else

                pgwet(i,j,k)=-1.e0

              end if

            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_prodctwg == DUMP_TARGET_prodctwg .and. .not. dump_done_prodctwg) then
  call dump_array_3d('pgwet_ref.bin', pgwet, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clir_ref.bin', clir, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clis_ref.bin', clis, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clig_ref.bin', clig, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsr_ref.bin', clsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsg_ref.bin', clsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsrn_ref.bin', clsrn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsgn_ref.bin', clsgn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_prodctwg = .true.
end if


call profile_stop(prof_id1, loop_len)

!!! -----

      end subroutine s_prodctwg

!-----7--------------------------------------------------------------7--

      end module m_prodctwg

!***********************************************************************
      module m_fallblk
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/27, 2000/11/17, 2001/04/15,
!                   2001/05/29, 2001/06/29, 2001/12/11, 2002/01/15,
!                   2002/04/02, 2002/12/02, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/10/31, 2003/11/28, 2003/12/12,
!                   2004/03/05, 2004/03/22, 2004/04/01, 2004/04/15,
!                   2004/05/07, 2004/05/31, 2004/06/10, 2004/08/01,
!                   2004/08/20, 2004/09/01, 2004/09/10, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/01/31, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/07/21,
!                   2006/09/30, 2007/05/14, 2007/10/19, 2007/11/26,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2009/11/05,
!                   2011/01/14, 2011/06/01, 2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the fall out.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkfall
      use m_comprofile
      use m_dump_kernel
      use m_comindx
      use m_commath
      use m_getiname
      use m_getrname
      use m_temparam
      use m_upwnp
      use m_upwqcg
      use m_upwqp

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: fallblk, s_fallblk

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface fallblk

        module procedure s_fallblk

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic min
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_fallblk(fpcphopt,fphaiopt,fpqcgopt,fpdz,             &
     &                     dtb,ni,nj,nk,jcb,rbr,rst,ucq,urq,            &
     &                     uiq,usq,ugq,uhq,ucn,urn,uin,usn,ugn,uhn,     &
     &                     qcf,qrf,qif,qsf,qgf,qhf,nccf,ncrf,ncif,ncsf, &
     &                     ncgf,nchf,qccf,qrcf,qicf,qscf,qgcf,qhcf,     &
     &                     prc,prr,pri,prs,prg,prh,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: ucq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud water

      real, intent(in) :: urq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water

      real, intent(in) :: uiq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud ice

      real, intent(in) :: usq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow

      real, intent(in) :: ugq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel

      real, intent(in) :: uhq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of hail

      real, intent(in) :: ucn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud water concentrations

      real, intent(in) :: urn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water concentrations

      real, intent(in) :: uin(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of cloud ice concentrations

      real, intent(in) :: usn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of snow concentrations

      real, intent(in) :: ugn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of graupel concentrations

      real, intent(in) :: uhn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of hail concentrations

! Input and output variables

      real, intent(inout) :: qcf(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixnig ratio at future

      real, intent(inout) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixnig ratio at future

      real, intent(inout) :: qif(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixnig ratio at future

      real, intent(inout) :: qsf(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixnig ratio at future

      real, intent(inout) :: qgf(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixnig ratio at future

      real, intent(inout) :: qhf(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixnig ratio at future

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

      real, intent(inout) :: nchf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of hail at future

      real, intent(inout) :: qccf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud water at future

      real, intent(inout) :: qrcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of rain water at future

      real, intent(inout) :: qicf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud ice at future

      real, intent(inout) :: qscf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of snow at future

      real, intent(inout) :: qgcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of graupel at future

      real, intent(inout) :: qhcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of hail at future

      real, intent(inout) :: prc(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for cloud water

      real, intent(inout) :: prr(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for rain

      real, intent(inout) :: pri(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for cloud ice

      real, intent(inout) :: prs(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for snow

      real, intent(inout) :: prg(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for graupel

      real, intent(inout) :: prh(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for hail

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution

      integer npstp    ! Number of steps for fall out time integration

      integer ipstp    ! Index of steps for fall out time integration

      real dz          ! Grid distance in z direction

      real dtp         ! Time steps interval of fall out integration

      real dtpc        ! Time steps interval of fall out integration
                       ! for cloud water

      real dtpr        ! Time steps interval of fall out integration
                       ! for rain water

      real dtpi        ! Time steps interval of fall out integration
                       ! for cloud ice

      real dtps        ! Time steps interval of fall out integration
                       ! for snow

      real dtpg        ! Time steps interval of fall out integration
                       ! for graupel

      real dtph        ! Time steps interval of fall out integration
                       ! for hail

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real dzjcb       ! dz x jcb


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_fallblk = 0
      integer, parameter :: DUMP_TARGET_fallblk = 360
      logical, save :: dump_done_fallblk = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getrname(fpdz,dz)

! -----

!! Perform the fall out.

! Initialize the processed variables.

      dtpc=dtb
      dtpr=dtb
      dtpi=dtb
      dtps=dtb
      dtpg=dtb
      dtph=dtb

! -----

! Get the minimum time interval.

!@llm start meta_info ----------------------------------------------------
! Location: fallblk.f90 :: subroutine s_fallblk
! Summary : Calculates minimum time step for precipitation fallout based on
!           CFL condition with terminal velocities of hydrometeors.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constant eps from commath via temparam.
!   - Uses OpenMP reduction(min:) for dtpc, dtpr, dtpi, dtps, dtpg, dtph.
!   - Uses intrinsic min() - GPU compatible.
!   - Multiple code paths based on flqcqi_opt, haiopt.
!   - Note: After parallel region, calls upwqp, upwnp, upwqcg in a loop.
! Next:
!   - Reduction requires GPU-compatible reduction pattern.
!   - Use atomicMin or warp-level reduction for CFL check.
!   - Consider computing reductions in separate kernel before main loop.
!   - The upwqp/upwnp/upwqcg calls need separate GPU porting.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 2.581s (0.09%)
!   - AvgTime: 7.169ms
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('fallblk.f90', 's_fallblk', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_fallblk = dump_call_count_fallblk + 1
if (dump_call_count_fallblk == DUMP_TARGET_fallblk .and. .not. dump_done_fallblk) then
  call dump_init('fallblk')
  call dump_scalar_i('fpcphopt', fpcphopt)
  call dump_scalar_i('fphaiopt', fphaiopt)
  call dump_scalar_i('fpqcgopt', fpqcgopt)
  call dump_scalar_i('fpdz', fpdz)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('dtb', dtb)
  call dump_array_3d('jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ucq.bin', ucq, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('urq.bin', urq, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('uiq.bin', uiq, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('usq.bin', usq, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ugq.bin', ugq, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('uhq.bin', uhq, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ucn.bin', ucn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('urn.bin', urn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('uin.bin', uin, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('usn.bin', usn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ugn.bin', ugn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('uhn.bin', uhn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qcf_in.bin', qcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qrf_in.bin', qrf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qif_in.bin', qif, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qsf_in.bin', qsf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qgf_in.bin', qgf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qhf_in.bin', qhf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nccf_in.bin', nccf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncrf_in.bin', ncrf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncif_in.bin', ncif, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncsf_in.bin', ncsf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncgf_in.bin', ncgf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nchf_in.bin', nchf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qccf_in.bin', qccf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qrcf_in.bin', qrcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qicf_in.bin', qicf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qscf_in.bin', qscf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qgcf_in.bin', qgcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qhcf_in.bin', qhcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('prc_in.bin', prc, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('prr_in.bin', prr, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('pri_in.bin', pri, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('prs_in.bin', prs, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('prg_in.bin', prg, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('prh_in.bin', prh, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('tmp1_in.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

      if(flqcqi_opt.eq.0) then

        if(haiopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,dzjcb) reduction(min: dtpr,dtps,dtpg)

            do j=1,nj-1
            do i=1,ni-1
              dzjcb=dz*jcb(i,j,k)

              dtpr=min(dzjcb/(urq(i,j,k)+eps),dtpr)
              dtps=min(dzjcb/(usq(i,j,k)+eps),dtps)
              dtpg=min(dzjcb/(ugq(i,j,k)+eps),dtpg)

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,dzjcb) reduction(min: dtpr,dtps,dtpg,dtph)

            do j=1,nj-1
            do i=1,ni-1
              dzjcb=dz*jcb(i,j,k)

              dtpr=min(dzjcb/(urq(i,j,k)+eps),dtpr)
              dtps=min(dzjcb/(usq(i,j,k)+eps),dtps)
              dtpg=min(dzjcb/(ugq(i,j,k)+eps),dtpg)
              dtph=min(dzjcb/(uhq(i,j,k)+eps),dtph)

            end do
            end do

!$omp end do

          end do

        end if

      else

        if(haiopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,dzjcb) reduction(min: dtpc,dtpr,dtpi,dtps,dtpg)

            do j=1,nj-1
            do i=1,ni-1
              dzjcb=dz*jcb(i,j,k)

              dtpc=min(dzjcb/(ucq(i,j,k)+eps),dtpc)
              dtpr=min(dzjcb/(urq(i,j,k)+eps),dtpr)
              dtpi=min(dzjcb/(uiq(i,j,k)+eps),dtpi)
              dtps=min(dzjcb/(usq(i,j,k)+eps),dtps)
              dtpg=min(dzjcb/(ugq(i,j,k)+eps),dtpg)

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,dzjcb)                           &
!$omp&   reduction(min: dtpc,dtpr,dtpi,dtps,dtpg,dtph)

            do j=1,nj-1
            do i=1,ni-1
              dzjcb=dz*jcb(i,j,k)

              dtpc=min(dzjcb/(ucq(i,j,k)+eps),dtpc)
              dtpr=min(dzjcb/(urq(i,j,k)+eps),dtpr)
              dtpi=min(dzjcb/(uiq(i,j,k)+eps),dtpi)
              dtps=min(dzjcb/(usq(i,j,k)+eps),dtps)
              dtpg=min(dzjcb/(ugq(i,j,k)+eps),dtpg)
              dtph=min(dzjcb/(uhq(i,j,k)+eps),dtph)

            end do
            end do

!$omp end do

          end do

        end if

      end if

!$omp end parallel

! Dump output data at target call
if (dump_call_count_fallblk == DUMP_TARGET_fallblk .and. .not. dump_done_fallblk) then
  call dump_array_3d('qcf_ref.bin', qcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qrf_ref.bin', qrf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qif_ref.bin', qif, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qsf_ref.bin', qsf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qgf_ref.bin', qgf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qhf_ref.bin', qhf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nccf_ref.bin', nccf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncrf_ref.bin', ncrf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncif_ref.bin', ncif, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncsf_ref.bin', ncsf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncgf_ref.bin', ncgf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nchf_ref.bin', nchf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qccf_ref.bin', qccf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qrcf_ref.bin', qrcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qicf_ref.bin', qicf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qscf_ref.bin', qscf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qgcf_ref.bin', qgcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qhcf_ref.bin', qhcf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('prc_ref.bin', prc, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('prr_ref.bin', prr, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('pri_ref.bin', pri, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('prs_ref.bin', prs, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('prg_ref.bin', prg, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('prh_ref.bin', prh, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('tmp1_ref.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_fallblk = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

! Get the time interval and number of steps.

      if(flqcqi_opt.eq.0) then

        if(haiopt.eq.0) then
          dtp=min(dtpr,dtps,dtpg)
        else
          dtp=min(dtpr,dtps,dtpg,dtph)
        end if

      else

        if(haiopt.eq.0) then
          dtp=min(dtpc,dtpr,dtpi,dtps,dtpg)
        else
          dtp=min(dtpc,dtpr,dtpi,dtps,dtpg,dtph)
        end if

      end if

      call chkfall(dtp)

      if(dtp.lt.dtb) then
        npstp=int((dtb+.001e0)/dtp)+1
        dtp=dtb/real(npstp)

      else
        npstp=1
        dtp=dtb

      end if

! -----

! Calculate the sedimentation and precipitation.

      if(flqcqi_opt.ne.0) then

        do ipstp=1,npstp

          call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,ucq,qcf,prc,  &
     &               tmp1)

          call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,uiq,qif,pri,  &
     &               tmp1)

          if(abs(cphopt).eq.4) then

            call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,ucn,nccf,tmp1)

          end if

          call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,uin,ncif,tmp1)

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,ucq,qccf,tmp1)

            end if

            call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,uiq,qicf,tmp1)

          end if

        end do

      end if

      do ipstp=1,npstp

        call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,urq,qrf,prr,    &
     &             tmp1)

        call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,usq,qsf,prs,    &
     &             tmp1)

        call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,ugq,qgf,prg,    &
     &             tmp1)

        if(haiopt.eq.1) then

          call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,uhq,qhf,prh,  &
     &               tmp1)

        end if

        if(abs(cphopt).eq.4) then

          call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,urn,ncrf,tmp1)

        end if

        if(abs(cphopt).ge.3) then

          call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,usn,ncsf,tmp1)
          call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,ugn,ncgf,tmp1)

          if(haiopt.eq.1) then

            call upwnp(iddziv,dtp,ni,nj,nk,rbr,rst,uhn,nchf,tmp1)

          end if

        end if

        if(cphopt.lt.0) then

          if(qcgopt.eq.2) then

            call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,urq,qrcf,tmp1)

          end if

          call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,usq,qscf,tmp1)
          call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,ugq,qgcf,tmp1)

          if(haiopt.eq.1) then

            call upwqcg(iddziv,dtp,ni,nj,nk,rbr,rst,uhq,qhcf,tmp1)

          end if

        end if

      end do

! -----

!! -----

      end subroutine s_fallblk

!-----7--------------------------------------------------------------7--

      end module m_fallblk

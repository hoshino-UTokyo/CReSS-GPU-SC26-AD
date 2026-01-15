!***********************************************************************
      module m_heatsfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/12/11, 2002/01/15, 2002/04/02, 2002/07/03,
!                   2003/01/04, 2003/01/20, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/10/31, 2003/12/12, 2004/04/15,
!                   2004/08/01, 2004/08/20, 2004/09/10, 2005/01/31,
!                   2007/05/21, 2007/07/30, 2007/10/19, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/02/27, 2011/09/22,
!                   2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the sensible and latent heat on the surface.

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

      public :: heatsfc, s_heatsfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface heatsfc

        module procedure s_heatsfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_heatsfc(fmois,ni,nj,nk,nund,t,qv,qvsfc,ct,cq,        &
     &                     land,kai,tund,tice,hs,le)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qvsfc(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio on surface

      real, intent(in) :: ct(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface heat flux

      real, intent(in) :: cq(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface moisture flux

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(in) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature

      real, intent(in) :: tice(0:ni+1,0:nj+1)
                       ! Mixed ice surface temperature

! Output variables

      real, intent(out) :: hs(0:ni+1,0:nj+1)
                       ! Sensible heat

      real, intent(out) :: le(0:ni+1,0:nj+1)
                       ! Latent heat

! Internal shared variable

      real cwmci       ! cw - ci

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real tmsfc       ! Temperature on mixed ice surface

      real lva         ! Latent heat of evaporation at lowest plane
      real lsa         ! Latent heat of sublimation at lowest plane


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_heatsfc = 0
      integer, parameter :: DUMP_TARGET_heatsfc = 361
      logical, save :: dump_done_heatsfc = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variable.

      cwmci=cw-ci

! -----

! Calculate the sensible and latent heat on the surface.

!@llm start meta_info ----------------------------------------------------
! Location: heatsfc.f90 :: s_heatsfc
! Summary : Compute sensible and latent heat fluxes on surface based on
!           land type, temperature, and moisture conditions.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls; uses intrinsic exp/log
!   - Complex branching based on fmois (dry/moist) and land type
!   - Reads from t, qv, qvsfc, ct, cq, kai, tund, tice
!   - Writes to hs, le arrays (output)
!   - No sync constructs
! Next:
!   - Convert to OpenACC with collapse(2) on j,i loops
!   - Branch logic based on land type may cause GPU thread divergence
!   - Consider separating dry/moist cases into different kernels
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.058s (0.00%)
!   - AvgTime: 0.161ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('heatsfc.f90', 's_heatsfc', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_heatsfc = dump_call_count_heatsfc + 1
if (dump_call_count_heatsfc == DUMP_TARGET_heatsfc .and. .not. dump_done_heatsfc) then
  call dump_init('heatsfc')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('nund', nund)
  call dump_scalar_r('cp', cp)
  call dump_scalar_r('t0', t0)
  call dump_scalar_r('lv0', lv0)
  call dump_scalar_r('lf0', lf0)
  call dump_array_3d('t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('qvsfc.bin', qvsfc, 0, ni+1, 0, nj+1)
  call dump_array_2d('ct.bin', ct, 0, ni+1, 0, nj+1)
  call dump_array_2d('cq.bin', cq, 0, ni+1, 0, nj+1)
  call dump_array_2d_int('land.bin', land, 0, ni+1, 0, nj+1)
  call dump_array_2d('kai.bin', kai, 0, ni+1, 0, nj+1)
  call dump_array_3d('tund.bin', tund, 0, ni+1, 0, nj+1, 1, nund)
  call dump_array_2d('tice.bin', tice, 0, ni+1, 0, nj+1)
  call dump_scalar_r('cwmci', cwmci)
  call dump_scalar_c('fmois', fmois)
end if

!$omp parallel default(shared)

      if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j,tmsfc)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).eq.1) then

            tmsfc=kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1)

            hs(i,j)=cp*ct(i,j)*(tmsfc-t(i,j,2))

          else

            hs(i,j)=cp*ct(i,j)*(tund(i,j,1)-t(i,j,2))

          end if

          le(i,j)=0.e0

        end do
        end do

!$omp end do

      else if(fmois(1:5).eq.'moist') then

!$omp do schedule(runtime) private(i,j,tmsfc,lva,lsa)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).lt.0) then

            lva=lv0*exp((.167e0+3.67e-4*t(i,j,2))*log(t0/t(i,j,2)))

            hs(i,j)=cp*ct(i,j)*(tund(i,j,1)-t(i,j,2))
            le(i,j)=lva*cq(i,j)*(qvsfc(i,j)-qv(i,j,2))

          else if(land(i,j).eq.1) then

            tmsfc=kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1)

            lva=lv0*exp((.167e0+3.67e-4*t(i,j,2))*log(t0/t(i,j,2)))

            hs(i,j)=cp*ct(i,j)*(tmsfc-t(i,j,2))
            le(i,j)=lva*cq(i,j)*(qvsfc(i,j)-qv(i,j,2))

          else

            if(t(i,j,2).gt.tlow) then

              lva=lv0*exp((.167e0+3.67e-4*t(i,j,2))*log(t0/t(i,j,2)))

              hs(i,j)=cp*ct(i,j)*(tund(i,j,1)-t(i,j,2))
              le(i,j)=lva*cq(i,j)*(qvsfc(i,j)-qv(i,j,2))

            else

              lsa=lv0*exp((.167e0+3.67e-4*t(i,j,2))*log(t0/t(i,j,2)))   &
     &          +(lf0+cwmci*(t(i,j,2)-t0))

              hs(i,j)=cp*ct(i,j)*(tund(i,j,1)-t(i,j,2))
              le(i,j)=lsa*cq(i,j)*(qvsfc(i,j)-qv(i,j,2))

            end if

          end if

        end do
        end do

!$omp end do

      end if

!$omp end parallel

! Dump output data at target call
if (dump_call_count_heatsfc == DUMP_TARGET_heatsfc .and. .not. dump_done_heatsfc) then
  call dump_array_2d('hs_ref.bin', hs, 0, ni+1, 0, nj+1)
  call dump_array_2d('le_ref.bin', le, 0, ni+1, 0, nj+1)
  call dump_finalize()
  dump_done_heatsfc = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_heatsfc

!-----7--------------------------------------------------------------7--

      end module m_heatsfc

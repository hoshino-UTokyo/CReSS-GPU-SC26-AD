!***********************************************************************
      module m_sadjstbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/03/06
!     Modification: 2006/08/08, 2006/09/30, 2007/03/29, 2007/10/19,
!                   2008/01/11, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/01/30, 2009/02/27, 2010/02/01, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the saturation adjustment for water in the case the
!     generated total water is more than critical value.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: sadjstbw, s_sadjstbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sadjstbw

        module procedure s_sadjstbw

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
      subroutine s_sadjstbw(fpthresq,ni,nj,nk,ptbr,pi,p,w,t,ptp,qv,qwtr,&
     &                      ptptmp,qvtmp,qwtmp)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpthresq
                       ! Formal parameter of unique index of thresq

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure [dyn/cm^2]

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio [g/g]

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk)
                       ! Total water mixing ratio [g/g]

! Input and output variables

      real, intent(inout) :: ptptmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary alternative
                       ! potential temperature perturbation

      real, intent(inout) :: qvtmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary alternative
                       ! water vapor mixing ratio [g/g]

      real, intent(inout) :: qwtmp(0:ni+1,0:nj+1,1:nk)
                       ! Temporary alternative
                       ! total water mixing ratio [g/g]

! Internal shared variables

      real thresq      ! Minimum threshold value of mixing ratio

      real es01        ! 10.0 x es0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real esw         ! Saturation vapor pressure for water
      real qvsw        ! Saturation mixing ratio for water

      real lvcpi       ! Latent heat of evaporation / (cp x pi)

      real dqw         ! Variation of total water mixing ratio

      real ttmp        ! Temporary alternative air temperature

      real a           ! Temporary variable
      real b           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpthresq,thresq)

! -----

! Set the common used variable.

      es01=10.e0*es0

! -----

! Perform the saturation adjustment.

!@llm start meta_info ----------------------------------------------------
! Location: sadjstbw.f90 :: s_sadjstbw
! Summary : Perform saturation adjustment for water when total water
!           exceeds critical value, computing new qv/qw/ptp from thermodynamics
! GPU diff: Medium
! Findings:
!   - No omp_get_thread usage
!   - No external function calls (only intrinsic exp, log)
!   - Complex conditional branching within loops
!   - Writes to ptptmp, qvtmp, qwtmp arrays
!   - Uses module variables from m_comphy (es0, t0, epsva, lv0, cp, qccrit)
!   - No synchronization constructs
! Next:
!   - Convert to OpenACC with parallel loop collapse(3)
!   - Ensure m_comphy module constants are accessible on device
!   - Consider branch divergence impact on GPU performance
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('sadjstbw.f90', 's_sadjstbw', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(2)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j,esw,qvsw,lvcpi,dqw,ttmp,a,b)

        do j=1,nj-1
        do i=1,ni-1

          if(w(i,j,k).gt.0.e0) then

            if(qwtr(i,j,k).gt.thresq) then

              if(qwtr(i,j,k-1).le.thresq.and.w(i,j,k).gt.0.e0) then

                if(qv(i,j,k).gt.thresq) then
                  qwtmp(i,j,k)=-11.e0
                else
                  qwtmp(i,j,k)=-1.e0
                end if

              else

                qwtmp(i,j,k)=-1.e0

              end if

            else

              a=1.e0/(t(i,j,k)-35.86e0)
              b=a*(t(i,j,k)-t0)

              esw=es01*exp(17.269e0*b)

              qvsw=epsva*esw/(p(i,j,k)-esw)

              if(qv(i,j,k).gt.qvsw) then

                lvcpi=lv0*exp((.167e0+3.67e-4*t(i,j,k))                 &
     &            *log(t0/t(i,j,k)))/(cp*pi(i,j,k))

                dqw=(qv(i,j,k)-qvsw)                                    &
     &            /(1.e0+17.269e0*a*(1.e0-b)*qvsw*lvcpi*pi(i,j,k))

                ptptmp(i,j,k)=ptp(i,j,k)+dqw*lvcpi

                qvtmp(i,j,k)=qv(i,j,k)-dqw
                qwtmp(i,j,k)=qwtr(i,j,k)+dqw

              else

                qwtmp(i,j,k)=-1.e0

              end if

            end if

          else

            qwtmp(i,j,k)=-1.e0

          end if

          if(qwtmp(i,j,k).gt.thresq) then

            ttmp=(ptbr(i,j,k)+ptptmp(i,j,k))*pi(i,j,k)

            a=1.e0/(ttmp-35.86e0)
            b=a*(ttmp-t0)

            esw=es01*exp(17.269e0*b)

            qvsw=epsva*esw/(p(i,j,k)-esw)

            lvcpi=lv0*exp((.167e0+3.67e-4*ttmp)*log(t0/ttmp))           &
     &        /(cp*pi(i,j,k))

            dqw=(qvtmp(i,j,k)-qvsw)                                     &
     &        /(1.e0+17.269e0*a*(1.e0-b)*qvsw*lvcpi*pi(i,j,k))

            qwtmp(i,j,k)=qwtmp(i,j,k)+dqw

            if(qwtmp(i,j,k).gt.qccrit) then

              ptptmp(i,j,k)=ptptmp(i,j,k)+dqw*lvcpi

              qvtmp(i,j,k)=qvtmp(i,j,k)-dqw

            else

              qwtmp(i,j,k)=-1.e0

            end if

          else

            qwtmp(i,j,k)=-1.e0

          end if

        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_sadjstbw

!-----7--------------------------------------------------------------7--

      end module m_sadjstbw

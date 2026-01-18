!***********************************************************************
      module m_buoywse
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 1999/11/19, 2000/01/17, 2000/07/05, 2001/06/29,
!                   2001/11/20, 2001/12/11, 2002/04/02, 2002/08/15,
!                   2003/01/04, 2003/04/30, 2003/05/19, 2003/11/05,
!                   2003/12/26, 2004/04/15, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the buoyancy in the small time steps integration for the
!     horizontally explicit and vertically explicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: buoywse, s_buoywse

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface buoywse

        module procedure s_buoywse

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_buoywse(fpgwmopt,ni,nj,nk,ptbr,rst,rcsq,pp,ptp,      &
     &                     wsml,wb8s)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

! Input and output variable

      real, intent(inout) :: wsml(0:ni+1,0:nj+1,1:nk)
                       ! Acoustic mood in w equation

! Internal shared variables

      integer gwmopt   ! Option for gravity wave mode integration

      real g05n        ! - 0.5 x g

      real, intent(inout) :: wb8s(0:ni+1,0:nj+1,1:nk)
                       ! 0.5 x buoyancy value
                       ! in small time steps at scalar points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpgwmopt,gwmopt)

! -----

! Set the common used variable.

      g05n=-.5e0*g

! -----

! Calculate the buoyancy in the small time steps.

!@llm start meta_info ----------------------------------------------------
! Location: buoywse.f90 :: subroutine s_buoywse
! Summary : Calculates buoyancy forcing for small time step integration
!           in horizontally/vertically explicit method.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constant g from comphy.
!   - No synchronization constructs.
!   - Conditional on gwmopt for buoyancy formulation selection.
!   - Two stages: (1) compute wb8s, (2) accumulate to wsml.
!   - All grid points are independent within each stage.
! Next:
!   - Direct OpenACC kernels for each loop nest.
!   - Consider fusing the two stages if wb8s is temporary.
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('buoywse.f90', 's_buoywse', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) &
     & * int((nj-2)-(2)+1,8) &
     & * int((ni-2)-(2)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      if(gwmopt.eq.0) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wb8s(i,j,k)=pp(i,j,k)*rst(i,j,k)/rcsq(i,j,k)*g05n
          end do
          end do

!$omp end do

        end do

      else

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wb8s(i,j,k)=rst(i,j,k)                                      &
     &        *(pp(i,j,k)/rcsq(i,j,k)-ptp(i,j,k)/ptbr(i,j,k))*g05n
          end do
          end do

!$omp end do

        end do

      end if

      do k=3,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wsml(i,j,k)=wsml(i,j,k)+(wb8s(i,j,k-1)+wb8s(i,j,k))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_buoywse

!-----7--------------------------------------------------------------7--

      end module m_buoywse

!***********************************************************************
      module m_getta3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the air temperature.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing
      use m_comprofile

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getta3d, s_getta3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getta3d

        module procedure s_getta3d

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
      subroutine s_getta3d(ni,nj,nk,ptbr,pi,ptp,t)
!***********************************************************************

! Input variables

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

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

! Output variable

      real, intent(out) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Calculate the air temperature.

!@llm start meta_info ----------------------------------------------------
! Location: getta3d.f90 :: s_getta3d
! Summary : Compute air temperature from base state potential temperature,
!           perturbation, and Exner function at all grid points.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Simple element-wise arithmetic: t = (ptbr + ptp) * pi
!   - No global writes, only output array t is modified
!   - No synchronization constructs other than implicit barrier at end do
! Next:
!   - Direct translation to OpenACC or OpenACC with collapsed loops
!   - Consider loop collapse for k,j,i dimensions
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('getta3d.f90', 's_getta3d', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          t(i,j,k)=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_getta3d

!-----7--------------------------------------------------------------7--

      end module m_getta3d

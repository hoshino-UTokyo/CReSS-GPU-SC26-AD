!***********************************************************************
      module m_move2d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 2000/01/17, 2001/01/15, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2004/01/09, 2004/07/10, 2007/01/31,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the grid moving velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname
      use m_comprofile

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: move2d, s_move2d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface move2d

        module procedure s_move2d

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
      subroutine s_move2d(fpumove,fpvmove,nlev,u1d,v1d)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpumove
                       ! Formal parameter of unique index of umove

      integer, intent(in) :: fpvmove
                       ! Formal parameter of unique index of vmove

      integer, intent(in) :: nlev
                       ! Horizontally averaged vertical dimension

! Input and output variables

      real, intent(inout) :: u1d(0:nlev)
                       ! Horizontally averaged x components of velocity

      real, intent(inout) :: v1d(0:nlev)
                       ! Horizontally averaged y components of velocity

! Internal shared variables

      real umove       ! x components of grid moving velocity
      real vmove       ! y components of grid moving velocity

! Internal private variable

      integer kl       ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpumove,umove)
      call getrname(fpvmove,vmove)

! -----

! Add the relative velocity to the u1d and the v1d.

!@llm start meta_info ----------------------------------------------------
! Location: move2d.f90 :: s_move2d
! Summary : Subtract grid moving velocity (umove, vmove) from horizontally
!           averaged velocity profiles u1d and v1d
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple subtraction operation on 1D arrays
!   - Writes to u1d, v1d arrays (in-place modification)
!   - Single worksharing construct with 1D loop
!   - No synchronization constructs besides implicit barriers at !$omp end do
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Small array size (nlev) may not benefit significantly from GPU offloading
!   - Consider keeping on CPU if nlev is small
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('move2d.f90', 's_move2d', &
   & 'OMP section 1')
end if
loop_len = int((nlev)-(0)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared)

!$omp do schedule(runtime) private(kl)

      do kl=0,nlev
        u1d(kl)=u1d(kl)-umove
        v1d(kl)=v1d(kl)-vmove
      end do

!$omp end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_move2d

!-----7--------------------------------------------------------------7--

      end module m_move2d

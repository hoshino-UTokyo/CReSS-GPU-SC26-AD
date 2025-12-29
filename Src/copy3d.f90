!***********************************************************************
      module m_copy3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/08/23, 1999/10/12,
!                   2000/01/17, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2006/01/10, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     copy 3 dimensional invar to outvar.

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

      public :: copy3d, s_copy3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface copy3d

        module procedure s_copy3d

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
      subroutine s_copy3d(imin,imax,jmin,jmax,kmin,kmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: invar(imin:imax,jmin:jmax,kmin:kmax)
                       ! Coping variable

! Output variable

      real, intent(out) :: outvar(imin:imax,jmin:jmax,kmin:kmax)
                       ! Copied variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Copy the invar to the outvar.

!@llm start meta_info ----------------------------------------------------
! Location: copy3d.f90 :: s_copy3d
! Summary : Simple 3D array copy from invar to outvar.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls
!   - Trivial memory copy operation
!   - Independent element-wise operations
!   - Outer k loop is sequential in current OpenMP structure
! Next:
!   - Collapse all three loops (k,j,i) for better GPU occupancy
!   - Consider using device-to-device memcpy for efficiency
!   - May be better to keep data resident on GPU and avoid copy calls
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('copy3d.f90', 's_copy3d', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(kmin)+1,8) &
     & * int((jmax)-(jmin)+1,8) &
     & * int((imax)-(imin)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=kmin,kmax

!$omp do schedule(runtime) private(i,j)

        do j=jmin,jmax
        do i=imin,imax
          outvar(i,j,k)=invar(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_copy3d

!-----7--------------------------------------------------------------7--

      end module m_copy3d

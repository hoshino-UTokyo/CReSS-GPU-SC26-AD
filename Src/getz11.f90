!***********************************************************************
      module m_getz11
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/25
!     Modification: 1999/04/06, 1999/05/10, 1999/11/01, 2000/01/17,
!                   2002/04/02, 2003/04/30, 2003/05/19, 2005/08/05,
!                   2007/01/31, 2007/06/27, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the zeta coordinates.

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

      public :: getz11, s_getz11

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getz11

        module procedure s_getz11

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getz11(fpdz,fpzsfc,nk,z)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: fpzsfc
                       ! Formal parameter of unique index of zsfc

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Output variable

      real, intent(out) :: z(1:nk)
                       ! zeta coordinates

! Internal shared variables

      real dz          ! Grid distance in z direction

      real zsfc        ! Sea surface terrain height

      real zsfc11      ! zsfc + 11.0

! Internal private variable

      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpdz,dz)
      call getrname(fpzsfc,zsfc)

! -----

! Set the common used variable.

      zsfc11=zsfc+11.e0

! -----

! Calculate the zeta coordinates.

!@llm start meta_info ----------------------------------------------------
! Location: getz11.f90 :: s_getz11
! Summary : Calculate 1D zeta vertical coordinates with 11m offset
!           from sea surface height for surface layer reference.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Uses intrinsic real function (GPU compatible)
!   - Simple 1D loop with arithmetic: z = zsfc11 + (k-2)*dz
!   - No global writes, only output array z is modified
! Next:
!   - 1D array with nk elements (typically small, <100)
!   - May not benefit from GPU offload due to small size
!   - If needed, use OpenACC with single team
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('getz11.f90', 's_getz11', &
   & 'OMP section 1')
end if
loop_len = int((nk)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=1,nk
        z(k)=zsfc11+real(k-2)*dz
      end do

!$omp end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_getz11

!-----7--------------------------------------------------------------7--

      end module m_getz11

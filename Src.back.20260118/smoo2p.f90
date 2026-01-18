!***********************************************************************
      module m_smoo2p
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/09/16
!     Modification: 1999/09/30, 1999/10/12, 1999/11/01, 1999/11/19,
!                   1999/11/24, 2000/01/17, 2001/04/15, 2001/06/06,
!                   2001/11/20, 2002/04/02, 2003/01/04, 2003/04/30,
!                   2003/05/19, 2003/07/15, 2003/11/28, 2004/03/05,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the 2nd order smoothing for the pressure.

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

      public :: smoo2p, s_smoo2p

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface smoo2p

        module procedure s_smoo2p

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
      subroutine s_smoo2p(fpsmhcoe,fpsmvcoe,ni,nj,nk,pp,pfrc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsmhcoe
                       ! Formal parameter of unique index of smhcoe

      integer, intent(in) :: fpsmvcoe
                       ! Formal parameter of unique index of smvcoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

! Input and output variable

      real, intent(inout) :: pfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

! Internal shared variables

      real smhcoe      ! Horizontal smoothig coefficient
      real smvcoe      ! Vertical smoothig coefficient

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpsmhcoe,smhcoe)
      call getrname(fpsmvcoe,smvcoe)

! -----

! Calculate the 2nd order pressure numerical smoothing.

!@llm start meta_info ----------------------------------------------------
! Location: smoo2p.f90 :: s_smoo2p
! Summary : Applies 2nd order numerical smoothing to pressure field using
!           horizontal and vertical smoothing coefficients in a 7-point stencil.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple stencil operation with neighbor access (i+/-1, j+/-1, k+/-1)
!   - All loops independent with private i,j,k and local temporary a
!   - No synchronization constructs
!   - Read-only access to pp array, update to pfrc array
! Next:
!   - Straightforward GPU port with 3D kernel
!   - Use OpenACC/OpenACC with collapse(3)
!   - Good memory access pattern for GPU (regular stencil)
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('smoo2p.f90', 's_smoo2p', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) &
     & * int((nj-2)-(2)+1,8) &
     & * int((ni-2)-(2)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-2
          a=2.e0*pp(i,j,k)

          pfrc(i,j,k)=pfrc(i,j,k)                                       &
     &      +(smvcoe*((pp(i,j,k+1)+pp(i,j,k-1))-a)                      &
     &      +smhcoe*(((pp(i+1,j,k)+pp(i-1,j,k))-a)                      &
     &      +((pp(i,j+1,k)+pp(i,j-1,k))-a)))

        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_smoo2p

!-----7--------------------------------------------------------------7--

      end module m_smoo2p

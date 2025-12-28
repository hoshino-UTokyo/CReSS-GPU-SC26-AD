!***********************************************************************
      module m_smoo2qv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/09/16
!     Modification: 1999/09/30, 1999/10/12, 1999/11/01, 1999/11/19,
!                   1999/11/24, 2000/01/17, 2001/04/15, 2001/06/06,
!                   2001/11/20, 2002/04/02, 2003/01/04, 2003/04/30,
!                   2003/05/19, 2003/07/15, 2003/11/28, 2003/12/12,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the 2nd order smoothing for the water vapor mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: smoo2qv, s_smoo2qv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface smoo2qv

        module procedure s_smoo2qv

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
      subroutine s_smoo2qv(fpsmhcoe,fpsmvcoe,                           &
     &                     ni,nj,nk,qvbr,rbr,qv,qvfrc,rbrqv)
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

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

! Input and output variable

      real, intent(inout) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in
                       ! water vapor mixing ratio equation

! Internal shared variables

      real smhcoe      ! Horizontal smoothig coefficient
      real smvcoe      ! Vertical smoothig coefficient

      real, intent(inout) :: rbrqv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x (qv - qvbr)

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpsmhcoe,smhcoe)
      call getrname(fpsmvcoe,smvcoe)

! -----

! Calculate the 2nd order smoothing.

!@llm start meta_info ----------------------------------------------------
! Location: smoo2qv.f90 :: s_smoo2qv
! Summary : Applies 2nd order numerical smoothing to water vapor mixing ratio
!           using density-weighted perturbation and 7-point stencil.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Two-phase computation: first computes rbrqv, then applies stencil
!   - Simple stencil operation with neighbor access (i+/-1, j+/-1, k+/-1)
!   - All loops independent with private i,j,k and local temporary a
!   - No synchronization constructs (implicit barrier between phases)
!   - Intermediate array rbrqv used between computation phases
! Next:
!   - GPU port needs to respect phase ordering (compute rbrqv first)
!   - Use OpenACC/OpenACC with collapse(3) for each phase
!   - Consider explicit barrier or separate kernels for two phases
!@llm end meta_info ------------------------------------------------------
!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          rbrqv(i,j,k)=rbr(i,j,k)*(qv(i,j,k)-qvbr(i,j,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-2
          a=2.e0*rbrqv(i,j,k)

          qvfrc(i,j,k)=qvfrc(i,j,k)                                     &
     &      +(smvcoe*((rbrqv(i,j,k+1)+rbrqv(i,j,k-1))-a)                &
     &      +smhcoe*(((rbrqv(i+1,j,k)+rbrqv(i-1,j,k))-a)                &
     &      +((rbrqv(i,j+1,k)+rbrqv(i,j-1,k))-a)))

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_smoo2qv

!-----7--------------------------------------------------------------7--

      end module m_smoo2qv

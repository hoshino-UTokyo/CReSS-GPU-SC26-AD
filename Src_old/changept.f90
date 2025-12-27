!***********************************************************************
      module m_changept
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/20
!     Modification: 2010/12/17, 2010/12/21

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     change potential temperature perturbation.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: changept, s_changept

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface changept

        module procedure s_changept

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
      subroutine s_changept(dtb,ni,nj,nk,htrsd,htrsu,htrld,htrlu,ptp)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: htrsd(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by downward short wave radiation

      real, intent(in) :: htrsu(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by upward short wave radiation

      real, intent(in) :: htrld(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by downward long wave radiation

      real, intent(in) :: htrlu(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by upward long wave radiation

! Input and output variable

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      nkm1=nk-1
      nkm2=nk-2

! -----

! Change potential temperature perturbation.

      do k=2,nk-2
      do j=1,nj-1
      do i=1,ni-1
        ptp(i,j,k)=ptp(i,j,k)                                           &
     &    +((htrsd(i,j,k)+htrsu(i,j,k))+(htrld(i,j,k)+htrlu(i,j,k)))*dtb
      end do
      end do
      end do

      do j=1,nj-1
      do i=1,ni-1
        ptp(i,j,1)=ptp(i,j,1)                                           &
     &    +((htrsd(i,j,2)+htrsu(i,j,2))+(htrld(i,j,2)+htrlu(i,j,2)))*dtb

        ptp(i,j,nkm1)=ptp(i,j,nkm1)+((htrsd(i,j,nkm2)+htrsu(i,j,nkm2))  &
     &    +(htrld(i,j,nkm2)+htrlu(i,j,nkm2)))*dtb

      end do
      end do

! -----

      end subroutine s_changept

!-----7--------------------------------------------------------------7--

      end module m_changept

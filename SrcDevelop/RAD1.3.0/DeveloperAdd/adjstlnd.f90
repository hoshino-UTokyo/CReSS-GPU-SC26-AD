!***********************************************************************
      module m_adjstlnd
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/11/17
!     Modification: 2010/12/17

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     adjust the land use categories to mstanx radiation scheme.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: adjstlnd, s_adjstlnd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstlnd

        module procedure s_adjstlnd

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
      subroutine s_adjstlnd(i,j,aland,ni,nj,land)
!***********************************************************************

! Input variables

      integer, intent(in) :: i
                       ! Array index in x direction

      integer, intent(in) :: j
                       ! Array index in y direction

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

! Input and output variable

      real, intent(out) :: aland
                       ! Adjusted real land use category

!-----7--------------------------------------------------------------7--

! Adjust the land use categories to mstanx radiation scheme.

      if(land(i,j).lt.0) then
        aland=1.1e0
      else

        if(land(i,j).lt.10) then

          if(land(i,j).lt.5) then
            aland=7.1e0
          else
            aland=6.1e0
          end if

        else

          if(land(i,j).eq.13) then
            aland=5.1e0
          else
            aland=4.1e0
          end if

        end if

      end if

! -----

      end subroutine s_adjstlnd

!-----7--------------------------------------------------------------7--

      end module m_adjstlnd

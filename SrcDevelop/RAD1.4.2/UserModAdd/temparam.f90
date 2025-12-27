!***********************************************************************
      module m_temparam
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2008/10/10
!     Modification: 2009/02/27, 2010/05/17, 2010/09/22, 2011/01/14,
!                   2011/06/01, 2011/09/22, 2011/11/10, 2013/02/11,
!                   2013/10/08

!     Author      : Hasegawa Koichi
!     Modification: 2023/02/07

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the constants for the exclusive use.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      public

! Exceptional access control

!     none

!-----7--------------------------------------------------------------7--

! Module variables

      integer, parameter :: elemnt_opt=0
                       ! Option for
                       ! estimeting elementary process of cloud physics

      integer, parameter :: evapor_opt=1
                       ! Option for evaporation rate

      integer, parameter :: nuc2nd_opt=1
                       ! Option for secondary nucleation rate

      integer, parameter :: flout_opt=1
!     integer, parameter :: flout_opt=0
                       ! Option to output falling velocity

!ORIG integer, parameter :: flqcqi_opt=22
!ORIG                  ! Option for additional falling term of 
!ORIG                  ! cloud_ice and cloud_water.
!ORIG                  ! [a] no fall of qx
!ORIG                  ! [b] constant terminal velocity
!ORIG                  ! [c] terminal velocity is depend on qx, Nx
!ORIG                  ! flqcqi_opt= 0: cloud_ice:[a], cloud_water:[a]
!ORIG                  ! flqcqi_opt= 1: cloud_ice:[a], cloud_water:[b]
!ORIG                  ! flqcqi_opt= 2: cloud_ice:[a], cloud_water:[c]
!ORIG                  ! flqcqi_opt=10: cloud_ice:[b], cloud_water:[a]
!ORIG                  ! flqcqi_opt=11: cloud_ice:[b], cloud_water:[b]
!ORIG                  ! flqcqi_opt=12: cloud_ice:[b], cloud_water:[c]
!ORIG                  ! flqcqi_opt=20: cloud_ice:[c], cloud_water:[a]
!ORIG                  ! flqcqi_opt=21: cloud_ice:[c], cloud_water:[b]
!ORIG                  ! flqcqi_opt=22: cloud_ice:[c], cloud_water:[c]
!ORIG       
!ORIG real, parameter :: ucqcst=.1e0
!ORIG                  ! Constant terminal velocity of
!ORIG                  ! cloud water mixing ratio

!ORIG real, parameter :: ucncst=.1e0
!ORIG                  ! Constant terminal velocity of
!ORIG                  ! cloud water concentrations

!ORIG real, parameter :: uiqcst=.1e0
!ORIG                  ! Constant terminal velocity of
!ORIG                  ! cloud ice mixing ratio

!ORIG real, parameter :: uincst=.1e0
!ORIG                  ! Constant terminal velocity of
!ORIG                  ! cloud ice concentrations

      real, parameter :: qpmin=2.778e-4
                       ! Lower limit of precipitation mixing ratio

!ORIG real, parameter :: qpmin=2.778e-5
!ORIG                  ! Lower limit of precipitation mixing ratio

!ORIG real, parameter :: qpmin=2.778e-6
!ORIG                  ! Lower limit of precipitation mixing ratio

!ORIG real, parameter :: qpmin=0.e0
!ORIG                  ! Lower limit of precipitation mixing ratio

      real, parameter :: adjqv=1.e0
                       ! Adjustment coefficient
                       ! of water vapor mixing ratio to radar

      real, parameter :: qvtop=14400.e0
                       ! Highest z physical coordinates of
                       ! water vapor nudging

      real, parameter :: rhqp=.9e0
                       ! Relative humidity under rain fall

! Module procedure

!     none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

!     none

!-----7--------------------------------------------------------------7--

      end module m_temparam

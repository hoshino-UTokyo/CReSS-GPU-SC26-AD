!***********************************************************************
      module m_commstrn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/08
!     Modification: 2010/12/17

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the parameters for mstranx radiation scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      public

! Exceptional access control

      private :: r8, selected_real_kind

!-----7--------------------------------------------------------------7--

! Module variables

      integer, parameter :: kln=100
                       ! Maximum number of layers

      integer, parameter :: kln1=kln+1
                       ! Maximum number of boundaries

      integer, parameter :: kmol=7
                       ! Number of gaserous spicies

      integer, parameter :: kpcl=11
                       ! Number of aerosol spicies

      integer, parameter :: kpclc=3
                       ! Aerosol distribution data

      integer, parameter :: kprg=5
                       ! Surface parameter

      integer, parameter :: kvar=6
                       ! Optical variables

      integer, parameter :: kbnd=29
                       ! Number of bands

      integer, parameter :: kbnd1=kbnd+1
                       ! Number of bands + 1

      integer, parameter :: kp=26
                       ! Number of pressure

      integer, parameter :: kt=3
                       ! Number of temperature

      integer, parameter :: kch=10
                       ! Number of subintervals

      integer, parameter :: ksfc=7
                       ! Number of surface spicies

      integer, parameter :: kmax2=6
                       ! Number of optical valuables

      integer, parameter :: kaplk=5
                       ! Number of dimension for planck function

      integer, parameter :: kplnk=3
                       ! Number of order for plank function

      integer, parameter :: kabs=3
                       ! Maximum number of absorption in each band

      integer, parameter :: kcfc=28
                       ! Number of CFC species

      integer, parameter :: kflg=7
                       ! Number of flags about bands

      integer, parameter :: iatm=1
                       ! Number of reading profile

      real(kind=r8), parameter :: pstd=1013.25e0
      real(kind=r8), parameter :: tstd=273.15e0
      real(kind=r8), parameter :: avg=6.0221367e23

      real(kind=r8), parameter :: rgas=8.31451e0
      real(kind=r8), parameter :: grav=9.80665e0
      real(kind=r8), parameter :: airm=29.e0

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

      end module m_commstrn

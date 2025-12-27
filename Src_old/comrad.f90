!***********************************************************************
      module m_comrad
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/08
!     Modification: 2010/12/17, 2010/12/21

!     Author      : Hasegawa Koichi
!     Modification: 2014/06/10, 2014/07/01, 2015/07/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the parameter and the array and variable
!     for radiation scheme.

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

      integer, parameter :: iug=81
                       ! Device unit number to read PARAG.${cbnd}

      integer, parameter :: iup=82
                       ! Device unit number to read PARAPC.${cbnd}

      integer, parameter :: iuv=83
                       ! Device unit number to read VARDATA.RM${cbnd}

      integer, parameter :: iud=84
                       ! Device unit number to read DATA.${ca}

      integer, parameter :: iuo=85
                       ! Device unit number to output data

      integer, parameter :: iuf=86
                       ! NO use

      integer, parameter :: iuh=87
                       ! NO use

      real(kind=r8) ams
                       ! Solar Zenith angle

      real(kind=r8) gtmp
                       ! Ground temperature

      real(kind=r8), allocatable, save :: prg(:)
                       ! Parameter packet for surface reflection

      real(kind=r8), allocatable, save :: zl(:)
                       ! z physical coordinates at CReSS layers

      real(kind=r8), allocatable, save :: pl(:)
                       ! Pressure at CReSS layers

      real(kind=r8), allocatable, save :: tl(:)
                       ! Temperature at CReSS layers

      real(kind=r8), allocatable, save :: pb(:)
                       ! Pressure at layer boundaries

      real(kind=r8), allocatable, save :: tb(:)
                       ! Temperature at layer boundaries

      real(kind=r8), allocatable, save :: zalt_data(:)
                       ! external data (altitude[km])

      real(kind=r8), allocatable, save :: patm_data(:)
                       ! external data (pressure[atm])

      real(kind=r8), allocatable, save :: tklv_data(:)
                       ! external data (temperature[K])

      real(kind=r8), allocatable, save :: cpcl_data(:,:,:)
                       ! external data (particulates concentrations)

      real(kind=r8), allocatable, save :: cgas_data(:,:)
                       ! external data (gas concentrations)

      real(kind=r8), allocatable, save :: cpcl(:,:,:)
                       ! Parameter packet for particulates

      real(kind=r8), allocatable, save :: gdcfrc(:)
                       ! Cloud cover late

      real(kind=r8), allocatable, save :: cgas(:,:)
                       ! Gas concentrations

      real(kind=r8), allocatable, save :: ccfc(:)
                       ! CFCs concentrations

      real(kind=r8), allocatable, save :: fd(:,:)
                       ! Downward flux at sublayer interface

      real(kind=r8), allocatable, save :: fu(:,:)
                       ! Upward flux at sublayer interface

      real, allocatable, save :: htrsd(:,:,:)
                       ! Heating rate by downward short wave radiation

      real, allocatable, save :: htrsu(:,:,:)
                       ! Heating rate by upward short wave radiation

      real, allocatable, save :: htrld(:,:,:)
                       ! Heating rate by downward long wave radiation

      real, allocatable, save :: htrlu(:,:,:)
                       ! Heating rate by upward long wave radiation

      real, allocatable, save :: rsnet(:,:)
                       ! Net short wave radiation (+:upward)

      real, allocatable, save :: rlnet(:,:)
                       ! Net long wave radiation (+:upward)

      integer, save :: iradcnt
                       ! Counter to run the radiation scheme

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

      end module m_comrad

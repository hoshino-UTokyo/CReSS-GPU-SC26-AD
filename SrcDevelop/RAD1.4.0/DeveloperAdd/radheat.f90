!***********************************************************************
      module m_radheat
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/20
!     Modification: 2010/12/17, 2010/12/21, 2010/12/27

!     Author      : Hasegawa Koichi
!     Modification: 2013/11/29, 2014/06/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get heating rate.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_comphy
      use m_commstrn

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: radheat, s_radheat

!-----7--------------------------------------------------------------7--

! Module variable

      real(kind=r8), parameter :: rr=8.31451e0
      real(kind=r8), parameter :: g0=9.80665e0
      real(kind=r8), parameter :: ambase=28.8e-3

      real(kind=r8), parameter ::                                       &
     &  am(1:kmol)=(/18.e-3,44.e-3,48.e-3,44.e-3,16.e-3,28.e-3,32.e-3/)

! Module procedure

      interface radheat

        module procedure s_radheat

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_radheat(i,j,ni,nj,nk,pb,cgas,                        &
     &                     fd,fu,htrsd,htrsu,htrld,htrlu,rsnet,rlnet)
!***********************************************************************

! Input variables

      integer, intent(in) :: i
                       ! Current array index in x direction

      integer, intent(in) :: j
                       ! Current array index in y direction

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real(kind=r8), intent(in) :: pb(1:nk-2)
                       ! Pressure at layer boundaries

      real(kind=r8), intent(in) :: cgas(1:nk-3,1:kmol)
                       ! Gas concentrations

      real(kind=r8), intent(in) :: fd(1:nk-2,1:2)
                       ! Downward flux at sublayer interface

      real(kind=r8), intent(in) :: fu(1:nk-2,1:2)
                       ! Upward flux at sublayer interface

! Output variables

      real, intent(out) :: htrsd(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by downward short wave radiation

      real, intent(out) :: htrsu(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by upward short wave radiation

      real, intent(out) :: htrld(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by downward long wave radiation

      real, intent(out) :: htrlu(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by upward long wave radiation

      real, intent(out) :: rsnet(0:ni+1,0:nj+1)
                       ! Net short wave radiation (+:upward)

      real, intent(out) :: rlnet(0:ni+1,0:nj+1)
                       ! Net long wave radiation (+:upward)

! Internal private variables

      integer k        ! Array index in z direction

      real(kind=r8) amair
                       ! Temporary variable

      real(kind=r8) pc
                       ! Pressure at grid center

      real(kind=r8) rddvcp
                       ! rd / cp

      real(kind=r8) cnv
                       ! Converter

      real(kind=r8) a  ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize variables

      rddvcp = real(rd,kind=r8)/real(cp,kind=r8)

! -----

! Get heating rate.

      do k=1,nk-3

        amair=ambase+1.e-6                                              &
     &    *(cgas(k,1)*(am(1)-ambase)+cgas(k,2)*(am(2)-ambase)           &
     &     +cgas(k,3)*(am(3)-ambase)+cgas(k,4)*(am(4)-ambase)           &
     &     +cgas(k,5)*(am(5)-ambase)+cgas(k,6)*(am(6)-ambase))

        a=-g0*amair/(3.5e2*rr*(pb(k)-pb(k+1)))

        !! for potential temperature -->

        pc  = 0.5d0*(pb(k)+pb(k+1))*1.d2           !! [hPa]=>[Pa]
        cnv = exp(rddvcp*log(real(p0,kind=r8)/pc))

        a = a * cnv

        !! <--

        htrsd(i,j,nk-k-1)=a*(fd(k,2)-fd(k+1,2))
        htrsu(i,j,nk-k-1)=a*(fu(k+1,2)-fu(k,2))

        htrld(i,j,nk-k-1)=a*(fd(k,1)-fd(k+1,1))
        htrlu(i,j,nk-k-1)=a*(fu(k+1,1)-fu(k,1))

      end do

      !! copy

      htrsd(i,j,nk-2) = htrsd(i,j,nk-3)
      htrsu(i,j,nk-2) = htrsu(i,j,nk-3)
      htrld(i,j,nk-2) = htrld(i,j,nk-3)
      htrlu(i,j,nk-2) = htrlu(i,j,nk-3)

! -----

! Get net flux to surface (+:upward) for calculating
! surface temperature.

      rlnet(i,j) = real(fu(nk-2,1)-fd(nk-2,1))
      rsnet(i,j) = real(fu(nk-2,2)-fd(nk-2,2))

! -----

      end subroutine s_radheat

!-----7--------------------------------------------------------------7--

      end module m_radheat

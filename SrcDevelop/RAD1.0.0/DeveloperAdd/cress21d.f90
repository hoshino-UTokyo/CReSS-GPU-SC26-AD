!***********************************************************************
      module m_cress21d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/20
!     Modification: 2010/12/17, 2013/02/13, 2013/10/08

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the one dimentional variables for mstranx radiation scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commstrn
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: cress21d, s_cress21d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface cress21d

        module procedure s_cress21d

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
      subroutine s_cress21d(fpcphopt,i,j,ni,nj,nk,nqw,nqi,zph,rbr,      &
     &                      p,t,qv,qwtr,qice,cpcl_data,gdcfrc_data,     &
     &                      cgas_data,zl,pl,tl,pb,tb,cpcl,gdcfrc,cgas)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

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

      integer, intent(in) :: nqw
                       ! Number of water hydrometeor array

      integer, intent(in) :: nqi
                       ! Number of ice hydrometeor array

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Temperature

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

      real(kind=r8), intent(in) :: cpcl_data(1:kln,1:kpcl,1:kpclc)
                       ! Parameter packet for particulates in sublayers

      real(kind=r8), intent(in) :: gdcfrc_data(1:kln)
                       ! Cloud cover late in sublayers

      real(kind=r8), intent(in) :: cgas_data(1:kln,1:kmol)
                       ! Gas concentrations in sublayers

! Output variables

      real(kind=r8), intent(out) :: zl(1:nk-3)
                       ! z physical coordinates at CReSS layers

      real(kind=r8), intent(out) :: pl(1:nk-3)
                       ! Pressure at CReSS layers

      real(kind=r8), intent(out) :: tl(1:nk-3)
                       ! Temperature at CReSS layers

      real(kind=r8), intent(out) :: pb(1:nk-2)
                       ! Pressure at layer boundaries

      real(kind=r8), intent(out) :: tb(1:nk-2)
                       ! Temperature at layer boundaries

      real(kind=r8), intent(out) :: cpcl(1:nk-3,1:kpcl,1:kpclc)
                       ! Parameter packet for particulates

      real(kind=r8), intent(out) :: gdcfrc(1:nk-3)
                       ! Cloud cover late

      real(kind=r8), intent(out) :: cgas(1:nk-3,1:kmol)
                       ! Gas concentrations

! Internal private variable

      integer cphopt   ! Option for cloud micro physics

      integer k        ! Array index in z direction

      integer kl       ! Index in layers

      integer ipcl     ! Index for kpcl
      integer ipclc    ! Index for kpclc

      integer imol     ! Index for kmol

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpcphopt,cphopt)

! -----

! Interpolate sublayers to CReSS layers.

      do k=2,nk-2
        zl(nk-k-1)=.0005e0*(zph(i,j,k)+zph(i,j,k+1))

        pl(nk-k-1)=1.e-2*p(i,j,k)
        tl(nk-k-1)=t(i,j,k)

      end do

      do k=2,nk-1
        pb(nk-k)=.5e-2*(p(i,j,k-1)+p(i,j,k))
        tb(nk-k)=.5e0*(t(i,j,k-1)+t(i,j,k))
      end do

      do k=1,nk-3
        kl=100-int(zl(k))

        if(kl.lt.1) then
          kl=1
        end if

        if(kl.gt.100) then
          kl=100
        end if

        do ipclc=1,kpclc
        do ipcl=1,kpcl
          cpcl(k,ipcl,ipclc)=cpcl_data(kl,ipcl,ipclc)
        end do
        end do

        gdcfrc(k)=gdcfrc_data(kl)

        do imol=1,kmol
          cgas(k,imol)=cgas_data(kl,imol)
        end do

      end do

      do k=2,nk-2
        cgas(nk-k-1,1)=1.e6*rbr(i,j,k)*qv(i,j,k)
      end do

      if(abs(cphopt).ge.1) then

        do k=2,nk-2
          cpcl(nk-k-1,1,1)=1.e6*rbr(i,j,k)*qwtr(i,j,k,1)
        end do

      end if

      if(abs(cphopt).ge.2) then

        do k=2,nk-2
          cpcl(nk-k-1,2,1)=1.e6*rbr(i,j,k)*(qice(i,j,k,1)+qice(i,j,k,2))
        end do

      end if

! -----

      end subroutine s_cress21d

!-----7--------------------------------------------------------------7--

      end module m_cress21d

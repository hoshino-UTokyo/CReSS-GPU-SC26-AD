!***********************************************************************
      module m_rrtmdrv
!***********************************************************************

!     Author      : Hasegawa Koichi
!     Date        : 2013/09/18
!     Modification: 2013/11/29

!     Author      : Oda Naotaka
!     Modification: 2013/11/12

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     run RRTM radiation scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy, only : cp, epsav, p0, rd, sun0
      use m_comrrtm
      use parkind, only : im => kind_im, rb => kind_rb

      use m_chkerr
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_setcst2d
      use m_setcst3d

      use rrtmg_lw_init
      use rrtmg_lw_rad
      use rrtmg_sw_init
      use rrtmg_sw_rad

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rrtmdrv, s_rrtmdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rrtmdrv

        module procedure s_rrtmdrv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rrtmdrv(ni,nj,nk,nqw,nqi,nund,                       &
     &                     htrsd,htrsu,htrld,htrlu,                     &
     &                     p,t,qv,tund,qwtr,qice,zph,albe,coseta)
!***********************************************************************

! Use parameter

      use parrrsw, only: nbndsw, naerec
      use parrrtm, only: nbndlw

! Input variables

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

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Temperature

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mass mixing ratio

      real, intent(in) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Ground temperature

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: albe(0:ni+1,0:nj+1)
                       ! Albedo

      real, intent(in) :: coseta(0:ni+1,0:nj+1)
                       ! Cosine of solar zenith angle

! Output variable

      real, intent(out) :: htrsd(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by downward short wave
                       ! radiation [K/s]

      real, intent(out) :: htrsu(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by upward short wave radiation
                       ! radiation [K/s]

      real, intent(out) :: htrld(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by downward long wave radiation
                       ! radiation [K/s]

      real, intent(out) :: htrlu(0:ni+1,0:nj+1,1:nk)
                       ! Heating rate by upward long wave radiation
                       ! radiation [K/s]

! Internal private variables

      integer(kind=im) ncol
                       ! Number of horizontal columns

      integer(kind=im) nlay
                       ! Number of model layers

      integer(kind=im) icld
                       ! Clear/cloud and cloud overlap flag [0,1,2,3]

      integer(kind=im) idrv
                       ! Flag for claculation of dFdT [0,1]

      integer(kind=im) inflglw
                       ! Flag [0,1,2] for cloud property method

      integer(kind=im) iceflglw
                       ! Flag [0,1,2,3] for ice cloud properties

      integer(kind=im) liqflglw
                       ! Flag [0,1] for liquid cloud properties

      integer(kind=im) inflgsw
                       ! Flag [0,1,2] for cloud property method

      integer(kind=im) iceflgsw
                       ! Flag [0,1,2,3] for ice cloud properties

      integer(kind=im) liqflgsw
                       ! Flag [0,1] for liquid cloud properties

      integer(kind=im) dyofyr
                       ! Day of the year (used to get Earth/
                       ! Sun distance if adjflx not provided)

      real(kind=rb) adjes
                       ! Flux adjustment for Earth/Sun distance

      real(kind=rb) tv 
                       ! virtual temperature

      real(kind=rb) rhod
                       ! density of dry air

      real(kind=rb) dz
                       ! grid interval in vertical

      real(kind=rb) cpdair
                       ! specific heat of dry air at constant pressure

      real(kind=rb) scon
                       ! sun light constant

      real rddvcp      ! rd / cp

! Internal private variables

      logical :: isw = .true.
      logical :: ilw = .true.
                       ! Flag for radiation scheme

      logical :: isw0 = .false.
      logical :: ilw0 = .false.
                       ! Flag for initial setting

      real :: cnv      ! converter

      integer stat     ! Runtime status
      integer cstat    ! Runtime status at current allocate statement

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction
      integer m        ! Array index

!-----7--------------------------------------------------------------7--

! Initialize.

      ncol = (ni-1) * (nj-1)
      nlay = nk-3

      cpdair = real(cp,kind=rb)
      scon   = real(sun0,kind=rb)

      rddvcp = rd/cp

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,htrsd)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,htrsu)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,htrld)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,htrlu)

! -----

! Allocate buffer for longwave/shortwave scheme

      stat = 0

      !## longwave/shortwave ##!

      if ( .not. allocated(play) ) then

        allocate(play(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(plev(1:ncol,1:nlay+1),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tlay(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tlev(1:ncol,1:nlay+1),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tsfc(1:ncol),stat=cstat)
        stat = stat + abs(cstat)
        allocate(h2ovmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(o3vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(co2vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(ch4vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(n2ovmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(o2vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(cldfr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(cicewp(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(cliqwp(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(reice(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(reliq(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tmp1(1:ncol,1:nlay+1),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tmp2(1:ncol,1:nlay+1),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tmp3(1:ncol,1:nlay+1),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tmp4(1:ncol,1:nlay+1),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tmp5(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)

      end if

      !## shortwave ##!

      if ( .not. allocated(asdir) .and. isw ) then

        isw0 = .true.

        allocate(asdir(1:ncol),stat=cstat)
        stat = stat + abs(cstat)
        allocate(asdif(1:ncol),stat=cstat)
        stat = stat + abs(cstat)
        allocate(aldir(1:ncol),stat=cstat)
        stat = stat + abs(cstat)
        allocate(aldif(1:ncol),stat=cstat)
        stat = stat + abs(cstat)
        allocate(coszen(1:ncol),stat=cstat)
        stat = stat + abs(cstat)
        allocate(taucldsw(1:nbndsw,1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(ssacld(1:nbndsw,1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(asmcld(1:nbndsw,1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(fsfcld(1:nbndsw,1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tauaersw(1:ncol,1:nlay,nbndsw),stat=cstat)
        stat = stat + abs(cstat)
        allocate(ssaaer(1:ncol,1:nlay,1:nbndsw),stat=cstat)
        stat = stat + abs(cstat)
        allocate(asmaer(1:ncol,1:nlay,1:nbndsw),stat=cstat)
        stat = stat + abs(cstat)
        allocate(ecaer(1:ncol,1:nlay,1:naerec),stat=cstat)
        stat = stat + abs(cstat)
        allocate(swhr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)

      else

        isw0 = .false.

      end if

      !## longwave ##!

      if ( .not. allocated(cfc11vmr) .and. ilw ) then

        ilw0 = .true.

        allocate(cfc11vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(cfc12vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(cfc22vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(ccl4vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(emis(1:ncol,1:nbndlw),stat=cstat)
        stat = stat + abs(cstat)
        allocate(taucld(1:nbndlw,1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(tauaer(1:ncol,1:nlay,1:nbndlw),stat=cstat)
        stat = stat + abs(cstat)
        allocate(lwhr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(duflx_dt(1:ncol,1:nlay+1),stat=cstat)
        stat = stat + abs(cstat)
        allocate(duflxc_dt(1:ncol,1:nlay+1),stat=cstat)
        stat = stat + abs(cstat)

      else

        ilw0 = .false.

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat<0) then

        if(mype.eq.-stat-1) then

          call destroy('rrtmdrv ',7,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('rrtmdrv ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

! Input buffer

      do k=2,nk-2
      do j=1,nj-1
      do i=1,ni-1

        m = i + (j-1)*(ni-1)

        !## pressure[hPa] at scalar point ##!

        play(m,k-1) = p(i,j,k) * 1.d-2

        !## temperature[K] at scalar point ##!

        tlay(m,k-1) = t(i,j,k)

        !## water vapor mixing ratio [m3/m3]                   ##!
        !## convert mass ratio [kg/kg] to volume ratio [m3/m3] ##!

        h2ovmr(m,k-1)   = qv(i,j,k) * 1.607793d0
        o3vmr(m,k-1)    = 0.e0
        co2vmr(m,k-1)   = 0.e0
        ch4vmr(m,k-1)   = 0.e0
        n2ovmr(m,k-1)   = 0.e0
        o2vmr(m,k-1)    = 0.e0

!ORG    o3vmr  = o3vmr  * 0.603428d0
!ORG    co2vmr = co2vmr * 0.658114d0
!ORG    ch4vmr = ch4vmr * 1.805423d0
!ORG    n2ovmr = n2ovmr * 0.658090d0
!ORG    o2vmr  = o2vmr  * 0.905140d0

        !## cloud fraction ##!

        if ( qwtr(i,j,k,1)+qice(i,j,k,1)+qice(i,j,k,2)>=5.d-5 ) then
          cldfr(m,k-1) = 1.d0
        else
          cldfr(m,k-1) = 0.d0
        endif

      end do
      end do
      end do

      !## variables at half level (w grid) ##!

      do j=1,nj-1
      do i=1,ni-1
        m = i + (j-1)*(ni-1)
        plev(m,1) = 0.5d0 * (p(i,j,1)+p(i,j,2)) * 1.d-2
        tlev(m,1) = 0.5d0 * (t(i,j,1)+t(i,j,2))
        tsfc(m)   = tund(i,j,1)
      end do
      end do

!oda  do k=2,nlay-1
!oda  do i=1,ncol
      do k=2,nlay
      do i=1,ncol
        plev(i,k) = 0.5d0 * (play(i,k)+play(i,k-1))
        tlev(i,k) = 0.5d0 * (tlay(i,k)+tlay(i,k-1))
      end do
      end do

!oda  ! extrapolation
!oda  do i=1,ncol
!oda    plev(i,nlay) = plev(i,nlay-1)+plev(i,nlay-1)-plev(i,nlay-2)
!oda    tlev(i,nlay) = tlev(i,nlay-1)+tlev(i,nlay-1)-tlev(i,nlay-2)
!oda  end do

      ! copy
      do i=1,ncol
        plev(i,nlay+1) = play(i,nlay)
        tlev(i,nlay+1) = tlay(i,nlay)
      end do

      !## cloud ice/liquid particle effective size (um) ##!

      call setcst2d(1,ncol,1,nlay,0.d0,reice)
      call setcst2d(1,ncol,1,nlay,0.d0,reliq)

      !## cloud ice/liquid water path [g/m2] ##!

      do k=2,nk-2
      do j=1,nj-1
      do i=1,ni-1

        m = (j-1)*(ni-1)+i

        tv   = t(i,j,k)*(1.d0+epsav*qv(i,j,k))/(1.d0+qv(i,j,k))
        rhod = p(i,j,k)/(tv*rd)/(1.d0+qv(i,j,k))
        dz   = zph(i,j,k+1)-zph(i,j,k)

        cliqwp(m,k-1) = dz * rhod * qwtr(i,j,k,1) * 1.d3
        cicewp(m,k-1) = dz * rhod * (qice(i,j,k,1)+qice(i,j,k,2)) * 1.d3

      end do
      end do
      end do

      do k=nlay-1,1,-1
      do i=1,ncol
        cliqwp(i,k) = cliqwp(i,k) + cliqwp(i,k+1)
        cicewp(i,k) = cicewp(i,k) + cicewp(i,k+1)
      end do
      end do

! -----

! RRTM shortwave scheme without McICA [v3.8]

      if(mype==root) write(*,*)

      if( isw ) then

        icld     = 0
        inflgsw  = 1
        iceflgsw = 1
        liqflgsw = 1

        dyofyr   = -1
        adjes    = 1.d0

        !## albedo, zenith angle ##!

        do j=1,nj-1
        do i=1,ni-1
          m = (j-1)*(ni-1)+i
          asdir(m)  = albe(i,j)
          aldir(m)  = albe(i,j)
          asdif(m)  = albe(i,j)
          aldif(m)  = albe(i,j)
          coszen(m) = coseta(i,j)
        end do
        end do

        !## variables about cloud ##!

        call setcst3d(1,nbndsw,1,ncol,1,nlay,0.d0,taucldsw)
        call setcst3d(1,nbndsw,1,ncol,1,nlay,0.d0,ssacld)
        call setcst3d(1,nbndsw,1,ncol,1,nlay,0.d0,asmcld)
        call setcst3d(1,nbndsw,1,ncol,1,nlay,0.d0,fsfcld)

        !## variables about aerosol ##!

        call setcst3d(1,ncol,1,nlay,1,nbndsw,0.d0,tauaersw)
        call setcst3d(1,ncol,1,nlay,1,nbndsw,0.d0,ssaaer)
        call setcst3d(1,ncol,1,nlay,1,nbndsw,0.d0,asmaer)
        call setcst3d(1,ncol,1,nlay,1,naerec,0.d0,asmaer)

        !## temporary ##!

        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp1)
        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp2)
        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp3)
        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp4)
        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp5)

        !## RRTM-shortwave scheme ##!

         if( mype==root ) then
           write(*,*)' ### run RRTM(sw) radiation scheme ###'
         end if

        if( isw0 ) then
          call rrtmg_sw_ini(cpdair)
          isw0 = .false.
        end if

        call rrtmg_sw(ncol    ,nlay    ,icld    ,                       &
     &                play    ,plev    ,tlay    ,tlev    ,tsfc  ,       &
     &                h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr,o2vmr, &
     &                asdir   ,asdif   ,aldir   ,aldif   ,              &
     &                coszen  ,adjes   ,dyofyr  ,scon    ,              &
     &                inflgsw ,iceflgsw,liqflgsw,cldfr   ,              &
     &                taucldsw,ssacld  ,asmcld  ,fsfcld  ,              &
     &                cicewp  ,cliqwp  ,reice   ,reliq   ,              &
     &                tauaersw,ssaaer  ,asmaer  ,ecaer   ,              &
     &                tmp1    ,tmp2    ,swhr    ,tmp3    ,tmp4  ,tmp5)

      end if

! -----

! RRTM longwave scheme without McICA [v4.85]

      if( ilw ) then

        icld     = 1
        idrv     = 0
        inflglw  = 1
        iceflglw = 1
        liqflglw = 1

        !## CFC ##!

        call setcst2d(1,ncol,1,nlay,0.d0,cfc11vmr)
        call setcst2d(1,ncol,1,nlay,0.d0,cfc12vmr)
        call setcst2d(1,ncol,1,nlay,0.d0,cfc22vmr)
        call setcst2d(1,ncol,1,nlay,0.d0,ccl4vmr)

!ORG    cfc11vmr = cfc11vmr * 0.210852d0
!ORG    cfc12vmr = cfc12vmr * 0.239546d0
!ORG    cfc22vmr = cfc22vmr * 0.33498788002d0
!ORG    ccl4vmr  = ccl4vmr  * 0.18831254103d0

        !## surface emissivity ##!

        do k=1,nbndlw
        do i=1,ncol
          emis(i,k) = 1.d0   !! as black-body radiation
        end do
        end do

        !## optical depth ##!

        call setcst3d(1,nbndlw,1,ncol,1,nlay,0.d0,taucld)
        call setcst3d(1,ncol,1,nlay,1,nbndlw,0.d0,tauaer)

        !## temporary ##!

        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp1)
        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp2)
        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp3)
        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp4)
        call setcst2d(1,ncol,1,nlay+1,0.d0,tmp5)

        !## RRTM-longwave scheme ##!

         if( mype==root ) then
           write(*,*)' ### run RRTM(lw) radiation scheme ###'
         end if

        if( ilw0 ) then
          call rrtmg_lw_ini(cpdair)
          ilw0 = .false.
        end if

        call rrtmg_lw(ncol    ,nlay    ,icld    ,idrv    ,              &
     &                play    ,plev    ,tlay    ,tlev    ,tsfc,         &
     &                h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr,o2vmr, &
     &                cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis  ,       &
     &                inflglw ,iceflglw,liqflglw,cldfr   ,              &
     &                taucld  ,cicewp  ,cliqwp  ,reice   ,reliq ,       &
     &                tauaer  ,                                         &
     &                tmp1    ,tmp2    ,lwhr    ,tmp3    ,tmp4  , tmp5, &
     &                duflx_dt,duflxc_dt)

      end if

      if(mype==root) write(*,*)

! -----

! Convert heating rate of temperature to rate of potential temperature.

      do k=2,nk-2
      do j=1,nj-1
      do i=1,ni-1

        cnv = exp(rddvcp*log(p0/p(i,j,k)))

        htrsd(i,j,k) = 0.e0
        htrld(i,j,k) = 0.e0

        htrsu(i,j,k) = cnv * real(swhr((j-1)*(ni-1)+i,k-1)/86400.d0)
        htrlu(i,j,k) = cnv * real(lwhr((j-1)*(ni-1)+i,k-1)/86400.d0)

      end do
      end do
      end do

! -----

      end subroutine s_rrtmdrv

!-----7--------------------------------------------------------------7--

      end module m_rrtmdrv

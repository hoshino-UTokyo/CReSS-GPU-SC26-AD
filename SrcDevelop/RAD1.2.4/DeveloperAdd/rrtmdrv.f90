!***********************************************************************
      module m_rrtmdrv
!***********************************************************************

!     Author      : Hasegawa Koichi
!     Date        : 2013/09/18
!     Modification: 2013/11/29, 2014/06/10, 2014/07/01, 2014/08/25,
!                   2015/06/08

!     Author      : Oda Naotaka
!     Modification: 2013/11/12, 2014/07/17

!     Author      : Mayumi Yoshioka
!     Modification: 2015/04/08

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     run RRTM radiation scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy, only : cp, epsav, p0, rd, rv, sun0
      use m_comrrtm
      use parkind, only : im => kind_im, rb => kind_rb

      use m_chkerr
      use m_comkind
      use m_commpi
      use m_commstrn
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_getrname
      use m_setcst2d
      use m_setcst3d
      use m_var8w8s
      use m_vint13

      use rrtmg_lw_init
      use rrtmg_lw_rad
      use rrtmg_lw_rad_mcica !oda
      use rrtmg_sw_init
      use rrtmg_sw_rad
      use rrtmg_sw_rad_mcica !oda

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
      subroutine s_rrtmdrv(radon,fpcphopt,fpradobj,fpmcaopt,fpovpopt,   &
     &                     fpcldlim,fpreliqum,fpreiceum,                &
     &                     icnt,ni,nj,nk,nqw,nqi,nund,                  &
     &                     htrsd,htrsu,htrld,htrlu,rsnet,rlnet,         &
     &                     p,t,qv,tund,qwtr,qice,zph,zph8s,             &
     &                     albe,coseta,cgas,                            &
     &                     nln,zalt_data,cgas_data,ccfc)
!***********************************************************************

! Use parameter

      use parrrsw, only: nbndsw, naerec
      use parrrtm, only: nbndlw

! Input variables

      character(len=6), intent(in) :: radon
                       ! Control flag of radiation scheme

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpradobj
                       ! Formal parameter of unique index of radobj

      integer, intent(in) :: fpmcaopt
                       ! Formal parameter of unique index of mcaopt

      integer, intent(in) :: fpovpopt
                       ! Formal parameter of unique index of ovpopt

      integer, intent(in) :: fpcldlim
                       ! Formal parameter of unique index of cldlim

      integer, intent(in) :: fpreliqum
                       ! Formal parameter of unique index of reliqum

      integer, intent(in) :: fpreiceum
                       ! Formal parameter of unique index of reiceum

      integer, intent(in) :: icnt
                       ! Counter

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

      integer, intent(in) :: nln
                       ! Layers number of external data

      real(kind=r8), intent(in) :: zalt_data(1:nln)
                       ! external data (altitude[km])

      real(kind=r8), intent(in) :: cgas_data(1:nln,1:kmol)
                       ! external data (gas concentrations)

      real(kind=r8), intent(in) :: ccfc(1:kcfc)
                       ! external data (CFCs concentrations)

! Input and output variables

      real, intent(inout) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(inout) :: cgas(0:ni+1,0:nj+1,1:nk)
                       ! gas conc.

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

      real, intent(out) :: rsnet(0:ni+1,0:nj+1)
                       ! Net short wave radiation (+:upward)

      real, intent(out) :: rlnet(0:ni+1,0:nj+1)
                       ! Net long wave radiation (+:upward)

! Internal private variables

      integer(kind=im) ncol
                       ! Number of horizontal columns

      integer(kind=im) nlay
                       ! Number of model layers

      integer(kind=im) imca
                       ! McICA scheme flag [0,1]
                       ! 0 = not use McICA, 1 = use McICA 

      integer(kind=im) icld
                       ! Clear/cloud and cloud overlap flag [0,1,2,3]
                       ! 0 = ignore cloud effect, 1 = random overlap
                       ! 2 = max/random,          3 = maximum overlap

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

      real(kind=rb) cldliq
                       ! total amount of liquid cloud

      real(kind=rb) cldice
                       ! total amount of ice cloud

      real(kind=rb) dz
                       ! grid interval in vertical

      real(kind=rb) cpdair
                       ! specific heat of dry air at constant pressure

      real(kind=rb) scon
                       ! sun light constant

      real rddvcp      ! rd / cp

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer radobj   ! Option for object substance in radiation
                       ! scheme
      integer mcaopt   ! Option for McICA method only at RRTM scheme
      integer ovpopt   ! Option for overlap method only at RRTM scheme

      real cldlim(1:2) ! Lower limit of mass mixing ratio [kg/kg]
                       ! to consider as cloudy grid
      real reliqum     ! Cloud water effective particle radius [um]
      real reiceum     ! Cloud ice effective particle radius [um]

! Internal private variables

      logical :: isw   ! Flag for radiation scheme
      logical :: ilw

      logical :: isw0 = .false.
      logical :: ilw0 = .false.
                       ! Flag for initial setting

      integer nsnd     ! External data number

      real dsnd(0:nln) ! Temporary
      real zsnd(0:nln)

      real(kind=rb) fct_qv
      real(kind=rb) fct_qc
      real(kind=rb) fct_qr
      real(kind=rb) fct_qi
      real(kind=rb) fct_qs
      real(kind=rb) fct_qg
      real(kind=rb) fct_other
                       ! factor

      real(kind=rb) uflx
                       ! upward net flux 

      real(kind=rb) dflx
                       ! downward net flux 

      character(len=108) ctmp1
      character(len=108) ctmp2
                       ! temporary

      real cnv         ! converter

      integer stat     ! Runtime status
      integer cstat    ! Runtime status at current allocate statement

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction
      integer m        ! Array index
      integer imol     ! Array index

!-----7--------------------------------------------------------------7--

! Initialize.

      if( trim(radon).eq.'lsw' ) then
        isw = .true.
        ilw = .true.
      else if( trim(radon).eq.'sw' ) then
        isw = .true.
        ilw = .false.
      else if( trim(radon).eq.'lw' ) then
        isw = .false.
        ilw = .true.
      else
        return
      end if
      
      ncol = (ni-1) * (nj-1)
      nlay = nk-3

      cpdair = real(cp,kind=rb)
      scon   = real(sun0,kind=rb)
      rddvcp = rd/cp

      call getiname(fpcphopt,cphopt)
      call getiname(fpradobj,radobj)
      call getiname(fpmcaopt,mcaopt)
      call getiname(fpovpopt,ovpopt)
      call getrname(fpcldlim,cldlim(1))
      call getrname(fpcldlim+1,cldlim(2))
      call getrname(fpreliqum,reliqum)
      call getrname(fpreiceum,reiceum)

      imca = mcaopt
      icld = ovpopt

! -----

! For message

      if( imca==0 ) then

        ctmp1 = 'without McICA'
        ctmp2 = ''
 
      else if( imca==1 ) then

        ctmp1 = 'with McICA'

        if( icld==0 ) then
          ctmp2 = '(clear sky)'
!!!       ctmp2 = '(clear sky, so ignore cloud effect)'
        else if( icld==1 ) then
          ctmp2 = '(random cloud overlap)'
        else if( icld==2 ) then
          ctmp2 = '(max-random cloud overlap)'
        else if( icld==3 ) then
          ctmp2 = '(maximum cloud overlap)'
        end if 

      end if

! -----

! Set factor

      if( radobj<=0 ) then
        fct_other = 0.d0
      else
        fct_other = 1.d0
      end if

      if( abs(radobj)==0 ) then
        fct_qv = 0.d0
        fct_qc = 0.d0
        fct_qr = 0.d0
      else if( abs(radobj)<10 ) then
        fct_qv = 1.d0
        fct_qc = 0.d0
        fct_qr = 0.d0
      else if( abs(radobj)<20 ) then
        fct_qv = 1.d0
        fct_qc = 1.d0
        fct_qr = 0.d0
      else
        fct_qv = 1.d0
        fct_qc = 1.d0
        fct_qr = 1.d0
      end if
 
      if( abs(radobj)<10 .or. mod(abs(radobj),10)==0 ) then
        fct_qi = 0.d0
        fct_qs = 0.d0
        fct_qg = 0.d0
      else if( mod(abs(radobj),10)==1) then
        fct_qi = 1.d0
        fct_qs = 0.d0
        fct_qg = 0.d0
      else if( mod(abs(radobj),10)==2) then
        fct_qi = 1.d0
        fct_qs = 1.d0
        fct_qg = 0.d0
      else
        fct_qi = 1.d0
        fct_qs = 1.d0
        fct_qg = 1.d0
      end if

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
        allocate(cfc11vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(cfc12vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(cfc22vmr(1:ncol,1:nlay),stat=cstat)
        stat = stat + abs(cstat)
        allocate(ccl4vmr(1:ncol,1:nlay),stat=cstat)
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
        allocate(tmp5(1:ncol,1:nlay+1),stat=cstat)
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

      if ( .not. allocated(emis) .and. ilw ) then

        ilw0 = .true.

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

!$omp parallel default(shared)
!$omp do schedule(runtime) private(i,j,k,m)

      do k=2,nk-2
      do j=1,nj-1
      do i=1,ni-1

        m = i + (j-1)*(ni-1)

        !## pressure[hPa] at scalar point ##!

        play(m,k-1) = p(i,j,k) * 1.d-2

        !## temperature[K] at scalar point ##!

        tlay(m,k-1) = t(i,j,k)

        !## water vapor conc.[m3/m3]                       ##!
        !##   : convert mass ratio [kg/kg] to volume ratio ##!

        h2ovmr(m,k-1) = dble(qv(i,j,k)) * 1.607793d0 * fct_qv

      end do
      end do
      end do

!$omp end do
!$omp end parallel

      !## gas conc.[m3/m3] from the external file only once ##!
      !##   : convert unit [ppmV] to unit [m3/m3]           ##!

      if(icnt==0) then

        call var8w8s(ni,nj,nk,zph,zph8s)

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,k,m,imol,cgas,nsnd,dsnd,zsnd)

        do imol=2,kmol

          nsnd=0
          do k=1,nln
            if(zalt_data(k)>30.d0) exit
            nsnd=nsnd+1
            zsnd(nsnd) = real(zalt_data(k))*1.e3
            dsnd(nsnd) = real(cgas_data(k,imol))
          end do

          zsnd(0) = min(2.e0*zsnd(1)-zsnd(2),0.e0)
          dsnd(0) = dsnd(1)

!!!       !! make a better modification than NEC by hasegawa  (2015/06/08)
!!!       !! modify for (detected heap, heap=nan) by NEC (2015/04/08)
!!!
!!!       nsnd=1
!!!       do k=1,nln
!!!         if(zalt_data(k)>30.d0) exit
!!!         zsnd(nsnd) = real(zalt_data(k))*1.e3
!!!         dsnd(nsnd) = real(cgas_data(k,imol))
!!!         nsnd=nsnd+1
!!!       end do
!!!
!!!       zsnd(0) = min(2.e0*zsnd(1)-zsnd(2),0.e0)
!!!       dsnd(0) = dsnd(1)
!!!
!!! NEC2015/04/08
!!!       zsnd(nsnd) = 0.0e0
!!!       dsnd(nsnd) = 0.0e0
!!! NEC2015/04/08

          call s_vint13('xx',0,ni+1,0,nj+1,2,nk-2,zph8s(0,0,2),         &
     &                  cgas(0,0,2),nsnd,zsnd,dsnd)

          do k=2,nk-2
          do j=1,nj-1
          do i=1,ni-1

            m = i + (j-1)*(ni-1)
            cgas(i,j,k) = cgas(i,j,k) * fct_other

            if(imol==2) co2vmr(m,k-1) = dble(cgas(i,j,k)) * 0.658114d-6
            if(imol==3) o3vmr(m,k-1)  = dble(cgas(i,j,k)) * 0.603428d-6
            if(imol==4) n2ovmr(m,k-1) = dble(cgas(i,j,k)) * 0.658090d-6
            if(imol==6) ch4vmr(m,k-1) = dble(cgas(i,j,k)) * 1.805423d-6
            if(imol==7) o2vmr(m,k-1)  = dble(cgas(i,j,k)) * 0.905140d-6

          end do
          end do
          end do

        end do

!$omp end do

!$omp do schedule(runtime) private(i,j,k,m)

        do k=2,nk-2
        do j=1,nj-1
        do i=1,ni-1
          m = i + (j-1)*(ni-1)
          cfc11vmr(m,k-1) = dble(ccfc( 1)) * 0.210852d-6 * fct_other
          cfc12vmr(m,k-1) = dble(ccfc( 2)) * 0.239546d-6 * fct_other
          cfc22vmr(m,k-1) = dble(ccfc( 9)) * 0.334987880d-6 * fct_other
          ccl4vmr(m,k-1)  = dble(ccfc(24)) * 0.188312541d-6 * fct_other
        end do
        end do
        end do

!$omp end do

!$omp end parallel

      end if

      !## variables at half level (w grid) ##!

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,m)

      do j=1,nj-1
      do i=1,ni-1
        m = i + (j-1)*(ni-1)
        plev(m,1) = 0.5d0 * (p(i,j,1)+p(i,j,2)) * 1.d-2
        tlev(m,1) = 0.5d0 * (t(i,j,1)+t(i,j,2))
        tsfc(m)   = tund(i,j,1)
!DBG    tsfc(m)   = 260.94d0   !! debug for mstrnX-CAOS
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,k)

      do k=2,nlay
      do i=1,ncol
        plev(i,k) = 0.5d0 * (play(i,k)+play(i,k-1))
        tlev(i,k) = 0.5d0 * (tlay(i,k)+tlay(i,k-1))
      end do
      end do

!$omp end do

!oda  !! extrapolation
!oda  do i=1,ncol
!oda    plev(i,nlay) = plev(i,nlay-1)+plev(i,nlay-1)-plev(i,nlay-2)
!oda    tlev(i,nlay) = tlev(i,nlay-1)+tlev(i,nlay-1)-tlev(i,nlay-2)
!oda  end do

      !! copy
!$omp do schedule(runtime) private(i)
      do i=1,ncol
        plev(i,nlay+1) = play(i,nlay)
        tlev(i,nlay+1) = tlay(i,nlay)
      end do
!$omp end do

!$omp end parallel

      !## cloud ice/liquid particle effective size (um) ##!

      call setcst2d(1,ncol,1,nlay,dble(reiceum),reice)
      call setcst2d(1,ncol,1,nlay,dble(reliqum),reliq)

      !## cloud ice/liquid water path [g/m2], cloud fraction [-] ##!

!$omp parallel default(shared)

      if(abs(cphopt)==1) then

        !! warm rain

!$omp do schedule(runtime) private(i,j,k,m,cldliq,cldice,dz,rhod,tv)

        do k=2,nk-2
        do j=1,nj-1
        do i=1,ni-1

          m = i + (j-1)*(ni-1)

          tv   = t(i,j,k)*(1.d0+epsav*qv(i,j,k))/(1.d0+qv(i,j,k))
          rhod = p(i,j,k)/(tv*rd)/(1.d0+qv(i,j,k))
          dz   = zph(i,j,k+1)-zph(i,j,k)

          cldliq = dble(qwtr(i,j,k,1))*fct_qc                           &
     &           + dble(qwtr(i,j,k,2))*fct_qr
          cldice = 0.d0

          if( imca==1 .and. icld>0 .and. cldlim(1)<cldlim(2) ) then

            if ( (cldliq+cldice)>=dble(cldlim(2)) ) then
              cldfr(m,k-1) = 1.d0
            else if ( (cldliq+cldice)<dble(cldlim(1)) ) then
              cldfr(m,k-1) = 0.d0
            else
              cldfr(m,k-1) = dble((cldliq+cldice)-cldlim(1))            &
     &                     / dble(cldlim(2)-cldlim(1))
            end if

          else

            if ( (cldliq+cldice)>=dble(cldlim(1)) ) then
              cldfr(m,k-1) = 1.d0
            else
              cldfr(m,k-1) = 0.d0
            end if

          end if

          cliqwp(m,k-1) = dz * rhod * cldliq * 1.d3
          cicewp(m,k-1) = dz * rhod * cldice * 1.d3

        end do
        end do
        end do

!$omp end do

      else if(abs(cphopt)>=2) then

        !! cold rain

!$omp do schedule(runtime) private(i,j,k,m,cldliq,cldice,dz,rhod,tv)

        do k=2,nk-2
        do j=1,nj-1
        do i=1,ni-1

          m = i + (j-1)*(ni-1)

          tv   = t(i,j,k)*(1.d0+epsav*qv(i,j,k))/(1.d0+qv(i,j,k))
          rhod = p(i,j,k)/(tv*rd)/(1.d0+qv(i,j,k))
          dz   = zph(i,j,k+1)-zph(i,j,k)

          cldliq = dble(qwtr(i,j,k,1))*fct_qc                           &
     &           + dble(qwtr(i,j,k,2))*fct_qr
          cldice = dble(qice(i,j,k,1))*fct_qi                           &
     &           + dble(qice(i,j,k,2))*fct_qs                           &
     &           + dble(qice(i,j,k,3))*fct_qg

!DBG      !!! ** debug for mstrnX-CAOS -->
!DBG
!DBG      cnv = 0.5e0*(zph(i,j,k+1)+zph(i,j,k))
!DBG
!DBG      !! low cloud
!DBG      if(cnv>0.0.and.cnv<=2000.0) then
!DBG!!!     cldliq = 6.7919d-03 * (rd/rv) * fct_qc   !! qc[ppmV]
!DBG        cldliq = 5.2556d-06 * fct_qc             !! qc[kg/kg]
!DBG        cldice = 0.d0
!DBG      else
!DBG        cldliq = 0.d0
!DBG        cldice = 0.d0
!DBG      end if
!DBG
!DBG      !! high cloud
!DBG      if(cnv>9000.0.and.cnv<10000.0) then
!DBG        cldliq = 0.d0
!DBG!!!     cldice = 1.3979d-02 * (rd/rv) * fct_qi   !! qi[ppmV]
!DBG        cldice = 1.0817d-05 * fct_qi             !! qi[kg/kg]
!DBG      else
!DBG        cldliq = 0.d0
!DBG        cldice = 0.d0
!DBG      end if
!DBG
!DBG      !!! <-- debug for mstrnX-CAOS **

          if( imca==1 .and. icld>0 .and. cldlim(1)<cldlim(2) ) then

            if ( (cldliq+cldice)>=dble(cldlim(2)) ) then
              cldfr(m,k-1) = 1.d0
            else if ( (cldliq+cldice)<dble(cldlim(1)) ) then
              cldfr(m,k-1) = 0.d0
            else
              cldfr(m,k-1) = dble((cldliq+cldice)-cldlim(1))            &
     &                     / dble(cldlim(2)-cldlim(1))
            end if

          else

            if ( (cldliq+cldice)>=dble(cldlim(1)) ) then
              cldfr(m,k-1) = 1.d0
            else
              cldfr(m,k-1) = 0.d0
            end if

          end if

          cliqwp(m,k-1) = dz * rhod * cldliq * 1.d3
          cicewp(m,k-1) = dz * rhod * cldice * 1.d3

        end do
        end do
        end do

!$omp end do

      end if

      !! integrate to calculate water path

!$omp do schedule(runtime) private(i,k)

      do k=nlay-1,1,-1
      do i=1,ncol
        cliqwp(i,k) = cliqwp(i,k) + cliqwp(i,k+1)
        cicewp(i,k) = cicewp(i,k) + cicewp(i,k+1)
      end do
      end do

!$omp end do

!$omp end parallel

! -----

! RRTM shortwave scheme without McICA [v3.8]

      if(mype==root) write(*,*)

      if( isw ) then

        inflgsw  = 2
        iceflgsw = 3
        liqflgsw = 1

        dyofyr = -1
        adjes  = 1.d0

        !## albedo, zenith angle ##!

!$omp parallel default(shared)
!$omp do schedule(runtime) private(i,j,m)

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

!$omp end do
!$omp end parallel

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
          write(*,*)' ### run RRTM(sw) radiation scheme ',              &
     &              trim(ctmp1),' ',trim(ctmp2),' ###'
        end if

        if( isw0 ) then
          call rrtmg_sw_ini(cpdair)
          isw0 = .false.
        end if

        if( imca==0 ) then

          icld = 1   !! '1' means consider cloud effect

          call rrtmg_sw(ncol  ,nlay    ,icld    ,                       &
     &                play    ,plev    ,tlay    ,tlev    ,tsfc  ,       &
     &                h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr,o2vmr, &
     &                asdir   ,asdif   ,aldir   ,aldif   ,              &
     &                coszen  ,adjes   ,dyofyr  ,scon    ,              &
     &                inflgsw ,iceflgsw,liqflgsw,cldfr   ,              &
     &                taucldsw,ssacld  ,asmcld  ,fsfcld  ,              &
     &                cicewp  ,cliqwp  ,reice   ,reliq   ,              &
     &                tauaersw,ssaaer  ,asmaer  ,ecaer   ,              &
     &                tmp1    ,tmp2    ,swhr    ,tmp3    ,tmp4  ,tmp5)

        else

          call rrtmg_sw_mcica(ncol     ,nlay    ,icld    ,              &
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


        !## total shortwave net flux (W/m2,+:upward) ##!

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,dflx,uflx)

        do j=1,nj-1
        do i=1,ni-1
          uflx = tmp1((j-1)*(ni-1)+i,1)
          dflx = tmp2((j-1)*(ni-1)+i,1)
          rsnet(i,j) = real(uflx-dflx)
        end do
        end do

!$omp end do

        !## convert heating rate of temperature to ##!
        !## rate of potential temperature.         ##!

!$omp do schedule(runtime) private(i,j,k,cnv)

        do k=2,nk-2
        do j=1,nj-1
        do i=1,ni-1
          cnv = exp(rddvcp*log(p0/p(i,j,k)))
          htrsd(i,j,k) = 0.e0
          htrsu(i,j,k) = cnv * real(swhr((j-1)*(ni-1)+i,k-1)/86400.d0)
        end do
        end do
        end do

!$omp end do

        !! copy

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          htrsu(i,j,nk-2) = htrsu(i,j,nk-3)
        end do
        end do

!$omp end do

!$omp end parallel

      end if

! -----

! RRTM longwave scheme without McICA [v4.85]

      if( ilw ) then

        idrv     = 0
        inflglw  = 2
        iceflglw = 3
        liqflglw = 1

        !## surface emissivity ##!

!$omp parallel default(shared)
!$omp do schedule(runtime) private(i,k)

        do k=1,nbndlw
        do i=1,ncol
          emis(i,k) = 1.d0   !! as black-body radiation
        end do
        end do

!$omp end do
!$omp end parallel

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
          write(*,*)' ### run RRTM(lw) radiation scheme ',              &
     &              trim(ctmp1),' ',trim(ctmp2),' ###'
        end if

        if( ilw0 ) then
          call rrtmg_lw_ini(cpdair)
          ilw0 = .false.
        end if

        if( imca==0 ) then

          icld = 1   !! '1' means consider cloud effect

          call rrtmg_lw(ncol  ,nlay    ,icld    ,idrv    ,              &
     &                play    ,plev    ,tlay    ,tlev    ,tsfc,         &
     &                h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr,o2vmr, &
     &                cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis  ,       &
     &                inflglw ,iceflglw,liqflglw,cldfr   ,              &
     &                taucld  ,cicewp  ,cliqwp  ,reice   ,reliq ,       &
     &                tauaer  ,                                         &
     &                tmp1    ,tmp2    ,lwhr    ,tmp3    ,tmp4  , tmp5, &
     &                duflx_dt,duflxc_dt)

        else

          call rrtmg_lw_mcica(ncol  ,nlay    ,icld    ,idrv    ,        &
     &                play    ,plev    ,tlay    ,tlev    ,tsfc,         &
     &                h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr,o2vmr, &
     &                cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis  ,       &
     &                inflglw ,iceflglw,liqflglw,cldfr   ,              &
     &                taucld  ,cicewp  ,cliqwp  ,reice   ,reliq ,       &
     &                tauaer  ,                                         &
     &                tmp1    ,tmp2    ,lwhr    ,tmp3    ,tmp4  , tmp5, &
     &                duflx_dt,duflxc_dt)

        end if

        !## total longwave net flux (W/m2,+:upward) ##!


        !## total longwave net flux (W/m2,+:upward) ##!

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,dflx,uflx)

        do j=1,nj-1
        do i=1,ni-1
          uflx = tmp1((j-1)*(ni-1)+i,1)
          dflx = tmp2((j-1)*(ni-1)+i,1)
          rlnet(i,j) = real(uflx-dflx)
        end do
        end do

!$omp end do

        !## convert heating rate of temperature to ##!
        !## rate of potential temperature.         ##!

!$omp do schedule(runtime) private(i,j,k,cnv)

        do k=2,nk-2
        do j=1,nj-1
        do i=1,ni-1
          cnv = exp(rddvcp*log(p0/p(i,j,k)))
          htrld(i,j,k) = 0.e0
          htrlu(i,j,k) = cnv * real(lwhr((j-1)*(ni-1)+i,k-1)/86400.d0)
        end do
        end do
        end do

!$omp end do

        !! copy

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          htrlu(i,j,nk-2) = htrlu(i,j,nk-3)
        end do
        end do

!$omp end do

!$omp end parallel

      end if

      if(mype==root) write(*,*)

! -----

      end subroutine s_rrtmdrv

!-----7--------------------------------------------------------------7--

      end module m_rrtmdrv

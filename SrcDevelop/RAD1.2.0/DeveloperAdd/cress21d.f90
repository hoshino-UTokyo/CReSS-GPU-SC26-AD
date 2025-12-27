!***********************************************************************
      module m_cress21d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/20
!     Modification: 2010/12/17, 2013/02/13, 2013/10/08

!     Author      : Hasegawa Koichi
!     Modification: 2013/12/03, 2014/06/10, 2014/07/01

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the one dimentional variables for mstranx radiation scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commstrn
      use m_comphy
      use m_getiname
      use m_getrname
      use m_vint11

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
      subroutine s_cress21d(fpcphopt,fpradobj,fpcldlim,fpreliqum,       &
     &                      fpreiceum,icnt,i,j,ni,nj,nk,nqw,nqi,        &
     &                      zph,rbr,p,t,qv,qwtr,qice,                   &
     &                      nln,zalt_data,cpcl_data,cgas_data,          &
     &                      zl,pl,tl,pb,tb,cpcl,gdcfrc,cgas,ccfc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpradobj
                       ! Formal parameter of unique index of radobj

      integer, intent(in) :: fpcldlim
                       ! Formal parameter of unique index of cldlim

      integer, intent(in) :: fpreliqum
                       ! Formal parameter of unique index of reliqum

      integer, intent(in) :: fpreiceum
                       ! Formal parameter of unique index of reiceum

      integer, intent(in) :: icnt
                       ! Counter

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

      integer, intent(in) :: nln
                       ! Layers number of external data

      real(kind=r8), intent(in) :: zalt_data(1:nln)
                       ! external data (altitude[km])

      real(kind=r8), intent(in) :: cpcl_data(1:nln,1:kpcl,1:kpclc)
                       ! external data (particulates concentrations)

      real(kind=r8), intent(in) :: cgas_data(1:nln,1:kmol)
                       ! external data (gas concentrations)

! Input and output variables

      real(kind=r8), intent(inout) :: ccfc(1:kcfc)
                       ! CFCs concentrations

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
                       ! 1 : concentrations in standard-state [ppmV]
                       ! 2 : mode radius [cm]
                       ! 3 : undefined

      real(kind=r8), intent(out) :: gdcfrc(1:nk-3)
                       ! Cloud cover late

      real(kind=r8), intent(out) :: cgas(1:nk-3,1:kmol)
                       ! Gas concentrations [ppmV]

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer radobj   ! Option for object substance in radiation
                       ! scheme
      real cldlim      ! Lower limit of mass mixing ratio [kg/kg]
                       ! to consider as cloudy grid
      real reliqum     ! Cloud water effective particle radius [um]
      real reiceum     ! Cloud ice effective particle radius [um]

! Internal private variable

      integer nsnd     ! External data number

      real d1d(1:nk-3) ! Temporary
      real z1d(1:nk-3)
      real dsnd(1:nln)
      real zsnd(1:nln)

      real(kind=r8) cnv
                       ! Factor to convert vapor/water/ice hydrometeor
                       ! mass mixing ratio to volume ratio.

      real(kind=r8) cldliq
                       ! total amount of liquid cloud

      real(kind=r8) cldice
                       ! total amount of ice cloud

      real(kind=r8) fct_qv
      real(kind=r8) fct_qc
      real(kind=r8) fct_qr
      real(kind=r8) fct_qi
      real(kind=r8) fct_qs
      real(kind=r8) fct_qg
      real(kind=r8) fct_other
                       ! Factor

      integer k        ! Array index in z direction

      integer kl       ! Index in layers

      integer ipcl     ! Index for kpcl
      integer ipclc    ! Index for kpclc

      integer imol     ! Index for kmol

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpcphopt,cphopt)
      call getiname(fpradobj,radobj)
      call getrname(fpcldlim,cldlim)
      call getrname(fpreliqum,reliqum)
      call getrname(fpreiceum,reiceum)

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

! Interpolate sublayers to CReSS layers.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=2,nk-2
        zl(nk-k-1)=.0005e0*(zph(i,j,k)+zph(i,j,k+1))
        pl(nk-k-1)=1.e-2*p(i,j,k)
        tl(nk-k-1)=t(i,j,k)
      end do

!$omp end do

!$omp do schedule(runtime) private(k)

      do k=2,nk-1
        pb(nk-k)=.5e-2*(p(i,j,k-1)+p(i,j,k))
        tb(nk-k)=.5e0*(t(i,j,k-1)+t(i,j,k))
      end do

!$omp end do

! -----

! Interpolate external data (gas & ptcl) to CReSS layers only once.

      if(icnt==0) then

!$omp do schedule(runtime) private(k)

        do k=2,nk-2
          z1d(k-1) = .0005e0*(zph(i,j,k)+zph(i,j,k+1))
        end do

!$omp end do

!$omp do schedule(runtime) private(k,imol,nsnd,dsnd,d1d,zsnd)

        do imol=2,kmol

          nsnd=1
          do k=1,nln
            if(zalt_data(k)>30.d0) exit
            zsnd(nsnd) = real(zalt_data(k))
            dsnd(nsnd) = real(cgas_data(k,imol))
            nsnd=nsnd+1
          end do

          call vint11(nk-3,z1d,d1d,nsnd,zsnd,dsnd)

          do k=1,nk-3
            cgas(k,imol) = dble(d1d(nk-2-k)) * fct_other
          end do

        end do

!$omp end do

        do ipclc=1,2

!$omp do schedule(runtime) private(k,ipcl,nsnd,dsnd,d1d,zsnd)

          do ipcl=3,kpcl

            nsnd=1
            do k=1,nln
              if(zalt_data(k)>30.d0) exit
              zsnd(nsnd) = real(zalt_data(k))
              dsnd(nsnd) = real(cpcl_data(k,ipcl,ipclc))
              nsnd=nsnd+1
            end do

            call vint11(nk-3,z1d,d1d,nsnd,zsnd,dsnd)

            do k=1,nk-3
              cpcl(k,ipcl,ipclc) = dble(d1d(nk-2-k)) * fct_other
            end do

          end do

!$omp end do

        end do

!$omp do schedule(runtime) private(imol)

        do imol=1,kcfc
          ccfc(imol) = ccfc(imol) * fct_other
        end do

!$omp end do

      end if

! -----

! Convert mass mixing ratio [kg/kg] of water vapor
! to volume mixing ratio [ppmV].

      !! @ dtrn32 comments -->
      !! cgas  R(nln,ngas)   gas concentration in the sublayer (ppmv)
      !!  1: h2o, 2: co2, 3: o3, 4: n2o, 5: co, 6: ch4, 7: o2

      cnv = dble(rv/rd) * 1.d6

!$omp do schedule(runtime) private(k)

      do k=2,nk-2
        cgas(nk-k-1,1) = dble(qv(i,j,k)) * cnv * fct_qv
      end do

!$omp end do

! -----

! Convert mass mixing ratio [kg/kg] of cloud
! to volume mixing ratio [ppmV].

      !! @ dtrn32 comments -->
      !!  cpcl R(nln,npcl,3) parameter packet for particulates
      !!  in sublayers. (l,i,j)
      !!   l: layer number (from top to bottom)
      !!   1: water cloud               2: ice cloud
      !!   3: dust-like                 4: soot
      !!   5: volcanic-ash              6: h2so4
      !!   7: rural                     8: sea salt
      !!   9: urban                    10: tropo.
      !!  11: yellow dust
      !! j=1: concentration in ppmv
      !!       = aerosol volume(cm3/cm3)*t*p0/t0/p*1.0e6
      !!       with t0=273, p0=1atm.
      !! j=2: mode radius(cm)
      !! j=3: undecided for future use.

      cnv = dble((p20/t0)/(rhow*rd)) * 1.d6

      !! cloud water/ice mixing ratio, cloud fraction

      if(abs(cphopt)==1) then

        !!! warm rain

!$omp do schedule(runtime) private(k,cldliq,cldice)

        do k=2,nk-2

          cldliq = dble(qwtr(i,j,k,1)) * fct_qc                         &
     &           + dble(qwtr(i,j,k,2)) * fct_qr
          cldice = 0.d0

          cpcl(nk-k-1,1,1) = cldliq * cnv
          cpcl(nk-k-1,2,1) = 0.d0
          cpcl(nk-k-1,1,2) = reliqum * 1.d-4   !! [cm] => [um]
          cpcl(nk-k-1,2,2) = 0.d0
          cpcl(nk-k-1,1,3) = 0.d0
          cpcl(nk-k-1,2,3) = 0.d0

          if ( (cldliq+cldice)>=dble(cldlim) ) then
            gdcfrc(nk-k-1) = 1.e0
          else
            gdcfrc(nk-k-1) = 0.e0
          end if

        end do

!$omp end do

      else if(abs(cphopt)>=2) then

        !!! cold rain

!$omp do schedule(runtime) private(k,cldliq,cldice)

        do k=2,nk-2

          cldliq = dble(qwtr(i,j,k,1)) * fct_qc                         &
     &           + dble(qwtr(i,j,k,2)) * fct_qr
          cldice = dble(qice(i,j,k,1)) * fct_qi                         &
     &           + dble(qice(i,j,k,2)) * fct_qs                         &
     &           + dble(qice(i,j,k,3)) * fct_qg

          ! ** debug for mstrnX-CAOS -->

          !! low cloud
!!!       if(zl(nk-k-1)>0.0 .and. zl(nk-k-1)<=2.0) then
!!!!!       cldliq = 6.7919d-03 * fct_qc / cnv   !! ppmV
!!!         cldliq = 5.2556d-06 * fct_qc         !! mix. of qc [kg/kg]
!!!         cldice = 0.d0
!!!       else
!!!         cldliq = 0.d0
!!!         cldice = 0.d0
!!!       end if
   
!!!       !! high cloud
!!!       if(zl(nk-k-1)>9.0 .and. zl(nk-k-1)<=10.0) then
!!!         cldliq = 0.d0
!!!!!       cldice = 1.3979d-02 * fct_qi / cnv   !! ppmV
!!!         cldice = 1.0817d-05 * fct_qi         !! mix. of qi [kg/kg]
!!!       else
!!!         cldliq = 0.d0
!!!         cldice = 0.d0
!!!       end if

          ! <-- debug for mstrnX-CAOS **

          cpcl(nk-k-1,1,1) = cldliq * cnv
          cpcl(nk-k-1,2,1) = cldice * cnv
          cpcl(nk-k-1,1,2) = reliqum * 1.d-4    !! [um] => [cm]
          cpcl(nk-k-1,2,2) = reiceum * 1.d-4
          cpcl(nk-k-1,1,3) = 0.d0
          cpcl(nk-k-1,2,3) = 0.d0

          if ( (cldliq+cldice)>=dble(cldlim) ) then
            gdcfrc(nk-k-1) = 1.e0
          else
            gdcfrc(nk-k-1) = 0.e0
          end if

        end do

!$omp end do

      end if

! -----

!$omp end parallel

!!! -----

      end subroutine s_cress21d

!-----7--------------------------------------------------------------7--

      end module m_cress21d

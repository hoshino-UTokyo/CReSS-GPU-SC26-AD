!***********************************************************************
      module m_comrrtm
!***********************************************************************

!     Author      : Oda Naotaka
!     Date        : 2013/11/12
!     Modification:

!     Author      : Hasegawa Koichi
!     Modification: 2013/11/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the parameter and the array and variable for RRTM
!     radiation scheme.

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

      real(kind=r8), allocatable, save :: play(:,:)
                      ! Layer pressures [hPa] (ncol,nlay)

      real(kind=r8), allocatable, save :: plev(:,:)
                      ! Interface pressures [hPa] (ncol,nlay+1)

      real(kind=r8), allocatable, save :: tlay(:,:)
                      ! Layer temperatures [K] (ncol,nlay)

      real(kind=r8), allocatable, save :: tlev(:,:)
                      ! Interface temperatures [K] (ncol,nlay+1)

      real(kind=r8), allocatable, save :: tsfc(:)
                      ! Surface temperature [K] (ncol)

      real(kind=r8), allocatable, save :: h2ovmr(:,:)
                      ! H2O volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: o3vmr(:,:)
                      ! O3 volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: co2vmr(:,:)
                      ! CO2 volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: ch4vmr(:,:)
                      ! Methane volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: n2ovmr(:,:)
                      ! Nitrous oxide volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: o2vmr(:,:)
                      ! Oxigen volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: cfc11vmr(:,:)
                      ! CFC11 volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: cfc12vmr(:,:)
                      ! CFC12 volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: cfc22vmr(:,:)
                      ! CFC22 volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: ccl4vmr(:,:)
                      ! CCL4 volume mixing ratio (ncol,nlay)

      real(kind=r8), allocatable, save :: emis(:,:)
                      ! Surface emissivity (ncol,nbndlw)

      real(kind=r8), allocatable, save :: cldfr(:,:)
                      ! Cloud fraction (ncol,nlay)

      real(kind=r8), allocatable, save :: taucld(:,:,:)
                      ! In-cloud optical depth (nbndlw,ncol,nlay)

      real(kind=r8), allocatable, save :: cicewp(:,:)
                      ! Cloud ice water path (g/m2)
                      ! (ncol,nlay)

      real(kind=r8), allocatable, save :: cliqwp(:,:)
                      ! Cloud liquid water path (g/m2)
                      ! (ncol,nlay)

      real(kind=r8), allocatable, save :: reice(:,:)
                      ! Cloud ice particle effective size (microns)
                      ! (ncol,nlay)

      real(kind=r8), allocatable, save :: reliq(:,:)
                      ! Cloud water drop effective radius (microns)
                      ! (ncol,nlay)

      real(kind=r8), allocatable, save :: tauaer(:,:,:)
                      ! Aerosol optical depth
                      ! at mid-point of LW spectral bands
                      ! (ncol,nlay,nbndlw)

!ORG  real(kind=r8), allocatable, save :: uflx(:,:)
!ORG                  ! Total sky longwave upward flux [W/m2]
!ORG                  ! (ncol,nlay+1)

!ORG  real(kind=r8), allocatable, save :: dflx(:,:)
!ORG                  ! Total sky longwave downward [W/m2]
!ORG                  ! (ncol,nlay+1)

      real(kind=r8), allocatable, save :: lwhr(:,:)
                      ! Total sky longwave radiative heating rate [K/d]
                      ! (ncol,nlay)

!ORG  real(kind=r8), allocatable, save :: uflxc(:,:)
!ORG                  ! Clear sky longwave upward flux [W/m2]
!ORG                  ! (ncol,nlay+1)

!ORG  real(kind=r8), allocatable, save :: dflxc(:,:)
!ORG                  ! Clear sky longwave downward flux [W/m2]
!ORG                  ! (ncol,nlay+1)

!ORG  real(kind=r8), allocatable, save :: hrc(:,:)
!ORG                  ! Clear sky longwave radiative heating rate [K/d]
!ORG                  ! (ncol,nlay)

      real(kind=r8), allocatable, save :: duflx_dt(:,:) ! ncol,nlay
                      ! Change in upward longwave flux
                      ! with respect to surface temperatuere [w/m2/K]
                      ! (ncol,nlay)

      real(kind=r8), allocatable, save :: duflxc_dt(:,:)
                      ! Change in clear sky upward longwave flux
                      ! with respect to surface temperature [W/m2/K]
                      ! (ncol,nlay)

! short wave

      real(kind=r8), allocatable, save :: asdir(:)
                      ! UV/vis surface albedo direct rad (ncol)

      real(kind=r8), allocatable, save :: asdif(:)
                      ! UV/vis surface albedo: diffuse rad (ncol)

      real(kind=r8), allocatable, save :: aldir(:)
                      ! Near-IR surface albedo direct rad (ncol)

      real(kind=r8), allocatable, save :: aldif(:)
                      ! Near-IR surface albedo: diffuse rad (ncol)

      real(kind=r8), allocatable, save :: coszen(:)
                      ! Cosine of solar zenith angle (ncol)

      real(kind=r8), allocatable, save :: taucldsw(:,:,:)
                      ! In-cloud optical depth (nbndsw,ncol,nlay)

      real(kind=r8), allocatable, save :: ssacld(:,:,:)
                      ! In-cloud single scattering albedo
                      ! (nbndsw,ncol,nlay)

      real(kind=r8), allocatable, save :: asmcld(:,:,:)
                      ! In-cloud asymmetry parameter
                      ! (nbndsw,ncol,nlay)

      real(kind=r8), allocatable, save :: fsfcld(:,:,:)
                      ! In-cloud forward scattering fraction
                      ! (nbndsw,ncol,nlay)

      real(kind=r8), allocatable, save :: tauaersw(:,:,:)
                      ! Aerosol optical depth [iaer=10 only]
                      ! (ncol,nlay,nbndsw)

      real(kind=r8), allocatable, save :: ssaaer(:,:,:)
                      ! Aerosol single scattering albedo [iaer=10 only]
                      ! (ncol,nlay,nbndsw)

      real(kind=r8), allocatable, save :: asmaer(:,:,:)
                      ! Aerosol asymmetry parameter [iaer=10 only]
                      ! (ncol,nlay,nbndsw)

      real(kind=r8), allocatable, save :: ecaer(:,:,:)
                      ! Aerosol optical depth at 0.55 micron
                      ! [iaer=6 only] (ncol,nlay,naerec)

!ORG  real(kind=r8), allocatable, save :: swuflx(:,:)
!ORG                  ! Total sky shortwave upward flux (W/m2)
!ORG                  ! (ncol,nlay+1)

!ORG  real(kind=r8), allocatable, save :: swdflx(:,:)
!ORG                  ! Total sky shortwave downward flux (W/m2)
!ORG                  ! (ncol,nlay+1)

      real(kind=r8), allocatable, save :: swhr(:,:)
                      ! Total sky shortwave radiative heating rate (K/d)
                      ! (ncol,nlay)

!ORG  real(kind=r8), allocatable, save :: swuflxc(:,:)
!ORG                  ! Clear sky shortwave upward flux (W/m2)
!ORG                  ! (ncol,nlay+1)

!ORG  real(kind=r8), allocatable, save :: swdflxc(:,:)
!ORG                  ! Clear sky shortwave downward flux (W/m2)
!ORG                  ! (ncol,nlay+1)

!ORG  real(kind=r8), allocatable, save :: swhrc(:,:)
!ORG                  ! Clear sky shortwave radiative heating rate (K/d)
!ORG                  ! (ncol,nlay)

! Temporary

      real(kind=r8), allocatable, save :: tmp1(:,:)   ! (ncol,nlay+1)
      real(kind=r8), allocatable, save :: tmp2(:,:)
      real(kind=r8), allocatable, save :: tmp3(:,:)
      real(kind=r8), allocatable, save :: tmp4(:,:)
      real(kind=r8), allocatable, save :: tmp5(:,:)   ! (ncol,nlay)

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

      end module m_comrrtm

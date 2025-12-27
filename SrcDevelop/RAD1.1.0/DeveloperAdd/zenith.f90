!***********************************************************************
      module m_zenith
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/20
!     Modification: 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculte the zenith angle.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comdays
      use m_commath

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: zenith, s_zenith

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface zenith

        module procedure s_zenith

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic sin
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_zenith(cdate,ni,nj,lat,lon,coseta)
!***********************************************************************

! Input variables

      character(len=12), intent(in) :: cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: lat(0:ni+1,0:nj+1)
                       ! Latitude

      real, intent(in) :: lon(0:ni+1,0:nj+1)
                       ! Longitude

! Output variable

      real, intent(out) :: coseta(0:ni+1,0:nj+1)
                       ! cos(Zenith angle)

! Internal shared variables

      integer cyr      ! Year of current forecast date
      integer cmo      ! Month of current forecast date
      integer cdy      ! Day of current forecast date
      integer chr      ! Hour of current forecast date
      integer cmn      ! Minite of current forecast date

      real jday        ! Number of elapse of days from start of year

      real eqt         ! Equation of local time

      real phs         ! Solar angle

      real rchr        ! real(chr)
      real rcmn        ! real(cmn) / 60.0

      real sinphs      ! sin(phs)
      real cosphs      ! cos(phs)

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real tlc         ! Local time

!-----7--------------------------------------------------------------7--

! Read out the integer variables from the input current forecast date
! with Gregorian calendar, yyyymmddhhmm.

      read(cdate(1:12),'(i4.4,4i2.2)') cyr,cmo,cdy,chr,cmn

! -----

! Calculate the solar angle.

      if(mod(cyr,400).eq.0                                              &
     &  .or.(mod(cyr,4).eq.0.and.mod(cyr,100).ne.0)) then

        jday=2.e0*cc*i366*real(elaitc(cmo-1)+cdy-1)

      else

        jday=2.e0*cc*i365*real(ela(cmo-1)+cdy-1)

      end if

      eqt=.000075e0+.001868e0*cos(jday)-.032077e0*sin(jday)             &
     &  -.014615e0*cos(2.e0*jday)-.040849e0*sin(2.e0*jday)

      phs=.006918e0-.399912e0*cos(jday)+.070257e0*sin(jday)             &
     &  -.006758e0*cos(2.e0*jday)+.000907e0*sin(2.e0*jday)              &
     &  -.002697e0*cos(3.e0*jday)+.001480e0*sin(3.e0*jday)

! -----

! Set the common used variables.

      rchr=real(chr)
      rcmn=oned60*real(cmn)

      sinphs=sin(phs)
      cosphs=cos(phs)

! -----

! Calculte the zenith angle.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,tlc)

      do j=1,nj-1
      do i=1,ni-1

        tlc=rchr+rcmn+oned15*lon(i,j)

        coseta(i,j)=sinphs*sin(lat(i,j)*d2r)                            &
     &    +cosphs*cos(lat(i,j)*d2r)*cos(eqt+15.e0*(tlc-12.e0)*d2r)

      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_zenith

!-----7--------------------------------------------------------------7--

      end module m_zenith

!     path:      $Source:
!     /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 22:22:22 $

!  --------------------------------------------------------------------------
! |
! |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc.
! (AER).  |
! |  This software may be used, copied, or redistributed as long as it
! is    |
! |  not sold and this copyright notice is reproduced on each copy made.
! |
! |  This model is provided as is without any express or implied
! warranties. |
! |                       (http://www.rtweb.aer.com/)
! |
! |
! |
!  --------------------------------------------------------------------------

! **************************************************************************
!      subroutine sw_kgbnn
! **************************************************************************
!  RRTM Shortwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original by J.Delamere, Atmospheric & Environmental Research.
!  Reformatted for F90: JJMorcrette, ECMWF
!  Further F90 and GCM revisions:  MJIacono, AER, July 2002
!
!  This file contains 14 subroutines that include the 
!  absorption coefficients and other data for each of the 14 shortwave
!  spectral bands used in RRTM_SW.  Here, the data are defined for 16
!  g-points, or sub-intervals, per band.  These data are combined and
!  weighted using a mapping procedure in routine RRTMG_SW_INIT to reduce
!  the total number of g-points from 224 to 112 for use in the GCM.
! **************************************************************************
!!prev      subroutine sw_kgb20
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for  sfluxrefo(:), absch4o(:), rayl
!!
      subroutine sw_kgb200
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg20, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
!!prev                            absch4o, rayl
      use rrsw_kg20, only : sfluxrefo, &
                            absch4o, rayl

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:) = (/ &
        & 9.34081_rb , 8.93720_rb    , 8.19346_rb    , 7.39196_rb    , &
        & 6.12127_rb , 5.23956_rb    , 4.24941_rb    , 3.20013_rb    , &
        & 2.16047_rb , 0.234509_rb   , 0.194593_rb   , 0.151512_rb   , &
        & 0.110315_rb, 7.09959e-02_rb, 2.70573e-02_rb, 3.36042e-03_rb /)
  
      absch4o(:) = (/ &  
        & 1.01381e-03_rb,6.33692e-03_rb,1.94185e-02_rb,4.83210e-02_rb, &
        & 2.36574e-03_rb,6.61973e-04_rb,5.64552e-04_rb,2.83183e-04_rb, &
        & 7.43623e-05_rb,8.90159e-07_rb,6.98728e-07_rb,6.51832e-08_rb, &
        & 2.96619e-08_rb,         0._rb,         0._rb,         0._rb /)

! Rayleigh extinction coefficient at v = 5670 cm-1.
      rayl = 4.12e-09_rb

      end subroutine sw_kgb200

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
!!prev      subroutine sw_kgb23
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for  sfluxrefo(:), raylo(:)
!!
      subroutine sw_kgb230
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg23, only : kao, selfrefo, forrefo, sfluxrefo, raylo
      use rrsw_kg23, only : sfluxrefo, raylo

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:) = (/ &
        & 53.2101_rb , 51.4143_rb, 49.3348_rb, 45.4612_rb    , &
        & 40.8294_rb , 35.1801_rb, 28.6947_rb, 21.5751_rb    , &
        & 14.6388_rb , 1.59111_rb, 1.31860_rb, 1.04018_rb    , &
        & 0.762140_rb,0.484214_rb,0.182275_rb, 2.54948e-02_rb /)

! Rayleigh extinction coefficient at all v 
      raylo(:) = (/ &
        & 5.94837e-08_rb,5.70593e-08_rb,6.27845e-08_rb,5.56602e-08_rb, &
        & 5.25571e-08_rb,4.73388e-08_rb,4.17466e-08_rb,3.98097e-08_rb, &
        & 4.00786e-08_rb,3.67478e-08_rb,3.45186e-08_rb,3.46156e-08_rb, &
        & 3.32155e-08_rb,3.23642e-08_rb,2.72590e-08_rb,2.96813e-08_rb /)

      end subroutine sw_kgb230

!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 22:22:22 $

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
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
!!prev      subroutine sw_kgb16
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for sfluxrefo(:), rayl
!!
      subroutine sw_kgb160
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg16, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg16, only : sfluxrefo, rayl

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:) = (/ &
        &  1.92269_rb    , 1.72844_rb    , 1.64326_rb    , 1.58451_rb     &
        &, 1.44031_rb    , 1.25108_rb    , 1.02724_rb    , 0.776759_rb    &
        &, 0.534444_rb   , 5.87755e-02_rb, 4.86706e-02_rb, 3.87989e-02_rb &
        &, 2.84532e-02_rb, 1.82431e-02_rb, 6.92320e-03_rb, 9.70770e-04_rb /)

! Rayleigh extinction coefficient at v = 2925 cm-1.
      rayl = 2.91e-10_rb

      end subroutine sw_kgb160

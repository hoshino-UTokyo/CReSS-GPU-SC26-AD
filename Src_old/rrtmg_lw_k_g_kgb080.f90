!     path:      $Source:
!     /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 21:04:30 $
!
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
!      subroutine lw_kgbnn
! **************************************************************************
!  RRTM Longwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original version:   E. J. Mlawer, et al.
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90 formatting:  Michael J. Iacono; June 2006
!
!  This file contains 16 subroutines that include the 
!  absorption coefficients and other data for each of the 16 longwave
!  spectral bands used in RRTM.  Here, the data are defined for 16
!  g-points, or sub-intervals, per band.  These data are combined and
!  weighted using a mapping procedure in routine RRTMG_LW_INIT to reduce
!  the total number of g-points from 256 to 140 for use in the GCM.
! **************************************************************************
!!prev      subroutine lw_kgb08
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:,:), fracrefbo(:,:), cfc12o(:), cfc22adjo(:)
!!
      subroutine lw_kgb080
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg08, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
!!prev                            kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
!!prev                            cfc12o, cfc22adjo
      use rrlw_kg08, only : fracrefao, fracrefbo, &
                            cfc12o, cfc22adjo

      implicit none
      save

! Planck fraction mapping level : P=473.4280 mb, T = 259.83 K
      fracrefao(:) = (/ &
        1.6004e-01_rb,1.5437e-01_rb,1.4502e-01_rb,1.3084e-01_rb,1.1523e-01_rb,9.7743e-02_rb, &
        8.0376e-02_rb,6.0261e-02_rb,4.1111e-02_rb,4.4772e-03_rb,3.6511e-03_rb,2.9154e-03_rb, &
        2.1184e-03_rb,1.3048e-03_rb,4.6637e-04_rb,6.5624e-05_rb/)

! Planck fraction mapping level : P=95.5835 mb, T= 215.7 K
      fracrefbo(:) = (/ &
        1.4987e-01_rb,1.4665e-01_rb,1.4154e-01_rb,1.3200e-01_rb,1.1902e-01_rb,1.0352e-01_rb, &
        8.4939e-02_rb,6.4105e-02_rb,4.3190e-02_rb,4.5129e-03_rb,3.7656e-03_rb,2.8733e-03_rb, &
        2.0947e-03_rb,1.3201e-03_rb,5.1832e-04_rb,7.7473e-05_rb/)

! Minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

      cfc12o(:) = (/ &
        85.4027_rb, 89.4696_rb, 74.0959_rb, 67.7480_rb, &
        61.2444_rb, 59.9073_rb, 60.8296_rb, 63.0998_rb, &
        59.6110_rb, 64.0735_rb, 57.2622_rb, 58.9721_rb, &
        43.5505_rb, 26.1192_rb, 32.7023_rb, 32.8667_rb/)
! Original CFC22 is multiplied by 1.485 to account for the 780-850 cm-1 
! and 1290-1335 cm-1 bands.
      cfc22adjo(:) = (/ &
        135.335_rb, 89.6642_rb, 76.2375_rb, 65.9748_rb, &
        63.1164_rb, 60.2935_rb, 64.0299_rb, 75.4264_rb, &
        51.3018_rb, 7.07911_rb, 5.86928_rb, 0.398693_rb, &
        2.82885_rb, 9.12751_rb, 6.28271_rb, 0._rb/)

      end subroutine lw_kgb080

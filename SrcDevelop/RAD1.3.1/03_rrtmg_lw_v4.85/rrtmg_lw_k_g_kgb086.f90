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
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb086
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg08, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
!!prev                            kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
!!prev                            cfc12o, cfc22adjo

      use rrlw_kg08, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &4.8166e-07_rb,3.7500e-07_rb,4.8978e-07_rb,5.9624e-07_rb,6.3742e-07_rb,7.5551e-07_rb, &
     &7.7706e-07_rb,6.8681e-07_rb,7.5212e-07_rb,8.0956e-07_rb,7.8117e-07_rb,7.4835e-07_rb, &
     &9.4118e-07_rb,1.2585e-06_rb,1.4976e-06_rb,1.4976e-06_rb/)
      forrefo(2,:) = (/ &
     &3.1320e-07_rb,4.0764e-07_rb,4.7468e-07_rb,5.9976e-07_rb,7.3324e-07_rb,8.1488e-07_rb, &
     &7.6442e-07_rb,8.2007e-07_rb,7.7721e-07_rb,7.6377e-07_rb,8.0327e-07_rb,7.1881e-07_rb, &
     &8.2148e-07_rb,1.0203e-06_rb,1.5033e-06_rb,1.5032e-06_rb/)
      forrefo(3,:) = (/ &
     &4.1831e-07_rb,5.5043e-07_rb,5.7783e-07_rb,6.1294e-07_rb,6.3396e-07_rb,6.2292e-07_rb, &
     &6.1719e-07_rb,6.4183e-07_rb,7.6180e-07_rb,9.5477e-07_rb,9.5901e-07_rb,1.0207e-06_rb, &
     &1.0387e-06_rb,1.1305e-06_rb,1.3602e-06_rb,1.5063e-06_rb/)
      forrefo(4,:) = (/ &
     &8.5878e-07_rb,6.0921e-07_rb,5.5773e-07_rb,5.3374e-07_rb,5.0495e-07_rb,4.9844e-07_rb, &
     &5.1536e-07_rb,5.2908e-07_rb,4.7977e-07_rb,5.3177e-07_rb,4.9266e-07_rb,4.5403e-07_rb, &
     &3.9695e-07_rb,3.4792e-07_rb,3.4912e-07_rb,3.4102e-07_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 3.16029e-02_rb, 2.74633e-02_rb, 2.38660e-02_rb, 2.07399e-02_rb, 1.80232e-02_rb, &
     & 1.56624e-02_rb, 1.36108e-02_rb, 1.18280e-02_rb, 1.02787e-02_rb, 8.93231e-03_rb/)
      selfrefo(:, 2) = (/ &
     & 3.10422e-02_rb, 2.71312e-02_rb, 2.37130e-02_rb, 2.07254e-02_rb, 1.81142e-02_rb, &
     & 1.58320e-02_rb, 1.38374e-02_rb, 1.20940e-02_rb, 1.05703e-02_rb, 9.23854e-03_rb/)
      selfrefo(:, 3) = (/ &
     & 3.08657e-02_rb, 2.69431e-02_rb, 2.35190e-02_rb, 2.05301e-02_rb, 1.79210e-02_rb, &
     & 1.56435e-02_rb, 1.36554e-02_rb, 1.19200e-02_rb, 1.04051e-02_rb, 9.08279e-03_rb/)
      selfrefo(:, 4) = (/ &
     & 3.02668e-02_rb, 2.64686e-02_rb, 2.31470e-02_rb, 2.02422e-02_rb, 1.77020e-02_rb, &
     & 1.54806e-02_rb, 1.35379e-02_rb, 1.18390e-02_rb, 1.03533e-02_rb, 9.05406e-03_rb/)
      selfrefo(:, 5) = (/ &
     & 2.98317e-02_rb, 2.61491e-02_rb, 2.29210e-02_rb, 2.00914e-02_rb, 1.76112e-02_rb, &
     & 1.54371e-02_rb, 1.35314e-02_rb, 1.18610e-02_rb, 1.03968e-02_rb, 9.11332e-03_rb/)
      selfrefo(:, 6) = (/ &
     & 2.95545e-02_rb, 2.59083e-02_rb, 2.27120e-02_rb, 1.99100e-02_rb, 1.74537e-02_rb, &
     & 1.53004e-02_rb, 1.34128e-02_rb, 1.17580e-02_rb, 1.03074e-02_rb, 9.03576e-03_rb/)
      selfrefo(:, 7) = (/ &
     & 2.97352e-02_rb, 2.60320e-02_rb, 2.27900e-02_rb, 1.99517e-02_rb, 1.74670e-02_rb, &
     & 1.52916e-02_rb, 1.33872e-02_rb, 1.17200e-02_rb, 1.02604e-02_rb, 8.98258e-03_rb/)
      selfrefo(:, 8) = (/ &
     & 2.96543e-02_rb, 2.59760e-02_rb, 2.27540e-02_rb, 1.99316e-02_rb, 1.74593e-02_rb, &
     & 1.52937e-02_rb, 1.33967e-02_rb, 1.17350e-02_rb, 1.02794e-02_rb, 9.00437e-03_rb/)
      selfrefo(:, 9) = (/ &
     & 2.97998e-02_rb, 2.60786e-02_rb, 2.28220e-02_rb, 1.99721e-02_rb, 1.74781e-02_rb, &
     & 1.52955e-02_rb, 1.33855e-02_rb, 1.17140e-02_rb, 1.02512e-02_rb, 8.97110e-03_rb/)
      selfrefo(:,10) = (/ &
     & 2.98826e-02_rb, 2.61096e-02_rb, 2.28130e-02_rb, 1.99326e-02_rb, 1.74159e-02_rb, &
     & 1.52170e-02_rb, 1.32957e-02_rb, 1.16170e-02_rb, 1.01502e-02_rb, 8.86867e-03_rb/)
      selfrefo(:,11) = (/ &
     & 2.94710e-02_rb, 2.58147e-02_rb, 2.26120e-02_rb, 1.98066e-02_rb, 1.73493e-02_rb, &
     & 1.51969e-02_rb, 1.33115e-02_rb, 1.16600e-02_rb, 1.02134e-02_rb, 8.94628e-03_rb/)
      selfrefo(:,12) = (/ &
     & 2.96297e-02_rb, 2.59544e-02_rb, 2.27350e-02_rb, 1.99149e-02_rb, 1.74446e-02_rb, &
     & 1.52808e-02_rb, 1.33853e-02_rb, 1.17250e-02_rb, 1.02706e-02_rb, 8.99663e-03_rb/)
      selfrefo(:,13) = (/ &
     & 2.96272e-02_rb, 2.59013e-02_rb, 2.26440e-02_rb, 1.97963e-02_rb, 1.73067e-02_rb, &
     & 1.51302e-02_rb, 1.32275e-02_rb, 1.15640e-02_rb, 1.01097e-02_rb, 8.83833e-03_rb/)
      selfrefo(:,14) = (/ &
     & 2.89906e-02_rb, 2.53971e-02_rb, 2.22490e-02_rb, 1.94911e-02_rb, 1.70751e-02_rb, &
     & 1.49585e-02_rb, 1.31044e-02_rb, 1.14800e-02_rb, 1.00570e-02_rb, 8.81038e-03_rb/)
      selfrefo(:,15) = (/ &
     & 2.80884e-02_rb, 2.46987e-02_rb, 2.17180e-02_rb, 1.90970e-02_rb, 1.67924e-02_rb, &
     & 1.47659e-02_rb, 1.29839e-02_rb, 1.14170e-02_rb, 1.00392e-02_rb, 8.82765e-03_rb/)
      selfrefo(:,16) = (/ &
     & 2.80884e-02_rb, 2.46987e-02_rb, 2.17180e-02_rb, 1.90970e-02_rb, 1.67924e-02_rb, &
     & 1.47659e-02_rb, 1.29839e-02_rb, 1.14170e-02_rb, 1.00392e-02_rb, 8.82765e-03_rb/)

      end subroutine lw_kgb086

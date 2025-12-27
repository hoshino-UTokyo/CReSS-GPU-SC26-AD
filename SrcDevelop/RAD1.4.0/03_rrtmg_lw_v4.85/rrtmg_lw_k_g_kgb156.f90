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
!!prev      subroutine lw_kgb15
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb156
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg15, only : fracrefao, kao, kao_mn2, selfrefo, forrefo
      use rrlw_kg15, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.1755e-06_rb,6.5398e-07_rb,4.3915e-07_rb,3.0753e-07_rb,1.9677e-07_rb,1.4362e-07_rb, &
     &9.4598e-08_rb,1.1848e-07_rb,1.4280e-07_rb,1.5821e-07_rb,1.5816e-07_rb,1.5769e-07_rb, &
     &1.5844e-07_rb,1.6016e-07_rb,1.6232e-07_rb,1.6320e-07_rb/)
      forrefo(2,:) = (/ &
     &1.0703e-06_rb,6.2783e-07_rb,4.7122e-07_rb,2.6300e-07_rb,1.8538e-07_rb,1.5076e-07_rb, &
     &1.9474e-07_rb,2.9543e-07_rb,2.0093e-07_rb,1.5819e-07_rb,1.5826e-07_rb,1.5737e-07_rb, &
     &1.5751e-07_rb,1.5910e-07_rb,1.6181e-07_rb,1.6320e-07_rb/)
      forrefo(3,:) = (/ &
     &1.0470e-06_rb,5.8184e-07_rb,4.8218e-07_rb,2.7771e-07_rb,1.9036e-07_rb,1.5737e-07_rb, &
     &1.8633e-07_rb,2.5754e-07_rb,4.0647e-07_rb,1.5839e-07_rb,1.5914e-07_rb,1.5788e-07_rb, &
     &1.5731e-07_rb,1.5836e-07_rb,1.6103e-07_rb,1.6320e-07_rb/)
      forrefo(4,:) = (/ &
     &1.3891e-06_rb,5.4901e-07_rb,2.8850e-07_rb,1.9176e-07_rb,1.4549e-07_rb,1.3603e-07_rb, &
     &1.7472e-07_rb,2.9796e-07_rb,3.2452e-07_rb,2.5231e-07_rb,2.8195e-07_rb,1.5527e-07_rb, &
     &1.5507e-07_rb,1.5442e-07_rb,1.5275e-07_rb,1.6057e-07_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 1.73980e-03_rb, 1.41928e-03_rb, 1.15780e-03_rb, 9.44496e-04_rb, 7.70490e-04_rb, &
     & 6.28541e-04_rb, 5.12744e-04_rb, 4.18280e-04_rb, 3.41219e-04_rb, 2.78356e-04_rb/)
      selfrefo(:, 2) = (/ &
     & 1.84082e-03_rb, 1.50228e-03_rb, 1.22600e-03_rb, 1.00053e-03_rb, 8.16525e-04_rb, &
     & 6.66359e-04_rb, 5.43811e-04_rb, 4.43800e-04_rb, 3.62182e-04_rb, 2.95574e-04_rb/)
      selfrefo(:, 3) = (/ &
     & 1.92957e-03_rb, 1.57727e-03_rb, 1.28930e-03_rb, 1.05390e-03_rb, 8.61484e-04_rb, &
     & 7.04197e-04_rb, 5.75627e-04_rb, 4.70530e-04_rb, 3.84622e-04_rb, 3.14399e-04_rb/)
      selfrefo(:, 4) = (/ &
     & 2.12958e-03_rb, 1.73572e-03_rb, 1.41470e-03_rb, 1.15305e-03_rb, 9.39798e-04_rb, &
     & 7.65984e-04_rb, 6.24317e-04_rb, 5.08850e-04_rb, 4.14739e-04_rb, 3.38034e-04_rb/)
      selfrefo(:, 5) = (/ &
     & 2.30636e-03_rb, 1.88401e-03_rb, 1.53900e-03_rb, 1.25717e-03_rb, 1.02695e-03_rb, &
     & 8.38891e-04_rb, 6.85270e-04_rb, 5.59780e-04_rb, 4.57270e-04_rb, 3.73533e-04_rb/)
      selfrefo(:, 6) = (/ &
     & 2.47824e-03_rb, 2.03278e-03_rb, 1.66740e-03_rb, 1.36769e-03_rb, 1.12185e-03_rb, &
     & 9.20206e-04_rb, 7.54803e-04_rb, 6.19130e-04_rb, 5.07844e-04_rb, 4.16561e-04_rb/)
      selfrefo(:, 7) = (/ &
     & 2.54196e-03_rb, 2.10768e-03_rb, 1.74760e-03_rb, 1.44904e-03_rb, 1.20148e-03_rb, &
     & 9.96215e-04_rb, 8.26019e-04_rb, 6.84900e-04_rb, 5.67890e-04_rb, 4.70870e-04_rb/)
      selfrefo(:, 8) = (/ &
     & 2.52650e-03_rb, 2.11773e-03_rb, 1.77510e-03_rb, 1.48790e-03_rb, 1.24717e-03_rb, &
     & 1.04539e-03_rb, 8.76251e-04_rb, 7.34480e-04_rb, 6.15646e-04_rb, 5.16039e-04_rb/)
      selfrefo(:, 9) = (/ &
     & 2.82351e-03_rb, 2.34652e-03_rb, 1.95010e-03_rb, 1.62065e-03_rb, 1.34686e-03_rb, &
     & 1.11933e-03_rb, 9.30232e-04_rb, 7.73080e-04_rb, 6.42477e-04_rb, 5.33939e-04_rb/)
      selfrefo(:,10) = (/ &
     & 2.98189e-03_rb, 2.46741e-03_rb, 2.04170e-03_rb, 1.68944e-03_rb, 1.39795e-03_rb, &
     & 1.15676e-03_rb, 9.57176e-04_rb, 7.92030e-04_rb, 6.55377e-04_rb, 5.42302e-04_rb/)
      selfrefo(:,11) = (/ &
     & 2.98239e-03_rb, 2.46774e-03_rb, 2.04190e-03_rb, 1.68954e-03_rb, 1.39799e-03_rb, &
     & 1.15675e-03_rb, 9.57137e-04_rb, 7.91970e-04_rb, 6.55305e-04_rb, 5.42224e-04_rb/)
      selfrefo(:,12) = (/ &
     & 2.97833e-03_rb, 2.46461e-03_rb, 2.03950e-03_rb, 1.68772e-03_rb, 1.39661e-03_rb, &
     & 1.15571e-03_rb, 9.56370e-04_rb, 7.91410e-04_rb, 6.54903e-04_rb, 5.41942e-04_rb/)
      selfrefo(:,13) = (/ &
     & 2.97779e-03_rb, 2.46463e-03_rb, 2.03990e-03_rb, 1.68836e-03_rb, 1.39741e-03_rb, &
     & 1.15659e-03_rb, 9.57278e-04_rb, 7.92310e-04_rb, 6.55771e-04_rb, 5.42762e-04_rb/)
      selfrefo(:,14) = (/ &
     & 2.98326e-03_rb, 2.46943e-03_rb, 2.04410e-03_rb, 1.69203e-03_rb, 1.40060e-03_rb, &
     & 1.15936e-03_rb, 9.59673e-04_rb, 7.94380e-04_rb, 6.57557e-04_rb, 5.44301e-04_rb/)
      selfrefo(:,15) = (/ &
     & 2.99407e-03_rb, 2.47825e-03_rb, 2.05130e-03_rb, 1.69790e-03_rb, 1.40539e-03_rb, &
     & 1.16327e-03_rb, 9.62862e-04_rb, 7.96980e-04_rb, 6.59676e-04_rb, 5.46028e-04_rb/)
      selfrefo(:,16) = (/ &
     & 3.00005e-03_rb, 2.48296e-03_rb, 2.05500e-03_rb, 1.70080e-03_rb, 1.40765e-03_rb, &
     & 1.16503e-03_rb, 9.64224e-04_rb, 7.98030e-04_rb, 6.60481e-04_rb, 5.46641e-04_rb/)

      end subroutine lw_kgb156

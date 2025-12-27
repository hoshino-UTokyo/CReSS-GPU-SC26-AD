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
!!prev      subroutine lw_kgb14
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb144
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg14, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg14, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &2.7075e-06_rb,2.2609e-06_rb,1.5633e-06_rb,8.7484e-07_rb,5.5470e-07_rb,4.8456e-07_rb, &
     &4.7463e-07_rb,4.6154e-07_rb,4.4425e-07_rb,4.2960e-07_rb,4.2626e-07_rb,4.1715e-07_rb, &
     &4.2607e-07_rb,3.6616e-07_rb,2.6366e-07_rb,2.6029e-07_rb/)
      forrefo(2,:) = (/ &
     &2.6759e-06_rb,2.2237e-06_rb,1.4466e-06_rb,9.3032e-07_rb,6.4927e-07_rb,5.4809e-07_rb, &
     &4.9504e-07_rb,4.6305e-07_rb,4.4873e-07_rb,4.2146e-07_rb,4.2176e-07_rb,4.2812e-07_rb, &
     &4.0529e-07_rb,4.0969e-07_rb,2.9442e-07_rb,2.6821e-07_rb/)
      forrefo(3,:) = (/ &
     &2.6608e-06_rb,2.1140e-06_rb,1.4838e-06_rb,9.2083e-07_rb,6.3350e-07_rb,5.7195e-07_rb, &
     &6.2253e-07_rb,5.1783e-07_rb,4.4749e-07_rb,4.3261e-07_rb,4.2553e-07_rb,4.2175e-07_rb, &
     &4.1085e-07_rb,4.0358e-07_rb,3.5340e-07_rb,2.7191e-07_rb/)
      forrefo(4,:) = (/ &
     &2.6412e-06_rb,1.9814e-06_rb,1.2672e-06_rb,8.1129e-07_rb,7.1447e-07_rb,7.5026e-07_rb, &
     &7.4386e-07_rb,7.2759e-07_rb,7.3583e-07_rb,7.6493e-07_rb,8.8959e-07_rb,7.5534e-07_rb, &
     &5.3734e-07_rb,4.5572e-07_rb,4.1676e-07_rb,3.6198e-07_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 4.67262e-03_rb, 3.95211e-03_rb, 3.34270e-03_rb, 2.82726e-03_rb, 2.39130e-03_rb, &
     & 2.02256e-03_rb, 1.71069e-03_rb, 1.44690e-03_rb, 1.22379e-03_rb, 1.03508e-03_rb/)
      selfrefo(:, 2) = (/ &
     & 4.42593e-03_rb, 3.73338e-03_rb, 3.14920e-03_rb, 2.65643e-03_rb, 2.24076e-03_rb, &
     & 1.89014e-03_rb, 1.59438e-03_rb, 1.34490e-03_rb, 1.13446e-03_rb, 9.56943e-04_rb/)
      selfrefo(:, 3) = (/ &
     & 3.96072e-03_rb, 3.33789e-03_rb, 2.81300e-03_rb, 2.37065e-03_rb, 1.99786e-03_rb, &
     & 1.68369e-03_rb, 1.41893e-03_rb, 1.19580e-03_rb, 1.00776e-03_rb, 8.49286e-04_rb/)
      selfrefo(:, 4) = (/ &
     & 3.71833e-03_rb, 3.10030e-03_rb, 2.58500e-03_rb, 2.15535e-03_rb, 1.79711e-03_rb, &
     & 1.49841e-03_rb, 1.24936e-03_rb, 1.04170e-03_rb, 8.68558e-04_rb, 7.24195e-04_rb/)
      selfrefo(:, 5) = (/ &
     & 3.55755e-03_rb, 2.95355e-03_rb, 2.45210e-03_rb, 2.03578e-03_rb, 1.69015e-03_rb, &
     & 1.40320e-03_rb, 1.16497e-03_rb, 9.67180e-04_rb, 8.02973e-04_rb, 6.66646e-04_rb/)
      selfrefo(:, 6) = (/ &
     & 3.47601e-03_rb, 2.88628e-03_rb, 2.39660e-03_rb, 1.99000e-03_rb, 1.65238e-03_rb, &
     & 1.37204e-03_rb, 1.13927e-03_rb, 9.45980e-04_rb, 7.85487e-04_rb, 6.52224e-04_rb/)
      selfrefo(:, 7) = (/ &
     & 3.44479e-03_rb, 2.86224e-03_rb, 2.37820e-03_rb, 1.97602e-03_rb, 1.64185e-03_rb, &
     & 1.36420e-03_rb, 1.13350e-03_rb, 9.41810e-04_rb, 7.82539e-04_rb, 6.50204e-04_rb/)
      selfrefo(:, 8) = (/ &
     & 3.40154e-03_rb, 2.82953e-03_rb, 2.35370e-03_rb, 1.95789e-03_rb, 1.62864e-03_rb, &
     & 1.35476e-03_rb, 1.12694e-03_rb, 9.37430e-04_rb, 7.79788e-04_rb, 6.48655e-04_rb/)
      selfrefo(:, 9) = (/ &
     & 3.39380e-03_rb, 2.82288e-03_rb, 2.34800e-03_rb, 1.95301e-03_rb, 1.62446e-03_rb, &
     & 1.35119e-03_rb, 1.12389e-03_rb, 9.34820e-04_rb, 7.77560e-04_rb, 6.46755e-04_rb/)
      selfrefo(:,10) = (/ &
     & 3.37185e-03_rb, 2.80654e-03_rb, 2.33600e-03_rb, 1.94435e-03_rb, 1.61837e-03_rb, &
     & 1.34704e-03_rb, 1.12120e-03_rb, 9.33220e-04_rb, 7.76759e-04_rb, 6.46530e-04_rb/)
      selfrefo(:,11) = (/ &
     & 3.37924e-03_rb, 2.81172e-03_rb, 2.33950e-03_rb, 1.94659e-03_rb, 1.61967e-03_rb, &
     & 1.34765e-03_rb, 1.12132e-03_rb, 9.33000e-04_rb, 7.76306e-04_rb, 6.45930e-04_rb/)
      selfrefo(:,12) = (/ &
     & 3.39658e-03_rb, 2.82289e-03_rb, 2.34610e-03_rb, 1.94984e-03_rb, 1.62051e-03_rb, &
     & 1.34680e-03_rb, 1.11933e-03_rb, 9.30270e-04_rb, 7.73146e-04_rb, 6.42561e-04_rb/)
      selfrefo(:,13) = (/ &
     & 3.36070e-03_rb, 2.79913e-03_rb, 2.33140e-03_rb, 1.94183e-03_rb, 1.61735e-03_rb, &
     & 1.34709e-03_rb, 1.12199e-03_rb, 9.34510e-04_rb, 7.78354e-04_rb, 6.48292e-04_rb/)
      selfrefo(:,14) = (/ &
     & 3.40428e-03_rb, 2.81994e-03_rb, 2.33590e-03_rb, 1.93495e-03_rb, 1.60282e-03_rb, &
     & 1.32770e-03_rb, 1.09980e-03_rb, 9.11020e-04_rb, 7.54645e-04_rb, 6.25111e-04_rb/)
      selfrefo(:,15) = (/ &
     & 3.27075e-03_rb, 2.70783e-03_rb, 2.24180e-03_rb, 1.85597e-03_rb, 1.53655e-03_rb, &
     & 1.27210e-03_rb, 1.05317e-03_rb, 8.71910e-04_rb, 7.21849e-04_rb, 5.97615e-04_rb/)
      selfrefo(:,16) = (/ &
     & 3.23123e-03_rb, 2.67891e-03_rb, 2.22100e-03_rb, 1.84136e-03_rb, 1.52661e-03_rb, &
     & 1.26567e-03_rb, 1.04932e-03_rb, 8.69960e-04_rb, 7.21256e-04_rb, 5.97970e-04_rb/)

      end subroutine lw_kgb144

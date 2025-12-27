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
!!prev      subroutine lw_kgb04
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb049
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg04, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg04, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &3.3839e-04_rb,2.4739e-04_rb,2.2846e-04_rb,2.3376e-04_rb,2.2622e-04_rb,2.3188e-04_rb, &
     &2.2990e-04_rb,2.2532e-04_rb,2.1233e-04_rb,2.0593e-04_rb,2.0716e-04_rb,2.0809e-04_rb, &
     &2.0889e-04_rb,2.0932e-04_rb,2.0944e-04_rb,2.0945e-04_rb/)
      forrefo(2,:) = (/ &
     &3.4391e-04_rb,2.6022e-04_rb,2.3449e-04_rb,2.4544e-04_rb,2.3831e-04_rb,2.3014e-04_rb, &
     &2.3729e-04_rb,2.2726e-04_rb,2.1892e-04_rb,1.9223e-04_rb,2.1291e-04_rb,2.1406e-04_rb, &
     &2.1491e-04_rb,2.1548e-04_rb,2.1562e-04_rb,2.1567e-04_rb/)
      forrefo(3,:) = (/ &
     &3.4219e-04_rb,2.7334e-04_rb,2.3727e-04_rb,2.4515e-04_rb,2.5272e-04_rb,2.4212e-04_rb, &
     &2.3824e-04_rb,2.3615e-04_rb,2.2724e-04_rb,2.2381e-04_rb,1.9634e-04_rb,2.1625e-04_rb, &
     &2.1963e-04_rb,2.2032e-04_rb,2.2057e-04_rb,2.2058e-04_rb/)
      forrefo(4,:) = (/ &
     &3.1684e-04_rb,2.4823e-04_rb,2.4890e-04_rb,2.4577e-04_rb,2.4106e-04_rb,2.4353e-04_rb, &
     &2.4038e-04_rb,2.3932e-04_rb,2.3604e-04_rb,2.3773e-04_rb,2.4243e-04_rb,2.2597e-04_rb, &
     &2.2879e-04_rb,2.2440e-04_rb,2.1104e-04_rb,2.1460e-04_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 2.62922e-01_rb, 2.29106e-01_rb, 1.99640e-01_rb, 1.73964e-01_rb, 1.51589e-01_rb, &
     & 1.32093e-01_rb, 1.15104e-01_rb, 1.00300e-01_rb, 8.74000e-02_rb, 7.61592e-02_rb/)
      selfrefo(:, 2) = (/ &
     & 2.45448e-01_rb, 2.13212e-01_rb, 1.85210e-01_rb, 1.60886e-01_rb, 1.39756e-01_rb, &
     & 1.21401e-01_rb, 1.05457e-01_rb, 9.16070e-02_rb, 7.95759e-02_rb, 6.91249e-02_rb/)
      selfrefo(:, 3) = (/ &
     & 2.41595e-01_rb, 2.09697e-01_rb, 1.82010e-01_rb, 1.57979e-01_rb, 1.37121e-01_rb, &
     & 1.19016e-01_rb, 1.03302e-01_rb, 8.96630e-02_rb, 7.78246e-02_rb, 6.75492e-02_rb/)
      selfrefo(:, 4) = (/ &
     & 2.44818e-01_rb, 2.12172e-01_rb, 1.83880e-01_rb, 1.59360e-01_rb, 1.38110e-01_rb, &
     & 1.19694e-01_rb, 1.03733e-01_rb, 8.99010e-02_rb, 7.79131e-02_rb, 6.75238e-02_rb/)
      selfrefo(:, 5) = (/ &
     & 2.43458e-01_rb, 2.10983e-01_rb, 1.82840e-01_rb, 1.58451e-01_rb, 1.37315e-01_rb, &
     & 1.18998e-01_rb, 1.03125e-01_rb, 8.93690e-02_rb, 7.74480e-02_rb, 6.71171e-02_rb/)
      selfrefo(:, 6) = (/ &
     & 2.40186e-01_rb, 2.08745e-01_rb, 1.81420e-01_rb, 1.57672e-01_rb, 1.37032e-01_rb, &
     & 1.19095e-01_rb, 1.03505e-01_rb, 8.99560e-02_rb, 7.81806e-02_rb, 6.79467e-02_rb/)
      selfrefo(:, 7) = (/ &
     & 2.42752e-01_rb, 2.10579e-01_rb, 1.82670e-01_rb, 1.58460e-01_rb, 1.37459e-01_rb, &
     & 1.19240e-01_rb, 1.03437e-01_rb, 8.97280e-02_rb, 7.78359e-02_rb, 6.75200e-02_rb/)
      selfrefo(:, 8) = (/ &
     & 2.39620e-01_rb, 2.08166e-01_rb, 1.80840e-01_rb, 1.57101e-01_rb, 1.36479e-01_rb, &
     & 1.18563e-01_rb, 1.03000e-01_rb, 8.94790e-02_rb, 7.77332e-02_rb, 6.75292e-02_rb/)
      selfrefo(:, 9) = (/ &
     & 2.38856e-01_rb, 2.07166e-01_rb, 1.79680e-01_rb, 1.55841e-01_rb, 1.35165e-01_rb, &
     & 1.17232e-01_rb, 1.01678e-01_rb, 8.81880e-02_rb, 7.64877e-02_rb, 6.63397e-02_rb/)
      selfrefo(:,10) = (/ &
     & 2.29821e-01_rb, 2.00586e-01_rb, 1.75070e-01_rb, 1.52800e-01_rb, 1.33363e-01_rb, &
     & 1.16398e-01_rb, 1.01591e-01_rb, 8.86680e-02_rb, 7.73887e-02_rb, 6.75443e-02_rb/)
      selfrefo(:,11) = (/ &
     & 2.39945e-01_rb, 2.08186e-01_rb, 1.80630e-01_rb, 1.56722e-01_rb, 1.35978e-01_rb, &
     & 1.17980e-01_rb, 1.02364e-01_rb, 8.88150e-02_rb, 7.70594e-02_rb, 6.68598e-02_rb/)
      selfrefo(:,12) = (/ &
     & 2.40271e-01_rb, 2.08465e-01_rb, 1.80870e-01_rb, 1.56927e-01_rb, 1.36154e-01_rb, &
     & 1.18131e-01_rb, 1.02494e-01_rb, 8.89260e-02_rb, 7.71545e-02_rb, 6.69412e-02_rb/)
      selfrefo(:,13) = (/ &
     & 2.40503e-01_rb, 2.08670e-01_rb, 1.81050e-01_rb, 1.57086e-01_rb, 1.36294e-01_rb, &
     & 1.18254e-01_rb, 1.02602e-01_rb, 8.90210e-02_rb, 7.72380e-02_rb, 6.70147e-02_rb/)
      selfrefo(:,14) = (/ &
     & 2.40670e-01_rb, 2.08811e-01_rb, 1.81170e-01_rb, 1.57188e-01_rb, 1.36380e-01_rb, &
     & 1.18327e-01_rb, 1.02663e-01_rb, 8.90730e-02_rb, 7.72819e-02_rb, 6.70517e-02_rb/)
      selfrefo(:,15) = (/ &
     & 2.40711e-01_rb, 2.08846e-01_rb, 1.81200e-01_rb, 1.57213e-01_rb, 1.36402e-01_rb, &
     & 1.18346e-01_rb, 1.02679e-01_rb, 8.90870e-02_rb, 7.72939e-02_rb, 6.70621e-02_rb/)
      selfrefo(:,16) = (/ &
     & 2.40727e-01_rb, 2.08859e-01_rb, 1.81210e-01_rb, 1.57221e-01_rb, 1.36408e-01_rb, &
     & 1.18350e-01_rb, 1.02682e-01_rb, 8.90890e-02_rb, 7.72952e-02_rb, 6.70627e-02_rb/)

      end subroutine lw_kgb049

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
!!prev      subroutine lw_kgb02
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:), forrefo(:)
!!
      subroutine lw_kgb024
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg02, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg02, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &2.8549e-03_rb,4.8281e-03_rb,6.2570e-03_rb,8.2731e-03_rb,7.9056e-03_rb,7.7840e-03_rb, &
     &1.0115e-02_rb,9.6599e-03_rb,1.0153e-02_rb,1.0921e-02_rb,1.2408e-02_rb,1.3496e-02_rb, &
     &1.5059e-02_rb,1.4636e-02_rb,1.6483e-02_rb,1.2394e-02_rb/)
      forrefo(2,:) = (/ &
     &3.0036e-03_rb,5.1093e-03_rb,5.7317e-03_rb,9.2246e-03_rb,8.9829e-03_rb,8.6477e-03_rb, &
     &1.1448e-02_rb,1.0391e-02_rb,1.0211e-02_rb,1.2921e-02_rb,1.2726e-02_rb,1.2426e-02_rb, &
     &1.4609e-02_rb,1.5783e-02_rb,1.6617e-02_rb,1.6858e-02_rb/)
      forrefo(3,:) = (/ &
     &3.0771e-03_rb,5.1206e-03_rb,5.8426e-03_rb,9.5727e-03_rb,1.0338e-02_rb,9.3737e-03_rb, &
     &1.2805e-02_rb,1.1272e-02_rb,1.1353e-02_rb,1.1837e-02_rb,1.1550e-02_rb,1.3020e-02_rb, &
     &1.3536e-02_rb,1.6226e-02_rb,1.6039e-02_rb,2.2578e-02_rb/)
      forrefo(4,:) = (/ &
     &3.3072e-03_rb,5.0240e-03_rb,6.8474e-03_rb,8.2736e-03_rb,8.6151e-03_rb,8.6762e-03_rb, &
     &1.1476e-02_rb,1.0246e-02_rb,1.0819e-02_rb,1.0640e-02_rb,1.0545e-02_rb,1.0533e-02_rb, &
     &1.0496e-02_rb,1.0142e-02_rb,9.7979e-03_rb,1.5255e-02_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 7.25695e-01_rb, 6.53591e-01_rb, 5.88650e-01_rb, 5.30162e-01_rb, 4.77485e-01_rb, &
     & 4.30042e-01_rb, 3.87313e-01_rb, 3.48830e-01_rb, 3.14170e-01_rb, 2.82954e-01_rb/)
      selfrefo(:, 2) = (/ &
     & 9.61996e-01_rb, 8.77853e-01_rb, 8.01070e-01_rb, 7.31003e-01_rb, 6.67064e-01_rb, &
     & 6.08718e-01_rb, 5.55476e-01_rb, 5.06890e-01_rb, 4.62554e-01_rb, 4.22096e-01_rb/)
      selfrefo(:, 3) = (/ &
     & 9.72584e-01_rb, 9.02658e-01_rb, 8.37760e-01_rb, 7.77527e-01_rb, 7.21626e-01_rb, &
     & 6.69743e-01_rb, 6.21591e-01_rb, 5.76900e-01_rb, 5.35423e-01_rb, 4.96927e-01_rb/)
      selfrefo(:, 4) = (/ &
     & 1.24790e+00_rb, 1.14353e+00_rb, 1.04790e+00_rb, 9.60263e-01_rb, 8.79956e-01_rb, &
     & 8.06364e-01_rb, 7.38927e-01_rb, 6.77130e-01_rb, 6.20501e-01_rb, 5.68608e-01_rb/)
      selfrefo(:, 5) = (/ &
     & 1.23574e+00_rb, 1.12928e+00_rb, 1.03200e+00_rb, 9.43096e-01_rb, 8.61851e-01_rb, &
     & 7.87605e-01_rb, 7.19755e-01_rb, 6.57750e-01_rb, 6.01087e-01_rb, 5.49305e-01_rb/)
      selfrefo(:, 6) = (/ &
     & 1.20921e+00_rb, 1.10660e+00_rb, 1.01270e+00_rb, 9.26766e-01_rb, 8.48124e-01_rb, &
     & 7.76155e-01_rb, 7.10293e-01_rb, 6.50020e-01_rb, 5.94861e-01_rb, 5.44384e-01_rb/)
      selfrefo(:, 7) = (/ &
     & 1.38112e+00_rb, 1.26727e+00_rb, 1.16280e+00_rb, 1.06694e+00_rb, 9.78990e-01_rb, &
     & 8.98287e-01_rb, 8.24236e-01_rb, 7.56290e-01_rb, 6.93945e-01_rb, 6.36739e-01_rb/)
      selfrefo(:, 8) = (/ &
     & 1.30321e+00_rb, 1.20127e+00_rb, 1.10730e+00_rb, 1.02068e+00_rb, 9.40840e-01_rb, &
     & 8.67243e-01_rb, 7.99403e-01_rb, 7.36870e-01_rb, 6.79229e-01_rb, 6.26096e-01_rb/)
      selfrefo(:, 9) = (/ &
     & 1.26713e+00_rb, 1.17927e+00_rb, 1.09750e+00_rb, 1.02140e+00_rb, 9.50575e-01_rb, &
     & 8.84662e-01_rb, 8.23319e-01_rb, 7.66230e-01_rb, 7.13099e-01_rb, 6.63653e-01_rb/)
      selfrefo(:,10) = (/ &
     & 1.49824e+00_rb, 1.37053e+00_rb, 1.25370e+00_rb, 1.14683e+00_rb, 1.04908e+00_rb, &
     & 9.59651e-01_rb, 8.77849e-01_rb, 8.03020e-01_rb, 7.34569e-01_rb, 6.71954e-01_rb/)
      selfrefo(:,11) = (/ &
     & 1.44786e+00_rb, 1.34594e+00_rb, 1.25120e+00_rb, 1.16313e+00_rb, 1.08125e+00_rb, &
     & 1.00514e+00_rb, 9.34392e-01_rb, 8.68620e-01_rb, 8.07477e-01_rb, 7.50639e-01_rb/)
      selfrefo(:,12) = (/ &
     & 1.38460e+00_rb, 1.30437e+00_rb, 1.22880e+00_rb, 1.15760e+00_rb, 1.09053e+00_rb, &
     & 1.02735e+00_rb, 9.67825e-01_rb, 9.11750e-01_rb, 8.58924e-01_rb, 8.09159e-01_rb/)
      selfrefo(:,13) = (/ &
     & 1.51953e+00_rb, 1.42822e+00_rb, 1.34240e+00_rb, 1.26173e+00_rb, 1.18592e+00_rb, &
     & 1.11465e+00_rb, 1.04768e+00_rb, 9.84720e-01_rb, 9.25548e-01_rb, 8.69932e-01_rb/)
      selfrefo(:,14) = (/ &
     & 1.62608e+00_rb, 1.51021e+00_rb, 1.40260e+00_rb, 1.30266e+00_rb, 1.20983e+00_rb, &
     & 1.12363e+00_rb, 1.04356e+00_rb, 9.69200e-01_rb, 9.00138e-01_rb, 8.35998e-01_rb/)
      selfrefo(:,15) = (/ &
     & 1.65383e+00_rb, 1.54808e+00_rb, 1.44910e+00_rb, 1.35644e+00_rb, 1.26971e+00_rb, &
     & 1.18853e+00_rb, 1.11254e+00_rb, 1.04140e+00_rb, 9.74813e-01_rb, 9.12484e-01_rb/)
      selfrefo(:,16) = (/ &
     & 1.78105e+00_rb, 1.61421e+00_rb, 1.46300e+00_rb, 1.32595e+00_rb, 1.20174e+00_rb, &
     & 1.08917e+00_rb, 9.87141e-01_rb, 8.94670e-01_rb, 8.10861e-01_rb, 7.34904e-01_rb/)

      end subroutine lw_kgb024

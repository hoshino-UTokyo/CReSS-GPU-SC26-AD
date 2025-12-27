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
!!prev      subroutine lw_kgb09
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb099
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg09, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
!!prev                            kbo_mn2o, selfrefo, forrefo
      use rrlw_kg09, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &7.5352e-06_rb,2.9812e-05_rb,1.4497e-04_rb,4.4006e-04_rb,1.0492e-03_rb,1.9676e-03_rb, &
     &1.9989e-03_rb,1.9099e-03_rb,2.2121e-03_rb,2.4491e-03_rb,2.9573e-03_rb,2.6344e-03_rb, &
     &3.0629e-03_rb,3.3547e-03_rb,5.0643e-03_rb,5.0642e-03_rb/)
      forrefo(2,:) = (/ &
     &6.6070e-06_rb,4.8618e-05_rb,3.1112e-04_rb,8.4235e-04_rb,1.4179e-03_rb,1.4315e-03_rb, &
     &1.4685e-03_rb,1.6554e-03_rb,2.1171e-03_rb,2.3545e-03_rb,2.5165e-03_rb,2.7680e-03_rb, &
     &2.6985e-03_rb,3.5345e-03_rb,4.2924e-03_rb,5.0712e-03_rb/)
      forrefo(3,:) = (/ &
     &6.5962e-06_rb,7.2595e-04_rb,1.3429e-03_rb,1.1675e-03_rb,9.8384e-04_rb,8.8787e-04_rb, &
     &8.7557e-04_rb,8.0589e-04_rb,7.7024e-04_rb,8.7518e-04_rb,9.5213e-04_rb,9.0849e-04_rb, &
     &1.2596e-03_rb,2.5106e-03_rb,3.9471e-03_rb,5.0742e-03_rb/)
      forrefo(4,:) = (/ &
     &3.6217e-04_rb,1.0709e-03_rb,1.0628e-03_rb,8.5640e-04_rb,8.9332e-04_rb,8.3372e-04_rb, &
     &7.8539e-04_rb,8.2828e-04_rb,8.3329e-04_rb,8.5118e-04_rb,8.2878e-04_rb,6.8570e-04_rb, &
     &6.3815e-04_rb,8.0648e-04_rb,2.3236e-03_rb,4.0321e-03_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 2.83453e-02_rb, 2.51439e-02_rb, 2.23040e-02_rb, 1.97849e-02_rb, 1.75503e-02_rb, &
     & 1.55681e-02_rb, 1.38097e-02_rb, 1.22500e-02_rb, 1.08664e-02_rb, 9.63912e-03_rb/)
      selfrefo(:, 2) = (/ &
     & 3.05185e-02_rb, 2.72374e-02_rb, 2.43090e-02_rb, 2.16955e-02_rb, 1.93629e-02_rb, &
     & 1.72811e-02_rb, 1.54232e-02_rb, 1.37650e-02_rb, 1.22851e-02_rb, 1.09643e-02_rb/)
      selfrefo(:, 3) = (/ &
     & 4.23833e-02_rb, 3.76250e-02_rb, 3.34010e-02_rb, 2.96512e-02_rb, 2.63223e-02_rb, &
     & 2.33672e-02_rb, 2.07439e-02_rb, 1.84150e-02_rb, 1.63476e-02_rb, 1.45123e-02_rb/)
      selfrefo(:, 4) = (/ &
     & 5.76481e-02_rb, 5.13686e-02_rb, 4.57730e-02_rb, 4.07870e-02_rb, 3.63441e-02_rb, &
     & 3.23851e-02_rb, 2.88574e-02_rb, 2.57140e-02_rb, 2.29130e-02_rb, 2.04171e-02_rb/)
      selfrefo(:, 5) = (/ &
     & 6.92255e-02_rb, 6.33521e-02_rb, 5.79770e-02_rb, 5.30580e-02_rb, 4.85563e-02_rb, &
     & 4.44365e-02_rb, 4.06663e-02_rb, 3.72160e-02_rb, 3.40584e-02_rb, 3.11687e-02_rb/)
      selfrefo(:, 6) = (/ &
     & 6.07694e-02_rb, 5.94182e-02_rb, 5.80970e-02_rb, 5.68052e-02_rb, 5.55422e-02_rb, &
     & 5.43072e-02_rb, 5.30997e-02_rb, 5.19190e-02_rb, 5.07646e-02_rb, 4.96358e-02_rb/)
      selfrefo(:, 7) = (/ &
     & 6.23749e-02_rb, 6.07744e-02_rb, 5.92150e-02_rb, 5.76956e-02_rb, 5.62152e-02_rb, &
     & 5.47728e-02_rb, 5.33674e-02_rb, 5.19980e-02_rb, 5.06638e-02_rb, 4.93638e-02_rb/)
      selfrefo(:, 8) = (/ &
     & 6.90744e-02_rb, 6.61811e-02_rb, 6.34090e-02_rb, 6.07530e-02_rb, 5.82083e-02_rb, &
     & 5.57702e-02_rb, 5.34342e-02_rb, 5.11960e-02_rb, 4.90516e-02_rb, 4.69970e-02_rb/)
      selfrefo(:, 9) = (/ &
     & 8.08992e-02_rb, 7.68876e-02_rb, 7.30750e-02_rb, 6.94514e-02_rb, 6.60075e-02_rb, &
     & 6.27344e-02_rb, 5.96236e-02_rb, 5.66670e-02_rb, 5.38570e-02_rb, 5.11864e-02_rb/)
      selfrefo(:,10) = (/ &
     & 8.70197e-02_rb, 8.27485e-02_rb, 7.86870e-02_rb, 7.48248e-02_rb, 7.11522e-02_rb, &
     & 6.76599e-02_rb, 6.43389e-02_rb, 6.11810e-02_rb, 5.81781e-02_rb, 5.53225e-02_rb/)
      selfrefo(:,11) = (/ &
     & 8.84776e-02_rb, 8.54262e-02_rb, 8.24800e-02_rb, 7.96354e-02_rb, 7.68890e-02_rb, &
     & 7.42373e-02_rb, 7.16770e-02_rb, 6.92050e-02_rb, 6.68183e-02_rb, 6.45139e-02_rb/)
      selfrefo(:,12) = (/ &
     & 9.82552e-02_rb, 9.25696e-02_rb, 8.72130e-02_rb, 8.21664e-02_rb, 7.74118e-02_rb, &
     & 7.29323e-02_rb, 6.87121e-02_rb, 6.47360e-02_rb, 6.09900e-02_rb, 5.74608e-02_rb/)
      selfrefo(:,13) = (/ &
     & 9.32447e-02_rb, 8.96818e-02_rb, 8.62550e-02_rb, 8.29592e-02_rb, 7.97893e-02_rb, &
     & 7.67405e-02_rb, 7.38082e-02_rb, 7.09880e-02_rb, 6.82755e-02_rb, 6.56667e-02_rb/)
      selfrefo(:,14) = (/ &
     & 1.15363e-01_rb, 1.08593e-01_rb, 1.02220e-01_rb, 9.62210e-02_rb, 9.05741e-02_rb, &
     & 8.52585e-02_rb, 8.02549e-02_rb, 7.55450e-02_rb, 7.11115e-02_rb, 6.69382e-02_rb/)
      selfrefo(:,15) = (/ &
     & 1.23179e-01_rb, 1.19247e-01_rb, 1.15440e-01_rb, 1.11755e-01_rb, 1.08187e-01_rb, &
     & 1.04734e-01_rb, 1.01391e-01_rb, 9.81540e-02_rb, 9.50207e-02_rb, 9.19875e-02_rb/)
      selfrefo(:,16) = (/ &
     & 1.44104e-01_rb, 1.36412e-01_rb, 1.29130e-01_rb, 1.22237e-01_rb, 1.15712e-01_rb, &
     & 1.09535e-01_rb, 1.03688e-01_rb, 9.81530e-02_rb, 9.29135e-02_rb, 8.79537e-02_rb/)

      end subroutine lw_kgb099

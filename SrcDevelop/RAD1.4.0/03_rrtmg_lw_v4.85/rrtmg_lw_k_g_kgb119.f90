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
!!prev      subroutine lw_kgb11
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb119
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg11, only : fracrefao, fracrefbo, kao, kbo, kao_mo2, &
!!prev                            kbo_mo2, selfrefo, forrefo
      use rrlw_kg11, only :  selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &2.8858e-02_rb,3.6879e-02_rb,4.0746e-02_rb,4.2561e-02_rb,4.2740e-02_rb,4.2707e-02_rb, &
     &4.4109e-02_rb,4.4540e-02_rb,4.5206e-02_rb,4.4679e-02_rb,4.5034e-02_rb,4.5364e-02_rb, &
     &4.6790e-02_rb,4.7857e-02_rb,4.8328e-02_rb,4.8084e-02_rb/)
      forrefo(2,:) = (/ &
     &2.7887e-02_rb,3.7376e-02_rb,4.0980e-02_rb,4.2986e-02_rb,4.3054e-02_rb,4.2975e-02_rb, &
     &4.3754e-02_rb,4.4352e-02_rb,4.4723e-02_rb,4.6236e-02_rb,4.5273e-02_rb,4.5360e-02_rb, &
     &4.5332e-02_rb,4.7587e-02_rb,4.7035e-02_rb,5.0267e-02_rb/)
      forrefo(3,:) = (/ &
     &2.5846e-02_rb,3.6753e-02_rb,4.2334e-02_rb,4.3806e-02_rb,4.3848e-02_rb,4.3215e-02_rb, &
     &4.3838e-02_rb,4.4278e-02_rb,4.4658e-02_rb,4.5403e-02_rb,4.5255e-02_rb,4.6347e-02_rb, &
     &4.4722e-02_rb,4.6612e-02_rb,4.6836e-02_rb,4.8720e-02_rb/)
      forrefo(4,:) = (/ &
     &2.8955e-02_rb,3.7608e-02_rb,4.1989e-02_rb,4.4919e-02_rb,4.2803e-02_rb,4.2842e-02_rb, &
     &4.2632e-02_rb,4.1056e-02_rb,4.0086e-02_rb,4.1401e-02_rb,4.2746e-02_rb,4.2142e-02_rb, &
     &4.1871e-02_rb,4.3917e-02_rb,4.5462e-02_rb,4.8359e-02_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 5.96496e-01_rb, 5.49171e-01_rb, 5.05600e-01_rb, 4.65486e-01_rb, 4.28555e-01_rb, &
     & 3.94554e-01_rb, 3.63250e-01_rb, 3.34430e-01_rb, 3.07897e-01_rb, 2.83468e-01_rb/)
      selfrefo(:, 2) = (/ &
     & 7.46455e-01_rb, 6.82459e-01_rb, 6.23950e-01_rb, 5.70457e-01_rb, 5.21550e-01_rb, &
     & 4.76836e-01_rb, 4.35956e-01_rb, 3.98580e-01_rb, 3.64409e-01_rb, 3.33167e-01_rb/)
      selfrefo(:, 3) = (/ &
     & 7.86805e-01_rb, 7.21186e-01_rb, 6.61040e-01_rb, 6.05910e-01_rb, 5.55378e-01_rb, &
     & 5.09059e-01_rb, 4.66605e-01_rb, 4.27690e-01_rb, 3.92021e-01_rb, 3.59327e-01_rb/)
      selfrefo(:, 4) = (/ &
     & 8.11740e-01_rb, 7.44359e-01_rb, 6.82570e-01_rb, 6.25910e-01_rb, 5.73954e-01_rb, &
     & 5.26311e-01_rb, 4.82622e-01_rb, 4.42560e-01_rb, 4.05823e-01_rb, 3.72136e-01_rb/)
      selfrefo(:, 5) = (/ &
     & 8.14870e-01_rb, 7.47200e-01_rb, 6.85150e-01_rb, 6.28253e-01_rb, 5.76081e-01_rb, &
     & 5.28241e-01_rb, 4.84374e-01_rb, 4.44150e-01_rb, 4.07266e-01_rb, 3.73446e-01_rb/)
      selfrefo(:, 6) = (/ &
     & 8.10104e-01_rb, 7.43259e-01_rb, 6.81930e-01_rb, 6.25661e-01_rb, 5.74035e-01_rb, &
     & 5.26669e-01_rb, 4.83212e-01_rb, 4.43340e-01_rb, 4.06758e-01_rb, 3.73195e-01_rb/)
      selfrefo(:, 7) = (/ &
     & 8.13119e-01_rb, 7.48127e-01_rb, 6.88330e-01_rb, 6.33312e-01_rb, 5.82692e-01_rb, &
     & 5.36118e-01_rb, 4.93267e-01_rb, 4.53840e-01_rb, 4.17565e-01_rb, 3.84189e-01_rb/)
      selfrefo(:, 8) = (/ &
     & 8.26137e-01_rb, 7.58984e-01_rb, 6.97290e-01_rb, 6.40611e-01_rb, 5.88539e-01_rb, &
     & 5.40699e-01_rb, 4.96748e-01_rb, 4.56370e-01_rb, 4.19274e-01_rb, 3.85193e-01_rb/)
      selfrefo(:, 9) = (/ &
     & 8.30566e-01_rb, 7.63984e-01_rb, 7.02740e-01_rb, 6.46405e-01_rb, 5.94587e-01_rb, &
     & 5.46922e-01_rb, 5.03079e-01_rb, 4.62750e-01_rb, 4.25654e-01_rb, 3.91532e-01_rb/)
      selfrefo(:,10) = (/ &
     & 8.67471e-01_rb, 7.91575e-01_rb, 7.22320e-01_rb, 6.59124e-01_rb, 6.01457e-01_rb, &
     & 5.48835e-01_rb, 5.00817e-01_rb, 4.57000e-01_rb, 4.17017e-01_rb, 3.80532e-01_rb/)
      selfrefo(:,11) = (/ &
     & 8.51029e-01_rb, 7.79373e-01_rb, 7.13750e-01_rb, 6.53652e-01_rb, 5.98615e-01_rb, &
     & 5.48212e-01_rb, 5.02053e-01_rb, 4.59780e-01_rb, 4.21067e-01_rb, 3.85613e-01_rb/)
      selfrefo(:,12) = (/ &
     & 8.36772e-01_rb, 7.68751e-01_rb, 7.06260e-01_rb, 6.48848e-01_rb, 5.96104e-01_rb, &
     & 5.47647e-01_rb, 5.03129e-01_rb, 4.62230e-01_rb, 4.24655e-01_rb, 3.90136e-01_rb/)
      selfrefo(:,13) = (/ &
     & 8.36551e-01_rb, 7.71089e-01_rb, 7.10750e-01_rb, 6.55133e-01_rb, 6.03867e-01_rb, &
     & 5.56614e-01_rb, 5.13058e-01_rb, 4.72910e-01_rb, 4.35904e-01_rb, 4.01794e-01_rb/)
      selfrefo(:,14) = (/ &
     & 8.84307e-01_rb, 8.11175e-01_rb, 7.44090e-01_rb, 6.82553e-01_rb, 6.26106e-01_rb, &
     & 5.74326e-01_rb, 5.26829e-01_rb, 4.83260e-01_rb, 4.43294e-01_rb, 4.06633e-01_rb/)
      selfrefo(:,15) = (/ &
     & 8.90356e-01_rb, 8.19830e-01_rb, 7.54890e-01_rb, 6.95094e-01_rb, 6.40035e-01_rb, &
     & 5.89337e-01_rb, 5.42655e-01_rb, 4.99670e-01_rb, 4.60090e-01_rb, 4.23646e-01_rb/)
      selfrefo(:,16) = (/ &
     & 9.67549e-01_rb, 8.79393e-01_rb, 7.99270e-01_rb, 7.26447e-01_rb, 6.60259e-01_rb, &
     & 6.00101e-01_rb, 5.45425e-01_rb, 4.95730e-01_rb, 4.50563e-01_rb, 4.09511e-01_rb/)

      end subroutine lw_kgb119

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
!!prev      subroutine lw_kgb10
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for kao(:,:,5:8)
!!
      subroutine lw_kgb102
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg10, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg10, only : kao

      implicit none
      save

!     The array KAO contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  The first
!     index in the array, JT, which runs from 1 to 5, corresponds to 
!     different temperatures.  More specifically, JT = 3 means that the 
!     data are for the corresponding TREF for this  pressure level, 
!     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
!     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
!     index, JP, runs from 1 to 13 and refers to the corresponding 
!     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
!     The third index, IG, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.

      kao(:, 1, 5) = (/ &
     &2.5886e-01_rb,2.7100e-01_rb,2.8220e-01_rb,2.9245e-01_rb,3.0214e-01_rb/)
      kao(:, 2, 5) = (/ &
     &2.0909e-01_rb,2.1918e-01_rb,2.2861e-01_rb,2.3740e-01_rb,2.4610e-01_rb/)
      kao(:, 3, 5) = (/ &
     &1.6622e-01_rb,1.7497e-01_rb,1.8311e-01_rb,1.9088e-01_rb,1.9860e-01_rb/)
      kao(:, 4, 5) = (/ &
     &1.3301e-01_rb,1.4069e-01_rb,1.4792e-01_rb,1.5478e-01_rb,1.6132e-01_rb/)
      kao(:, 5, 5) = (/ &
     &1.0688e-01_rb,1.1355e-01_rb,1.1989e-01_rb,1.2584e-01_rb,1.3154e-01_rb/)
      kao(:, 6, 5) = (/ &
     &8.5419e-02_rb,9.1274e-02_rb,9.6848e-02_rb,1.0202e-01_rb,1.0704e-01_rb/)
      kao(:, 7, 5) = (/ &
     &6.8100e-02_rb,7.3262e-02_rb,7.8047e-02_rb,8.2594e-02_rb,8.7027e-02_rb/)
      kao(:, 8, 5) = (/ &
     &5.4252e-02_rb,5.8815e-02_rb,6.2989e-02_rb,6.6937e-02_rb,7.0783e-02_rb/)
      kao(:, 9, 5) = (/ &
     &4.3117e-02_rb,4.7224e-02_rb,5.0916e-02_rb,5.4328e-02_rb,5.7817e-02_rb/)
      kao(:,10, 5) = (/ &
     &2.9151e-02_rb,3.5667e-02_rb,4.0092e-02_rb,4.3888e-02_rb,4.6870e-02_rb/)
      kao(:,11, 5) = (/ &
     &2.3941e-02_rb,2.7450e-02_rb,2.9750e-02_rb,3.5348e-02_rb,3.9440e-02_rb/)
      kao(:,12, 5) = (/ &
     &2.0570e-02_rb,2.3238e-02_rb,2.5236e-02_rb,2.7327e-02_rb,3.2816e-02_rb/)
      kao(:,13, 5) = (/ &
     &1.7253e-02_rb,1.9755e-02_rb,2.1401e-02_rb,2.3591e-02_rb,2.8114e-02_rb/)
      kao(:, 1, 6) = (/ &
     &5.6993e-01_rb,6.0360e-01_rb,6.3489e-01_rb,6.6230e-01_rb,6.8659e-01_rb/)
      kao(:, 2, 6) = (/ &
     &4.6502e-01_rb,4.9401e-01_rb,5.2007e-01_rb,5.4304e-01_rb,5.6306e-01_rb/)
      kao(:, 3, 6) = (/ &
     &3.7271e-01_rb,3.9757e-01_rb,4.1967e-01_rb,4.3927e-01_rb,4.5631e-01_rb/)
      kao(:, 4, 6) = (/ &
     &2.9874e-01_rb,3.1990e-01_rb,3.3890e-01_rb,3.5594e-01_rb,3.7139e-01_rb/)
      kao(:, 5, 6) = (/ &
     &2.3859e-01_rb,2.5711e-01_rb,2.7375e-01_rb,2.8893e-01_rb,3.0280e-01_rb/)
      kao(:, 6, 6) = (/ &
     &1.8917e-01_rb,2.0549e-01_rb,2.2017e-01_rb,2.3369e-01_rb,2.4611e-01_rb/)
      kao(:, 7, 6) = (/ &
     &1.4938e-01_rb,1.6367e-01_rb,1.7681e-01_rb,1.8887e-01_rb,1.9989e-01_rb/)
      kao(:, 8, 6) = (/ &
     &1.1758e-01_rb,1.3011e-01_rb,1.4157e-01_rb,1.5226e-01_rb,1.6197e-01_rb/)
      kao(:, 9, 6) = (/ &
     &9.1609e-02_rb,1.0230e-01_rb,1.1226e-01_rb,1.2154e-01_rb,1.2990e-01_rb/)
      kao(:,10, 6) = (/ &
     &7.2629e-02_rb,8.2527e-02_rb,9.1786e-02_rb,1.0065e-01_rb,1.0832e-01_rb/)
      kao(:,11, 6) = (/ &
     &6.0173e-02_rb,6.9095e-02_rb,7.8339e-02_rb,8.6340e-02_rb,9.2996e-02_rb/)
      kao(:,12, 6) = (/ &
     &5.1074e-02_rb,5.9514e-02_rb,6.6969e-02_rb,7.3447e-02_rb,7.9277e-02_rb/)
      kao(:,13, 6) = (/ &
     &4.3826e-02_rb,5.0848e-02_rb,5.7125e-02_rb,6.2916e-02_rb,6.7915e-02_rb/)
      kao(:, 1, 7) = (/ &
     &1.4129e+00_rb,1.4750e+00_rb,1.5277e+00_rb,1.5727e+00_rb,1.6117e+00_rb/)
      kao(:, 2, 7) = (/ &
     &1.1749e+00_rb,1.2287e+00_rb,1.2750e+00_rb,1.3159e+00_rb,1.3523e+00_rb/)
      kao(:, 3, 7) = (/ &
     &9.5348e-01_rb,1.0008e+00_rb,1.0433e+00_rb,1.0822e+00_rb,1.1168e+00_rb/)
      kao(:, 4, 7) = (/ &
     &7.6995e-01_rb,8.1381e-01_rb,8.5421e-01_rb,8.9012e-01_rb,9.2111e-01_rb/)
      kao(:, 5, 7) = (/ &
     &6.2211e-01_rb,6.6255e-01_rb,6.9934e-01_rb,7.3149e-01_rb,7.5879e-01_rb/)
      kao(:, 6, 7) = (/ &
     &5.0056e-01_rb,5.3742e-01_rb,5.7038e-01_rb,5.9882e-01_rb,6.2305e-01_rb/)
      kao(:, 7, 7) = (/ &
     &4.0115e-01_rb,4.3352e-01_rb,4.6244e-01_rb,4.8765e-01_rb,5.0920e-01_rb/)
      kao(:, 8, 7) = (/ &
     &3.1902e-01_rb,3.4754e-01_rb,3.7305e-01_rb,3.9544e-01_rb,4.1482e-01_rb/)
      kao(:, 9, 7) = (/ &
     &2.5207e-01_rb,2.7718e-01_rb,2.9962e-01_rb,3.1933e-01_rb,3.3666e-01_rb/)
      kao(:,10, 7) = (/ &
     &1.9498e-01_rb,2.1588e-01_rb,2.3464e-01_rb,2.5127e-01_rb,2.6668e-01_rb/)
      kao(:,11, 7) = (/ &
     &1.6583e-01_rb,1.8332e-01_rb,1.9835e-01_rb,2.1191e-01_rb,2.2523e-01_rb/)
      kao(:,12, 7) = (/ &
     &1.3827e-01_rb,1.5202e-01_rb,1.6547e-01_rb,1.7854e-01_rb,1.8991e-01_rb/)
      kao(:,13, 7) = (/ &
     &1.1534e-01_rb,1.2733e-01_rb,1.3962e-01_rb,1.4987e-01_rb,1.6008e-01_rb/)
      kao(:, 1, 8) = (/ &
     &3.5920e+00_rb,3.7800e+00_rb,3.9524e+00_rb,4.1109e+00_rb,4.2529e+00_rb/)
      kao(:, 2, 8) = (/ &
     &3.2111e+00_rb,3.3886e+00_rb,3.5521e+00_rb,3.6948e+00_rb,3.8225e+00_rb/)
      kao(:, 3, 8) = (/ &
     &2.7787e+00_rb,2.9477e+00_rb,3.0955e+00_rb,3.2271e+00_rb,3.3442e+00_rb/)
      kao(:, 4, 8) = (/ &
     &2.3588e+00_rb,2.5119e+00_rb,2.6472e+00_rb,2.7686e+00_rb,2.8775e+00_rb/)
      kao(:, 5, 8) = (/ &
     &1.9705e+00_rb,2.1082e+00_rb,2.2326e+00_rb,2.3452e+00_rb,2.4470e+00_rb/)
      kao(:, 6, 8) = (/ &
     &1.6168e+00_rb,1.7418e+00_rb,1.8565e+00_rb,1.9616e+00_rb,2.0551e+00_rb/)
      kao(:, 7, 8) = (/ &
     &1.3089e+00_rb,1.4235e+00_rb,1.5288e+00_rb,1.6233e+00_rb,1.7079e+00_rb/)
      kao(:, 8, 8) = (/ &
     &1.0505e+00_rb,1.1539e+00_rb,1.2478e+00_rb,1.3319e+00_rb,1.4095e+00_rb/)
      kao(:, 9, 8) = (/ &
     &8.3583e-01_rb,9.2738e-01_rb,1.0101e+00_rb,1.0861e+00_rb,1.1567e+00_rb/)
      kao(:,10, 8) = (/ &
     &6.6552e-01_rb,7.4516e-01_rb,8.1795e-01_rb,8.8505e-01_rb,9.4596e-01_rb/)
      kao(:,11, 8) = (/ &
     &5.5048e-01_rb,6.1499e-01_rb,6.7571e-01_rb,7.3145e-01_rb,7.8064e-01_rb/)
      kao(:,12, 8) = (/ &
     &4.5804e-01_rb,5.1260e-01_rb,5.6222e-01_rb,6.0624e-01_rb,6.4709e-01_rb/)
      kao(:,13, 8) = (/ &
     &3.8513e-01_rb,4.3049e-01_rb,4.7109e-01_rb,5.0795e-01_rb,5.3980e-01_rb/)

      end subroutine lw_kgb102

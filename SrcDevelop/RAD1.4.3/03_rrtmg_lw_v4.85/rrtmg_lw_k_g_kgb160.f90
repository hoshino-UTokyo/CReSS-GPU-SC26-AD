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
!!prev      subroutine lw_kgb16
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:,:), fracrefbo(:)
!!
      subroutine lw_kgb160
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg16, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg16, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level: P = 387.6100 mbar, T = 250.17 K
      fracrefao(:, 1) = (/ &
     &  1.1593e-01_rb,2.3390e-01_rb,1.9120e-01_rb,1.3121e-01_rb,1.0590e-01_rb,8.4852e-02_rb, &
     &  6.4168e-02_rb,4.2537e-02_rb,2.3220e-02_rb,2.1767e-03_rb,1.8203e-03_rb,1.3724e-03_rb, &
     &  9.5452e-04_rb,5.5015e-04_rb,1.9348e-04_rb,2.7344e-05_rb/)
      fracrefao(:, 2) = (/ &
     &  2.8101e-01_rb,1.9773e-01_rb,1.4749e-01_rb,1.1399e-01_rb,8.8190e-02_rb,7.0531e-02_rb, &
     &  4.6356e-02_rb,3.0774e-02_rb,1.7332e-02_rb,2.0054e-03_rb,1.5950e-03_rb,1.2760e-03_rb, &
     &  9.5034e-04_rb,5.4992e-04_rb,1.9349e-04_rb,2.7309e-05_rb/)
      fracrefao(:, 3) = (/ &
     &  2.9054e-01_rb,2.1263e-01_rb,1.4133e-01_rb,1.1083e-01_rb,8.5107e-02_rb,6.5247e-02_rb, &
     &  4.4542e-02_rb,2.7205e-02_rb,1.6495e-02_rb,1.8453e-03_rb,1.5222e-03_rb,1.1884e-03_rb, &
     &  8.1094e-04_rb,4.9173e-04_rb,1.9344e-04_rb,2.7286e-05_rb/)
      fracrefao(:, 4) = (/ &
     &  2.9641e-01_rb,2.1738e-01_rb,1.4228e-01_rb,1.0830e-01_rb,8.2837e-02_rb,6.1359e-02_rb, &
     &  4.4683e-02_rb,2.5027e-02_rb,1.6057e-02_rb,1.7558e-03_rb,1.4193e-03_rb,1.0970e-03_rb, &
     &  7.8281e-04_rb,4.3260e-04_rb,1.4837e-04_rb,2.2958e-05_rb/)
      fracrefao(:, 5) = (/ &
     &  2.9553e-01_rb,2.2139e-01_rb,1.4816e-01_rb,1.0601e-01_rb,8.0048e-02_rb,6.0082e-02_rb, &
     &  4.3952e-02_rb,2.3788e-02_rb,1.5734e-02_rb,1.6586e-03_rb,1.3434e-03_rb,1.0281e-03_rb, &
     &  7.0256e-04_rb,4.2577e-04_rb,1.2803e-04_rb,1.3315e-05_rb/)
      fracrefao(:, 6) = (/ &
     &  2.9313e-01_rb,2.2476e-01_rb,1.5470e-01_rb,1.0322e-01_rb,7.8904e-02_rb,5.8175e-02_rb, &
     &  4.3097e-02_rb,2.3618e-02_rb,1.5385e-02_rb,1.5942e-03_rb,1.2702e-03_rb,9.5566e-04_rb, &
     &  6.5421e-04_rb,4.0165e-04_rb,1.2805e-04_rb,1.3355e-05_rb/)
      fracrefao(:, 7) = (/ &
     &  2.9069e-01_rb,2.2823e-01_rb,1.5995e-01_rb,1.0170e-01_rb,7.7287e-02_rb,5.6780e-02_rb, &
     &  4.1752e-02_rb,2.3899e-02_rb,1.4937e-02_rb,1.4916e-03_rb,1.1909e-03_rb,9.1307e-04_rb, &
     &  6.3518e-04_rb,3.9866e-04_rb,1.2805e-04_rb,1.3298e-05_rb/)
      fracrefao(:, 8) = (/ &
     &  2.8446e-01_rb,2.2651e-01_rb,1.7133e-01_rb,1.0299e-01_rb,7.4231e-02_rb,5.6031e-02_rb, &
     &  4.1368e-02_rb,2.4318e-02_rb,1.4135e-02_rb,1.4216e-03_rb,1.1465e-03_rb,8.9800e-04_rb, &
     &  6.3553e-04_rb,3.9536e-04_rb,1.2749e-04_rb,1.3298e-05_rb/)
      fracrefao(:, 9) = (/ &
     &  2.0568e-01_rb,2.5049e-01_rb,2.0568e-01_rb,1.1781e-01_rb,7.5579e-02_rb,5.8136e-02_rb, &
     &  4.2397e-02_rb,2.6544e-02_rb,1.3067e-02_rb,1.4061e-03_rb,1.1455e-03_rb,8.9408e-04_rb, &
     &  6.3652e-04_rb,3.9450e-04_rb,1.2841e-04_rb,1.3315e-05_rb/)

! Planck fraction mapping level : P=95.58350 mb, T = 215.70 K
      fracrefbo(:) = (/ &
     &  1.8111e-01_rb,2.2612e-01_rb,1.6226e-01_rb,1.1872e-01_rb,9.9048e-02_rb,8.0390e-02_rb, &
     &  6.1648e-02_rb,4.1704e-02_rb,2.2976e-02_rb,1.9263e-03_rb,1.4694e-03_rb,1.1498e-03_rb, &
     &  7.9906e-04_rb,4.8310e-04_rb,1.6188e-04_rb,2.2651e-05_rb/)

      end subroutine lw_kgb160

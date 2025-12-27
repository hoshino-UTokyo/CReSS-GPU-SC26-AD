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
!!          set parameters for fracrefao(:,:), fracrefbo(:,:)
!!
      subroutine lw_kgb040
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg04, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg04, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level : P = 142.5940 mbar, T = 215.70 K
      fracrefao(:, 1) = (/ &
     &   1.5572e-01_rb,1.4925e-01_rb,1.4107e-01_rb,1.3126e-01_rb,1.1791e-01_rb,1.0173e-01_rb, &
     &   8.2949e-02_rb,6.2393e-02_rb,4.2146e-02_rb,4.5907e-03_rb,3.7965e-03_rb,2.9744e-03_rb, &
     &   2.2074e-03_rb,1.4063e-03_rb,5.3012e-04_rb,7.4595e-05_rb/)
      fracrefao(:, 2) = (/ &
     &   1.5572e-01_rb,1.4925e-01_rb,1.4107e-01_rb,1.3126e-01_rb,1.1791e-01_rb,1.0173e-01_rb, &
     &   8.2949e-02_rb,6.2392e-02_rb,4.2146e-02_rb,4.5906e-03_rb,3.7965e-03_rb,2.9745e-03_rb, &
     &   2.2074e-03_rb,1.4063e-03_rb,5.3012e-04_rb,7.4595e-05_rb/)
      fracrefao(:, 3) = (/ &
     &   1.5572e-01_rb,1.4925e-01_rb,1.4107e-01_rb,1.3126e-01_rb,1.1791e-01_rb,1.0173e-01_rb, &
     &   8.2949e-02_rb,6.2393e-02_rb,4.2146e-02_rb,4.5907e-03_rb,3.7965e-03_rb,2.9745e-03_rb, &
     &   2.2074e-03_rb,1.4063e-03_rb,5.3012e-04_rb,7.4595e-05_rb/)
      fracrefao(:, 4) = (/ &
     &   1.5572e-01_rb,1.4925e-01_rb,1.4107e-01_rb,1.3126e-01_rb,1.1791e-01_rb,1.0173e-01_rb, &
     &   8.2949e-02_rb,6.2393e-02_rb,4.2146e-02_rb,4.5907e-03_rb,3.7964e-03_rb,2.9744e-03_rb, &
     &   2.2074e-03_rb,1.4063e-03_rb,5.3012e-04_rb,7.4595e-05_rb/)
      fracrefao(:, 5) = (/ &
     &   1.5572e-01_rb,1.4925e-01_rb,1.4107e-01_rb,1.3126e-01_rb,1.1791e-01_rb,1.0173e-01_rb, &
     &   8.2949e-02_rb,6.2393e-02_rb,4.2146e-02_rb,4.5907e-03_rb,3.7965e-03_rb,2.9744e-03_rb, &
     &   2.2074e-03_rb,1.4063e-03_rb,5.3012e-04_rb,7.4595e-05_rb/)
      fracrefao(:, 6) = (/ &
     &   1.5572e-01_rb,1.4925e-01_rb,1.4107e-01_rb,1.3126e-01_rb,1.1791e-01_rb,1.0173e-01_rb, &
     &   8.2949e-02_rb,6.2393e-02_rb,4.2146e-02_rb,4.5907e-03_rb,3.7965e-03_rb,2.9744e-03_rb, &
     &   2.2074e-03_rb,1.4063e-03_rb,5.3012e-04_rb,7.4595e-05_rb/)
      fracrefao(:, 7) = (/ &
     &   1.5572e-01_rb,1.4926e-01_rb,1.4107e-01_rb,1.3126e-01_rb,1.1791e-01_rb,1.0173e-01_rb, &
     &   8.2949e-02_rb,6.2393e-02_rb,4.2146e-02_rb,4.5908e-03_rb,3.7964e-03_rb,2.9745e-03_rb, &
     &   2.2074e-03_rb,1.4063e-03_rb,5.3012e-04_rb,7.4595e-05_rb/)
      fracrefao(:, 8) = (/ &
     &   1.5571e-01_rb,1.4926e-01_rb,1.4107e-01_rb,1.3125e-01_rb,1.1791e-01_rb,1.0173e-01_rb, &
     &   8.2949e-02_rb,6.2393e-02_rb,4.2146e-02_rb,4.5907e-03_rb,3.7964e-03_rb,2.9744e-03_rb, &
     &   2.2074e-03_rb,1.4063e-03_rb,5.3012e-04_rb,7.4595e-05_rb/)
      fracrefao(:, 9) = (/ &
     &   1.5952e-01_rb,1.5155e-01_rb,1.4217e-01_rb,1.3077e-01_rb,1.1667e-01_rb,1.0048e-01_rb, &
     &   8.1511e-02_rb,6.1076e-02_rb,4.1111e-02_rb,4.4432e-03_rb,3.6910e-03_rb,2.9076e-03_rb, &
     &   2.1329e-03_rb,1.3566e-03_rb,5.2235e-04_rb,7.9935e-05_rb/)

! Planck fraction mapping level : P = 95.58350 mb, T = 215.70 K
      fracrefbo(:, 1) = (/ &
     &   1.5558e-01_rb,1.4931e-01_rb,1.4104e-01_rb,1.3124e-01_rb,1.1793e-01_rb,1.0160e-01_rb, &
     &   8.3142e-02_rb,6.2403e-02_rb,4.2170e-02_rb,4.5935e-03_rb,3.7976e-03_rb,2.9986e-03_rb, &
     &   2.1890e-03_rb,1.4061e-03_rb,5.3005e-04_rb,7.4587e-05_rb/)
      fracrefbo(:, 2) = (/ &
     &   1.5558e-01_rb,1.4932e-01_rb,1.4104e-01_rb,1.3124e-01_rb,1.1792e-01_rb,1.0159e-01_rb, &
     &   8.3142e-02_rb,6.2403e-02_rb,4.2170e-02_rb,4.5935e-03_rb,3.7976e-03_rb,2.9986e-03_rb, &
     &   2.1890e-03_rb,1.4061e-03_rb,5.3005e-04_rb,7.4587e-05_rb/)
      fracrefbo(:, 3) = (/ &
     &   1.5558e-01_rb,1.4933e-01_rb,1.4103e-01_rb,1.3124e-01_rb,1.1792e-01_rb,1.0159e-01_rb, &
     &   8.3142e-02_rb,6.2403e-02_rb,4.2170e-02_rb,4.5935e-03_rb,3.7976e-03_rb,2.9986e-03_rb, &
     &   2.1890e-03_rb,1.4061e-03_rb,5.3005e-04_rb,7.4587e-05_rb/)
      fracrefbo(:, 4) = (/ &
     &   1.5569e-01_rb,1.4926e-01_rb,1.4102e-01_rb,1.3122e-01_rb,1.1791e-01_rb,1.0159e-01_rb, &
     &   8.3141e-02_rb,6.2403e-02_rb,4.2170e-02_rb,4.5935e-03_rb,3.7976e-03_rb,2.9986e-03_rb, &
     &   2.1890e-03_rb,1.4061e-03_rb,5.3005e-04_rb,7.4587e-05_rb/)
      fracrefbo(:, 5) = (/ &
     &   1.5947e-01_rb,1.5132e-01_rb,1.4195e-01_rb,1.3061e-01_rb,1.1680e-01_rb,1.0054e-01_rb, &
     &   8.1785e-02_rb,6.1212e-02_rb,4.1276e-02_rb,4.4424e-03_rb,3.6628e-03_rb,2.8943e-03_rb, &
     &   2.1134e-03_rb,1.3457e-03_rb,5.1024e-04_rb,7.3998e-05_rb/)

      end subroutine lw_kgb040

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
!!prev      subroutine lw_kgb12
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:,:)
!!
      subroutine lw_kgb120
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg12, only : fracrefao, kao, selfrefo, forrefo
      use rrlw_kg12, only : fracrefao

      implicit none
      save

! Planck fraction mapping level : P = 174.1640 mbar, T= 215.78 K
      fracrefao(:, 1) = (/ &
     &  1.3984e-01_rb,1.6809e-01_rb,1.8072e-01_rb,1.5400e-01_rb,1.2613e-01_rb,9.6959e-02_rb, &
     &  5.9713e-02_rb,3.8631e-02_rb,2.6937e-02_rb,3.1711e-03_rb,2.3458e-03_rb,1.4653e-03_rb, &
     &  1.0567e-03_rb,6.6504e-04_rb,2.4957e-04_rb,3.5172e-05_rb/)
      fracrefao(:, 2) = (/ &
     &  1.2745e-01_rb,1.6107e-01_rb,1.6568e-01_rb,1.5436e-01_rb,1.3183e-01_rb,1.0166e-01_rb, &
     &  6.4506e-02_rb,4.7756e-02_rb,3.4472e-02_rb,3.7189e-03_rb,2.9349e-03_rb,2.1469e-03_rb, &
     &  1.3746e-03_rb,7.1691e-04_rb,2.8057e-04_rb,5.6242e-05_rb/)
      fracrefao(:, 3) = (/ &
     &  1.2181e-01_rb,1.5404e-01_rb,1.6540e-01_rb,1.5255e-01_rb,1.3736e-01_rb,9.8856e-02_rb, &
     &  6.8927e-02_rb,5.1385e-02_rb,3.7046e-02_rb,4.0302e-03_rb,3.0949e-03_rb,2.3772e-03_rb, &
     &  1.6538e-03_rb,8.9641e-04_rb,4.6991e-04_rb,1.1251e-04_rb/)
      fracrefao(:, 4) = (/ &
     &  1.1794e-01_rb,1.4864e-01_rb,1.6316e-01_rb,1.5341e-01_rb,1.3986e-01_rb,9.6656e-02_rb, &
     &  7.2478e-02_rb,5.5061e-02_rb,3.8886e-02_rb,4.3398e-03_rb,3.3576e-03_rb,2.4891e-03_rb, &
     &  1.7674e-03_rb,1.0764e-03_rb,7.7689e-04_rb,1.1251e-04_rb/)
      fracrefao(:, 5) = (/ &
     &  1.1635e-01_rb,1.4342e-01_rb,1.5924e-01_rb,1.5670e-01_rb,1.3740e-01_rb,9.7087e-02_rb, &
     &  7.6250e-02_rb,5.7802e-02_rb,4.0808e-02_rb,4.4113e-03_rb,3.6035e-03_rb,2.6269e-03_rb, &
     &  1.7586e-03_rb,1.6498e-03_rb,7.7689e-04_rb,1.1251e-04_rb/)
      fracrefao(:, 6) = (/ &
     &  1.1497e-01_rb,1.3751e-01_rb,1.5587e-01_rb,1.5904e-01_rb,1.3140e-01_rb,1.0159e-01_rb, &
     &  7.9729e-02_rb,6.1475e-02_rb,4.2382e-02_rb,4.5291e-03_rb,3.8161e-03_rb,2.7683e-03_rb, &
     &  1.9899e-03_rb,2.0395e-03_rb,7.7720e-04_rb,1.1251e-04_rb/)
      fracrefao(:, 7) = (/ &
     &  1.1331e-01_rb,1.3015e-01_rb,1.5574e-01_rb,1.5489e-01_rb,1.2697e-01_rb,1.0746e-01_rb, &
     &  8.4777e-02_rb,6.5145e-02_rb,4.4293e-02_rb,4.7426e-03_rb,3.8383e-03_rb,2.9065e-03_rb, &
     &  2.8430e-03_rb,2.0401e-03_rb,7.7689e-04_rb,1.1251e-04_rb/)
      fracrefao(:, 8) = (/ &
     &  1.0993e-01_rb,1.2320e-01_rb,1.4893e-01_rb,1.4573e-01_rb,1.3174e-01_rb,1.1149e-01_rb, &
     &  9.3326e-02_rb,6.9942e-02_rb,4.6762e-02_rb,4.9309e-03_rb,3.8583e-03_rb,4.1889e-03_rb, &
     &  3.0415e-03_rb,2.0406e-03_rb,7.7720e-04_rb,1.1251e-04_rb/)
      fracrefao(:, 9) = (/ &
     &  1.2028e-01_rb,1.2091e-01_rb,1.3098e-01_rb,1.3442e-01_rb,1.3574e-01_rb,1.1739e-01_rb, &
     &  9.5343e-02_rb,7.0224e-02_rb,5.3456e-02_rb,6.0206e-03_rb,5.0758e-03_rb,4.1906e-03_rb, &
     &  3.0431e-03_rb,2.0400e-03_rb,7.7689e-04_rb,1.1251e-04_rb/)

      end subroutine lw_kgb120

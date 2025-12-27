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
!!prev      subroutine lw_kgb03
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:,:), fracrefbo(:,:)
!!
      subroutine lw_kgb030
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg03, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
!!prev                            kbo_mn2o, selfrefo, forrefo
      use rrlw_kg03, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level: P=212.7250 mbar, T = 223.06 K
      fracrefao(:, 1) = (/ &
     &   1.6251e-01_rb,1.5572e-01_rb,1.4557e-01_rb,1.3208e-01_rb,1.1582e-01_rb,9.6895e-02_rb, &
     &   7.8720e-02_rb,5.8462e-02_rb,3.9631e-02_rb,4.3001e-03_rb,3.5555e-03_rb,2.8101e-03_rb, &
     &   2.0547e-03_rb,1.3109e-03_rb,4.9403e-04_rb,6.9515e-05_rb/)
      fracrefao(:, 2) = (/ &
     &   1.6006e-01_rb,1.5576e-01_rb,1.4609e-01_rb,1.3276e-01_rb,1.1594e-01_rb,9.7336e-02_rb, &
     &   7.9035e-02_rb,5.8696e-02_rb,3.9723e-02_rb,4.3001e-03_rb,3.5555e-03_rb,2.8101e-03_rb, &
     &   2.0547e-03_rb,1.3109e-03_rb,4.9403e-04_rb,6.9515e-05_rb/)
      fracrefao(:, 3) = (/ &
     &   1.5952e-01_rb,1.5566e-01_rb,1.4590e-01_rb,1.3294e-01_rb,1.1599e-01_rb,9.7511e-02_rb, &
     &   7.9127e-02_rb,5.8888e-02_rb,3.9874e-02_rb,4.3001e-03_rb,3.5555e-03_rb,2.8102e-03_rb, &
     &   2.0547e-03_rb,1.3109e-03_rb,4.9403e-04_rb,6.9515e-05_rb/)
      fracrefao(:, 4) = (/ &
     &   1.5907e-01_rb,1.5541e-01_rb,1.4585e-01_rb,1.3316e-01_rb,1.1596e-01_rb,9.7647e-02_rb, &
     &   7.9243e-02_rb,5.9024e-02_rb,4.0028e-02_rb,4.3112e-03_rb,3.5555e-03_rb,2.8102e-03_rb, &
     &   2.0547e-03_rb,1.3109e-03_rb,4.9403e-04_rb,6.9515e-05_rb/)
      fracrefao(:, 5) = (/ &
     &   1.5862e-01_rb,1.5517e-01_rb,1.4588e-01_rb,1.3328e-01_rb,1.1585e-01_rb,9.7840e-02_rb, &
     &   7.9364e-02_rb,5.9174e-02_rb,4.0160e-02_rb,4.3403e-03_rb,3.5900e-03_rb,2.8102e-03_rb, &
     &   2.0547e-03_rb,1.3109e-03_rb,4.9403e-04_rb,6.9515e-05_rb/)
      fracrefao(:, 6) = (/ &
     &   1.5830e-01_rb,1.5490e-01_rb,1.4582e-01_rb,1.3331e-01_rb,1.1567e-01_rb,9.8079e-02_rb, &
     &   7.9510e-02_rb,5.9369e-02_rb,4.0326e-02_rb,4.3343e-03_rb,3.5908e-03_rb,2.8527e-03_rb, &
     &   2.0655e-03_rb,1.3109e-03_rb,4.9403e-04_rb,6.9515e-05_rb/)
      fracrefao(:, 7) = (/ &
     &   1.5789e-01_rb,1.5435e-01_rb,1.4595e-01_rb,1.3304e-01_rb,1.1566e-01_rb,9.8426e-02_rb, &
     &   7.9704e-02_rb,5.9618e-02_rb,4.0520e-02_rb,4.3812e-03_rb,3.6147e-03_rb,2.8395e-03_rb, &
     &   2.1301e-03_rb,1.3145e-03_rb,4.9403e-04_rb,6.9515e-05_rb/)
      fracrefao(:, 8) = (/ &
     &   1.5704e-01_rb,1.5398e-01_rb,1.4564e-01_rb,1.3222e-01_rb,1.1586e-01_rb,9.9230e-02_rb, &
     &   8.0011e-02_rb,6.0149e-02_rb,4.0790e-02_rb,4.4253e-03_rb,3.6534e-03_rb,2.9191e-03_rb, &
     &   2.1373e-03_rb,1.3558e-03_rb,5.1631e-04_rb,7.8794e-05_rb/)
      fracrefao(:, 9) = (/ &
     &   1.5270e-01_rb,1.5126e-01_rb,1.4264e-01_rb,1.3106e-01_rb,1.1740e-01_rb,1.0137e-01_rb, &
     &   8.3057e-02_rb,6.2282e-02_rb,4.2301e-02_rb,4.6486e-03_rb,3.8159e-03_rb,3.0472e-03_rb, &
     &   2.2870e-03_rb,1.4818e-03_rb,5.6773e-04_rb,7.8794e-05_rb/)

! Planck fraction mapping level: p = 95.8 mbar, t = 215.7 k
      fracrefbo(:, 1) = (/ &
     &   1.6413e-01_rb,1.5665e-01_rb,1.4606e-01_rb,1.3184e-01_rb,1.1517e-01_rb,9.6243e-02_rb, &
     &   7.7982e-02_rb,5.8165e-02_rb,3.9311e-02_rb,4.2586e-03_rb,3.5189e-03_rb,2.7793e-03_rb, &
     &   2.0376e-03_rb,1.2938e-03_rb,4.8853e-04_rb,6.8745e-05_rb/)
      fracrefbo(:, 2) = (/ &
     &   1.6254e-01_rb,1.5674e-01_rb,1.4652e-01_rb,1.3221e-01_rb,1.1535e-01_rb,9.6439e-02_rb, &
     &   7.8155e-02_rb,5.8254e-02_rb,3.9343e-02_rb,4.2586e-03_rb,3.5189e-03_rb,2.7793e-03_rb, &
     &   2.0376e-03_rb,1.2938e-03_rb,4.8853e-04_rb,6.8745e-05_rb/)
      fracrefbo(:, 3) = (/ &
     &   1.6177e-01_rb,1.5664e-01_rb,1.4669e-01_rb,1.3242e-01_rb,1.1541e-01_rb,9.6536e-02_rb, &
     &   7.8257e-02_rb,5.8387e-02_rb,3.9431e-02_rb,4.2587e-03_rb,3.5189e-03_rb,2.7793e-03_rb, &
     &   2.0376e-03_rb,1.2938e-03_rb,4.8853e-04_rb,6.8745e-05_rb/)
      fracrefbo(:, 4) = (/ &
     &   1.6077e-01_rb,1.5679e-01_rb,1.4648e-01_rb,1.3273e-01_rb,1.1546e-01_rb,9.6779e-02_rb, &
     &   7.8371e-02_rb,5.8546e-02_rb,3.9611e-02_rb,4.2772e-03_rb,3.5190e-03_rb,2.7793e-03_rb, &
     &   2.0376e-03_rb,1.2938e-03_rb,4.8853e-04_rb,6.8745e-05_rb/)
      fracrefbo(:, 5) = (/ &
     &   1.6067e-01_rb,1.5608e-01_rb,1.4247e-01_rb,1.2881e-01_rb,1.1449e-01_rb,9.8802e-02_rb, &
     &   8.0828e-02_rb,6.0977e-02_rb,4.1494e-02_rb,4.5116e-03_rb,3.7290e-03_rb,2.9460e-03_rb, &
     &   2.1948e-03_rb,1.3778e-03_rb,5.4552e-04_rb,7.9969e-05_rb/)

      end subroutine lw_kgb030

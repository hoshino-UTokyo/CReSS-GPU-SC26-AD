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
!!          set parameters for fracrefao(:,:), fracrefbo(:,:)
!!
      subroutine lw_kgb090
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg09, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
!!prev                            kbo_mn2o, selfrefo, forrefo
      use rrlw_kg09, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fractions mapping level : P=212.7250 mb, T = 223.06 K
      fracrefao(:, 1) = (/ &
     &  1.8129e-01_rb,1.6119e-01_rb,1.3308e-01_rb,1.2342e-01_rb,1.1259e-01_rb,9.7580e-02_rb, &
     &  7.9176e-02_rb,5.8541e-02_rb,3.9084e-02_rb,4.2419e-03_rb,3.4314e-03_rb,2.6935e-03_rb, &
     &  1.9404e-03_rb,1.2218e-03_rb,4.5263e-04_rb,6.0909e-05_rb/)
      fracrefao(:, 2) = (/ &
     &  1.9665e-01_rb,1.5640e-01_rb,1.3101e-01_rb,1.2153e-01_rb,1.1037e-01_rb,9.6043e-02_rb, &
     &  7.7856e-02_rb,5.7547e-02_rb,3.8670e-02_rb,4.1955e-03_rb,3.4104e-03_rb,2.6781e-03_rb, &
     &  1.9245e-03_rb,1.2093e-03_rb,4.4113e-04_rb,6.0913e-05_rb/)
      fracrefao(:, 3) = (/ &
     &  2.0273e-01_rb,1.5506e-01_rb,1.3044e-01_rb,1.2043e-01_rb,1.0952e-01_rb,9.5384e-02_rb, &
     &  7.7157e-02_rb,5.7176e-02_rb,3.8379e-02_rb,4.1584e-03_rb,3.3836e-03_rb,2.6412e-03_rb, &
     &  1.8865e-03_rb,1.1791e-03_rb,4.2094e-04_rb,4.7410e-05_rb/)
      fracrefao(:, 4) = (/ &
     &  2.0272e-01_rb,1.5963e-01_rb,1.2913e-01_rb,1.2060e-01_rb,1.0820e-01_rb,9.4685e-02_rb, &
     &  7.6544e-02_rb,5.6851e-02_rb,3.8155e-02_rb,4.0913e-03_rb,3.3442e-03_rb,2.6054e-03_rb, &
     &  1.8875e-03_rb,1.1263e-03_rb,3.7743e-04_rb,4.7410e-05_rb/)
      fracrefao(:, 5) = (/ &
     &  2.0280e-01_rb,1.6353e-01_rb,1.2910e-01_rb,1.1968e-01_rb,1.0725e-01_rb,9.4112e-02_rb, &
     &  7.5828e-02_rb,5.6526e-02_rb,3.7972e-02_rb,4.0205e-03_rb,3.3063e-03_rb,2.5681e-03_rb, &
     &  1.8386e-03_rb,1.0757e-03_rb,3.5301e-04_rb,4.7410e-05_rb/)
      fracrefao(:, 6) = (/ &
     &  2.0294e-01_rb,1.6840e-01_rb,1.2852e-01_rb,1.1813e-01_rb,1.0724e-01_rb,9.2946e-02_rb, &
     &  7.5029e-02_rb,5.6158e-02_rb,3.7744e-02_rb,3.9632e-03_rb,3.2434e-03_rb,2.5275e-03_rb, &
     &  1.7558e-03_rb,1.0080e-03_rb,3.5301e-04_rb,4.7410e-05_rb/)
      fracrefao(:, 7) = (/ &
     &  2.0313e-01_rb,1.7390e-01_rb,1.2864e-01_rb,1.1689e-01_rb,1.0601e-01_rb,9.1791e-02_rb, &
     &  7.4224e-02_rb,5.5500e-02_rb,3.7374e-02_rb,3.9214e-03_rb,3.1984e-03_rb,2.4162e-03_rb, &
     &  1.6394e-03_rb,9.7275e-04_rb,3.5299e-04_rb,4.7410e-05_rb/)
      fracrefao(:, 8) = (/ &
     &  2.0332e-01_rb,1.7800e-01_rb,1.3286e-01_rb,1.1555e-01_rb,1.0407e-01_rb,9.0475e-02_rb, &
     &  7.2452e-02_rb,5.4566e-02_rb,3.6677e-02_rb,3.7889e-03_rb,3.0351e-03_rb,2.2587e-03_rb, &
     &  1.5764e-03_rb,9.7270e-04_rb,3.5300e-04_rb,4.7410e-05_rb/)
      fracrefao(:, 9) = (/ &
     &  1.9624e-01_rb,1.6519e-01_rb,1.3663e-01_rb,1.1535e-01_rb,1.0719e-01_rb,9.4156e-02_rb, &
     &  7.6745e-02_rb,5.6987e-02_rb,3.8135e-02_rb,4.1626e-03_rb,3.4243e-03_rb,2.7116e-03_rb, &
     &  1.7095e-03_rb,9.7271e-04_rb,3.5299e-04_rb,4.7410e-05_rb/)

! Planck fraction mapping level : p=3.20e-2 mb, t = 197.92 k
      fracrefbo(:) = (/ &
     &  2.0914e-01_rb,1.5077e-01_rb,1.2878e-01_rb,1.1856e-01_rb,1.0695e-01_rb,9.3048e-02_rb, &
     &  7.7645e-02_rb,6.0785e-02_rb,4.0642e-02_rb,4.0499e-03_rb,3.3931e-03_rb,2.6363e-03_rb, &
     &  1.9151e-03_rb,1.1963e-03_rb,4.3471e-04_rb,5.1421e-05_rb/)

      end subroutine lw_kgb090

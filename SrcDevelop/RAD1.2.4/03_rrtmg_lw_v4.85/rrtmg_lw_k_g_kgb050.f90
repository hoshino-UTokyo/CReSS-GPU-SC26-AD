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
!!prev      subroutine lw_kgb05
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:,:), fracrefbo(:,:), ccl4o(:)
!!
      subroutine lw_kgb050
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg05, only : fracrefao, fracrefbo, kao, kbo, kao_mo3, &
!!prev                            selfrefo, forrefo, ccl4o
      use rrlw_kg05, only : fracrefao, fracrefbo, &
                            ccl4o

      implicit none
      save

! Planck fraction mapping level : P = 473.42 mb, T = 259.83
      fracrefao(:, 1) = (/ &
        1.4111e-01_rb,1.4222e-01_rb,1.3802e-01_rb,1.3101e-01_rb,1.2244e-01_rb,1.0691e-01_rb, &
        8.8703e-02_rb,6.7130e-02_rb,4.5509e-02_rb,4.9866e-03_rb,4.1214e-03_rb,3.2557e-03_rb, &
        2.3805e-03_rb,1.5450e-03_rb,5.8423e-04_rb,8.2275e-05_rb/)
      fracrefao(:, 2) = (/ &
        1.4152e-01_rb,1.4271e-01_rb,1.3784e-01_rb,1.3075e-01_rb,1.2215e-01_rb,1.0674e-01_rb, &
        8.8686e-02_rb,6.7135e-02_rb,4.5508e-02_rb,4.9866e-03_rb,4.1214e-03_rb,3.2558e-03_rb, &
        2.3805e-03_rb,1.5450e-03_rb,5.8423e-04_rb,8.2275e-05_rb/)
      fracrefao(:, 3) = (/ &
        1.4159e-01_rb,1.4300e-01_rb,1.3781e-01_rb,1.3094e-01_rb,1.2192e-01_rb,1.0661e-01_rb, &
        8.8529e-02_rb,6.7127e-02_rb,4.5511e-02_rb,4.9877e-03_rb,4.1214e-03_rb,3.2558e-03_rb, &
        2.3805e-03_rb,1.5450e-03_rb,5.8423e-04_rb,8.2275e-05_rb/)
      fracrefao(:, 4) = (/ &
        1.4162e-01_rb,1.4337e-01_rb,1.3774e-01_rb,1.3122e-01_rb,1.2172e-01_rb,1.0641e-01_rb, &
        8.8384e-02_rb,6.7056e-02_rb,4.5514e-02_rb,4.9880e-03_rb,4.1214e-03_rb,3.2557e-03_rb, &
        2.3805e-03_rb,1.5450e-03_rb,5.8423e-04_rb,8.2275e-05_rb/)
      fracrefao(:, 5) = (/ &
        1.4161e-01_rb,1.4370e-01_rb,1.3770e-01_rb,1.3143e-01_rb,1.2173e-01_rb,1.0613e-01_rb, &
        8.8357e-02_rb,6.6874e-02_rb,4.5509e-02_rb,4.9883e-03_rb,4.1214e-03_rb,3.2558e-03_rb, &
        2.3804e-03_rb,1.5450e-03_rb,5.8423e-04_rb,8.2275e-05_rb/)
      fracrefao(:, 6) = (/ &
        1.4154e-01_rb,1.4405e-01_rb,1.3771e-01_rb,1.3169e-01_rb,1.2166e-01_rb,1.0603e-01_rb, &
        8.8193e-02_rb,6.6705e-02_rb,4.5469e-02_rb,4.9902e-03_rb,4.1214e-03_rb,3.2558e-03_rb, &
        2.3804e-03_rb,1.5450e-03_rb,5.8423e-04_rb,8.2275e-05_rb/)
      fracrefao(:, 7) = (/ &
        1.4126e-01_rb,1.4440e-01_rb,1.3790e-01_rb,1.3214e-01_rb,1.2153e-01_rb,1.0603e-01_rb, &
        8.7908e-02_rb,6.6612e-02_rb,4.5269e-02_rb,4.9900e-03_rb,4.1256e-03_rb,3.2558e-03_rb, &
        2.3804e-03_rb,1.5451e-03_rb,5.8423e-04_rb,8.2275e-05_rb/)
      fracrefao(:, 8) = (/ &
        1.4076e-01_rb,1.4415e-01_rb,1.3885e-01_rb,1.3286e-01_rb,1.2147e-01_rb,1.0612e-01_rb, &
        8.7579e-02_rb,6.6280e-02_rb,4.4977e-02_rb,4.9782e-03_rb,4.1200e-03_rb,3.2620e-03_rb, &
        2.3820e-03_rb,1.5452e-03_rb,5.8423e-04_rb,8.2275e-05_rb/)
      fracrefao(:, 9) = (/ &
        1.4205e-01_rb,1.4496e-01_rb,1.4337e-01_rb,1.3504e-01_rb,1.2260e-01_rb,1.0428e-01_rb, &
        8.4946e-02_rb,6.3625e-02_rb,4.2951e-02_rb,4.7313e-03_rb,3.9157e-03_rb,3.0879e-03_rb, &
        2.2666e-03_rb,1.5193e-03_rb,5.7469e-04_rb,8.1674e-05_rb/)

! Planck fraction mapping level : P = 0.2369280 mbar, T = 253.60 K
      fracrefbo(:, 1) = (/ &
        1.4075e-01_rb,1.4196e-01_rb,1.3833e-01_rb,1.3345e-01_rb,1.2234e-01_rb,1.0718e-01_rb, &
        8.8004e-02_rb,6.6308e-02_rb,4.5028e-02_rb,4.9029e-03_rb,4.0377e-03_rb,3.1870e-03_rb, &
        2.3503e-03_rb,1.5146e-03_rb,5.7165e-04_rb,8.2371e-05_rb/)
      fracrefbo(:, 2) = (/ &
        1.4081e-01_rb,1.4225e-01_rb,1.3890e-01_rb,1.3410e-01_rb,1.2254e-01_rb,1.0680e-01_rb, &
        8.7391e-02_rb,6.5819e-02_rb,4.4725e-02_rb,4.9121e-03_rb,4.0420e-03_rb,3.1869e-03_rb, &
        2.3504e-03_rb,1.5146e-03_rb,5.7165e-04_rb,8.2371e-05_rb/)
      fracrefbo(:, 3) = (/ &
        1.4087e-01_rb,1.4227e-01_rb,1.3920e-01_rb,1.3395e-01_rb,1.2270e-01_rb,1.0694e-01_rb, &
        8.7229e-02_rb,6.5653e-02_rb,4.4554e-02_rb,4.8797e-03_rb,4.0460e-03_rb,3.1939e-03_rb, &
        2.3505e-03_rb,1.5146e-03_rb,5.7165e-04_rb,8.1910e-05_rb/)
      fracrefbo(:, 4) = (/ &
        1.4089e-01_rb,1.4238e-01_rb,1.3956e-01_rb,1.3379e-01_rb,1.2284e-01_rb,1.0688e-01_rb, &
        8.7192e-02_rb,6.5490e-02_rb,4.4390e-02_rb,4.8395e-03_rb,4.0173e-03_rb,3.2070e-03_rb, &
        2.3559e-03_rb,1.5146e-03_rb,5.7165e-04_rb,8.2371e-05_rb/)
      fracrefbo(:, 5) = (/ &
        1.4091e-01_rb,1.4417e-01_rb,1.4194e-01_rb,1.3457e-01_rb,1.2167e-01_rb,1.0551e-01_rb, &
        8.6450e-02_rb,6.4889e-02_rb,4.3584e-02_rb,4.7551e-03_rb,3.9509e-03_rb,3.1374e-03_rb, &
        2.3226e-03_rb,1.4942e-03_rb,5.7545e-04_rb,8.0887e-05_rb/)

! Minor gas mapping level :
!     lower - o3, p = 317.34 mbar, t = 240.77 k
!     lower - ccl4

      ccl4o(:) = (/ &
        26.1407_rb, 53.9776_rb, 63.8085_rb, 36.1701_rb, 15.4099_rb, 10.23116_rb, &
        4.82948_rb, 5.03836_rb, 1.75558_rb, 0._rb,      0._rb,      0._rb,       &
        0._rb,      0._rb,      0._rb,      0._rb/)

      end subroutine lw_kgb050

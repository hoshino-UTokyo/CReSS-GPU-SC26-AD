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
!!          set parameters for kao(:,:,1:4)
!!
      subroutine lw_kgb101
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

      kao(:, 1, 1) = (/ &
     &2.7213e-02_rb,2.9029e-02_rb,3.0838e-02_rb,3.2767e-02_rb,3.4630e-02_rb/)
      kao(:, 2, 1) = (/ &
     &2.1383e-02_rb,2.2832e-02_rb,2.4365e-02_rb,2.5925e-02_rb,2.7376e-02_rb/)
      kao(:, 3, 1) = (/ &
     &1.6478e-02_rb,1.7675e-02_rb,1.8942e-02_rb,2.0177e-02_rb,2.1374e-02_rb/)
      kao(:, 4, 1) = (/ &
     &1.2848e-02_rb,1.3809e-02_rb,1.4825e-02_rb,1.5852e-02_rb,1.6799e-02_rb/)
      kao(:, 5, 1) = (/ &
     &1.0029e-02_rb,1.0847e-02_rb,1.1686e-02_rb,1.2512e-02_rb,1.3297e-02_rb/)
      kao(:, 6, 1) = (/ &
     &7.8313e-03_rb,8.5460e-03_rb,9.2025e-03_rb,9.9079e-03_rb,1.0554e-02_rb/)
      kao(:, 7, 1) = (/ &
     &6.1234e-03_rb,6.6846e-03_rb,7.2818e-03_rb,7.8451e-03_rb,8.4144e-03_rb/)
      kao(:, 8, 1) = (/ &
     &4.8747e-03_rb,5.2881e-03_rb,5.7474e-03_rb,6.2355e-03_rb,6.7085e-03_rb/)
      kao(:, 9, 1) = (/ &
     &4.1059e-03_rb,4.5848e-03_rb,4.9152e-03_rb,5.2401e-03_rb,5.5908e-03_rb/)
      kao(:,10, 1) = (/ &
     &3.5412e-03_rb,4.0565e-03_rb,4.5689e-03_rb,5.1092e-03_rb,5.6716e-03_rb/)
      kao(:,11, 1) = (/ &
     &3.0492e-03_rb,3.6436e-03_rb,4.0799e-03_rb,4.5260e-03_rb,4.9802e-03_rb/)
      kao(:,12, 1) = (/ &
     &2.5821e-03_rb,3.0995e-03_rb,3.5069e-03_rb,3.8699e-03_rb,4.2575e-03_rb/)
      kao(:,13, 1) = (/ &
     &2.1558e-03_rb,2.5890e-03_rb,2.9127e-03_rb,3.2139e-03_rb,3.5455e-03_rb/)
      kao(:, 1, 2) = (/ &
     &5.2889e-02_rb,5.6315e-02_rb,5.9927e-02_rb,6.3408e-02_rb,6.6329e-02_rb/)
      kao(:, 2, 2) = (/ &
     &4.1932e-02_rb,4.4934e-02_rb,4.8030e-02_rb,5.0855e-02_rb,5.3372e-02_rb/)
      kao(:, 3, 2) = (/ &
     &3.2754e-02_rb,3.5198e-02_rb,3.7786e-02_rb,4.0294e-02_rb,4.2438e-02_rb/)
      kao(:, 4, 2) = (/ &
     &2.5838e-02_rb,2.7680e-02_rb,2.9873e-02_rb,3.1963e-02_rb,3.3931e-02_rb/)
      kao(:, 5, 2) = (/ &
     &2.0528e-02_rb,2.2079e-02_rb,2.3740e-02_rb,2.5501e-02_rb,2.7219e-02_rb/)
      kao(:, 6, 2) = (/ &
     &1.6350e-02_rb,1.7491e-02_rb,1.8902e-02_rb,2.0280e-02_rb,2.1774e-02_rb/)
      kao(:, 7, 2) = (/ &
     &1.2732e-02_rb,1.3953e-02_rb,1.5005e-02_rb,1.6145e-02_rb,1.7323e-02_rb/)
      kao(:, 8, 2) = (/ &
     &9.7464e-03_rb,1.1159e-02_rb,1.2150e-02_rb,1.2940e-02_rb,1.3856e-02_rb/)
      kao(:, 9, 2) = (/ &
     &7.5337e-03_rb,8.5370e-03_rb,9.6841e-03_rb,1.0825e-02_rb,1.1894e-02_rb/)
      kao(:,10, 2) = (/ &
     &7.1476e-03_rb,7.8468e-03_rb,8.6104e-03_rb,9.4234e-03_rb,1.0165e-02_rb/)
      kao(:,11, 2) = (/ &
     &7.4726e-03_rb,7.8619e-03_rb,8.4110e-03_rb,8.8097e-03_rb,9.3621e-03_rb/)
      kao(:,12, 2) = (/ &
     &6.8212e-03_rb,7.3104e-03_rb,7.6521e-03_rb,8.0570e-03_rb,8.3782e-03_rb/)
      kao(:,13, 2) = (/ &
     &5.7626e-03_rb,6.1115e-03_rb,6.3578e-03_rb,6.6833e-03_rb,6.9834e-03_rb/)
      kao(:, 1, 3) = (/ &
     &9.2909e-02_rb,9.6713e-02_rb,9.9436e-02_rb,1.0242e-01_rb,1.0613e-01_rb/)
      kao(:, 2, 3) = (/ &
     &7.4548e-02_rb,7.7785e-02_rb,8.0363e-02_rb,8.2840e-02_rb,8.6028e-02_rb/)
      kao(:, 3, 3) = (/ &
     &5.8714e-02_rb,6.1823e-02_rb,6.4284e-02_rb,6.6702e-02_rb,6.9231e-02_rb/)
      kao(:, 4, 3) = (/ &
     &4.6349e-02_rb,4.9440e-02_rb,5.1743e-02_rb,5.3890e-02_rb,5.5994e-02_rb/)
      kao(:, 5, 3) = (/ &
     &3.6507e-02_rb,3.9309e-02_rb,4.1637e-02_rb,4.3586e-02_rb,4.5365e-02_rb/)
      kao(:, 6, 3) = (/ &
     &2.8427e-02_rb,3.1177e-02_rb,3.3328e-02_rb,3.5180e-02_rb,3.6733e-02_rb/)
      kao(:, 7, 3) = (/ &
     &2.2397e-02_rb,2.4640e-02_rb,2.6638e-02_rb,2.8311e-02_rb,2.9747e-02_rb/)
      kao(:, 8, 3) = (/ &
     &1.7861e-02_rb,1.9252e-02_rb,2.1033e-02_rb,2.2648e-02_rb,2.3985e-02_rb/)
      kao(:, 9, 3) = (/ &
     &1.4398e-02_rb,1.5884e-02_rb,1.6962e-02_rb,1.7935e-02_rb,1.8840e-02_rb/)
      kao(:,10, 3) = (/ &
     &1.3336e-02_rb,1.5085e-02_rb,1.6848e-02_rb,1.8578e-02_rb,1.9024e-02_rb/)
      kao(:,11, 3) = (/ &
     &1.0996e-02_rb,1.2680e-02_rb,1.4233e-02_rb,1.6184e-02_rb,1.7835e-02_rb/)
      kao(:,12, 3) = (/ &
     &9.6066e-03_rb,1.0608e-02_rb,1.2101e-02_rb,1.3713e-02_rb,1.5338e-02_rb/)
      kao(:,13, 3) = (/ &
     &8.0007e-03_rb,8.9344e-03_rb,1.0260e-02_rb,1.1651e-02_rb,1.2914e-02_rb/)
      kao(:, 1, 4) = (/ &
     &1.4098e-01_rb,1.4735e-01_rb,1.5390e-01_rb,1.6007e-01_rb,1.6623e-01_rb/)
      kao(:, 2, 4) = (/ &
     &1.1373e-01_rb,1.1930e-01_rb,1.2491e-01_rb,1.3057e-01_rb,1.3593e-01_rb/)
      kao(:, 3, 4) = (/ &
     &9.0856e-02_rb,9.5503e-02_rb,1.0012e-01_rb,1.0479e-01_rb,1.0971e-01_rb/)
      kao(:, 4, 4) = (/ &
     &7.2695e-02_rb,7.6746e-02_rb,8.0662e-02_rb,8.4622e-02_rb,8.8881e-02_rb/)
      kao(:, 5, 4) = (/ &
     &5.8402e-02_rb,6.1962e-02_rb,6.5265e-02_rb,6.8694e-02_rb,7.2371e-02_rb/)
      kao(:, 6, 4) = (/ &
     &4.7100e-02_rb,5.0100e-02_rb,5.2955e-02_rb,5.5897e-02_rb,5.9012e-02_rb/)
      kao(:, 7, 4) = (/ &
     &3.7714e-02_rb,4.0303e-02_rb,4.2858e-02_rb,4.5396e-02_rb,4.8032e-02_rb/)
      kao(:, 8, 4) = (/ &
     &2.9938e-02_rb,3.2418e-02_rb,3.4599e-02_rb,3.6801e-02_rb,3.9055e-02_rb/)
      kao(:, 9, 4) = (/ &
     &2.3570e-02_rb,2.5475e-02_rb,2.7580e-02_rb,2.9598e-02_rb,3.1549e-02_rb/)
      kao(:,10, 4) = (/ &
     &2.3291e-02_rb,2.1658e-02_rb,2.1418e-02_rb,2.1490e-02_rb,2.3536e-02_rb/)
      kao(:,11, 4) = (/ &
     &2.1808e-02_rb,2.2150e-02_rb,2.2833e-02_rb,2.0354e-02_rb,1.9688e-02_rb/)
      kao(:,12, 4) = (/ &
     &1.9246e-02_rb,1.9438e-02_rb,1.9839e-02_rb,2.0164e-02_rb,1.7716e-02_rb/)
      kao(:,13, 4) = (/ &
     &1.6236e-02_rb,1.6164e-02_rb,1.6551e-02_rb,1.6579e-02_rb,1.4849e-02_rb/)

      end subroutine lw_kgb101

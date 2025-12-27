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
!!          set parameters for kao(:,:,9:12)
!!
      subroutine lw_kgb103
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

      kao(:, 1, 9) = (/ &
     &1.0443e+01_rb,1.1048e+01_rb,1.1589e+01_rb,1.2078e+01_rb,1.2523e+01_rb/)
      kao(:, 2, 9) = (/ &
     &1.0305e+01_rb,1.0938e+01_rb,1.1513e+01_rb,1.2045e+01_rb,1.2531e+01_rb/)
      kao(:, 3, 9) = (/ &
     &9.8576e+00_rb,1.0529e+01_rb,1.1150e+01_rb,1.1712e+01_rb,1.2221e+01_rb/)
      kao(:, 4, 9) = (/ &
     &9.2371e+00_rb,9.9233e+00_rb,1.0551e+01_rb,1.1124e+01_rb,1.1647e+01_rb/)
      kao(:, 5, 9) = (/ &
     &8.4700e+00_rb,9.1602e+00_rb,9.7951e+00_rb,1.0380e+01_rb,1.0911e+01_rb/)
      kao(:, 6, 9) = (/ &
     &7.5880e+00_rb,8.2763e+00_rb,8.9154e+00_rb,9.5007e+00_rb,1.0024e+01_rb/)
      kao(:, 7, 9) = (/ &
     &6.6790e+00_rb,7.3503e+00_rb,7.9748e+00_rb,8.5357e+00_rb,9.0412e+00_rb/)
      kao(:, 8, 9) = (/ &
     &5.7712e+00_rb,6.4120e+00_rb,7.0009e+00_rb,7.5325e+00_rb,8.0088e+00_rb/)
      kao(:, 9, 9) = (/ &
     &4.8989e+00_rb,5.4909e+00_rb,6.0328e+00_rb,6.5241e+00_rb,6.9680e+00_rb/)
      kao(:,10, 9) = (/ &
     &4.1164e+00_rb,4.6491e+00_rb,5.1402e+00_rb,5.5886e+00_rb,6.0036e+00_rb/)
      kao(:,11, 9) = (/ &
     &3.6151e+00_rb,4.0753e+00_rb,4.5022e+00_rb,4.8994e+00_rb,5.2546e+00_rb/)
      kao(:,12, 9) = (/ &
     &3.1350e+00_rb,3.5329e+00_rb,3.9073e+00_rb,4.2435e+00_rb,4.5415e+00_rb/)
      kao(:,13, 9) = (/ &
     &2.6929e+00_rb,3.0365e+00_rb,3.3480e+00_rb,3.6295e+00_rb,3.8842e+00_rb/)
      kao(:, 1,10) = (/ &
     &1.9924e+01_rb,2.0919e+01_rb,2.2000e+01_rb,2.2972e+01_rb,2.3773e+01_rb/)
      kao(:, 2,10) = (/ &
     &2.0689e+01_rb,2.1909e+01_rb,2.3034e+01_rb,2.3987e+01_rb,2.4788e+01_rb/)
      kao(:, 3,10) = (/ &
     &2.0776e+01_rb,2.2053e+01_rb,2.3211e+01_rb,2.4276e+01_rb,2.5292e+01_rb/)
      kao(:, 4,10) = (/ &
     &2.0236e+01_rb,2.1839e+01_rb,2.3315e+01_rb,2.4662e+01_rb,2.5879e+01_rb/)
      kao(:, 5,10) = (/ &
     &1.9987e+01_rb,2.1808e+01_rb,2.3390e+01_rb,2.4727e+01_rb,2.5958e+01_rb/)
      kao(:, 6,10) = (/ &
     &1.9523e+01_rb,2.1285e+01_rb,2.2798e+01_rb,2.4193e+01_rb,2.5643e+01_rb/)
      kao(:, 7,10) = (/ &
     &1.8326e+01_rb,2.0057e+01_rb,2.1638e+01_rb,2.3310e+01_rb,2.4996e+01_rb/)
      kao(:, 8,10) = (/ &
     &1.6812e+01_rb,1.8592e+01_rb,2.0423e+01_rb,2.2320e+01_rb,2.4092e+01_rb/)
      kao(:, 9,10) = (/ &
     &1.5126e+01_rb,1.7059e+01_rb,1.9069e+01_rb,2.1023e+01_rb,2.2713e+01_rb/)
      kao(:,10,10) = (/ &
     &1.3603e+01_rb,1.5643e+01_rb,1.7659e+01_rb,1.9476e+01_rb,2.1017e+01_rb/)
      kao(:,11,10) = (/ &
     &1.3136e+01_rb,1.5065e+01_rb,1.6788e+01_rb,1.8298e+01_rb,1.9745e+01_rb/)
      kao(:,12,10) = (/ &
     &1.2291e+01_rb,1.3925e+01_rb,1.5389e+01_rb,1.6819e+01_rb,1.8254e+01_rb/)
      kao(:,13,10) = (/ &
     &1.1103e+01_rb,1.2502e+01_rb,1.3897e+01_rb,1.5277e+01_rb,1.6528e+01_rb/)
      kao(:, 1,11) = (/ &
     &2.4296e+01_rb,2.5293e+01_rb,2.6167e+01_rb,2.7129e+01_rb,2.8181e+01_rb/)
      kao(:, 2,11) = (/ &
     &2.5960e+01_rb,2.7094e+01_rb,2.8248e+01_rb,2.9459e+01_rb,3.0485e+01_rb/)
      kao(:, 3,11) = (/ &
     &2.6865e+01_rb,2.8307e+01_rb,2.9751e+01_rb,3.1102e+01_rb,3.2279e+01_rb/)
      kao(:, 4,11) = (/ &
     &2.7354e+01_rb,2.8995e+01_rb,3.0628e+01_rb,3.2100e+01_rb,3.3421e+01_rb/)
      kao(:, 5,11) = (/ &
     &2.6969e+01_rb,2.8892e+01_rb,3.0767e+01_rb,3.2621e+01_rb,3.4273e+01_rb/)
      kao(:, 6,11) = (/ &
     &2.6094e+01_rb,2.8469e+01_rb,3.0801e+01_rb,3.2980e+01_rb,3.4919e+01_rb/)
      kao(:, 7,11) = (/ &
     &2.5366e+01_rb,2.8149e+01_rb,3.0864e+01_rb,3.3154e+01_rb,3.5069e+01_rb/)
      kao(:, 8,11) = (/ &
     &2.4570e+01_rb,2.7568e+01_rb,3.0172e+01_rb,3.2348e+01_rb,3.4412e+01_rb/)
      kao(:, 9,11) = (/ &
     &2.3141e+01_rb,2.6027e+01_rb,2.8564e+01_rb,3.0892e+01_rb,3.3292e+01_rb/)
      kao(:,10,11) = (/ &
     &2.1375e+01_rb,2.4163e+01_rb,2.6748e+01_rb,2.9381e+01_rb,3.2113e+01_rb/)
      kao(:,11,11) = (/ &
     &2.0563e+01_rb,2.3343e+01_rb,2.6200e+01_rb,2.9082e+01_rb,3.1737e+01_rb/)
      kao(:,12,11) = (/ &
     &1.9759e+01_rb,2.2669e+01_rb,2.5592e+01_rb,2.8225e+01_rb,3.0484e+01_rb/)
      kao(:,13,11) = (/ &
     &1.9005e+01_rb,2.1823e+01_rb,2.4335e+01_rb,2.6560e+01_rb,2.8649e+01_rb/)
      kao(:, 1,12) = (/ &
     &2.8554e+01_rb,3.0109e+01_rb,3.1534e+01_rb,3.2820e+01_rb,3.3912e+01_rb/)
      kao(:, 2,12) = (/ &
     &3.1883e+01_rb,3.3547e+01_rb,3.5069e+01_rb,3.6425e+01_rb,3.7803e+01_rb/)
      kao(:, 3,12) = (/ &
     &3.5025e+01_rb,3.6824e+01_rb,3.8424e+01_rb,3.9924e+01_rb,4.1396e+01_rb/)
      kao(:, 4,12) = (/ &
     &3.7112e+01_rb,3.9231e+01_rb,4.1147e+01_rb,4.2932e+01_rb,4.4646e+01_rb/)
      kao(:, 5,12) = (/ &
     &3.8447e+01_rb,4.0861e+01_rb,4.3125e+01_rb,4.5199e+01_rb,4.7297e+01_rb/)
      kao(:, 6,12) = (/ &
     &3.8994e+01_rb,4.1741e+01_rb,4.4237e+01_rb,4.6637e+01_rb,4.8930e+01_rb/)
      kao(:, 7,12) = (/ &
     &3.8693e+01_rb,4.1674e+01_rb,4.4392e+01_rb,4.7287e+01_rb,5.0041e+01_rb/)
      kao(:, 8,12) = (/ &
     &3.7231e+01_rb,4.0534e+01_rb,4.3968e+01_rb,4.7486e+01_rb,5.0745e+01_rb/)
      kao(:, 9,12) = (/ &
     &3.5314e+01_rb,3.9308e+01_rb,4.3458e+01_rb,4.7418e+01_rb,5.1164e+01_rb/)
      kao(:,10,12) = (/ &
     &3.3674e+01_rb,3.8340e+01_rb,4.2937e+01_rb,4.7301e+01_rb,5.1085e+01_rb/)
      kao(:,11,12) = (/ &
     &3.4306e+01_rb,3.9228e+01_rb,4.3781e+01_rb,4.7771e+01_rb,5.1424e+01_rb/)
      kao(:,12,12) = (/ &
     &3.4432e+01_rb,3.9064e+01_rb,4.3274e+01_rb,4.7292e+01_rb,5.1289e+01_rb/)
      kao(:,13,12) = (/ &
     &3.3504e+01_rb,3.7950e+01_rb,4.2337e+01_rb,4.6726e+01_rb,5.1022e+01_rb/)

      end subroutine lw_kgb103

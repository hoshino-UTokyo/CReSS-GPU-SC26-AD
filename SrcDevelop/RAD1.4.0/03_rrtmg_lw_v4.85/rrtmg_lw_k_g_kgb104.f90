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
!!          set parameters for kao(:,:,13:16)
!!
      subroutine lw_kgb104
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

      kao(:, 1,13) = (/ &
     &3.3083e+01_rb,3.5231e+01_rb,3.7169e+01_rb,3.8902e+01_rb,4.0496e+01_rb/)
      kao(:, 2,13) = (/ &
     &3.7816e+01_rb,4.0283e+01_rb,4.2419e+01_rb,4.4328e+01_rb,4.6116e+01_rb/)
      kao(:, 3,13) = (/ &
     &4.2916e+01_rb,4.5668e+01_rb,4.8087e+01_rb,5.0367e+01_rb,5.2460e+01_rb/)
      kao(:, 4,13) = (/ &
     &4.7941e+01_rb,5.0937e+01_rb,5.3613e+01_rb,5.6174e+01_rb,5.8534e+01_rb/)
      kao(:, 5,13) = (/ &
     &5.2708e+01_rb,5.6040e+01_rb,5.9163e+01_rb,6.2033e+01_rb,6.4593e+01_rb/)
      kao(:, 6,13) = (/ &
     &5.6670e+01_rb,6.0546e+01_rb,6.4072e+01_rb,6.7245e+01_rb,7.0176e+01_rb/)
      kao(:, 7,13) = (/ &
     &5.9263e+01_rb,6.3946e+01_rb,6.8113e+01_rb,7.1753e+01_rb,7.5271e+01_rb/)
      kao(:, 8,13) = (/ &
     &6.0990e+01_rb,6.6317e+01_rb,7.1009e+01_rb,7.5310e+01_rb,7.9458e+01_rb/)
      kao(:, 9,13) = (/ &
     &6.1500e+01_rb,6.7310e+01_rb,7.2473e+01_rb,7.7444e+01_rb,8.1955e+01_rb/)
      kao(:,10,13) = (/ &
     &6.1095e+01_rb,6.7188e+01_rb,7.2950e+01_rb,7.8366e+01_rb,8.3528e+01_rb/)
      kao(:,11,13) = (/ &
     &6.2204e+01_rb,6.8688e+01_rb,7.5122e+01_rb,8.1365e+01_rb,8.7577e+01_rb/)
      kao(:,12,13) = (/ &
     &6.2610e+01_rb,7.0078e+01_rb,7.7313e+01_rb,8.4495e+01_rb,9.1274e+01_rb/)
      kao(:,13,13) = (/ &
     &6.3252e+01_rb,7.1343e+01_rb,7.9418e+01_rb,8.6956e+01_rb,9.3893e+01_rb/)
      kao(:, 1,14) = (/ &
     &4.3344e+01_rb,4.5304e+01_rb,4.7200e+01_rb,4.8977e+01_rb,5.0634e+01_rb/)
      kao(:, 2,14) = (/ &
     &4.8672e+01_rb,5.1219e+01_rb,5.3748e+01_rb,5.6207e+01_rb,5.8584e+01_rb/)
      kao(:, 3,14) = (/ &
     &5.3835e+01_rb,5.7385e+01_rb,6.0811e+01_rb,6.4094e+01_rb,6.7183e+01_rb/)
      kao(:, 4,14) = (/ &
     &5.9051e+01_rb,6.3463e+01_rb,6.7847e+01_rb,7.2047e+01_rb,7.6005e+01_rb/)
      kao(:, 5,14) = (/ &
     &6.4910e+01_rb,7.0372e+01_rb,7.5734e+01_rb,8.0867e+01_rb,8.5466e+01_rb/)
      kao(:, 6,14) = (/ &
     &7.1177e+01_rb,7.7830e+01_rb,8.4454e+01_rb,9.0502e+01_rb,9.5851e+01_rb/)
      kao(:, 7,14) = (/ &
     &7.8645e+01_rb,8.6401e+01_rb,9.3877e+01_rb,1.0067e+02_rb,1.0680e+02_rb/)
      kao(:, 8,14) = (/ &
     &8.6837e+01_rb,9.5972e+01_rb,1.0455e+02_rb,1.1218e+02_rb,1.1889e+02_rb/)
      kao(:, 9,14) = (/ &
     &9.4588e+01_rb,1.0523e+02_rb,1.1528e+02_rb,1.2405e+02_rb,1.3172e+02_rb/)
      kao(:,10,14) = (/ &
     &1.0185e+02_rb,1.1428e+02_rb,1.2585e+02_rb,1.3577e+02_rb,1.4474e+02_rb/)
      kao(:,11,14) = (/ &
     &1.1416e+02_rb,1.2817e+02_rb,1.4064e+02_rb,1.5135e+02_rb,1.6079e+02_rb/)
      kao(:,12,14) = (/ &
     &1.2658e+02_rb,1.4157e+02_rb,1.5496e+02_rb,1.6630e+02_rb,1.7621e+02_rb/)
      kao(:,13,14) = (/ &
     &1.3763e+02_rb,1.5352e+02_rb,1.6709e+02_rb,1.7925e+02_rb,1.9013e+02_rb/)
      kao(:, 1,15) = (/ &
     &5.4748e+01_rb,5.6924e+01_rb,5.8795e+01_rb,6.0375e+01_rb,6.1716e+01_rb/)
      kao(:, 2,15) = (/ &
     &6.5566e+01_rb,6.8287e+01_rb,7.0585e+01_rb,7.2551e+01_rb,7.4160e+01_rb/)
      kao(:, 3,15) = (/ &
     &7.7583e+01_rb,8.1027e+01_rb,8.3928e+01_rb,8.6358e+01_rb,8.8397e+01_rb/)
      kao(:, 4,15) = (/ &
     &8.9871e+01_rb,9.4173e+01_rb,9.7817e+01_rb,1.0084e+02_rb,1.0339e+02_rb/)
      kao(:, 5,15) = (/ &
     &1.0246e+02_rb,1.0776e+02_rb,1.1220e+02_rb,1.1592e+02_rb,1.1961e+02_rb/)
      kao(:, 6,15) = (/ &
     &1.1487e+02_rb,1.2133e+02_rb,1.2665e+02_rb,1.3188e+02_rb,1.3707e+02_rb/)
      kao(:, 7,15) = (/ &
     &1.2683e+02_rb,1.3446e+02_rb,1.4160e+02_rb,1.4872e+02_rb,1.5558e+02_rb/)
      kao(:, 8,15) = (/ &
     &1.3795e+02_rb,1.4756e+02_rb,1.5674e+02_rb,1.6606e+02_rb,1.7523e+02_rb/)
      kao(:, 9,15) = (/ &
     &1.4965e+02_rb,1.6170e+02_rb,1.7284e+02_rb,1.8422e+02_rb,1.9597e+02_rb/)
      kao(:,10,15) = (/ &
     &1.6376e+02_rb,1.7796e+02_rb,1.9116e+02_rb,2.0528e+02_rb,2.1939e+02_rb/)
      kao(:,11,15) = (/ &
     &1.8752e+02_rb,2.0335e+02_rb,2.1913e+02_rb,2.3662e+02_rb,2.5359e+02_rb/)
      kao(:,12,15) = (/ &
     &2.1430e+02_rb,2.3299e+02_rb,2.5241e+02_rb,2.7277e+02_rb,2.9286e+02_rb/)
      kao(:,13,15) = (/ &
     &2.4469e+02_rb,2.6735e+02_rb,2.9101e+02_rb,3.1438e+02_rb,3.3674e+02_rb/)
      kao(:, 1,16) = (/ &
     &5.6182e+01_rb,5.8534e+01_rb,6.0580e+01_rb,6.2336e+01_rb,6.3846e+01_rb/)
      kao(:, 2,16) = (/ &
     &6.7714e+01_rb,7.1253e+01_rb,7.3886e+01_rb,7.6167e+01_rb,7.8129e+01_rb/)
      kao(:, 3,16) = (/ &
     &8.2239e+01_rb,8.6274e+01_rb,8.9810e+01_rb,9.2865e+01_rb,9.5500e+01_rb/)
      kao(:, 4,16) = (/ &
     &9.7623e+01_rb,1.0299e+02_rb,1.0769e+02_rb,1.1178e+02_rb,1.1532e+02_rb/)
      kao(:, 5,16) = (/ &
     &1.1485e+02_rb,1.2193e+02_rb,1.2815e+02_rb,1.3359e+02_rb,1.3829e+02_rb/)
      kao(:, 6,16) = (/ &
     &1.3405e+02_rb,1.4336e+02_rb,1.5157e+02_rb,1.5875e+02_rb,1.6501e+02_rb/)
      kao(:, 7,16) = (/ &
     &1.5562e+02_rb,1.6780e+02_rb,1.7860e+02_rb,1.8805e+02_rb,1.9632e+02_rb/)
      kao(:, 8,16) = (/ &
     &1.7987e+02_rb,1.9571e+02_rb,2.0979e+02_rb,2.2214e+02_rb,2.3298e+02_rb/)
      kao(:, 9,16) = (/ &
     &2.0674e+02_rb,2.2717e+02_rb,2.4538e+02_rb,2.6145e+02_rb,2.7552e+02_rb/)
      kao(:,10,16) = (/ &
     &2.3822e+02_rb,2.6423e+02_rb,2.8739e+02_rb,3.0794e+02_rb,3.2592e+02_rb/)
      kao(:,11,16) = (/ &
     &2.8871e+02_rb,3.1996e+02_rb,3.4774e+02_rb,3.7226e+02_rb,3.9354e+02_rb/)
      kao(:,12,16) = (/ &
     &3.4916e+02_rb,3.8637e+02_rb,4.1932e+02_rb,4.4825e+02_rb,4.7321e+02_rb/)
      kao(:,13,16) = (/ &
     &4.2045e+02_rb,4.6429e+02_rb,5.0303e+02_rb,5.3671e+02_rb,5.6546e+02_rb/)

      end subroutine lw_kgb104

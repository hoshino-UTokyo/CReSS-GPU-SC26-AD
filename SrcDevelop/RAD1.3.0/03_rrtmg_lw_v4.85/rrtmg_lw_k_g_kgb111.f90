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
!!prev      subroutine lw_kgb11
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for kao(:,:,1:8)
!!
      subroutine lw_kgb111
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg11, only : fracrefao, fracrefbo, kao, kbo, kao_mo2, &
!!prev                            kbo_mo2, selfrefo, forrefo
      use rrlw_kg11, only : kao

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
     &4.9423e-02_rb,4.8938e-02_rb,4.8236e-02_rb,4.7630e-02_rb,4.7027e-02_rb/)
      kao(:, 2, 1) = (/ &
     &4.0264e-02_rb,3.9991e-02_rb,3.9414e-02_rb,3.8921e-02_rb,3.8455e-02_rb/)
      kao(:, 3, 1) = (/ &
     &3.2762e-02_rb,3.2612e-02_rb,3.2225e-02_rb,3.1842e-02_rb,3.1448e-02_rb/)
      kao(:, 4, 1) = (/ &
     &2.6929e-02_rb,2.6828e-02_rb,2.6636e-02_rb,2.6311e-02_rb,2.5979e-02_rb/)
      kao(:, 5, 1) = (/ &
     &2.2254e-02_rb,2.2190e-02_rb,2.2139e-02_rb,2.1844e-02_rb,2.1606e-02_rb/)
      kao(:, 6, 1) = (/ &
     &1.8428e-02_rb,1.8406e-02_rb,1.8346e-02_rb,1.8223e-02_rb,1.7986e-02_rb/)
      kao(:, 7, 1) = (/ &
     &1.5302e-02_rb,1.5279e-02_rb,1.5227e-02_rb,1.5173e-02_rb,1.5001e-02_rb/)
      kao(:, 8, 1) = (/ &
     &1.2917e-02_rb,1.2821e-02_rb,1.2762e-02_rb,1.2673e-02_rb,1.2614e-02_rb/)
      kao(:, 9, 1) = (/ &
     &1.4361e-02_rb,1.3355e-02_rb,1.2836e-02_rb,1.2424e-02_rb,1.2069e-02_rb/)
      kao(:,10, 1) = (/ &
     &1.9078e-02_rb,1.9077e-02_rb,1.9301e-02_rb,1.9494e-02_rb,1.9600e-02_rb/)
      kao(:,11, 1) = (/ &
     &1.6651e-02_rb,1.7279e-02_rb,1.7762e-02_rb,1.8212e-02_rb,1.8508e-02_rb/)
      kao(:,12, 1) = (/ &
     &1.4359e-02_rb,1.4912e-02_rb,1.5370e-02_rb,1.5766e-02_rb,1.6068e-02_rb/)
      kao(:,13, 1) = (/ &
     &1.2190e-02_rb,1.2693e-02_rb,1.3086e-02_rb,1.3355e-02_rb,1.3602e-02_rb/)
      kao(:, 1, 2) = (/ &
     &1.3874e-01_rb,1.3507e-01_rb,1.3188e-01_rb,1.2875e-01_rb,1.2573e-01_rb/)
      kao(:, 2, 2) = (/ &
     &1.1449e-01_rb,1.1143e-01_rb,1.0880e-01_rb,1.0626e-01_rb,1.0384e-01_rb/)
      kao(:, 3, 2) = (/ &
     &9.4610e-02_rb,9.2157e-02_rb,8.9899e-02_rb,8.7836e-02_rb,8.5813e-02_rb/)
      kao(:, 4, 2) = (/ &
     &7.8921e-02_rb,7.6917e-02_rb,7.5042e-02_rb,7.3283e-02_rb,7.1619e-02_rb/)
      kao(:, 5, 2) = (/ &
     &6.6025e-02_rb,6.4395e-02_rb,6.2788e-02_rb,6.1406e-02_rb,5.9938e-02_rb/)
      kao(:, 6, 2) = (/ &
     &5.5198e-02_rb,5.3979e-02_rb,5.2652e-02_rb,5.1421e-02_rb,5.0291e-02_rb/)
      kao(:, 7, 2) = (/ &
     &4.5982e-02_rb,4.5142e-02_rb,4.4142e-02_rb,4.3045e-02_rb,4.2121e-02_rb/)
      kao(:, 8, 2) = (/ &
     &3.8128e-02_rb,3.7632e-02_rb,3.6913e-02_rb,3.6127e-02_rb,3.5207e-02_rb/)
      kao(:, 9, 2) = (/ &
     &2.9988e-02_rb,3.0199e-02_rb,2.9924e-02_rb,2.9588e-02_rb,2.9106e-02_rb/)
      kao(:,10, 2) = (/ &
     &4.8557e-02_rb,4.4974e-02_rb,3.7588e-02_rb,3.1508e-02_rb,2.7381e-02_rb/)
      kao(:,11, 2) = (/ &
     &4.6113e-02_rb,4.4271e-02_rb,4.2375e-02_rb,4.0024e-02_rb,3.7665e-02_rb/)
      kao(:,12, 2) = (/ &
     &4.0305e-02_rb,3.9408e-02_rb,3.7504e-02_rb,3.5981e-02_rb,3.4390e-02_rb/)
      kao(:,13, 2) = (/ &
     &3.3779e-02_rb,3.2735e-02_rb,3.1178e-02_rb,2.9814e-02_rb,2.8271e-02_rb/)
      kao(:, 1, 3) = (/ &
     &2.4150e-01_rb,2.3616e-01_rb,2.3111e-01_rb,2.2673e-01_rb,2.2290e-01_rb/)
      kao(:, 2, 3) = (/ &
     &2.0007e-01_rb,1.9568e-01_rb,1.9176e-01_rb,1.8825e-01_rb,1.8518e-01_rb/)
      kao(:, 3, 3) = (/ &
     &1.6573e-01_rb,1.6221e-01_rb,1.5908e-01_rb,1.5616e-01_rb,1.5376e-01_rb/)
      kao(:, 4, 3) = (/ &
     &1.3876e-01_rb,1.3573e-01_rb,1.3298e-01_rb,1.3059e-01_rb,1.2856e-01_rb/)
      kao(:, 5, 3) = (/ &
     &1.1663e-01_rb,1.1399e-01_rb,1.1156e-01_rb,1.0951e-01_rb,1.0785e-01_rb/)
      kao(:, 6, 3) = (/ &
     &9.8259e-02_rb,9.5861e-02_rb,9.3802e-02_rb,9.1948e-02_rb,9.0492e-02_rb/)
      kao(:, 7, 3) = (/ &
     &8.2887e-02_rb,8.0616e-02_rb,7.8766e-02_rb,7.7206e-02_rb,7.5856e-02_rb/)
      kao(:, 8, 3) = (/ &
     &6.9730e-02_rb,6.7730e-02_rb,6.5964e-02_rb,6.4565e-02_rb,6.3437e-02_rb/)
      kao(:, 9, 3) = (/ &
     &5.8101e-02_rb,5.6452e-02_rb,5.5028e-02_rb,5.3674e-02_rb,5.2557e-02_rb/)
      kao(:,10, 3) = (/ &
     &3.7373e-02_rb,3.5041e-02_rb,3.8621e-02_rb,4.1558e-02_rb,4.3012e-02_rb/)
      kao(:,11, 3) = (/ &
     &5.3584e-02_rb,4.0890e-02_rb,3.3882e-02_rb,2.9034e-02_rb,2.6383e-02_rb/)
      kao(:,12, 3) = (/ &
     &5.8428e-02_rb,4.4689e-02_rb,3.5111e-02_rb,2.9172e-02_rb,2.4695e-02_rb/)
      kao(:,13, 3) = (/ &
     &4.7904e-02_rb,3.6965e-02_rb,2.8601e-02_rb,2.3919e-02_rb,2.0680e-02_rb/)
      kao(:, 1, 4) = (/ &
     &4.2548e-01_rb,4.1807e-01_rb,4.1094e-01_rb,4.0423e-01_rb,3.9806e-01_rb/)
      kao(:, 2, 4) = (/ &
     &3.5494e-01_rb,3.4863e-01_rb,3.4235e-01_rb,3.3657e-01_rb,3.3155e-01_rb/)
      kao(:, 3, 4) = (/ &
     &2.9648e-01_rb,2.9099e-01_rb,2.8545e-01_rb,2.8071e-01_rb,2.7648e-01_rb/)
      kao(:, 4, 4) = (/ &
     &2.4939e-01_rb,2.4483e-01_rb,2.4021e-01_rb,2.3626e-01_rb,2.3256e-01_rb/)
      kao(:, 5, 4) = (/ &
     &2.1025e-01_rb,2.0654e-01_rb,2.0285e-01_rb,1.9937e-01_rb,1.9615e-01_rb/)
      kao(:, 6, 4) = (/ &
     &1.7714e-01_rb,1.7417e-01_rb,1.7124e-01_rb,1.6826e-01_rb,1.6547e-01_rb/)
      kao(:, 7, 4) = (/ &
     &1.4884e-01_rb,1.4658e-01_rb,1.4422e-01_rb,1.4179e-01_rb,1.3940e-01_rb/)
      kao(:, 8, 4) = (/ &
     &1.2471e-01_rb,1.2286e-01_rb,1.2117e-01_rb,1.1917e-01_rb,1.1721e-01_rb/)
      kao(:, 9, 4) = (/ &
     &1.0337e-01_rb,1.0218e-01_rb,1.0077e-01_rb,9.9333e-02_rb,9.7869e-02_rb/)
      kao(:,10, 4) = (/ &
     &8.5198e-02_rb,8.3413e-02_rb,8.1099e-02_rb,7.8395e-02_rb,7.6293e-02_rb/)
      kao(:,11, 4) = (/ &
     &6.0820e-02_rb,6.7915e-02_rb,6.7626e-02_rb,6.6844e-02_rb,6.6371e-02_rb/)
      kao(:,12, 4) = (/ &
     &4.0000e-02_rb,4.8190e-02_rb,5.3574e-02_rb,5.5575e-02_rb,5.5151e-02_rb/)
      kao(:,13, 4) = (/ &
     &3.3587e-02_rb,4.0386e-02_rb,4.5135e-02_rb,4.6401e-02_rb,4.5805e-02_rb/)
      kao(:, 1, 5) = (/ &
     &7.7524e-01_rb,7.6536e-01_rb,7.5522e-01_rb,7.4431e-01_rb,7.3311e-01_rb/)
      kao(:, 2, 5) = (/ &
     &6.5185e-01_rb,6.4405e-01_rb,6.3568e-01_rb,6.2684e-01_rb,6.1765e-01_rb/)
      kao(:, 3, 5) = (/ &
     &5.4580e-01_rb,5.3942e-01_rb,5.3289e-01_rb,5.2546e-01_rb,5.1768e-01_rb/)
      kao(:, 4, 5) = (/ &
     &4.5945e-01_rb,4.5416e-01_rb,4.4868e-01_rb,4.4229e-01_rb,4.3592e-01_rb/)
      kao(:, 5, 5) = (/ &
     &3.8777e-01_rb,3.8350e-01_rb,3.7876e-01_rb,3.7352e-01_rb,3.6853e-01_rb/)
      kao(:, 6, 5) = (/ &
     &3.2750e-01_rb,3.2403e-01_rb,3.2009e-01_rb,3.1598e-01_rb,3.1216e-01_rb/)
      kao(:, 7, 5) = (/ &
     &2.7610e-01_rb,2.7333e-01_rb,2.7016e-01_rb,2.6699e-01_rb,2.6402e-01_rb/)
      kao(:, 8, 5) = (/ &
     &2.3237e-01_rb,2.3021e-01_rb,2.2764e-01_rb,2.2511e-01_rb,2.2263e-01_rb/)
      kao(:, 9, 5) = (/ &
     &1.9461e-01_rb,1.9312e-01_rb,1.9116e-01_rb,1.8911e-01_rb,1.8708e-01_rb/)
      kao(:,10, 5) = (/ &
     &1.4772e-01_rb,1.4997e-01_rb,1.5067e-01_rb,1.5164e-01_rb,1.5169e-01_rb/)
      kao(:,11, 5) = (/ &
     &1.1317e-01_rb,1.1322e-01_rb,1.1552e-01_rb,1.1824e-01_rb,1.1886e-01_rb/)
      kao(:,12, 5) = (/ &
     &9.7948e-02_rb,9.4310e-02_rb,9.3527e-02_rb,9.3094e-02_rb,9.5273e-02_rb/)
      kao(:,13, 5) = (/ &
     &8.1443e-02_rb,7.8686e-02_rb,7.8076e-02_rb,7.8298e-02_rb,8.0221e-02_rb/)
      kao(:, 1, 6) = (/ &
     &1.5697e+00_rb,1.5485e+00_rb,1.5270e+00_rb,1.5062e+00_rb,1.4845e+00_rb/)
      kao(:, 2, 6) = (/ &
     &1.3465e+00_rb,1.3286e+00_rb,1.3101e+00_rb,1.2912e+00_rb,1.2719e+00_rb/)
      kao(:, 3, 6) = (/ &
     &1.1474e+00_rb,1.1322e+00_rb,1.1161e+00_rb,1.0997e+00_rb,1.0837e+00_rb/)
      kao(:, 4, 6) = (/ &
     &9.7861e-01_rb,9.6591e-01_rb,9.5236e-01_rb,9.3871e-01_rb,9.2522e-01_rb/)
      kao(:, 5, 6) = (/ &
     &8.3345e-01_rb,8.2282e-01_rb,8.1147e-01_rb,8.0015e-01_rb,7.8882e-01_rb/)
      kao(:, 6, 6) = (/ &
     &7.0768e-01_rb,6.9919e-01_rb,6.8975e-01_rb,6.8018e-01_rb,6.7024e-01_rb/)
      kao(:, 7, 6) = (/ &
     &5.9916e-01_rb,5.9240e-01_rb,5.8444e-01_rb,5.7632e-01_rb,5.6762e-01_rb/)
      kao(:, 8, 6) = (/ &
     &5.0612e-01_rb,5.0051e-01_rb,4.9423e-01_rb,4.8731e-01_rb,4.7989e-01_rb/)
      kao(:, 9, 6) = (/ &
     &4.2622e-01_rb,4.2196e-01_rb,4.1690e-01_rb,4.1124e-01_rb,4.0525e-01_rb/)
      kao(:,10, 6) = (/ &
     &3.5041e-01_rb,3.4937e-01_rb,3.4686e-01_rb,3.4326e-01_rb,3.3900e-01_rb/)
      kao(:,11, 6) = (/ &
     &2.8138e-01_rb,2.8270e-01_rb,2.8314e-01_rb,2.8109e-01_rb,2.7863e-01_rb/)
      kao(:,12, 6) = (/ &
     &2.2569e-01_rb,2.3115e-01_rb,2.3249e-01_rb,2.3233e-01_rb,2.3100e-01_rb/)
      kao(:,13, 6) = (/ &
     &1.8817e-01_rb,1.9206e-01_rb,1.9341e-01_rb,1.9316e-01_rb,1.9216e-01_rb/)
      kao(:, 1, 7) = (/ &
     &3.4196e+00_rb,3.3801e+00_rb,3.3399e+00_rb,3.2946e+00_rb,3.2477e+00_rb/)
      kao(:, 2, 7) = (/ &
     &3.0295e+00_rb,2.9903e+00_rb,2.9469e+00_rb,2.9032e+00_rb,2.8599e+00_rb/)
      kao(:, 3, 7) = (/ &
     &2.6483e+00_rb,2.6096e+00_rb,2.5696e+00_rb,2.5316e+00_rb,2.4936e+00_rb/)
      kao(:, 4, 7) = (/ &
     &2.3044e+00_rb,2.2705e+00_rb,2.2374e+00_rb,2.2055e+00_rb,2.1724e+00_rb/)
      kao(:, 5, 7) = (/ &
     &1.9970e+00_rb,1.9688e+00_rb,1.9416e+00_rb,1.9143e+00_rb,1.8852e+00_rb/)
      kao(:, 6, 7) = (/ &
     &1.7224e+00_rb,1.6992e+00_rb,1.6760e+00_rb,1.6526e+00_rb,1.6276e+00_rb/)
      kao(:, 7, 7) = (/ &
     &1.4770e+00_rb,1.4579e+00_rb,1.4384e+00_rb,1.4185e+00_rb,1.3976e+00_rb/)
      kao(:, 8, 7) = (/ &
     &1.2593e+00_rb,1.2433e+00_rb,1.2274e+00_rb,1.2112e+00_rb,1.1939e+00_rb/)
      kao(:, 9, 7) = (/ &
     &1.0697e+00_rb,1.0563e+00_rb,1.0434e+00_rb,1.0295e+00_rb,1.0150e+00_rb/)
      kao(:,10, 7) = (/ &
     &9.0213e-01_rb,8.9139e-01_rb,8.8119e-01_rb,8.6954e-01_rb,8.5716e-01_rb/)
      kao(:,11, 7) = (/ &
     &7.5061e-01_rb,7.4446e-01_rb,7.3597e-01_rb,7.2638e-01_rb,7.1615e-01_rb/)
      kao(:,12, 7) = (/ &
     &6.2457e-01_rb,6.1924e-01_rb,6.1219e-01_rb,6.0465e-01_rb,5.9633e-01_rb/)
      kao(:,13, 7) = (/ &
     &5.2204e-01_rb,5.1758e-01_rb,5.1182e-01_rb,5.0530e-01_rb,4.9787e-01_rb/)
      kao(:, 1, 8) = (/ &
     &7.9795e+00_rb,7.8457e+00_rb,7.7132e+00_rb,7.5890e+00_rb,7.4645e+00_rb/)
      kao(:, 2, 8) = (/ &
     &7.5076e+00_rb,7.3733e+00_rb,7.2519e+00_rb,7.1316e+00_rb,7.0099e+00_rb/)
      kao(:, 3, 8) = (/ &
     &6.9340e+00_rb,6.8148e+00_rb,6.7036e+00_rb,6.5875e+00_rb,6.4766e+00_rb/)
      kao(:, 4, 8) = (/ &
     &6.3414e+00_rb,6.2366e+00_rb,6.1298e+00_rb,6.0215e+00_rb,5.9183e+00_rb/)
      kao(:, 5, 8) = (/ &
     &5.7373e+00_rb,5.6423e+00_rb,5.5412e+00_rb,5.4426e+00_rb,5.3434e+00_rb/)
      kao(:, 6, 8) = (/ &
     &5.1373e+00_rb,5.0482e+00_rb,4.9577e+00_rb,4.8621e+00_rb,4.7659e+00_rb/)
      kao(:, 7, 8) = (/ &
     &4.5492e+00_rb,4.4684e+00_rb,4.3821e+00_rb,4.2928e+00_rb,4.2067e+00_rb/)
      kao(:, 8, 8) = (/ &
     &3.9884e+00_rb,3.9140e+00_rb,3.8349e+00_rb,3.7563e+00_rb,3.6821e+00_rb/)
      kao(:, 9, 8) = (/ &
     &3.4656e+00_rb,3.3995e+00_rb,3.3300e+00_rb,3.2636e+00_rb,3.2010e+00_rb/)
      kao(:,10, 8) = (/ &
     &2.9841e+00_rb,2.9267e+00_rb,2.8680e+00_rb,2.8124e+00_rb,2.7592e+00_rb/)
      kao(:,11, 8) = (/ &
     &2.5294e+00_rb,2.4795e+00_rb,2.4321e+00_rb,2.3862e+00_rb,2.3418e+00_rb/)
      kao(:,12, 8) = (/ &
     &2.1315e+00_rb,2.0918e+00_rb,2.0527e+00_rb,2.0146e+00_rb,1.9776e+00_rb/)
      kao(:,13, 8) = (/ &
     &1.7932e+00_rb,1.7598e+00_rb,1.7272e+00_rb,1.6957e+00_rb,1.6643e+00_rb/)

      end subroutine lw_kgb111

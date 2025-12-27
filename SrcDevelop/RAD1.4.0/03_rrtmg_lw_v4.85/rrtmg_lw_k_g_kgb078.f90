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
!!prev      subroutine lw_kgb07
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for kbo_mco2(:,:)
!!
      subroutine lw_kgb078
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg07, only : fracrefao, fracrefbo, kao, kbo, kao_mco2, &
!!prev                            kbo_mco2, selfrefo, forrefo
      use rrlw_kg07, only : kbo_mco2

      implicit none
      save

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

      kbo_mco2(:, 1) = (/ &
     & 3.72069e-06_rb, 4.81866e-06_rb, 6.24064e-06_rb, 8.08226e-06_rb, 1.04673e-05_rb, &
     & 1.35562e-05_rb, 1.75567e-05_rb, 2.27376e-05_rb, 2.94475e-05_rb, 3.81375e-05_rb, &
     & 4.93918e-05_rb, 6.39674e-05_rb, 8.28441e-05_rb, 1.07291e-04_rb, 1.38953e-04_rb, &
     & 1.79958e-04_rb, 2.33064e-04_rb, 3.01840e-04_rb, 3.90913e-04_rb/)
      kbo_mco2(:, 2) = (/ &
     & 8.14357e-06_rb, 1.06031e-05_rb, 1.38056e-05_rb, 1.79752e-05_rb, 2.34041e-05_rb, &
     & 3.04728e-05_rb, 3.96763e-05_rb, 5.16596e-05_rb, 6.72622e-05_rb, 8.75770e-05_rb, &
     & 1.14027e-04_rb, 1.48467e-04_rb, 1.93307e-04_rb, 2.51691e-04_rb, 3.27708e-04_rb, &
     & 4.26685e-04_rb, 5.55555e-04_rb, 7.23346e-04_rb, 9.41814e-04_rb/)
      kbo_mco2(:, 3) = (/ &
     & 1.09367e-05_rb, 1.42063e-05_rb, 1.84533e-05_rb, 2.39701e-05_rb, 3.11362e-05_rb, &
     & 4.04446e-05_rb, 5.25358e-05_rb, 6.82417e-05_rb, 8.86432e-05_rb, 1.15144e-04_rb, &
     & 1.49567e-04_rb, 1.94281e-04_rb, 2.52363e-04_rb, 3.27809e-04_rb, 4.25810e-04_rb, &
     & 5.53109e-04_rb, 7.18466e-04_rb, 9.33256e-04_rb, 1.21226e-03_rb/)
      kbo_mco2(:, 4) = (/ &
     & 1.76192e-05_rb, 2.27752e-05_rb, 2.94401e-05_rb, 3.80553e-05_rb, 4.91916e-05_rb, &
     & 6.35867e-05_rb, 8.21944e-05_rb, 1.06247e-04_rb, 1.37339e-04_rb, 1.77529e-04_rb, &
     & 2.29480e-04_rb, 2.96635e-04_rb, 3.83440e-04_rb, 4.95648e-04_rb, 6.40691e-04_rb, &
     & 8.28180e-04_rb, 1.07054e-03_rb, 1.38381e-03_rb, 1.78876e-03_rb/)
      kbo_mco2(:, 5) = (/ &
     & 3.72142e-05_rb, 4.78603e-05_rb, 6.15520e-05_rb, 7.91605e-05_rb, 1.01806e-04_rb, &
     & 1.30931e-04_rb, 1.68387e-04_rb, 2.16558e-04_rb, 2.78510e-04_rb, 3.58185e-04_rb, &
     & 4.60653e-04_rb, 5.92435e-04_rb, 7.61915e-04_rb, 9.79881e-04_rb, 1.26020e-03_rb, &
     & 1.62071e-03_rb, 2.08436e-03_rb, 2.68064e-03_rb, 3.44751e-03_rb/)
      kbo_mco2(:, 6) = (/ &
     & 7.74131e-05_rb, 9.98876e-05_rb, 1.28887e-04_rb, 1.66305e-04_rb, 2.14587e-04_rb, &
     & 2.76886e-04_rb, 3.57272e-04_rb, 4.60994e-04_rb, 5.94831e-04_rb, 7.67521e-04_rb, &
     & 9.90348e-04_rb, 1.27787e-03_rb, 1.64886e-03_rb, 2.12755e-03_rb, 2.74522e-03_rb, &
     & 3.54221e-03_rb, 4.57059e-03_rb, 5.89752e-03_rb, 7.60968e-03_rb/)
      kbo_mco2(:, 7) = (/ &
     & 1.32294e-04_rb, 1.70977e-04_rb, 2.20973e-04_rb, 2.85587e-04_rb, 3.69095e-04_rb, &
     & 4.77022e-04_rb, 6.16507e-04_rb, 7.96779e-04_rb, 1.02976e-03_rb, 1.33088e-03_rb, &
     & 1.72004e-03_rb, 2.22299e-03_rb, 2.87301e-03_rb, 3.71310e-03_rb, 4.79884e-03_rb, &
     & 6.20207e-03_rb, 8.01561e-03_rb, 1.03594e-02_rb, 1.33886e-02_rb/)
      kbo_mco2(:, 8) = (/ &
     & 3.59868e-05_rb, 4.63611e-05_rb, 5.97261e-05_rb, 7.69439e-05_rb, 9.91253e-05_rb, &
     & 1.27701e-04_rb, 1.64515e-04_rb, 2.11941e-04_rb, 2.73040e-04_rb, 3.51752e-04_rb, &
     & 4.53155e-04_rb, 5.83790e-04_rb, 7.52085e-04_rb, 9.68897e-04_rb, 1.24821e-03_rb, &
     & 1.60804e-03_rb, 2.07161e-03_rb, 2.66882e-03_rb, 3.43818e-03_rb/)
      kbo_mco2(:, 9) = (/ &
     & 5.09543e-05_rb, 6.60510e-05_rb, 8.56205e-05_rb, 1.10988e-04_rb, 1.43872e-04_rb, &
     & 1.86498e-04_rb, 2.41753e-04_rb, 3.13380e-04_rb, 4.06228e-04_rb, 5.26585e-04_rb, &
     & 6.82601e-04_rb, 8.84842e-04_rb, 1.14700e-03_rb, 1.48684e-03_rb, 1.92735e-03_rb, &
     & 2.49839e-03_rb, 3.23861e-03_rb, 4.19814e-03_rb, 5.44196e-03_rb/)
      kbo_mco2(:,10) = (/ &
     & 2.08253e-05_rb, 2.64900e-05_rb, 3.36954e-05_rb, 4.28609e-05_rb, 5.45194e-05_rb, &
     & 6.93491e-05_rb, 8.82125e-05_rb, 1.12207e-04_rb, 1.42728e-04_rb, 1.81551e-04_rb, &
     & 2.30935e-04_rb, 2.93751e-04_rb, 3.73653e-04_rb, 4.75290e-04_rb, 6.04572e-04_rb, &
     & 7.69021e-04_rb, 9.78201e-04_rb, 1.24428e-03_rb, 1.58273e-03_rb/)
      kbo_mco2(:,11) = (/ &
     & 2.08953e-05_rb, 2.65543e-05_rb, 3.37459e-05_rb, 4.28852e-05_rb, 5.44996e-05_rb, &
     & 6.92595e-05_rb, 8.80169e-05_rb, 1.11854e-04_rb, 1.42147e-04_rb, 1.80644e-04_rb, &
     & 2.29568e-04_rb, 2.91741e-04_rb, 3.70752e-04_rb, 4.71161e-04_rb, 5.98764e-04_rb, &
     & 7.60925e-04_rb, 9.67005e-04_rb, 1.22889e-03_rb, 1.56171e-03_rb/)
      kbo_mco2(:,12) = (/ &
     & 2.65295e-05_rb, 3.36318e-05_rb, 4.26356e-05_rb, 5.40498e-05_rb, 6.85198e-05_rb, &
     & 8.68636e-05_rb, 1.10118e-04_rb, 1.39599e-04_rb, 1.76972e-04_rb, 2.24350e-04_rb, &
     & 2.84412e-04_rb, 3.60553e-04_rb, 4.57079e-04_rb, 5.79446e-04_rb, 7.34572e-04_rb, &
     & 9.31230e-04_rb, 1.18053e-03_rb, 1.49658e-03_rb, 1.89724e-03_rb/)
      kbo_mco2(:,13) = (/ &
     & 3.45358e-05_rb, 4.36743e-05_rb, 5.52309e-05_rb, 6.98455e-05_rb, 8.83273e-05_rb, &
     & 1.11700e-04_rb, 1.41256e-04_rb, 1.78634e-04_rb, 2.25902e-04_rb, 2.85678e-04_rb, &
     & 3.61271e-04_rb, 4.56867e-04_rb, 5.77758e-04_rb, 7.30639e-04_rb, 9.23973e-04_rb, &
     & 1.16847e-03_rb, 1.47765e-03_rb, 1.86865e-03_rb, 2.36311e-03_rb/)
      kbo_mco2(:,14) = (/ &
     & 3.99721e-05_rb, 5.12343e-05_rb, 6.56698e-05_rb, 8.41725e-05_rb, 1.07888e-04_rb, &
     & 1.38286e-04_rb, 1.77249e-04_rb, 2.27190e-04_rb, 2.91201e-04_rb, 3.73248e-04_rb, &
     & 4.78412e-04_rb, 6.13207e-04_rb, 7.85980e-04_rb, 1.00743e-03_rb, 1.29128e-03_rb, &
     & 1.65510e-03_rb, 2.12144e-03_rb, 2.71916e-03_rb, 3.48529e-03_rb/)
      kbo_mco2(:,15) = (/ &
     & 8.51533e-06_rb, 1.23021e-05_rb, 1.77730e-05_rb, 2.56767e-05_rb, 3.70953e-05_rb, &
     & 5.35918e-05_rb, 7.74243e-05_rb, 1.11855e-04_rb, 1.61598e-04_rb, 2.33461e-04_rb, &
     & 3.37283e-04_rb, 4.87275e-04_rb, 7.03968e-04_rb, 1.01703e-03_rb, 1.46930e-03_rb, &
     & 2.12271e-03_rb, 3.06670e-03_rb, 4.43047e-03_rb, 6.40072e-03_rb/)
      kbo_mco2(:,16) = (/ &
     & 2.93050e-06_rb, 3.65298e-06_rb, 4.55358e-06_rb, 5.67622e-06_rb, 7.07564e-06_rb, &
     & 8.82006e-06_rb, 1.09945e-05_rb, 1.37051e-05_rb, 1.70840e-05_rb, 2.12959e-05_rb, &
     & 2.65461e-05_rb, 3.30908e-05_rb, 4.12490e-05_rb, 5.14185e-05_rb, 6.40952e-05_rb, &
     & 7.98972e-05_rb, 9.95951e-05_rb, 1.24149e-04_rb, 1.54757e-04_rb/)

      end subroutine lw_kgb078

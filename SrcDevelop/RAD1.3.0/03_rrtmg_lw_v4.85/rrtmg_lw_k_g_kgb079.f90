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
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb079
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg07, only : fracrefao, fracrefbo, kao, kbo, kao_mco2, &
!!prev                            kbo_mco2, selfrefo, forrefo
      use rrlw_kg07, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296_rb,260_rb,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &2.0677e-07_rb,2.0363e-07_rb,2.0583e-07_rb,2.0547e-07_rb,2.0267e-07_rb,2.0154e-07_rb, &
     &2.0190e-07_rb,2.0103e-07_rb,1.9869e-07_rb,1.9663e-07_rb,1.9701e-07_rb,2.0103e-07_rb, &
     &2.0527e-07_rb,2.0206e-07_rb,2.0364e-07_rb,2.0364e-07_rb/)
      forrefo(2,:) = (/ &
     &2.2427e-07_rb,2.1489e-07_rb,2.0453e-07_rb,1.9710e-07_rb,1.9650e-07_rb,1.9738e-07_rb, &
     &1.9767e-07_rb,1.9769e-07_rb,1.9940e-07_rb,1.9846e-07_rb,1.9898e-07_rb,1.9853e-07_rb, &
     &2.0000e-07_rb,2.0517e-07_rb,2.0482e-07_rb,2.0482e-07_rb/)
      forrefo(3,:) = (/ &
     &2.2672e-07_rb,2.1706e-07_rb,2.0571e-07_rb,1.9747e-07_rb,1.9706e-07_rb,1.9698e-07_rb, &
     &1.9781e-07_rb,1.9774e-07_rb,1.9724e-07_rb,1.9714e-07_rb,1.9751e-07_rb,1.9758e-07_rb, &
     &1.9840e-07_rb,1.9968e-07_rb,1.9931e-07_rb,1.9880e-07_rb/)
      forrefo(4,:) = (/ &
     &2.2191e-07_rb,2.0899e-07_rb,2.0265e-07_rb,2.0101e-07_rb,2.0034e-07_rb,2.0021e-07_rb, &
     &1.9987e-07_rb,1.9978e-07_rb,1.9902e-07_rb,1.9742e-07_rb,1.9672e-07_rb,1.9615e-07_rb, &
     &1.9576e-07_rb,1.9540e-07_rb,1.9588e-07_rb,1.9590e-07_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 5.18832e-02_rb, 4.28690e-02_rb, 3.54210e-02_rb, 2.92670e-02_rb, 2.41822e-02_rb, &
     & 1.99808e-02_rb, 1.65093e-02_rb, 1.36410e-02_rb, 1.12710e-02_rb, 9.31280e-03_rb/)
      selfrefo(:, 2) = (/ &
     & 4.36030e-02_rb, 3.78379e-02_rb, 3.28350e-02_rb, 2.84936e-02_rb, 2.47262e-02_rb, &
     & 2.14569e-02_rb, 1.86199e-02_rb, 1.61580e-02_rb, 1.40216e-02_rb, 1.21677e-02_rb/)
      selfrefo(:, 3) = (/ &
     & 4.26492e-02_rb, 3.71443e-02_rb, 3.23500e-02_rb, 2.81745e-02_rb, 2.45379e-02_rb, &
     & 2.13707e-02_rb, 1.86124e-02_rb, 1.62100e-02_rb, 1.41177e-02_rb, 1.22955e-02_rb/)
      selfrefo(:, 4) = (/ &
     & 4.03591e-02_rb, 3.54614e-02_rb, 3.11580e-02_rb, 2.73769e-02_rb, 2.40546e-02_rb, &
     & 2.11355e-02_rb, 1.85706e-02_rb, 1.63170e-02_rb, 1.43369e-02_rb, 1.25970e-02_rb/)
      selfrefo(:, 5) = (/ &
     & 3.94512e-02_rb, 3.46232e-02_rb, 3.03860e-02_rb, 2.66674e-02_rb, 2.34038e-02_rb, &
     & 2.05397e-02_rb, 1.80260e-02_rb, 1.58200e-02_rb, 1.38839e-02_rb, 1.21848e-02_rb/)
      selfrefo(:, 6) = (/ &
     & 3.90567e-02_rb, 3.40694e-02_rb, 2.97190e-02_rb, 2.59241e-02_rb, 2.26138e-02_rb, &
     & 1.97261e-02_rb, 1.72072e-02_rb, 1.50100e-02_rb, 1.30933e-02_rb, 1.14214e-02_rb/)
      selfrefo(:, 7) = (/ &
     & 3.85397e-02_rb, 3.36462e-02_rb, 2.93740e-02_rb, 2.56443e-02_rb, 2.23881e-02_rb, &
     & 1.95454e-02_rb, 1.70636e-02_rb, 1.48970e-02_rb, 1.30055e-02_rb, 1.13541e-02_rb/)
      selfrefo(:, 8) = (/ &
     & 3.79692e-02_rb, 3.31360e-02_rb, 2.89180e-02_rb, 2.52369e-02_rb, 2.20245e-02_rb, &
     & 1.92209e-02_rb, 1.67742e-02_rb, 1.46390e-02_rb, 1.27756e-02_rb, 1.11493e-02_rb/)
      selfrefo(:, 9) = (/ &
     & 3.68819e-02_rb, 3.22827e-02_rb, 2.82570e-02_rb, 2.47333e-02_rb, 2.16490e-02_rb, &
     & 1.89494e-02_rb, 1.65863e-02_rb, 1.45180e-02_rb, 1.27076e-02_rb, 1.11229e-02_rb/)
      selfrefo(:,10) = (/ &
     & 3.65157e-02_rb, 3.20121e-02_rb, 2.80640e-02_rb, 2.46028e-02_rb, 2.15685e-02_rb, &
     & 1.89084e-02_rb, 1.65764e-02_rb, 1.45320e-02_rb, 1.27397e-02_rb, 1.11685e-02_rb/)
      selfrefo(:,11) = (/ &
     & 3.59917e-02_rb, 3.16727e-02_rb, 2.78720e-02_rb, 2.45274e-02_rb, 2.15841e-02_rb, &
     & 1.89940e-02_rb, 1.67148e-02_rb, 1.47090e-02_rb, 1.29439e-02_rb, 1.13907e-02_rb/)
      selfrefo(:,12) = (/ &
     & 3.66963e-02_rb, 3.20483e-02_rb, 2.79890e-02_rb, 2.44439e-02_rb, 2.13478e-02_rb, &
     & 1.86438e-02_rb, 1.62824e-02_rb, 1.42200e-02_rb, 1.24189e-02_rb, 1.08459e-02_rb/)
      selfrefo(:,13) = (/ &
     & 3.66422e-02_rb, 3.19026e-02_rb, 2.77760e-02_rb, 2.41832e-02_rb, 2.10551e-02_rb, &
     & 1.83317e-02_rb, 1.59605e-02_rb, 1.38960e-02_rb, 1.20986e-02_rb, 1.05336e-02_rb/)
      selfrefo(:,14) = (/ &
     & 3.81260e-02_rb, 3.29322e-02_rb, 2.84460e-02_rb, 2.45709e-02_rb, 2.12237e-02_rb, &
     & 1.83325e-02_rb, 1.58352e-02_rb, 1.36780e-02_rb, 1.18147e-02_rb, 1.02052e-02_rb/)
      selfrefo(:,15) = (/ &
     & 3.51264e-02_rb, 3.05081e-02_rb, 2.64970e-02_rb, 2.30133e-02_rb, 1.99876e-02_rb, &
     & 1.73597e-02_rb, 1.50773e-02_rb, 1.30950e-02_rb, 1.13733e-02_rb, 9.87800e-03_rb/)
      selfrefo(:,16) = (/ &
     & 3.51264e-02_rb, 3.05081e-02_rb, 2.64970e-02_rb, 2.30133e-02_rb, 1.99876e-02_rb, &
     & 1.73597e-02_rb, 1.50773e-02_rb, 1.30950e-02_rb, 1.13733e-02_rb, 9.87800e-03_rb/)

      end subroutine lw_kgb079

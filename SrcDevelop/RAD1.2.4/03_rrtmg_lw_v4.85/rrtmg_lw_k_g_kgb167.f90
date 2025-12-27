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
!!prev      subroutine lw_kgb16
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb167
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg16, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg16, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &5.1629e-06_rb,7.7578e-06_rb,1.9043e-05_rb,1.4802e-04_rb,2.2980e-04_rb,2.8057e-04_rb, &
     &3.2824e-04_rb,3.4913e-04_rb,3.6515e-04_rb,3.8271e-04_rb,3.7499e-04_rb,3.6966e-04_rb, &
     &3.7424e-04_rb,3.8884e-04_rb,3.7117e-04_rb,4.3710e-04_rb/)
      forrefo(2,:) = (/ &
     &5.0804e-06_rb,1.3466e-05_rb,7.2606e-05_rb,1.6940e-04_rb,2.1022e-04_rb,2.5900e-04_rb, &
     &2.9106e-04_rb,3.2261e-04_rb,3.2066e-04_rb,3.5421e-04_rb,3.7128e-04_rb,3.8144e-04_rb, &
     &3.7854e-04_rb,3.8347e-04_rb,3.8921e-04_rb,3.7339e-04_rb/)
      forrefo(3,:) = (/ &
     &5.4797e-05_rb,1.0026e-04_rb,1.2422e-04_rb,1.6386e-04_rb,1.8378e-04_rb,1.9616e-04_rb, &
     &2.0711e-04_rb,2.2492e-04_rb,2.5240e-04_rb,2.6187e-04_rb,2.6058e-04_rb,2.4892e-04_rb, &
     &2.6526e-04_rb,3.2105e-04_rb,3.6903e-04_rb,3.7213e-04_rb/)
      forrefo(4,:) = (/ &
     &4.2782e-05_rb,1.4775e-04_rb,1.4588e-04_rb,1.6964e-04_rb,1.6667e-04_rb,1.7192e-04_rb, &
     &1.9057e-04_rb,2.0180e-04_rb,2.1177e-04_rb,2.2326e-04_rb,2.3801e-04_rb,2.9308e-04_rb, &
     &3.1130e-04_rb,3.1829e-04_rb,3.5035e-04_rb,3.7782e-04_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 1.27793e-03_rb, 1.05944e-03_rb, 8.78300e-04_rb, 7.28133e-04_rb, 6.03641e-04_rb, &
     & 5.00434e-04_rb, 4.14873e-04_rb, 3.43940e-04_rb, 2.85135e-04_rb, 2.36384e-04_rb/)
      selfrefo(:, 2) = (/ &
     & 1.42785e-03_rb, 1.17602e-03_rb, 9.68600e-04_rb, 7.97765e-04_rb, 6.57060e-04_rb, &
     & 5.41172e-04_rb, 4.45724e-04_rb, 3.67110e-04_rb, 3.02361e-04_rb, 2.49033e-04_rb/)
      selfrefo(:, 3) = (/ &
     & 2.94095e-03_rb, 2.27102e-03_rb, 1.75370e-03_rb, 1.35422e-03_rb, 1.04574e-03_rb, &
     & 8.07525e-04_rb, 6.23577e-04_rb, 4.81530e-04_rb, 3.71841e-04_rb, 2.87138e-04_rb/)
      selfrefo(:, 4) = (/ &
     & 3.94894e-03_rb, 3.48184e-03_rb, 3.07000e-03_rb, 2.70687e-03_rb, 2.38669e-03_rb, &
     & 2.10439e-03_rb, 1.85547e-03_rb, 1.63600e-03_rb, 1.44249e-03_rb, 1.27187e-03_rb/)
      selfrefo(:, 5) = (/ &
     & 4.19971e-03_rb, 3.86333e-03_rb, 3.55390e-03_rb, 3.26925e-03_rb, 3.00740e-03_rb, &
     & 2.76652e-03_rb, 2.54494e-03_rb, 2.34110e-03_rb, 2.15359e-03_rb, 1.98110e-03_rb/)
      selfrefo(:, 6) = (/ &
     & 4.95922e-03_rb, 4.57134e-03_rb, 4.21380e-03_rb, 3.88422e-03_rb, 3.58042e-03_rb, &
     & 3.30038e-03_rb, 3.04225e-03_rb, 2.80430e-03_rb, 2.58496e-03_rb, 2.38278e-03_rb/)
      selfrefo(:, 7) = (/ &
     & 5.27379e-03_rb, 4.91005e-03_rb, 4.57140e-03_rb, 4.25611e-03_rb, 3.96256e-03_rb, &
     & 3.68925e-03_rb, 3.43480e-03_rb, 3.19790e-03_rb, 2.97734e-03_rb, 2.77199e-03_rb/)
      selfrefo(:, 8) = (/ &
     & 5.75341e-03_rb, 5.31533e-03_rb, 4.91060e-03_rb, 4.53669e-03_rb, 4.19126e-03_rb, &
     & 3.87212e-03_rb, 3.57729e-03_rb, 3.30490e-03_rb, 3.05325e-03_rb, 2.82077e-03_rb/)
      selfrefo(:, 9) = (/ &
     & 5.49849e-03_rb, 5.14295e-03_rb, 4.81040e-03_rb, 4.49935e-03_rb, 4.20842e-03_rb, &
     & 3.93629e-03_rb, 3.68177e-03_rb, 3.44370e-03_rb, 3.22102e-03_rb, 3.01275e-03_rb/)
      selfrefo(:,10) = (/ &
     & 6.04962e-03_rb, 5.60945e-03_rb, 5.20130e-03_rb, 4.82285e-03_rb, 4.47194e-03_rb, &
     & 4.14656e-03_rb, 3.84485e-03_rb, 3.56510e-03_rb, 3.30570e-03_rb, 3.06518e-03_rb/)
      selfrefo(:,11) = (/ &
     & 6.40108e-03_rb, 5.87551e-03_rb, 5.39310e-03_rb, 4.95029e-03_rb, 4.54385e-03_rb, &
     & 4.17077e-03_rb, 3.82833e-03_rb, 3.51400e-03_rb, 3.22548e-03_rb, 2.96065e-03_rb/)
      selfrefo(:,12) = (/ &
     & 6.77938e-03_rb, 6.15713e-03_rb, 5.59200e-03_rb, 5.07874e-03_rb, 4.61259e-03_rb, &
     & 4.18922e-03_rb, 3.80472e-03_rb, 3.45550e-03_rb, 3.13834e-03_rb, 2.85029e-03_rb/)
      selfrefo(:,13) = (/ &
     & 6.90020e-03_rb, 6.26766e-03_rb, 5.69310e-03_rb, 5.17121e-03_rb, 4.69717e-03_rb, &
     & 4.26658e-03_rb, 3.87546e-03_rb, 3.52020e-03_rb, 3.19750e-03_rb, 2.90439e-03_rb/)
      selfrefo(:,14) = (/ &
     & 6.92759e-03_rb, 6.32882e-03_rb, 5.78180e-03_rb, 5.28206e-03_rb, 4.82552e-03_rb, &
     & 4.40843e-03_rb, 4.02740e-03_rb, 3.67930e-03_rb, 3.36129e-03_rb, 3.07076e-03_rb/)
      selfrefo(:,15) = (/ &
     & 7.54539e-03_rb, 6.81161e-03_rb, 6.14920e-03_rb, 5.55120e-03_rb, 5.01136e-03_rb, &
     & 4.52402e-03_rb, 4.08407e-03_rb, 3.68690e-03_rb, 3.32836e-03_rb, 3.00468e-03_rb/)
      selfrefo(:,16) = (/ &
     & 7.62039e-03_rb, 7.10834e-03_rb, 6.63070e-03_rb, 6.18515e-03_rb, 5.76955e-03_rb, &
     & 5.38186e-03_rb, 5.02023e-03_rb, 4.68290e-03_rb, 4.36823e-03_rb, 4.07471e-03_rb/)

      end subroutine lw_kgb167

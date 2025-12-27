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
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb109
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg10, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg10, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.0515e-02_rb,1.4860e-02_rb,1.7181e-02_rb,1.6642e-02_rb,1.6644e-02_rb,1.5649e-02_rb, &
     &1.7734e-02_rb,1.7521e-02_rb,1.7868e-02_rb,1.8400e-02_rb,1.9361e-02_rb,2.1487e-02_rb, &
     &2.0192e-02_rb,1.6545e-02_rb,2.0922e-02_rb,2.0922e-02_rb/)
      forrefo(2,:) = (/ &
     &1.0423e-02_rb,1.4593e-02_rb,1.6329e-02_rb,1.7071e-02_rb,1.7252e-02_rb,1.6188e-02_rb, &
     &1.7752e-02_rb,1.7913e-02_rb,1.7551e-02_rb,1.8203e-02_rb,1.7946e-02_rb,1.9828e-02_rb, &
     &2.1566e-02_rb,1.9707e-02_rb,2.0944e-02_rb,2.0944e-02_rb/)
      forrefo(3,:) = (/ &
     &9.2770e-03_rb,1.2818e-02_rb,1.7181e-02_rb,1.7858e-02_rb,1.7888e-02_rb,1.7121e-02_rb, &
     &1.8116e-02_rb,1.8230e-02_rb,1.7719e-02_rb,1.7833e-02_rb,1.8438e-02_rb,1.7995e-02_rb, &
     &2.0895e-02_rb,2.1525e-02_rb,2.0517e-02_rb,2.0954e-02_rb/)
      forrefo(4,:) = (/ &
     &8.3290e-03_rb,1.3483e-02_rb,1.5432e-02_rb,2.0793e-02_rb,1.8404e-02_rb,1.7470e-02_rb, &
     &1.7253e-02_rb,1.7132e-02_rb,1.7119e-02_rb,1.7376e-02_rb,1.7030e-02_rb,1.6847e-02_rb, &
     &1.5562e-02_rb,1.6836e-02_rb,1.8746e-02_rb,2.1233e-02_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 2.41120e-01_rb, 2.27071e-01_rb, 2.13840e-01_rb, 2.01380e-01_rb, 1.89646e-01_rb, &
     & 1.78596e-01_rb, 1.68190e-01_rb, 1.58390e-01_rb, 1.49161e-01_rb, 1.40470e-01_rb/)
      selfrefo(:, 2) = (/ &
     & 3.11156e-01_rb, 2.92249e-01_rb, 2.74490e-01_rb, 2.57810e-01_rb, 2.42144e-01_rb, &
     & 2.27430e-01_rb, 2.13610e-01_rb, 2.00630e-01_rb, 1.88439e-01_rb, 1.76988e-01_rb/)
      selfrefo(:, 3) = (/ &
     & 3.37148e-01_rb, 3.17767e-01_rb, 2.99500e-01_rb, 2.82283e-01_rb, 2.66056e-01_rb, &
     & 2.50762e-01_rb, 2.36347e-01_rb, 2.22760e-01_rb, 2.09955e-01_rb, 1.97885e-01_rb/)
      selfrefo(:, 4) = (/ &
     & 3.57139e-01_rb, 3.32763e-01_rb, 3.10050e-01_rb, 2.88888e-01_rb, 2.69170e-01_rb, &
     & 2.50798e-01_rb, 2.33680e-01_rb, 2.17730e-01_rb, 2.02869e-01_rb, 1.89022e-01_rb/)
      selfrefo(:, 5) = (/ &
     & 3.60626e-01_rb, 3.35433e-01_rb, 3.12000e-01_rb, 2.90204e-01_rb, 2.69931e-01_rb, &
     & 2.51074e-01_rb, 2.33534e-01_rb, 2.17220e-01_rb, 2.02045e-01_rb, 1.87931e-01_rb/)
      selfrefo(:, 6) = (/ &
     & 3.42420e-01_rb, 3.18795e-01_rb, 2.96800e-01_rb, 2.76323e-01_rb, 2.57258e-01_rb, &
     & 2.39509e-01_rb, 2.22985e-01_rb, 2.07600e-01_rb, 1.93277e-01_rb, 1.79942e-01_rb/)
      selfrefo(:, 7) = (/ &
     & 3.65491e-01_rb, 3.41599e-01_rb, 3.19270e-01_rb, 2.98400e-01_rb, 2.78895e-01_rb, &
     & 2.60664e-01_rb, 2.43625e-01_rb, 2.27700e-01_rb, 2.12816e-01_rb, 1.98905e-01_rb/)
      selfrefo(:, 8) = (/ &
     & 3.70354e-01_rb, 3.45005e-01_rb, 3.21390e-01_rb, 2.99392e-01_rb, 2.78899e-01_rb, &
     & 2.59809e-01_rb, 2.42026e-01_rb, 2.25460e-01_rb, 2.10028e-01_rb, 1.95652e-01_rb/)
      selfrefo(:, 9) = (/ &
     & 3.60483e-01_rb, 3.37846e-01_rb, 3.16630e-01_rb, 2.96747e-01_rb, 2.78112e-01_rb, &
     & 2.60648e-01_rb, 2.44280e-01_rb, 2.28940e-01_rb, 2.14563e-01_rb, 2.01090e-01_rb/)
      selfrefo(:,10) = (/ &
     & 3.71845e-01_rb, 3.48164e-01_rb, 3.25990e-01_rb, 3.05229e-01_rb, 2.85790e-01_rb, &
     & 2.67588e-01_rb, 2.50547e-01_rb, 2.34590e-01_rb, 2.19650e-01_rb, 2.05661e-01_rb/)
      selfrefo(:,11) = (/ &
     & 3.60606e-01_rb, 3.40789e-01_rb, 3.22060e-01_rb, 3.04361e-01_rb, 2.87634e-01_rb, &
     & 2.71826e-01_rb, 2.56888e-01_rb, 2.42770e-01_rb, 2.29428e-01_rb, 2.16819e-01_rb/)
      selfrefo(:,12) = (/ &
     & 3.90046e-01_rb, 3.68879e-01_rb, 3.48860e-01_rb, 3.29928e-01_rb, 3.12023e-01_rb, &
     & 2.95089e-01_rb, 2.79075e-01_rb, 2.63930e-01_rb, 2.49607e-01_rb, 2.36061e-01_rb/)
      selfrefo(:,13) = (/ &
     & 4.38542e-01_rb, 4.05139e-01_rb, 3.74280e-01_rb, 3.45771e-01_rb, 3.19434e-01_rb, &
     & 2.95103e-01_rb, 2.72626e-01_rb, 2.51860e-01_rb, 2.32676e-01_rb, 2.14953e-01_rb/)
      selfrefo(:,14) = (/ &
     & 4.19448e-01_rb, 3.81920e-01_rb, 3.47750e-01_rb, 3.16637e-01_rb, 2.88307e-01_rb, &
     & 2.62513e-01_rb, 2.39026e-01_rb, 2.17640e-01_rb, 1.98168e-01_rb, 1.80438e-01_rb/)
      selfrefo(:,15) = (/ &
     & 4.20276e-01_rb, 3.92281e-01_rb, 3.66150e-01_rb, 3.41760e-01_rb, 3.18995e-01_rb, &
     & 2.97746e-01_rb, 2.77912e-01_rb, 2.59400e-01_rb, 2.42121e-01_rb, 2.25993e-01_rb/)
      selfrefo(:,16) = (/ &
     & 4.20276e-01_rb, 3.92281e-01_rb, 3.66150e-01_rb, 3.41760e-01_rb, 3.18995e-01_rb, &
     & 2.97746e-01_rb, 2.77912e-01_rb, 2.59400e-01_rb, 2.42121e-01_rb, 2.25993e-01_rb/)

      end subroutine lw_kgb109

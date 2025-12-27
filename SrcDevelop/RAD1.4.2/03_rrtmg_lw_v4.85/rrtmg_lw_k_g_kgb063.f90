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
!!prev      subroutine lw_kgb06
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb063
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg06, only : fracrefao, kao, kao_mco2, selfrefo, forrefo, &
!!prev                            cfc11adjo, cfc12o
      use rrlw_kg06, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &3.2710e-07_rb,5.2119e-07_rb,8.4740e-07_rb,1.6908e-06_rb,2.3433e-06_rb,4.4129e-06_rb, &
     &3.8930e-06_rb,2.3338e-06_rb,2.4115e-06_rb,2.4271e-06_rb,2.4836e-06_rb,2.6470e-06_rb, &
     &2.9559e-06_rb,2.3940e-06_rb,2.9711e-06_rb,2.9511e-06_rb/)
      forrefo(2,:) = (/ &
     &6.5125e-07_rb,1.2128e-06_rb,1.7249e-06_rb,2.7126e-06_rb,3.1780e-06_rb,2.1444e-06_rb, &
     &1.8265e-06_rb,1.7385e-06_rb,1.4574e-06_rb,1.6135e-06_rb,2.4966e-06_rb,2.8127e-06_rb, &
     &2.5229e-06_rb,2.3251e-06_rb,2.5353e-06_rb,3.0200e-06_rb/)
      forrefo(3,:) = (/ &
     &1.4969e-06_rb,1.8516e-06_rb,2.5791e-06_rb,2.7846e-06_rb,1.9789e-06_rb,1.6688e-06_rb, &
     &1.1037e-06_rb,9.9065e-07_rb,1.1557e-06_rb,7.0847e-07_rb,5.7758e-07_rb,4.0425e-07_rb, &
     &3.2427e-07_rb,3.2267e-07_rb,3.1444e-07_rb,2.6046e-07_rb/)
      forrefo(4,:) = (/ &
     &1.7567e-06_rb,1.6891e-06_rb,2.1003e-06_rb,2.0957e-06_rb,2.3664e-06_rb,2.1538e-06_rb, &
     &1.5275e-06_rb,1.0487e-06_rb,8.7390e-07_rb,7.9360e-07_rb,7.7778e-07_rb,8.1445e-07_rb, &
     &8.2121e-07_rb,5.4395e-07_rb,3.1273e-07_rb,3.1848e-07_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 7.73921e-02_rb, 6.45225e-02_rb, 5.37930e-02_rb, 4.48477e-02_rb, 3.73900e-02_rb, &
     & 3.11723e-02_rb, 2.59887e-02_rb, 2.16670e-02_rb, 1.80640e-02_rb, 1.50601e-02_rb/)
      selfrefo(:, 2) = (/ &
     & 8.47756e-02_rb, 7.10616e-02_rb, 5.95660e-02_rb, 4.99301e-02_rb, 4.18529e-02_rb, &
     & 3.50824e-02_rb, 2.94072e-02_rb, 2.46500e-02_rb, 2.06624e-02_rb, 1.73199e-02_rb/)
      selfrefo(:, 3) = (/ &
     & 8.84829e-02_rb, 7.46093e-02_rb, 6.29110e-02_rb, 5.30469e-02_rb, 4.47295e-02_rb, &
     & 3.77161e-02_rb, 3.18025e-02_rb, 2.68160e-02_rb, 2.26114e-02_rb, 1.90661e-02_rb/)
      selfrefo(:, 4) = (/ &
     & 9.27003e-02_rb, 7.88864e-02_rb, 6.71310e-02_rb, 5.71273e-02_rb, 4.86144e-02_rb, &
     & 4.13700e-02_rb, 3.52052e-02_rb, 2.99590e-02_rb, 2.54946e-02_rb, 2.16955e-02_rb/)
      selfrefo(:, 5) = (/ &
     & 9.14315e-02_rb, 7.85661e-02_rb, 6.75110e-02_rb, 5.80115e-02_rb, 4.98487e-02_rb, &
     & 4.28344e-02_rb, 3.68072e-02_rb, 3.16280e-02_rb, 2.71776e-02_rb, 2.33534e-02_rb/)
      selfrefo(:, 6) = (/ &
     & 7.72984e-02_rb, 6.91044e-02_rb, 6.17790e-02_rb, 5.52301e-02_rb, 4.93755e-02_rb, &
     & 4.41414e-02_rb, 3.94622e-02_rb, 3.52790e-02_rb, 3.15392e-02_rb, 2.81959e-02_rb/)
      selfrefo(:, 7) = (/ &
     & 7.46998e-02_rb, 6.66597e-02_rb, 5.94850e-02_rb, 5.30825e-02_rb, 4.73691e-02_rb, &
     & 4.22707e-02_rb, 3.77210e-02_rb, 3.36610e-02_rb, 3.00380e-02_rb, 2.68049e-02_rb/)
      selfrefo(:, 8) = (/ &
     & 7.59386e-02_rb, 6.66263e-02_rb, 5.84560e-02_rb, 5.12876e-02_rb, 4.49982e-02_rb, &
     & 3.94801e-02_rb, 3.46387e-02_rb, 3.03910e-02_rb, 2.66642e-02_rb, 2.33944e-02_rb/)
      selfrefo(:, 9) = (/ &
     & 7.26921e-02_rb, 6.43261e-02_rb, 5.69230e-02_rb, 5.03719e-02_rb, 4.45747e-02_rb, &
     & 3.94447e-02_rb, 3.49051e-02_rb, 3.08880e-02_rb, 2.73332e-02_rb, 2.41875e-02_rb/)
      selfrefo(:,10) = (/ &
     & 7.43684e-02_rb, 6.58735e-02_rb, 5.83490e-02_rb, 5.16840e-02_rb, 4.57803e-02_rb, &
     & 4.05509e-02_rb, 3.59189e-02_rb, 3.18160e-02_rb, 2.81818e-02_rb, 2.49626e-02_rb/)
      selfrefo(:,11) = (/ &
     & 8.97599e-02_rb, 7.73727e-02_rb, 6.66950e-02_rb, 5.74908e-02_rb, 4.95569e-02_rb, &
     & 4.27179e-02_rb, 3.68227e-02_rb, 3.17410e-02_rb, 2.73606e-02_rb, 2.35848e-02_rb/)
      selfrefo(:,12) = (/ &
     & 9.12262e-02_rb, 7.84848e-02_rb, 6.75230e-02_rb, 5.80922e-02_rb, 4.99786e-02_rb, &
     & 4.29982e-02_rb, 3.69927e-02_rb, 3.18260e-02_rb, 2.73809e-02_rb, 2.35567e-02_rb/)
      selfrefo(:,13) = (/ &
     & 9.03254e-02_rb, 7.83291e-02_rb, 6.79260e-02_rb, 5.89046e-02_rb, 5.10813e-02_rb, &
     & 4.42970e-02_rb, 3.84139e-02_rb, 3.33120e-02_rb, 2.88877e-02_rb, 2.50511e-02_rb/)
      selfrefo(:,14) = (/ &
     & 9.22803e-02_rb, 7.94172e-02_rb, 6.83470e-02_rb, 5.88199e-02_rb, 5.06209e-02_rb, &
     & 4.35647e-02_rb, 3.74921e-02_rb, 3.22660e-02_rb, 2.77684e-02_rb, 2.38977e-02_rb/)
      selfrefo(:,15) = (/ &
     & 9.36819e-02_rb, 8.10810e-02_rb, 7.01750e-02_rb, 6.07359e-02_rb, 5.25665e-02_rb, &
     & 4.54959e-02_rb, 3.93764e-02_rb, 3.40800e-02_rb, 2.94960e-02_rb, 2.55286e-02_rb/)
      selfrefo(:,16) = (/ &
     & 1.00195e-01_rb, 8.58713e-02_rb, 7.35950e-02_rb, 6.30737e-02_rb, 5.40566e-02_rb, &
     & 4.63286e-02_rb, 3.97054e-02_rb, 3.40290e-02_rb, 2.91641e-02_rb, 2.49948e-02_rb/)

      end subroutine lw_kgb063

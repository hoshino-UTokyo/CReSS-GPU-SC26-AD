!     path:      $Source:
!     /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 22:22:22 $

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
!      subroutine sw_kgbnn
! **************************************************************************
!  RRTM Shortwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original by J.Delamere, Atmospheric & Environmental Research.
!  Reformatted for F90: JJMorcrette, ECMWF
!  Further F90 and GCM revisions:  MJIacono, AER, July 2002
!
!  This file contains 14 subroutines that include the 
!  absorption coefficients and other data for each of the 14 shortwave
!  spectral bands used in RRTM_SW.  Here, the data are defined for 16
!  g-points, or sub-intervals, per band.  These data are combined and
!  weighted using a mapping procedure in routine RRTMG_SW_INIT to reduce
!  the total number of g-points from 224 to 112 for use in the GCM.
! **************************************************************************
!!prev      subroutine sw_kgb24
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine sw_kgb247
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg24, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
!!prev                            raylao, raylbo, abso3ao, abso3bo
      use rrsw_kg24, only : selfrefo, forrefo

      implicit none
      save

!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.515619e-08_rb, 0.131078e-06_rb, 0.349038e-06_rb /)
      forrefo(:, 2) = (/ 0.329605e-07_rb, 0.430497e-06_rb, 0.458569e-05_rb /)
      forrefo(:, 3) = (/ 0.188244e-06_rb, 0.792931e-06_rb, 0.267176e-05_rb /)
      forrefo(:, 4) = (/ 0.611237e-06_rb, 0.798868e-06_rb, 0.411583e-06_rb /)
      forrefo(:, 5) = (/ 0.111903e-05_rb, 0.914895e-06_rb, 0.444828e-06_rb /)
      forrefo(:, 6) = (/ 0.235399e-05_rb, 0.269099e-05_rb, 0.739855e-06_rb /)
      forrefo(:, 7) = (/ 0.400131e-05_rb, 0.378135e-05_rb, 0.231265e-06_rb /)
      forrefo(:, 8) = (/ 0.464257e-05_rb, 0.371927e-05_rb, 0.460611e-06_rb /)
      forrefo(:, 9) = (/ 0.476792e-05_rb, 0.311841e-05_rb, 0.934811e-06_rb /)
      forrefo(:,10) = (/ 0.555683e-05_rb, 0.238129e-05_rb, 0.400334e-07_rb /)
      forrefo(:,11) = (/ 0.569068e-05_rb, 0.196039e-05_rb, 0.374476e-07_rb /)
      forrefo(:,12) = (/ 0.554154e-05_rb, 0.131724e-05_rb, 0.399720e-07_rb /)
      forrefo(:,13) = (/ 0.462684e-05_rb, 0.238826e-07_rb, 0.325793e-07_rb /)
      forrefo(:,14) = (/ 0.808644e-06_rb, 0.105126e-11_rb, 0.148691e-07_rb /)
      forrefo(:,15) = (/ 0.865024e-12_rb, 0.822434e-12_rb, 0.825756e-12_rb /)
      forrefo(:,16) = (/ 0.945747e-12_rb, 0.802065e-12_rb, 0.724732e-12_rb /)

!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
        &  0.637755e-05_rb, 0.403921e-05_rb, 0.255823e-05_rb, 0.162025e-05_rb, 0.102618e-05_rb, &
        &  0.649930e-06_rb, 0.411632e-06_rb, 0.260707e-06_rb, 0.165118e-06_rb, 0.104577e-06_rb /)
      selfrefo(:, 2) = (/ &
        &  0.180887e-04_rb, 0.108890e-04_rb, 0.655493e-05_rb, 0.394592e-05_rb, 0.237536e-05_rb, &
        &  0.142991e-05_rb, 0.860774e-06_rb, 0.518167e-06_rb, 0.311925e-06_rb, 0.187772e-06_rb /)
      selfrefo(:, 3) = (/ &
        &  0.212261e-04_rb, 0.150697e-04_rb, 0.106989e-04_rb, 0.759581e-05_rb, 0.539274e-05_rb, &
        &  0.382864e-05_rb, 0.271819e-05_rb, 0.192981e-05_rb, 0.137009e-05_rb, 0.972711e-06_rb /)
      selfrefo(:, 4) = (/ &
        &  0.132497e-04_rb, 0.118071e-04_rb, 0.105216e-04_rb, 0.937599e-05_rb, 0.835516e-05_rb, &
        &  0.744547e-05_rb, 0.663482e-05_rb, 0.591243e-05_rb, 0.526870e-05_rb, 0.469506e-05_rb /)
      selfrefo(:, 5) = (/ &
        &  0.124069e-04_rb, 0.120785e-04_rb, 0.117589e-04_rb, 0.114477e-04_rb, 0.111447e-04_rb, &
        &  0.108498e-04_rb, 0.105626e-04_rb, 0.102831e-04_rb, 0.100109e-04_rb, 0.974601e-05_rb /)
      selfrefo(:, 6) = (/ &
        &  0.411994e-04_rb, 0.372560e-04_rb, 0.336901e-04_rb, 0.304654e-04_rb, 0.275494e-04_rb, &
        &  0.249126e-04_rb, 0.225281e-04_rb, 0.203718e-04_rb, 0.184219e-04_rb, 0.166587e-04_rb /)
      selfrefo(:, 7) = (/ &
        &  0.537376e-04_rb, 0.501002e-04_rb, 0.467090e-04_rb, 0.435473e-04_rb, 0.405996e-04_rb, &
        &  0.378515e-04_rb, 0.352893e-04_rb, 0.329006e-04_rb, 0.306736e-04_rb, 0.285974e-04_rb /)
      selfrefo(:, 8) = (/ &
        &  0.494279e-04_rb, 0.475365e-04_rb, 0.457175e-04_rb, 0.439681e-04_rb, 0.422857e-04_rb, &
        &  0.406676e-04_rb, 0.391114e-04_rb, 0.376148e-04_rb, 0.361755e-04_rb, 0.347912e-04_rb /)
      selfrefo(:, 9) = (/ &
        &  0.377444e-04_rb, 0.378199e-04_rb, 0.378956e-04_rb, 0.379715e-04_rb, 0.380475e-04_rb, &
        &  0.381236e-04_rb, 0.381999e-04_rb, 0.382763e-04_rb, 0.383529e-04_rb, 0.384297e-04_rb /)
      selfrefo(:,10) = (/ &
        &  0.245916e-04_rb, 0.267183e-04_rb, 0.290289e-04_rb, 0.315394e-04_rb, 0.342669e-04_rb, &
        &  0.372304e-04_rb, 0.404501e-04_rb, 0.439483e-04_rb, 0.477490e-04_rb, 0.518784e-04_rb /)
      selfrefo(:,11) = (/ &
        &  0.186528e-04_rb, 0.211417e-04_rb, 0.239628e-04_rb, 0.271603e-04_rb, 0.307845e-04_rb, &
        &  0.348923e-04_rb, 0.395482e-04_rb, 0.448254e-04_rb, 0.508068e-04_rb, 0.575863e-04_rb /)
      selfrefo(:,12) = (/ &
        &  0.109896e-04_rb, 0.133794e-04_rb, 0.162890e-04_rb, 0.198312e-04_rb, 0.241438e-04_rb, &
        &  0.293942e-04_rb, 0.357864e-04_rb, 0.435686e-04_rb, 0.530432e-04_rb, 0.645781e-04_rb /)
      selfrefo(:,13) = (/ &
        &  0.183885e-06_rb, 0.391019e-06_rb, 0.831472e-06_rb, 0.176806e-05_rb, 0.375966e-05_rb, &
        &  0.799463e-05_rb, 0.170000e-04_rb, 0.361492e-04_rb, 0.768686e-04_rb, 0.163455e-03_rb /)
      selfrefo(:,14) = (/ &
        &  0.466057e-07_rb, 0.937419e-07_rb, 0.188551e-06_rb, 0.379248e-06_rb, 0.762813e-06_rb, &
        &  0.153431e-05_rb, 0.308608e-05_rb, 0.620729e-05_rb, 0.124852e-04_rb, 0.251126e-04_rb /)
      selfrefo(:,15) = (/ &
        &  0.248961e-06_rb, 0.216780e-06_rb, 0.188758e-06_rb, 0.164358e-06_rb, 0.143113e-06_rb, &
        &  0.124613e-06_rb, 0.108505e-06_rb, 0.944795e-07_rb, 0.822667e-07_rb, 0.716326e-07_rb /)
      selfrefo(:,16) = (/ &
        &  0.252246e-06_rb, 0.220335e-06_rb, 0.192462e-06_rb, 0.168114e-06_rb, 0.146847e-06_rb, &
        &  0.128270e-06_rb, 0.112043e-06_rb, 0.978688e-07_rb, 0.854878e-07_rb, 0.746731e-07_rb /)
     
      end subroutine sw_kgb247

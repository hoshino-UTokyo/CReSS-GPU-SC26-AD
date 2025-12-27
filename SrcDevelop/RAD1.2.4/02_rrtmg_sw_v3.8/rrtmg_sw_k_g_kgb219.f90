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
!!prev      subroutine sw_kgb21
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine sw_kgb219
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
      use rrsw_kg21, only : selfrefo, forrefo

      implicit none
      save

!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.110008e-06_rb, 0.630912e-06_rb, 0.363159e-05_rb, 0.616892e-05_rb /)
      forrefo(:, 2) = (/ 0.429709e-05_rb, 0.789174e-05_rb, 0.217416e-04_rb, 0.639393e-04_rb /)
      forrefo(:, 3) = (/ 0.436283e-04_rb, 0.526247e-04_rb, 0.116341e-03_rb, 0.205616e-03_rb /)
      forrefo(:, 4) = (/ 0.215627e-03_rb, 0.234522e-03_rb, 0.280497e-03_rb, 0.838668e-03_rb /)
      forrefo(:, 5) = (/ 0.529283e-03_rb, 0.620848e-03_rb, 0.935561e-03_rb, 0.171252e-02_rb /)
      forrefo(:, 6) = (/ 0.212267e-02_rb, 0.218564e-02_rb, 0.222227e-02_rb, 0.199650e-02_rb /)
      forrefo(:, 7) = (/ 0.291120e-02_rb, 0.281168e-02_rb, 0.259543e-02_rb, 0.210159e-02_rb /)
      forrefo(:, 8) = (/ 0.316249e-02_rb, 0.310695e-02_rb, 0.279501e-02_rb, 0.208076e-02_rb /)
      forrefo(:, 9) = (/ 0.354993e-02_rb, 0.336989e-02_rb, 0.298930e-02_rb, 0.180424e-02_rb /)
      forrefo(:,10) = (/ 0.397729e-02_rb, 0.367409e-02_rb, 0.328982e-02_rb, 0.177807e-02_rb /)
      forrefo(:,11) = (/ 0.408831e-02_rb, 0.398792e-02_rb, 0.352727e-02_rb, 0.192470e-02_rb /)
      forrefo(:,12) = (/ 0.433926e-02_rb, 0.420667e-02_rb, 0.383894e-02_rb, 0.220836e-02_rb /)
      forrefo(:,13) = (/ 0.436397e-02_rb, 0.433769e-02_rb, 0.425752e-02_rb, 0.237343e-02_rb /)
      forrefo(:,14) = (/ 0.440525e-02_rb, 0.449018e-02_rb, 0.451881e-02_rb, 0.269169e-02_rb /)
      forrefo(:,15) = (/ 0.491350e-02_rb, 0.481760e-02_rb, 0.475799e-02_rb, 0.362666e-02_rb /)
      forrefo(:,16) = (/ 0.561641e-02_rb, 0.524553e-02_rb, 0.512473e-02_rb, 0.493802e-02_rb /)

!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
        &  0.115887e-03_rb, 0.926537e-04_rb, 0.740783e-04_rb, 0.592270e-04_rb, 0.473530e-04_rb, &
        &  0.378596e-04_rb, 0.302694e-04_rb, 0.242010e-04_rb, 0.193491e-04_rb, 0.154700e-04_rb /)
      selfrefo(:, 2) = (/ &
        &  0.459557e-03_rb, 0.381962e-03_rb, 0.317469e-03_rb, 0.263866e-03_rb, 0.219313e-03_rb, &
        &  0.182283e-03_rb, 0.151505e-03_rb, 0.125924e-03_rb, 0.104662e-03_rb, 0.869904e-04_rb /)
      selfrefo(:, 3) = (/ &
        &  0.166821e-02_rb, 0.151103e-02_rb, 0.136866e-02_rb, 0.123970e-02_rb, 0.112290e-02_rb, &
        &  0.101710e-02_rb, 0.921266e-03_rb, 0.834463e-03_rb, 0.755839e-03_rb, 0.684623e-03_rb /)
      selfrefo(:, 4) = (/ &
        &  0.460175e-02_rb, 0.421372e-02_rb, 0.385842e-02_rb, 0.353307e-02_rb, 0.323516e-02_rb, &
        &  0.296236e-02_rb, 0.271257e-02_rb, 0.248385e-02_rb, 0.227440e-02_rb, 0.208262e-02_rb /)
      selfrefo(:, 5) = (/ &
        &  0.101589e-01_rb, 0.924742e-02_rb, 0.841772e-02_rb, 0.766247e-02_rb, 0.697497e-02_rb, &
        &  0.634917e-02_rb, 0.577951e-02_rb, 0.526096e-02_rb, 0.478893e-02_rb, 0.435926e-02_rb /)
      selfrefo(:, 6) = (/ &
        &  0.328043e-01_rb, 0.300853e-01_rb, 0.275917e-01_rb, 0.253048e-01_rb, 0.232075e-01_rb, &
        &  0.212839e-01_rb, 0.195198e-01_rb, 0.179020e-01_rb, 0.164182e-01_rb, 0.150574e-01_rb /)
      selfrefo(:, 7) = (/ &
        &  0.405936e-01_rb, 0.376032e-01_rb, 0.348331e-01_rb, 0.322671e-01_rb, 0.298901e-01_rb, &
        &  0.276883e-01_rb, 0.256486e-01_rb, 0.237591e-01_rb, 0.220089e-01_rb, 0.203876e-01_rb /)
      selfrefo(:, 8) = (/ &
        &  0.448362e-01_rb, 0.413811e-01_rb, 0.381923e-01_rb, 0.352492e-01_rb, 0.325329e-01_rb, &
        &  0.300259e-01_rb, 0.277121e-01_rb, 0.255766e-01_rb, 0.236056e-01_rb, 0.217866e-01_rb /)
      selfrefo(:, 9) = (/ &
        &  0.479741e-01_rb, 0.445389e-01_rb, 0.413497e-01_rb, 0.383889e-01_rb, 0.356400e-01_rb, &
        &  0.330880e-01_rb, 0.307188e-01_rb, 0.285191e-01_rb, 0.264770e-01_rb, 0.245812e-01_rb /)
      selfrefo(:,10) = (/ &
        &  0.519308e-01_rb, 0.484130e-01_rb, 0.451335e-01_rb, 0.420761e-01_rb, 0.392259e-01_rb, &
        &  0.365687e-01_rb, 0.340916e-01_rb, 0.317822e-01_rb, 0.296293e-01_rb, 0.276222e-01_rb /)
      selfrefo(:,11) = (/ &
        &  0.572039e-01_rb, 0.527780e-01_rb, 0.486945e-01_rb, 0.449270e-01_rb, 0.414510e-01_rb, &
        &  0.382439e-01_rb, 0.352849e-01_rb, 0.325549e-01_rb, 0.300361e-01_rb, 0.277122e-01_rb /)
      selfrefo(:,12) = (/ &
        &  0.601046e-01_rb, 0.554411e-01_rb, 0.511395e-01_rb, 0.471716e-01_rb, 0.435116e-01_rb, &
        &  0.401356e-01_rb, 0.370215e-01_rb, 0.341490e-01_rb, 0.314994e-01_rb, 0.290554e-01_rb /)
      selfrefo(:,13) = (/ &
        &  0.616595e-01_rb, 0.567145e-01_rb, 0.521662e-01_rb, 0.479826e-01_rb, 0.441346e-01_rb, &
        &  0.405951e-01_rb, 0.373395e-01_rb, 0.343450e-01_rb, 0.315906e-01_rb, 0.290571e-01_rb /)
      selfrefo(:,14) = (/ &
        &  0.647916e-01_rb, 0.592493e-01_rb, 0.541811e-01_rb, 0.495465e-01_rb, 0.453083e-01_rb, &
        &  0.414326e-01_rb, 0.378885e-01_rb, 0.346475e-01_rb, 0.316837e-01_rb, 0.289735e-01_rb /)
      selfrefo(:,15) = (/ &
        &  0.694231e-01_rb, 0.637703e-01_rb, 0.585777e-01_rb, 0.538079e-01_rb, 0.494265e-01_rb, &
        &  0.454019e-01_rb, 0.417050e-01_rb, 0.383091e-01_rb, 0.351897e-01_rb, 0.323244e-01_rb /)
      selfrefo(:,16) = (/ &
        &  0.761764e-01_rb, 0.701815e-01_rb, 0.646584e-01_rb, 0.595700e-01_rb, 0.548820e-01_rb, &
        &  0.505629e-01_rb, 0.465838e-01_rb, 0.429178e-01_rb, 0.395403e-01_rb, 0.364286e-01_rb /)
     
      end subroutine sw_kgb219

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
!!prev      subroutine sw_kgb17
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine sw_kgb179
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg17, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg17, only : selfrefo, forrefo

      implicit none
      save

!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.553258e-03_rb, 0.555486e-03_rb, 0.601339e-03_rb, 0.708280e-03_rb /)
      forrefo(:, 2) = (/ 0.158558e-02_rb, 0.162957e-02_rb, 0.204991e-02_rb, 0.475881e-02_rb /)
      forrefo(:, 3) = (/ 0.772542e-02_rb, 0.784562e-02_rb, 0.111979e-01_rb, 0.229016e-01_rb /)
      forrefo(:, 4) = (/ 0.255097e-01_rb, 0.256272e-01_rb, 0.270691e-01_rb, 0.259505e-01_rb /)
      forrefo(:, 5) = (/ 0.323263e-01_rb, 0.324495e-01_rb, 0.305535e-01_rb, 0.263993e-01_rb /)
      forrefo(:, 6) = (/ 0.346920e-01_rb, 0.348255e-01_rb, 0.323586e-01_rb, 0.276357e-01_rb /)
      forrefo(:, 7) = (/ 0.366509e-01_rb, 0.366412e-01_rb, 0.344434e-01_rb, 0.319223e-01_rb /)
      forrefo(:, 8) = (/ 0.378451e-01_rb, 0.375341e-01_rb, 0.374369e-01_rb, 0.320334e-01_rb /)
      forrefo(:, 9) = (/ 0.407348e-01_rb, 0.396203e-01_rb, 0.393988e-01_rb, 0.318343e-01_rb /)
      forrefo(:,10) = (/ 0.433035e-01_rb, 0.426488e-01_rb, 0.408085e-01_rb, 0.332749e-01_rb /)
      forrefo(:,11) = (/ 0.428254e-01_rb, 0.441151e-01_rb, 0.408887e-01_rb, 0.327077e-01_rb /)
      forrefo(:,12) = (/ 0.443226e-01_rb, 0.446690e-01_rb, 0.404676e-01_rb, 0.350492e-01_rb /)
      forrefo(:,13) = (/ 0.466103e-01_rb, 0.460809e-01_rb, 0.401286e-01_rb, 0.370427e-01_rb /)
      forrefo(:,14) = (/ 0.483928e-01_rb, 0.477284e-01_rb, 0.380684e-01_rb, 0.387940e-01_rb /)
      forrefo(:,15) = (/ 0.506987e-01_rb, 0.490016e-01_rb, 0.467069e-01_rb, 0.368998e-01_rb /)
      forrefo(:,16) = (/ 0.510836e-01_rb, 0.522771e-01_rb, 0.500130e-01_rb, 0.483406e-01_rb /)

!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
        &  0.160537e-01_rb, 0.149038e-01_rb, 0.138363e-01_rb, 0.128452e-01_rb, 0.119251e-01_rb, &
        &  0.110709e-01_rb, 0.102779e-01_rb, 0.954175e-02_rb, 0.885829e-02_rb, 0.822379e-02_rb /)
      selfrefo(:, 2) = (/ &
        &  0.365753e-01_rb, 0.342267e-01_rb, 0.320288e-01_rb, 0.299720e-01_rb, 0.280474e-01_rb, &
        &  0.262463e-01_rb, 0.245609e-01_rb, 0.229837e-01_rb, 0.215078e-01_rb, 0.201267e-01_rb /)
      selfrefo(:, 3) = (/ &
        &  0.127419e+00_rb, 0.118553e+00_rb, 0.110304e+00_rb, 0.102629e+00_rb, 0.954883e-01_rb, &
        &  0.888442e-01_rb, 0.826624e-01_rb, 0.769107e-01_rb, 0.715593e-01_rb, 0.665802e-01_rb /)
      selfrefo(:, 4) = (/ &
        &  0.378687e+00_rb, 0.348961e+00_rb, 0.321568e+00_rb, 0.296325e+00_rb, 0.273064e+00_rb, &
        &  0.251629e+00_rb, 0.231876e+00_rb, 0.213674e+00_rb, 0.196901e+00_rb, 0.181444e+00_rb /)
      selfrefo(:, 5) = (/ &
        &  0.472822e+00_rb, 0.435018e+00_rb, 0.400236e+00_rb, 0.368236e+00_rb, 0.338794e+00_rb, &
        &  0.311706e+00_rb, 0.286783e+00_rb, 0.263854e+00_rb, 0.242757e+00_rb, 0.223348e+00_rb /)
      selfrefo(:, 6) = (/ &
        &  0.505620e+00_rb, 0.465050e+00_rb, 0.427736e+00_rb, 0.393416e+00_rb, 0.361849e+00_rb, &
        &  0.332815e+00_rb, 0.306111e+00_rb, 0.281550e+00_rb, 0.258959e+00_rb, 0.238181e+00_rb /)
      selfrefo(:, 7) = (/ &
        &  0.530488e+00_rb, 0.487993e+00_rb, 0.448902e+00_rb, 0.412943e+00_rb, 0.379864e+00_rb, &
        &  0.349434e+00_rb, 0.321443e+00_rb, 0.295694e+00_rb, 0.272007e+00_rb, 0.250218e+00_rb /)
      selfrefo(:, 8) = (/ &
        &  0.540222e+00_rb, 0.497746e+00_rb, 0.458610e+00_rb, 0.422551e+00_rb, 0.389327e+00_rb, &
        &  0.358716e+00_rb, 0.330511e+00_rb, 0.304524e+00_rb, 0.280580e+00_rb, 0.258519e+00_rb /)
      selfrefo(:, 9) = (/ &
        &  0.565727e+00_rb, 0.522899e+00_rb, 0.483313e+00_rb, 0.446724e+00_rb, 0.412905e+00_rb, &
        &  0.381646e+00_rb, 0.352753e+00_rb, 0.326048e+00_rb, 0.301365e+00_rb, 0.278550e+00_rb /)
      selfrefo(:,10) = (/ &
        &  0.610122e+00_rb, 0.562337e+00_rb, 0.518295e+00_rb, 0.477702e+00_rb, 0.440289e+00_rb, &
        &  0.405806e+00_rb, 0.374023e+00_rb, 0.344730e+00_rb, 0.317730e+00_rb, 0.292846e+00_rb /)
      selfrefo(:,11) = (/ &
        &  0.645176e+00_rb, 0.588957e+00_rb, 0.537636e+00_rb, 0.490788e+00_rb, 0.448022e+00_rb, &
        &  0.408982e+00_rb, 0.373344e+00_rb, 0.340812e+00_rb, 0.311114e+00_rb, 0.284004e+00_rb /)
      selfrefo(:,12) = (/ &
        &  0.651737e+00_rb, 0.596547e+00_rb, 0.546031e+00_rb, 0.499792e+00_rb, 0.457469e+00_rb, &
        &  0.418730e+00_rb, 0.383272e+00_rb, 0.350816e+00_rb, 0.321108e+00_rb, 0.293916e+00_rb /)
      selfrefo(:,13) = (/ &
        &  0.661086e+00_rb, 0.607954e+00_rb, 0.559093e+00_rb, 0.514159e+00_rb, 0.472836e+00_rb, &
        &  0.434834e+00_rb, 0.399886e+00_rb, 0.367747e+00_rb, 0.338191e+00_rb, 0.311011e+00_rb /)
      selfrefo(:,14) = (/ &
        &  0.692554e+00_rb, 0.635574e+00_rb, 0.583282e+00_rb, 0.535293e+00_rb, 0.491251e+00_rb, &
        &  0.450834e+00_rb, 0.413741e+00_rb, 0.379701e+00_rb, 0.348461e+00_rb, 0.319791e+00_rb /)
      selfrefo(:,15) = (/ &
        &  0.714646e+00_rb, 0.657179e+00_rb, 0.604334e+00_rb, 0.555737e+00_rb, 0.511049e+00_rb, &
        &  0.469954e+00_rb, 0.432164e+00_rb, 0.397412e+00_rb, 0.365455e+00_rb, 0.336068e+00_rb /)
      selfrefo(:,16) = (/ &
        &  0.782126e+00_rb, 0.710682e+00_rb, 0.645764e+00_rb, 0.586776e+00_rb, 0.533177e+00_rb, &
        &  0.484473e+00_rb, 0.440219e+00_rb, 0.400007e+00_rb, 0.363468e+00_rb, 0.330266e+00_rb /)
     
      end subroutine sw_kgb179

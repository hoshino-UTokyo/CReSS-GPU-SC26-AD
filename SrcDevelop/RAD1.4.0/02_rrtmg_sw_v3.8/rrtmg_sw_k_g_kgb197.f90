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
!!prev      subroutine sw_kgb19
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine sw_kgb197
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
      use rrsw_kg19, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl

      implicit none
      save
  
!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.106275e-05_rb, 0.104185e-05_rb, 0.420154e-05_rb /)
      forrefo(:, 2) = (/ 0.154343e-05_rb, 0.653193e-05_rb, 0.174596e-04_rb /)
      forrefo(:, 3) = (/ 0.348917e-05_rb, 0.108420e-04_rb, 0.540849e-04_rb /)
      forrefo(:, 4) = (/ 0.145822e-04_rb, 0.156027e-04_rb, 0.881263e-04_rb /)
      forrefo(:, 5) = (/ 0.220204e-04_rb, 0.819892e-04_rb, 0.817937e-04_rb /)
      forrefo(:, 6) = (/ 0.447840e-04_rb, 0.121116e-03_rb, 0.932635e-04_rb /)
      forrefo(:, 7) = (/ 0.166516e-03_rb, 0.147640e-03_rb, 0.754029e-04_rb /)
      forrefo(:, 8) = (/ 0.234756e-03_rb, 0.145934e-03_rb, 0.771734e-04_rb /)
      forrefo(:, 9) = (/ 0.289207e-03_rb, 0.146768e-03_rb, 0.677806e-04_rb /)
      forrefo(:,10) = (/ 0.334959e-03_rb, 0.125513e-03_rb, 0.636648e-04_rb /)
      forrefo(:,11) = (/ 0.333755e-03_rb, 0.136575e-03_rb, 0.593651e-04_rb /)
      forrefo(:,12) = (/ 0.340042e-03_rb, 0.116259e-03_rb, 0.595192e-04_rb /)
      forrefo(:,13) = (/ 0.422470e-03_rb, 0.148691e-03_rb, 0.630266e-04_rb /)
      forrefo(:,14) = (/ 0.440655e-03_rb, 0.461917e-04_rb, 0.108222e-04_rb /)
      forrefo(:,15) = (/ 0.486207e-03_rb, 0.428458e-03_rb, 0.108086e-04_rb /)
      forrefo(:,16) = (/ 0.657463e-03_rb, 0.657446e-03_rb, 0.126190e-04_rb /)

!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).
     
      selfrefo(:, 1) = (/ &
        &  0.331728e-03_rb, 0.287480e-03_rb, 0.249135e-03_rb, 0.215904e-03_rb, 0.187106e-03_rb, &
        &  0.162149e-03_rb, 0.140520e-03_rb, 0.121777e-03_rb, 0.105534e-03_rb, 0.914573e-04_rb /)
      selfrefo(:, 2) = (/ &
        &  0.882628e-03_rb, 0.698914e-03_rb, 0.553439e-03_rb, 0.438244e-03_rb, 0.347026e-03_rb, &
        &  0.274795e-03_rb, 0.217598e-03_rb, 0.172306e-03_rb, 0.136442e-03_rb, 0.108042e-03_rb /)
      selfrefo(:, 3) = (/ &
        &  0.115461e-02_rb, 0.937203e-03_rb, 0.760730e-03_rb, 0.617486e-03_rb, 0.501215e-03_rb, &
        &  0.406837e-03_rb, 0.330231e-03_rb, 0.268049e-03_rb, 0.217576e-03_rb, 0.176607e-03_rb /)
      selfrefo(:, 4) = (/ &
        &  0.103450e-02_rb, 0.960268e-03_rb, 0.891360e-03_rb, 0.827397e-03_rb, 0.768024e-03_rb, &
        &  0.712911e-03_rb, 0.661754e-03_rb, 0.614267e-03_rb, 0.570188e-03_rb, 0.529272e-03_rb /)
      selfrefo(:, 5) = (/ &
        &  0.289040e-02_rb, 0.240129e-02_rb, 0.199495e-02_rb, 0.165737e-02_rb, 0.137692e-02_rb, &
        &  0.114392e-02_rb, 0.950351e-03_rb, 0.789535e-03_rb, 0.655933e-03_rb, 0.544938e-03_rb /)
      selfrefo(:, 6) = (/ &
        &  0.361772e-02_rb, 0.306611e-02_rb, 0.259861e-02_rb, 0.220239e-02_rb, 0.186659e-02_rb, &
        &  0.158198e-02_rb, 0.134077e-02_rb, 0.113634e-02_rb, 0.963078e-03_rb, 0.816234e-03_rb /)
      selfrefo(:, 7) = (/ &
        &  0.329878e-02_rb, 0.318245e-02_rb, 0.307021e-02_rb, 0.296194e-02_rb, 0.285749e-02_rb, &
        &  0.275671e-02_rb, 0.265950e-02_rb, 0.256571e-02_rb, 0.247522e-02_rb, 0.238793e-02_rb /)
      selfrefo(:, 8) = (/ &
        &  0.293562e-02_rb, 0.300077e-02_rb, 0.306737e-02_rb, 0.313544e-02_rb, 0.320503e-02_rb, &
        &  0.327615e-02_rb, 0.334886e-02_rb, 0.342318e-02_rb, 0.349915e-02_rb, 0.357680e-02_rb /)
      selfrefo(:, 9) = (/ &
        &  0.281453e-02_rb, 0.295894e-02_rb, 0.311076e-02_rb, 0.327038e-02_rb, 0.343818e-02_rb, &
        &  0.361459e-02_rb, 0.380006e-02_rb, 0.399504e-02_rb, 0.420002e-02_rb, 0.441553e-02_rb /)
      selfrefo(:,10) = (/ &
        &  0.239488e-02_rb, 0.262487e-02_rb, 0.287696e-02_rb, 0.315325e-02_rb, 0.345607e-02_rb, &
        &  0.378798e-02_rb, 0.415176e-02_rb, 0.455048e-02_rb, 0.498749e-02_rb, 0.546647e-02_rb /)
      selfrefo(:,11) = (/ &
        &  0.271001e-02_rb, 0.292235e-02_rb, 0.315134e-02_rb, 0.339826e-02_rb, 0.366453e-02_rb, &
        &  0.395167e-02_rb, 0.426131e-02_rb, 0.459521e-02_rb, 0.495527e-02_rb, 0.534354e-02_rb /)
      selfrefo(:,12) = (/ &
        &  0.206702e-02_rb, 0.232254e-02_rb, 0.260966e-02_rb, 0.293226e-02_rb, 0.329475e-02_rb, &
        &  0.370204e-02_rb, 0.415969e-02_rb, 0.467391e-02_rb, 0.525169e-02_rb, 0.590090e-02_rb /)
      selfrefo(:,13) = (/ &
        &  0.227023e-02_rb, 0.257331e-02_rb, 0.291685e-02_rb, 0.330626e-02_rb, 0.374766e-02_rb, &
        &  0.424799e-02_rb, 0.481511e-02_rb, 0.545794e-02_rb, 0.618660e-02_rb, 0.701253e-02_rb /)
      selfrefo(:,14) = (/ &
        &  0.851078e-03_rb, 0.111512e-02_rb, 0.146109e-02_rb, 0.191439e-02_rb, 0.250832e-02_rb, &
        &  0.328653e-02_rb, 0.430617e-02_rb, 0.564215e-02_rb, 0.739261e-02_rb, 0.968616e-02_rb /)
      selfrefo(:,15) = (/ &
        &  0.742711e-02_rb, 0.721347e-02_rb, 0.700598e-02_rb, 0.680446e-02_rb, 0.660873e-02_rb, &
        &  0.641863e-02_rb, 0.623400e-02_rb, 0.605468e-02_rb, 0.588052e-02_rb, 0.571137e-02_rb /)
      selfrefo(:,16) = (/ &
        &  0.107170e-01_rb, 0.101913e-01_rb, 0.969138e-02_rb, 0.921599e-02_rb, 0.876392e-02_rb, &
        &  0.833402e-02_rb, 0.792521e-02_rb, 0.753646e-02_rb, 0.716677e-02_rb, 0.681522e-02_rb /)
  
      end subroutine sw_kgb197

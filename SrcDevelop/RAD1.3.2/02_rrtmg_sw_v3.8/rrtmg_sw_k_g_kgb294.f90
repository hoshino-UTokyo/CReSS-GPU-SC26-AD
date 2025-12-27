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
!!prev      subroutine sw_kgb29
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for kao(:,:,:)
!!
      subroutine sw_kgb294
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg29, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
!!prev                            absh2oo, absco2o, rayl
      use rrsw_kg29, only : selfrefo, forrefo

      implicit none
      save

!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.299818e-05_rb, 0.209282e-05_rb, 0.988353e-04_rb, 0.632178e-03_rb /)
      forrefo(:, 2) = (/ 0.633648e-05_rb, 0.509214e-04_rb, 0.650535e-03_rb, 0.264019e-02_rb /)
      forrefo(:, 3) = (/ 0.636782e-04_rb, 0.136577e-03_rb, 0.166500e-02_rb, 0.750821e-02_rb /)
      forrefo(:, 4) = (/ 0.472314e-03_rb, 0.988296e-03_rb, 0.585751e-02_rb, 0.187352e-01_rb /)
      forrefo(:, 5) = (/ 0.558635e-02_rb, 0.856489e-02_rb, 0.157438e-01_rb, 0.181471e-01_rb /)
      forrefo(:, 6) = (/ 0.217395e-01_rb, 0.229156e-01_rb, 0.230125e-01_rb, 0.143821e-01_rb /)
      forrefo(:, 7) = (/ 0.277222e-01_rb, 0.299252e-01_rb, 0.208929e-01_rb, 0.826748e-02_rb /)
      forrefo(:, 8) = (/ 0.252119e-01_rb, 0.262911e-01_rb, 0.187663e-01_rb, 0.417110e-02_rb /)
      forrefo(:, 9) = (/ 0.304941e-01_rb, 0.175545e-01_rb, 0.971224e-02_rb, 0.142023e-02_rb /)
      forrefo(:,10) = (/ 0.327200e-01_rb, 0.215788e-01_rb, 0.346831e-02_rb, 0.157989e-02_rb /)
      forrefo(:,11) = (/ 0.324955e-01_rb, 0.228571e-01_rb, 0.171749e-02_rb, 0.226853e-02_rb /)
      forrefo(:,12) = (/ 0.326588e-01_rb, 0.198544e-01_rb, 0.532339e-06_rb, 0.279086e-02_rb /)
      forrefo(:,13) = (/ 0.345157e-01_rb, 0.168679e-01_rb, 0.505361e-06_rb, 0.276647e-02_rb /)
      forrefo(:,14) = (/ 0.448765e-01_rb, 0.123791e-02_rb, 0.488367e-06_rb, 0.122245e-02_rb /)
      forrefo(:,15) = (/ 0.486925e-01_rb, 0.464371e-06_rb, 0.464241e-06_rb, 0.753846e-06_rb /)
      forrefo(:,16) = (/ 0.530511e-01_rb, 0.376234e-06_rb, 0.409824e-06_rb, 0.470650e-06_rb /)

!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
        &  0.118069e+00_rb, 0.713523e-01_rb, 0.431199e-01_rb, 0.260584e-01_rb, 0.157477e-01_rb, &
        &  0.951675e-02_rb, 0.575121e-02_rb, 0.347560e-02_rb, 0.210039e-02_rb, 0.126932e-02_rb /)
      selfrefo(:, 2) = (/ &
        &  0.137081e-01_rb, 0.139046e-01_rb, 0.141040e-01_rb, 0.143061e-01_rb, 0.145112e-01_rb, &
        &  0.147193e-01_rb, 0.149303e-01_rb, 0.151443e-01_rb, 0.153614e-01_rb, 0.155816e-01_rb /)
      selfrefo(:, 3) = (/ &
        &  0.166575e-01_rb, 0.164916e-01_rb, 0.163273e-01_rb, 0.161647e-01_rb, 0.160037e-01_rb, &
        &  0.158443e-01_rb, 0.156864e-01_rb, 0.155302e-01_rb, 0.153755e-01_rb, 0.152224e-01_rb /)
      selfrefo(:, 4) = (/ &
        &  0.597379e-01_rb, 0.509517e-01_rb, 0.434579e-01_rb, 0.370662e-01_rb, 0.316145e-01_rb, &
        &  0.269647e-01_rb, 0.229988e-01_rb, 0.196162e-01_rb, 0.167311e-01_rb, 0.142703e-01_rb /)
      selfrefo(:, 5) = (/ &
        &  0.227517e+00_rb, 0.198401e+00_rb, 0.173011e+00_rb, 0.150870e+00_rb, 0.131563e+00_rb, &
        &  0.114726e+00_rb, 0.100044e+00_rb, 0.872415e-01_rb, 0.760769e-01_rb, 0.663411e-01_rb /)
      selfrefo(:, 6) = (/ &
        &  0.453235e+00_rb, 0.414848e+00_rb, 0.379712e+00_rb, 0.347552e+00_rb, 0.318116e+00_rb, &
        &  0.291173e+00_rb, 0.266512e+00_rb, 0.243940e+00_rb, 0.223279e+00_rb, 0.204368e+00_rb /)
      selfrefo(:, 7) = (/ &
        &  0.569263e+00_rb, 0.516415e+00_rb, 0.468473e+00_rb, 0.424982e+00_rb, 0.385528e+00_rb, &
        &  0.349737e+00_rb, 0.317269e+00_rb, 0.287815e+00_rb, 0.261095e+00_rb, 0.236856e+00_rb /)
      selfrefo(:, 8) = (/ &
        &  0.490314e+00_rb, 0.448042e+00_rb, 0.409413e+00_rb, 0.374116e+00_rb, 0.341861e+00_rb, &
        &  0.312387e+00_rb, 0.285455e+00_rb, 0.260844e+00_rb, 0.238355e+00_rb, 0.217805e+00_rb /)
      selfrefo(:, 9) = (/ &
        &  0.258162e+00_rb, 0.265085e+00_rb, 0.272193e+00_rb, 0.279493e+00_rb, 0.286988e+00_rb, &
        &  0.294684e+00_rb, 0.302586e+00_rb, 0.310701e+00_rb, 0.319033e+00_rb, 0.327588e+00_rb /)
      selfrefo(:,10) = (/ &
        &  0.332019e+00_rb, 0.331902e+00_rb, 0.331784e+00_rb, 0.331666e+00_rb, 0.331549e+00_rb, &
        &  0.331431e+00_rb, 0.331314e+00_rb, 0.331197e+00_rb, 0.331079e+00_rb, 0.330962e+00_rb /)
      selfrefo(:,11) = (/ &
        &  0.357523e+00_rb, 0.353154e+00_rb, 0.348839e+00_rb, 0.344576e+00_rb, 0.340366e+00_rb, &
        &  0.336207e+00_rb, 0.332099e+00_rb, 0.328041e+00_rb, 0.324032e+00_rb, 0.320073e+00_rb /)
      selfrefo(:,12) = (/ &
        &  0.294662e+00_rb, 0.299043e+00_rb, 0.303488e+00_rb, 0.308000e+00_rb, 0.312579e+00_rb, &
        &  0.317226e+00_rb, 0.321941e+00_rb, 0.326727e+00_rb, 0.331585e+00_rb, 0.336514e+00_rb /)
      selfrefo(:,13) = (/ &
        &  0.227445e+00_rb, 0.241545e+00_rb, 0.256519e+00_rb, 0.272422e+00_rb, 0.289311e+00_rb, &
        &  0.307247e+00_rb, 0.326294e+00_rb, 0.346523e+00_rb, 0.368005e+00_rb, 0.390820e+00_rb /)
      selfrefo(:,14) = (/ &
        &  0.616203e-02_rb, 0.113523e-01_rb, 0.209144e-01_rb, 0.385307e-01_rb, 0.709852e-01_rb, &
        &  0.130776e+00_rb, 0.240929e+00_rb, 0.443865e+00_rb, 0.817733e+00_rb, 0.150651e+01_rb /)
      selfrefo(:,15) = (/ &
        &  0.279552e-03_rb, 0.808472e-03_rb, 0.233812e-02_rb, 0.676192e-02_rb, 0.195557e-01_rb, &
        &  0.565555e-01_rb, 0.163560e+00_rb, 0.473020e+00_rb, 0.136799e+01_rb, 0.395626e+01_rb /)
      selfrefo(:,16) = (/ &
        &  0.261006e-03_rb, 0.771043e-03_rb, 0.227776e-02_rb, 0.672879e-02_rb, 0.198777e-01_rb, &
        &  0.587212e-01_rb, 0.173470e+00_rb, 0.512452e+00_rb, 0.151385e+01_rb, 0.447209e+01_rb /)
     
      end subroutine sw_kgb294

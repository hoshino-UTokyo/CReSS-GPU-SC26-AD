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
!!prev      subroutine sw_kgb18
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine sw_kgb187
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg18, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg18, only : selfrefo, forrefo

      implicit none
      save

!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.860560e-06_rb, 0.130439e-05_rb, 0.382378e-05_rb /)
      forrefo(:, 2) = (/ 0.817926e-06_rb, 0.158599e-05_rb, 0.658771e-04_rb /)
      forrefo(:, 3) = (/ 0.129369e-05_rb, 0.824406e-05_rb, 0.952778e-04_rb /)
      forrefo(:, 4) = (/ 0.438918e-05_rb, 0.375356e-04_rb, 0.119111e-03_rb /)
      forrefo(:, 5) = (/ 0.306057e-04_rb, 0.622798e-04_rb, 0.100740e-03_rb /)
      forrefo(:, 6) = (/ 0.891934e-04_rb, 0.856393e-04_rb, 0.635583e-04_rb /)
      forrefo(:, 7) = (/ 0.171959e-03_rb, 0.173431e-03_rb, 0.611721e-04_rb /)
      forrefo(:, 8) = (/ 0.357795e-03_rb, 0.247261e-03_rb, 0.488864e-04_rb /)
      forrefo(:, 9) = (/ 0.326623e-03_rb, 0.289471e-03_rb, 0.548834e-04_rb /)
      forrefo(:,10) = (/ 0.345103e-03_rb, 0.320898e-03_rb, 0.633214e-04_rb /)
      forrefo(:,11) = (/ 0.392567e-03_rb, 0.325153e-03_rb, 0.744479e-04_rb /)
      forrefo(:,12) = (/ 0.349277e-03_rb, 0.345610e-03_rb, 0.916479e-04_rb /)
      forrefo(:,13) = (/ 0.425161e-03_rb, 0.348452e-03_rb, 0.125788e-03_rb /)
      forrefo(:,14) = (/ 0.407594e-03_rb, 0.435836e-03_rb, 0.287583e-03_rb /)
      forrefo(:,15) = (/ 0.521605e-03_rb, 0.486596e-03_rb, 0.483511e-03_rb /)
      forrefo(:,16) = (/ 0.773790e-03_rb, 0.737247e-03_rb, 0.665939e-03_rb /)

!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).
           
      selfrefo(:, 1) = (/ &
        &  0.750370e-03_rb, 0.644938e-03_rb, 0.554321e-03_rb, 0.476436e-03_rb, 0.409494e-03_rb, &
        &  0.351957e-03_rb, 0.302505e-03_rb, 0.260002e-03_rb, 0.223470e-03_rb, 0.192071e-03_rb /)
      selfrefo(:, 2) = (/ &
        &  0.136135e-02_rb, 0.113187e-02_rb, 0.941076e-03_rb, 0.782440e-03_rb, 0.650546e-03_rb, &
        &  0.540885e-03_rb, 0.449709e-03_rb, 0.373902e-03_rb, 0.310874e-03_rb, 0.258471e-03_rb /)
      selfrefo(:, 3) = (/ &
        &  0.333950e-02_rb, 0.256391e-02_rb, 0.196845e-02_rb, 0.151129e-02_rb, 0.116030e-02_rb, &
        &  0.890824e-03_rb, 0.683934e-03_rb, 0.525093e-03_rb, 0.403143e-03_rb, 0.309515e-03_rb /)
      selfrefo(:, 4) = (/ &
        &  0.793392e-02_rb, 0.589865e-02_rb, 0.438548e-02_rb, 0.326048e-02_rb, 0.242408e-02_rb, &
        &  0.180223e-02_rb, 0.133991e-02_rb, 0.996186e-03_rb, 0.740636e-03_rb, 0.550642e-03_rb /)
      selfrefo(:, 5) = (/ &
        &  0.828169e-02_rb, 0.703139e-02_rb, 0.596984e-02_rb, 0.506856e-02_rb, 0.430335e-02_rb, &
        &  0.365366e-02_rb, 0.310206e-02_rb, 0.263374e-02_rb, 0.223612e-02_rb, 0.189852e-02_rb /)
      selfrefo(:, 6) = (/ &
        &  0.834190e-02_rb, 0.780225e-02_rb, 0.729750e-02_rb, 0.682541e-02_rb, 0.638386e-02_rb, &
        &  0.597087e-02_rb, 0.558460e-02_rb, 0.522332e-02_rb, 0.488541e-02_rb, 0.456936e-02_rb /)
      selfrefo(:, 7) = (/ &
        &  0.119082e-01_rb, 0.112566e-01_rb, 0.106406e-01_rb, 0.100583e-01_rb, 0.950785e-02_rb, &
        &  0.898755e-02_rb, 0.849571e-02_rb, 0.803080e-02_rb, 0.759132e-02_rb, 0.717590e-02_rb /)
      selfrefo(:, 8) = (/ &
        &  0.144004e-01_rb, 0.141762e-01_rb, 0.139554e-01_rb, 0.137381e-01_rb, 0.135241e-01_rb, &
        &  0.133135e-01_rb, 0.131062e-01_rb, 0.129021e-01_rb, 0.127011e-01_rb, 0.125033e-01_rb /)
      selfrefo(:, 9) = (/ &
        &  0.186171e-01_rb, 0.175281e-01_rb, 0.165027e-01_rb, 0.155373e-01_rb, 0.146284e-01_rb, &
        &  0.137726e-01_rb, 0.129670e-01_rb, 0.122084e-01_rb, 0.114942e-01_rb, 0.108218e-01_rb /)
      selfrefo(:,10) = (/ &
        &  0.209396e-01_rb, 0.195077e-01_rb, 0.181737e-01_rb, 0.169309e-01_rb, 0.157731e-01_rb, &
        &  0.146945e-01_rb, 0.136897e-01_rb, 0.127535e-01_rb, 0.118814e-01_rb, 0.110689e-01_rb /)
      selfrefo(:,11) = (/ &
        &  0.203661e-01_rb, 0.193311e-01_rb, 0.183487e-01_rb, 0.174163e-01_rb, 0.165312e-01_rb, &
        &  0.156911e-01_rb, 0.148937e-01_rb, 0.141368e-01_rb, 0.134184e-01_rb, 0.127365e-01_rb /)
      selfrefo(:,12) = (/ &
        &  0.226784e-01_rb, 0.210210e-01_rb, 0.194848e-01_rb, 0.180608e-01_rb, 0.167409e-01_rb, &
        &  0.155174e-01_rb, 0.143834e-01_rb, 0.133322e-01_rb, 0.123579e-01_rb, 0.114547e-01_rb /)
      selfrefo(:,13) = (/ &
        &  0.221773e-01_rb, 0.210306e-01_rb, 0.199433e-01_rb, 0.189122e-01_rb, 0.179344e-01_rb, &
        &  0.170071e-01_rb, 0.161278e-01_rb, 0.152939e-01_rb, 0.145032e-01_rb, 0.137533e-01_rb /)
      selfrefo(:,14) = (/ &
        &  0.275920e-01_rb, 0.252595e-01_rb, 0.231241e-01_rb, 0.211693e-01_rb, 0.193797e-01_rb, &
        &  0.177415e-01_rb, 0.162417e-01_rb, 0.148687e-01_rb, 0.136117e-01_rb, 0.124610e-01_rb /)
      selfrefo(:,15) = (/ &
        &  0.288687e-01_rb, 0.269968e-01_rb, 0.252462e-01_rb, 0.236092e-01_rb, 0.220783e-01_rb, &
        &  0.206466e-01_rb, 0.193078e-01_rb, 0.180559e-01_rb, 0.168851e-01_rb, 0.157902e-01_rb /)
      selfrefo(:,16) = (/ &
        &  0.371842e-01_rb, 0.347595e-01_rb, 0.324929e-01_rb, 0.303741e-01_rb, 0.283934e-01_rb, &
        &  0.265419e-01_rb, 0.248112e-01_rb, 0.231933e-01_rb, 0.216809e-01_rb, 0.202671e-01_rb /)
  
      end subroutine sw_kgb187

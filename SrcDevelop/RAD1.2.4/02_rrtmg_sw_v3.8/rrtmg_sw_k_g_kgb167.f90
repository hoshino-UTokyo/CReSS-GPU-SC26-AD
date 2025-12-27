!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 22:22:22 $

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
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
!!prev      subroutine sw_kgb16
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine sw_kgb167
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg16, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg16, only : selfrefo, forrefo

      implicit none
      save

!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.525585e-05_rb, 0.527618e-05_rb, 0.746929e-04_rb /)
      forrefo(:, 2) = (/ 0.794660e-05_rb, 0.136902e-04_rb, 0.849878e-04_rb /)
      forrefo(:, 3) = (/ 0.197099e-04_rb, 0.733094e-04_rb, 0.121687e-03_rb /)
      forrefo(:, 4) = (/ 0.148274e-03_rb, 0.169776e-03_rb, 0.164848e-03_rb /)
      forrefo(:, 5) = (/ 0.230296e-03_rb, 0.210384e-03_rb, 0.182028e-03_rb /)
      forrefo(:, 6) = (/ 0.280575e-03_rb, 0.259217e-03_rb, 0.196080e-03_rb /)
      forrefo(:, 7) = (/ 0.329034e-03_rb, 0.291575e-03_rb, 0.207044e-03_rb /)
      forrefo(:, 8) = (/ 0.349989e-03_rb, 0.323471e-03_rb, 0.225712e-03_rb /)
      forrefo(:, 9) = (/ 0.366097e-03_rb, 0.321519e-03_rb, 0.253150e-03_rb /)
      forrefo(:,10) = (/ 0.383589e-03_rb, 0.355314e-03_rb, 0.262555e-03_rb /)
      forrefo(:,11) = (/ 0.375933e-03_rb, 0.372443e-03_rb, 0.261313e-03_rb /)
      forrefo(:,12) = (/ 0.370652e-03_rb, 0.382366e-03_rb, 0.250070e-03_rb /)
      forrefo(:,13) = (/ 0.375092e-03_rb, 0.379542e-03_rb, 0.265794e-03_rb /)
      forrefo(:,14) = (/ 0.389705e-03_rb, 0.384274e-03_rb, 0.322135e-03_rb /)
      forrefo(:,15) = (/ 0.372084e-03_rb, 0.390422e-03_rb, 0.370035e-03_rb /)
      forrefo(:,16) = (/ 0.437802e-03_rb, 0.373406e-03_rb, 0.373222e-03_rb /)
      
!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).
!     -----------------------------------------------------------------
      
      selfrefo(:, 1) = (/ &
        &  0.126758e-02_rb, 0.105253e-02_rb, 0.873963e-03_rb, 0.725690e-03_rb, 0.602573e-03_rb, &
        &  0.500344e-03_rb, 0.415458e-03_rb, 0.344973e-03_rb, 0.286447e-03_rb, 0.237849e-03_rb /)
      selfrefo(:, 2) = (/ &
        &  0.144006e-02_rb, 0.118514e-02_rb, 0.975351e-03_rb, 0.802697e-03_rb, 0.660606e-03_rb, &
        &  0.543667e-03_rb, 0.447429e-03_rb, 0.368226e-03_rb, 0.303044e-03_rb, 0.249400e-03_rb /)
      selfrefo(:, 3) = (/ &
        &  0.294018e-02_rb, 0.227428e-02_rb, 0.175920e-02_rb, 0.136077e-02_rb, 0.105258e-02_rb, &
        &  0.814189e-03_rb, 0.629789e-03_rb, 0.487153e-03_rb, 0.376821e-03_rb, 0.291478e-03_rb /)
      selfrefo(:, 4) = (/ &
        &  0.395290e-02_rb, 0.348405e-02_rb, 0.307081e-02_rb, 0.270658e-02_rb, 0.238556e-02_rb, &
        &  0.210261e-02_rb, 0.185322e-02_rb, 0.163341e-02_rb, 0.143967e-02_rb, 0.126891e-02_rb /)
      selfrefo(:, 5) = (/ &
        &  0.419122e-02_rb, 0.385638e-02_rb, 0.354829e-02_rb, 0.326481e-02_rb, 0.300398e-02_rb, &
        &  0.276399e-02_rb, 0.254317e-02_rb, 0.234000e-02_rb, 0.215305e-02_rb, 0.198104e-02_rb /)
      selfrefo(:, 6) = (/ &
        &  0.495659e-02_rb, 0.456777e-02_rb, 0.420945e-02_rb, 0.387924e-02_rb, 0.357494e-02_rb, &
        &  0.329450e-02_rb, 0.303606e-02_rb, 0.279790e-02_rb, 0.257842e-02_rb, 0.237615e-02_rb /)
      selfrefo(:, 7) = (/ &
        &  0.526981e-02_rb, 0.490687e-02_rb, 0.456893e-02_rb, 0.425426e-02_rb, 0.396126e-02_rb, &
        &  0.368844e-02_rb, 0.343441e-02_rb, 0.319788e-02_rb, 0.297764e-02_rb, 0.277256e-02_rb /)
      selfrefo(:, 8) = (/ &
        &  0.575426e-02_rb, 0.531597e-02_rb, 0.491106e-02_rb, 0.453699e-02_rb, 0.419141e-02_rb, &
        &  0.387216e-02_rb, 0.357722e-02_rb, 0.330475e-02_rb, 0.305303e-02_rb, 0.282048e-02_rb /)
      selfrefo(:, 9) = (/ &
        &  0.549881e-02_rb, 0.514328e-02_rb, 0.481074e-02_rb, 0.449970e-02_rb, 0.420877e-02_rb, &
        &  0.393665e-02_rb, 0.368213e-02_rb, 0.344406e-02_rb, 0.322138e-02_rb, 0.301310e-02_rb /)
      selfrefo(:,10) = (/ &
        &  0.605357e-02_rb, 0.561246e-02_rb, 0.520349e-02_rb, 0.482432e-02_rb, 0.447278e-02_rb, &
        &  0.414686e-02_rb, 0.384469e-02_rb, 0.356453e-02_rb, 0.330479e-02_rb, 0.306398e-02_rb /)
      selfrefo(:,11) = (/ &
        &  0.640504e-02_rb, 0.587858e-02_rb, 0.539540e-02_rb, 0.495194e-02_rb, 0.454492e-02_rb, &
        &  0.417136e-02_rb, 0.382850e-02_rb, 0.351382e-02_rb, 0.322501e-02_rb, 0.295993e-02_rb /)
      selfrefo(:,12) = (/ &
        &  0.677803e-02_rb, 0.615625e-02_rb, 0.559152e-02_rb, 0.507859e-02_rb, 0.461271e-02_rb, &
        &  0.418957e-02_rb, 0.380524e-02_rb, 0.345617e-02_rb, 0.313913e-02_rb, 0.285116e-02_rb /)
      selfrefo(:,13) = (/ &
        &  0.690347e-02_rb, 0.627003e-02_rb, 0.569472e-02_rb, 0.517219e-02_rb, 0.469761e-02_rb, &
        &  0.426658e-02_rb, 0.387509e-02_rb, 0.351953e-02_rb, 0.319659e-02_rb, 0.290328e-02_rb /)
      selfrefo(:,14) = (/ &
        &  0.692680e-02_rb, 0.632795e-02_rb, 0.578087e-02_rb, 0.528109e-02_rb, 0.482452e-02_rb, &
        &  0.440742e-02_rb, 0.402638e-02_rb, 0.367828e-02_rb, 0.336028e-02_rb, 0.306977e-02_rb /)
      selfrefo(:,15) = (/ &
        &  0.754894e-02_rb, 0.681481e-02_rb, 0.615207e-02_rb, 0.555378e-02_rb, 0.501367e-02_rb, &
        &  0.452609e-02_rb, 0.408593e-02_rb, 0.368857e-02_rb, 0.332986e-02_rb, 0.300603e-02_rb /)
      selfrefo(:,16) = (/ &
        &  0.760689e-02_rb, 0.709755e-02_rb, 0.662232e-02_rb, 0.617891e-02_rb, 0.576519e-02_rb, &
        &  0.537917e-02_rb, 0.501899e-02_rb, 0.468293e-02_rb, 0.436938e-02_rb, 0.407682e-02_rb /)
     
      end subroutine sw_kgb167

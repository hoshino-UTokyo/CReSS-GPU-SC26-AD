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
!!prev      subroutine sw_kgb20
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for  sfluxrefo(:), absch4o(:), rayl
!!
      subroutine sw_kgb204
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg20, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
!!prev                            absch4o, rayl
      use rrsw_kg20, only : selfrefo, forrefo

      implicit none
      save

!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.214504e-06_rb, 0.460418e-06_rb, 0.357608e-05_rb, 0.192037e-05_rb /)
      forrefo(:, 2) = (/ 0.142576e-05_rb, 0.364463e-05_rb, 0.117033e-04_rb, 0.112085e-04_rb /)
      forrefo(:, 3) = (/ 0.101536e-04_rb, 0.124096e-04_rb, 0.509190e-04_rb, 0.565282e-04_rb /)
      forrefo(:, 4) = (/ 0.143394e-03_rb, 0.154700e-03_rb, 0.466498e-03_rb, 0.918829e-03_rb /)
      forrefo(:, 5) = (/ 0.251631e-02_rb, 0.241729e-02_rb, 0.240057e-02_rb, 0.350408e-02_rb /)
      forrefo(:, 6) = (/ 0.410309e-02_rb, 0.416851e-02_rb, 0.390925e-02_rb, 0.383694e-02_rb /)
      forrefo(:, 7) = (/ 0.445387e-02_rb, 0.448657e-02_rb, 0.432310e-02_rb, 0.370739e-02_rb /)
      forrefo(:, 8) = (/ 0.458150e-02_rb, 0.460014e-02_rb, 0.450245e-02_rb, 0.336718e-02_rb /)
      forrefo(:, 9) = (/ 0.465423e-02_rb, 0.465595e-02_rb, 0.467006e-02_rb, 0.368061e-02_rb /)
      forrefo(:,10) = (/ 0.493955e-02_rb, 0.490181e-02_rb, 0.481941e-02_rb, 0.367577e-02_rb /)
      forrefo(:,11) = (/ 0.511876e-02_rb, 0.490981e-02_rb, 0.493303e-02_rb, 0.357423e-02_rb /)
      forrefo(:,12) = (/ 0.509845e-02_rb, 0.511556e-02_rb, 0.504031e-02_rb, 0.355915e-02_rb /)
      forrefo(:,13) = (/ 0.523822e-02_rb, 0.530473e-02_rb, 0.523811e-02_rb, 0.414259e-02_rb /)
      forrefo(:,14) = (/ 0.551133e-02_rb, 0.535831e-02_rb, 0.546702e-02_rb, 0.473875e-02_rb /)
      forrefo(:,15) = (/ 0.609781e-02_rb, 0.589859e-02_rb, 0.561187e-02_rb, 0.528981e-02_rb /)
      forrefo(:,16) = (/ 0.644958e-02_rb, 0.631718e-02_rb, 0.625201e-02_rb, 0.600448e-02_rb /)

!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
        &  0.217058e-03_rb, 0.176391e-03_rb, 0.143342e-03_rb, 0.116486e-03_rb, 0.946614e-04_rb, &
        &  0.769257e-04_rb, 0.625131e-04_rb, 0.508007e-04_rb, 0.412828e-04_rb, 0.335481e-04_rb /)
      selfrefo(:, 2) = (/ &
        &  0.598055e-03_rb, 0.484805e-03_rb, 0.393000e-03_rb, 0.318580e-03_rb, 0.258252e-03_rb, &
        &  0.209348e-03_rb, 0.169705e-03_rb, 0.137569e-03_rb, 0.111518e-03_rb, 0.904008e-04_rb /)
      selfrefo(:, 3) = (/ &
        &  0.102691e-02_rb, 0.930281e-03_rb, 0.842740e-03_rb, 0.763437e-03_rb, 0.691596e-03_rb, &
        &  0.626516e-03_rb, 0.567560e-03_rb, 0.514152e-03_rb, 0.465769e-03_rb, 0.421940e-03_rb /)
      selfrefo(:, 4) = (/ &
        &  0.388569e-02_rb, 0.365098e-02_rb, 0.343045e-02_rb, 0.322324e-02_rb, 0.302854e-02_rb, &
        &  0.284561e-02_rb, 0.267372e-02_rb, 0.251222e-02_rb, 0.236047e-02_rb, 0.221789e-02_rb /)
      selfrefo(:, 5) = (/ &
        &  0.349845e-01_rb, 0.326678e-01_rb, 0.305045e-01_rb, 0.284845e-01_rb, 0.265982e-01_rb, &
        &  0.248369e-01_rb, 0.231921e-01_rb, 0.216563e-01_rb, 0.202222e-01_rb, 0.188831e-01_rb /)
      selfrefo(:, 6) = (/ &
        &  0.613705e-01_rb, 0.562676e-01_rb, 0.515890e-01_rb, 0.472994e-01_rb, 0.433665e-01_rb, &
        &  0.397606e-01_rb, 0.364545e-01_rb, 0.334233e-01_rb, 0.306442e-01_rb, 0.280961e-01_rb /)
      selfrefo(:, 7) = (/ &
        &  0.656981e-01_rb, 0.602660e-01_rb, 0.552830e-01_rb, 0.507120e-01_rb, 0.465190e-01_rb, &
        &  0.426726e-01_rb, 0.391443e-01_rb, 0.359077e-01_rb, 0.329387e-01_rb, 0.302153e-01_rb /)
      selfrefo(:, 8) = (/ &
        &  0.671782e-01_rb, 0.616461e-01_rb, 0.565695e-01_rb, 0.519110e-01_rb, 0.476361e-01_rb, &
        &  0.437132e-01_rb, 0.401134e-01_rb, 0.368100e-01_rb, 0.337787e-01_rb, 0.309970e-01_rb /)
      selfrefo(:, 9) = (/ &
        &  0.675902e-01_rb, 0.620888e-01_rb, 0.570351e-01_rb, 0.523928e-01_rb, 0.481284e-01_rb, &
        &  0.442110e-01_rb, 0.406125e-01_rb, 0.373069e-01_rb, 0.342703e-01_rb, 0.314809e-01_rb /)
      selfrefo(:,10) = (/ &
        &  0.708308e-01_rb, 0.651419e-01_rb, 0.599099e-01_rb, 0.550981e-01_rb, 0.506728e-01_rb, &
        &  0.466030e-01_rb, 0.428600e-01_rb, 0.394176e-01_rb, 0.362517e-01_rb, 0.333401e-01_rb /)
      selfrefo(:,11) = (/ &
        &  0.698445e-01_rb, 0.646584e-01_rb, 0.598573e-01_rb, 0.554128e-01_rb, 0.512982e-01_rb, &
        &  0.474892e-01_rb, 0.439630e-01_rb, 0.406986e-01_rb, 0.376766e-01_rb, 0.348791e-01_rb /)
      selfrefo(:,12) = (/ &
        &  0.743921e-01_rb, 0.682057e-01_rb, 0.625337e-01_rb, 0.573334e-01_rb, 0.525655e-01_rb, &
        &  0.481942e-01_rb, 0.441863e-01_rb, 0.405118e-01_rb, 0.371428e-01_rb, 0.340540e-01_rb /)
      selfrefo(:,13) = (/ &
        &  0.775758e-01_rb, 0.709818e-01_rb, 0.649484e-01_rb, 0.594277e-01_rb, 0.543764e-01_rb, &
        &  0.497544e-01_rb, 0.455253e-01_rb, 0.416556e-01_rb, 0.381149e-01_rb, 0.348751e-01_rb /)
      selfrefo(:,14) = (/ &
        &  0.776545e-01_rb, 0.714761e-01_rb, 0.657894e-01_rb, 0.605550e-01_rb, 0.557372e-01_rb, &
        &  0.513026e-01_rb, 0.472209e-01_rb, 0.434639e-01_rb, 0.400058e-01_rb, 0.368229e-01_rb /)
      selfrefo(:,15) = (/ &
        &  0.855675e-01_rb, 0.787337e-01_rb, 0.724456e-01_rb, 0.666598e-01_rb, 0.613360e-01_rb, &
        &  0.564374e-01_rb, 0.519301e-01_rb, 0.477827e-01_rb, 0.439666e-01_rb, 0.404552e-01_rb /)
      selfrefo(:,16) = (/ &
        &  0.934781e-01_rb, 0.855190e-01_rb, 0.782376e-01_rb, 0.715761e-01_rb, 0.654819e-01_rb, &
        &  0.599065e-01_rb, 0.548058e-01_rb, 0.501394e-01_rb, 0.458704e-01_rb, 0.419648e-01_rb /)
     
      end subroutine sw_kgb204

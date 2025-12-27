!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 21:04:30 $
!
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
!!prev      subroutine lw_kgb01
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb015
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg01, only : fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
!!prev                            selfrefo, forrefo
      use rrlw_kg01, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &3.6742e-02_rb,1.0664e-01_rb,2.6132e-01_rb,2.7906e-01_rb,2.8151e-01_rb,2.7465e-01_rb, &
     &2.8530e-01_rb,2.9123e-01_rb,3.0697e-01_rb,3.1801e-01_rb,3.2444e-01_rb,2.7746e-01_rb, &
     &3.1994e-01_rb,2.9750e-01_rb,2.1226e-01_rb,1.2847e-01_rb/)
      forrefo(2,:) = (/ &
     &4.0450e-02_rb,1.1085e-01_rb,2.9205e-01_rb,3.1934e-01_rb,3.1739e-01_rb,3.1450e-01_rb, &
     &3.2797e-01_rb,3.2223e-01_rb,3.3099e-01_rb,3.4800e-01_rb,3.4046e-01_rb,3.5700e-01_rb, &
     &3.8264e-01_rb,3.6679e-01_rb,3.3481e-01_rb,3.2113e-01_rb/)
      forrefo(3,:) = (/ &
     &4.6952e-02_rb,1.1999e-01_rb,3.1473e-01_rb,3.7015e-01_rb,3.6913e-01_rb,3.6352e-01_rb, &
     &3.7754e-01_rb,3.7402e-01_rb,3.7113e-01_rb,3.7720e-01_rb,3.8365e-01_rb,4.0876e-01_rb, &
     &4.2968e-01_rb,4.4186e-01_rb,4.3468e-01_rb,4.7083e-01_rb/)
      forrefo(4,:) = (/ &
     &7.0645e-02_rb,1.6618e-01_rb,2.8516e-01_rb,3.1819e-01_rb,3.0131e-01_rb,2.9552e-01_rb, &
     &2.8972e-01_rb,2.9348e-01_rb,2.8668e-01_rb,2.8483e-01_rb,2.8130e-01_rb,2.7757e-01_rb, &
     &2.9735e-01_rb,3.1684e-01_rb,3.0681e-01_rb,3.6778e-01_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 2.16803e+00_rb, 1.98236e+00_rb, 1.81260e+00_rb, 1.65737e+00_rb, 1.51544e+00_rb, &
     & 1.38567e+00_rb, 1.26700e+00_rb, 1.15850e+00_rb, 1.05929e+00_rb, 9.68576e-01_rb/)
      selfrefo(:, 2) = (/ &
     & 3.70149e+00_rb, 3.43145e+00_rb, 3.18110e+00_rb, 2.94902e+00_rb, 2.73387e+00_rb, &
     & 2.53441e+00_rb, 2.34951e+00_rb, 2.17810e+00_rb, 2.01919e+00_rb, 1.87188e+00_rb/)
      selfrefo(:, 3) = (/ &
     & 6.17433e+00_rb, 5.62207e+00_rb, 5.11920e+00_rb, 4.66131e+00_rb, 4.24438e+00_rb, &
     & 3.86474e+00_rb, 3.51906e+00_rb, 3.20430e+00_rb, 2.91769e+00_rb, 2.65672e+00_rb/)
      selfrefo(:, 4) = (/ &
     & 6.56459e+00_rb, 5.94787e+00_rb, 5.38910e+00_rb, 4.88282e+00_rb, 4.42410e+00_rb, &
     & 4.00848e+00_rb, 3.63190e+00_rb, 3.29070e+00_rb, 2.98155e+00_rb, 2.70145e+00_rb/)
      selfrefo(:, 5) = (/ &
     & 6.49581e+00_rb, 5.91114e+00_rb, 5.37910e+00_rb, 4.89494e+00_rb, 4.45436e+00_rb, &
     & 4.05344e+00_rb, 3.68860e+00_rb, 3.35660e+00_rb, 3.05448e+00_rb, 2.77956e+00_rb/)
      selfrefo(:, 6) = (/ &
     & 6.50189e+00_rb, 5.89381e+00_rb, 5.34260e+00_rb, 4.84294e+00_rb, 4.39001e+00_rb, &
     & 3.97944e+00_rb, 3.60727e+00_rb, 3.26990e+00_rb, 2.96409e+00_rb, 2.68687e+00_rb/)
      selfrefo(:, 7) = (/ &
     & 6.64768e+00_rb, 6.01719e+00_rb, 5.44650e+00_rb, 4.92993e+00_rb, 4.46236e+00_rb, &
     & 4.03914e+00_rb, 3.65605e+00_rb, 3.30930e+00_rb, 2.99543e+00_rb, 2.71134e+00_rb/)
      selfrefo(:, 8) = (/ &
     & 6.43744e+00_rb, 5.87166e+00_rb, 5.35560e+00_rb, 4.88490e+00_rb, 4.45557e+00_rb, &
     & 4.06397e+00_rb, 3.70679e+00_rb, 3.38100e+00_rb, 3.08384e+00_rb, 2.81281e+00_rb/)
      selfrefo(:, 9) = (/ &
     & 6.55466e+00_rb, 5.99777e+00_rb, 5.48820e+00_rb, 5.02192e+00_rb, 4.59525e+00_rb, &
     & 4.20484e+00_rb, 3.84759e+00_rb, 3.52070e+00_rb, 3.22158e+00_rb, 2.94787e+00_rb/)
      selfrefo(:,10) = (/ &
     & 6.84510e+00_rb, 6.26933e+00_rb, 5.74200e+00_rb, 5.25902e+00_rb, 4.81667e+00_rb, &
     & 4.41152e+00_rb, 4.04046e+00_rb, 3.70060e+00_rb, 3.38933e+00_rb, 3.10424e+00_rb/)
      selfrefo(:,11) = (/ &
     & 6.83128e+00_rb, 6.25536e+00_rb, 5.72800e+00_rb, 5.24510e+00_rb, 4.80291e+00_rb, &
     & 4.39799e+00_rb, 4.02722e+00_rb, 3.68770e+00_rb, 3.37681e+00_rb, 3.09212e+00_rb/)
      selfrefo(:,12) = (/ &
     & 7.35969e+00_rb, 6.61719e+00_rb, 5.94960e+00_rb, 5.34936e+00_rb, 4.80968e+00_rb, &
     & 4.32445e+00_rb, 3.88817e+00_rb, 3.49590e+00_rb, 3.14321e+00_rb, 2.82610e+00_rb/)
      selfrefo(:,13) = (/ &
     & 7.50064e+00_rb, 6.80749e+00_rb, 6.17840e+00_rb, 5.60744e+00_rb, 5.08925e+00_rb, &
     & 4.61894e+00_rb, 4.19210e+00_rb, 3.80470e+00_rb, 3.45310e+00_rb, 3.13399e+00_rb/)
      selfrefo(:,14) = (/ &
     & 7.40801e+00_rb, 6.71328e+00_rb, 6.08370e+00_rb, 5.51316e+00_rb, 4.99613e+00_rb, &
     & 4.52759e+00_rb, 4.10298e+00_rb, 3.71820e+00_rb, 3.36950e+00_rb, 3.05351e+00_rb/)
      selfrefo(:,15) = (/ &
     & 7.51895e+00_rb, 6.68846e+00_rb, 5.94970e+00_rb, 5.29254e+00_rb, 4.70796e+00_rb, &
     & 4.18795e+00_rb, 3.72538e+00_rb, 3.31390e+00_rb, 2.94787e+00_rb, 2.62227e+00_rb/)
      selfrefo(:,16) = (/ &
     & 7.84774e+00_rb, 6.80673e+00_rb, 5.90380e+00_rb, 5.12065e+00_rb, 4.44138e+00_rb, &
     & 3.85223e+00_rb, 3.34122e+00_rb, 2.89800e+00_rb, 2.51357e+00_rb, 2.18014e+00_rb/)

      end subroutine lw_kgb015

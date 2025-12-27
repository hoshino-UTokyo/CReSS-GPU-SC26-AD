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
!!prev      subroutine lw_kgb12
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb125
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg12, only : fracrefao, kao, selfrefo, forrefo
      use rrlw_kg12, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.4739e-04_rb,3.1686e-04_rb,8.5973e-04_rb,1.9039e-03_rb,3.1820e-03_rb,3.6596e-03_rb, &
     &3.8724e-03_rb,3.6785e-03_rb,3.7141e-03_rb,3.7646e-03_rb,4.2955e-03_rb,4.6343e-03_rb, &
     &5.0612e-03_rb,4.0227e-03_rb,4.2966e-03_rb,4.6622e-03_rb/)
      forrefo(2,:) = (/ &
     &1.9397e-04_rb,3.6322e-04_rb,8.9797e-04_rb,2.1001e-03_rb,3.0307e-03_rb,3.5563e-03_rb, &
     &3.8498e-03_rb,3.5741e-03_rb,3.5914e-03_rb,3.7658e-03_rb,3.8895e-03_rb,4.4072e-03_rb, &
     &4.7112e-03_rb,4.2230e-03_rb,4.2666e-03_rb,4.6634e-03_rb/)
      forrefo(3,:) = (/ &
     &3.1506e-04_rb,7.3687e-04_rb,1.9678e-03_rb,2.5531e-03_rb,2.8345e-03_rb,2.7809e-03_rb, &
     &2.9124e-03_rb,2.7125e-03_rb,2.6644e-03_rb,2.4907e-03_rb,2.7032e-03_rb,4.0967e-03_rb, &
     &4.1971e-03_rb,4.4507e-03_rb,4.2293e-03_rb,4.6633e-03_rb/)
      forrefo(4,:) = (/ &
     &8.8196e-04_rb,2.1125e-03_rb,2.8042e-03_rb,2.8891e-03_rb,2.4362e-03_rb,1.8733e-03_rb, &
     &1.4078e-03_rb,1.1987e-03_rb,1.2808e-03_rb,8.9050e-04_rb,9.4375e-04_rb,7.8351e-04_rb, &
     &1.0756e-03_rb,1.6586e-03_rb,1.7511e-03_rb,4.7803e-03_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 2.37879e-02_rb, 2.10719e-02_rb, 1.86660e-02_rb, 1.65348e-02_rb, 1.46469e-02_rb, &
     & 1.29746e-02_rb, 1.14932e-02_rb, 1.01810e-02_rb, 9.01858e-03_rb, 7.98888e-03_rb/)
      selfrefo(:, 2) = (/ &
     & 3.10625e-02_rb, 2.82664e-02_rb, 2.57220e-02_rb, 2.34066e-02_rb, 2.12997e-02_rb, &
     & 1.93824e-02_rb, 1.76377e-02_rb, 1.60500e-02_rb, 1.46053e-02_rb, 1.32906e-02_rb/)
      selfrefo(:, 3) = (/ &
     & 5.19103e-02_rb, 4.80004e-02_rb, 4.43850e-02_rb, 4.10419e-02_rb, 3.79506e-02_rb, &
     & 3.50922e-02_rb, 3.24491e-02_rb, 3.00050e-02_rb, 2.77450e-02_rb, 2.56553e-02_rb/)
      selfrefo(:, 4) = (/ &
     & 9.12444e-02_rb, 8.38675e-02_rb, 7.70870e-02_rb, 7.08547e-02_rb, 6.51263e-02_rb, &
     & 5.98610e-02_rb, 5.50214e-02_rb, 5.05730e-02_rb, 4.64843e-02_rb, 4.27262e-02_rb/)
      selfrefo(:, 5) = (/ &
     & 1.11323e-01_rb, 1.04217e-01_rb, 9.75650e-02_rb, 9.13376e-02_rb, 8.55076e-02_rb, &
     & 8.00498e-02_rb, 7.49403e-02_rb, 7.01570e-02_rb, 6.56790e-02_rb, 6.14868e-02_rb/)
      selfrefo(:, 6) = (/ &
     & 1.25301e-01_rb, 1.16877e-01_rb, 1.09020e-01_rb, 1.01691e-01_rb, 9.48543e-02_rb, &
     & 8.84774e-02_rb, 8.25293e-02_rb, 7.69810e-02_rb, 7.18057e-02_rb, 6.69784e-02_rb/)
      selfrefo(:, 7) = (/ &
     & 1.34063e-01_rb, 1.24662e-01_rb, 1.15920e-01_rb, 1.07791e-01_rb, 1.00232e-01_rb, &
     & 9.32035e-02_rb, 8.66676e-02_rb, 8.05900e-02_rb, 7.49386e-02_rb, 6.96836e-02_rb/)
      selfrefo(:, 8) = (/ &
     & 1.26997e-01_rb, 1.18306e-01_rb, 1.10210e-01_rb, 1.02668e-01_rb, 9.56417e-02_rb, &
     & 8.90964e-02_rb, 8.29991e-02_rb, 7.73190e-02_rb, 7.20276e-02_rb, 6.70984e-02_rb/)
      selfrefo(:, 9) = (/ &
     & 1.28823e-01_rb, 1.20235e-01_rb, 1.12220e-01_rb, 1.04739e-01_rb, 9.77569e-02_rb, &
     & 9.12402e-02_rb, 8.51579e-02_rb, 7.94810e-02_rb, 7.41826e-02_rb, 6.92374e-02_rb/)
      selfrefo(:,10) = (/ &
     & 1.35802e-01_rb, 1.25981e-01_rb, 1.16870e-01_rb, 1.08418e-01_rb, 1.00577e-01_rb, &
     & 9.33034e-02_rb, 8.65557e-02_rb, 8.02960e-02_rb, 7.44890e-02_rb, 6.91020e-02_rb/)
      selfrefo(:,11) = (/ &
     & 1.35475e-01_rb, 1.27572e-01_rb, 1.20130e-01_rb, 1.13122e-01_rb, 1.06523e-01_rb, &
     & 1.00309e-01_rb, 9.44573e-02_rb, 8.89470e-02_rb, 8.37582e-02_rb, 7.88721e-02_rb/)
      selfrefo(:,12) = (/ &
     & 1.51195e-01_rb, 1.41159e-01_rb, 1.31790e-01_rb, 1.23043e-01_rb, 1.14876e-01_rb, &
     & 1.07251e-01_rb, 1.00132e-01_rb, 9.34860e-02_rb, 8.72809e-02_rb, 8.14877e-02_rb/)
      selfrefo(:,13) = (/ &
     & 1.57538e-01_rb, 1.47974e-01_rb, 1.38990e-01_rb, 1.30552e-01_rb, 1.22626e-01_rb, &
     & 1.15181e-01_rb, 1.08188e-01_rb, 1.01620e-01_rb, 9.54505e-02_rb, 8.96556e-02_rb/)
      selfrefo(:,14) = (/ &
     & 1.53567e-01_rb, 1.41564e-01_rb, 1.30500e-01_rb, 1.20300e-01_rb, 1.10898e-01_rb, &
     & 1.02231e-01_rb, 9.42406e-02_rb, 8.68750e-02_rb, 8.00851e-02_rb, 7.38259e-02_rb/)
      selfrefo(:,15) = (/ &
     & 1.53687e-01_rb, 1.42981e-01_rb, 1.33020e-01_rb, 1.23753e-01_rb, 1.15132e-01_rb, &
     & 1.07112e-01_rb, 9.96500e-02_rb, 9.27080e-02_rb, 8.62496e-02_rb, 8.02412e-02_rb/)
      selfrefo(:,16) = (/ &
     & 1.65129e-01_rb, 1.53285e-01_rb, 1.42290e-01_rb, 1.32084e-01_rb, 1.22610e-01_rb, &
     & 1.13815e-01_rb, 1.05651e-01_rb, 9.80730e-02_rb, 9.10384e-02_rb, 8.45083e-02_rb/)

      end subroutine lw_kgb125

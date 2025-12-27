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
!!prev      subroutine lw_kgb13
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb138
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg13, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
!!prev                            kbo_mo3, selfrefo, forrefo
      use rrlw_kg13, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.6586e-05_rb,1.9995e-05_rb,1.8582e-05_rb,1.3988e-05_rb,1.3650e-05_rb,1.1079e-05_rb, &
     &9.5855e-06_rb,8.4062e-06_rb,1.3558e-05_rb,1.8620e-05_rb,2.2652e-05_rb,1.7883e-05_rb, &
     &2.6241e-05_rb,3.1171e-05_rb,3.9386e-05_rb,4.4415e-05_rb/)
      forrefo(2,:) = (/ &
     &2.0730e-05_rb,2.3258e-05_rb,2.1543e-05_rb,1.5660e-05_rb,9.7872e-06_rb,8.1078e-06_rb, &
     &7.0246e-06_rb,6.0428e-06_rb,4.8793e-06_rb,4.4937e-06_rb,4.7078e-06_rb,4.6898e-06_rb, &
     &6.9481e-06_rb,8.6269e-06_rb,3.1761e-06_rb,3.1440e-06_rb/)
      forrefo(3,:) = (/ &
     &1.5737e-05_rb,2.2501e-05_rb,2.3520e-05_rb,2.0288e-05_rb,1.2083e-05_rb,6.8256e-06_rb, &
     &6.0637e-06_rb,5.5434e-06_rb,4.3888e-06_rb,3.8435e-06_rb,3.8477e-06_rb,3.8314e-06_rb, &
     &3.8251e-06_rb,3.3637e-06_rb,3.1950e-06_rb,3.1440e-06_rb/)
      forrefo(4,:) = (/ &
     &1.1400e-05_rb,7.9751e-06_rb,8.8659e-06_rb,1.5884e-05_rb,1.9118e-05_rb,1.9429e-05_rb, &
     &2.0532e-05_rb,2.2155e-05_rb,2.3894e-05_rb,2.2984e-05_rb,2.3731e-05_rb,2.4538e-05_rb, &
     &2.6697e-05_rb,1.9329e-05_rb,3.3306e-06_rb,3.2018e-06_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 9.62275e-03_rb, 8.29909e-03_rb, 7.15750e-03_rb, 6.17294e-03_rb, 5.32382e-03_rb, &
     & 4.59150e-03_rb, 3.95991e-03_rb, 3.41520e-03_rb, 2.94542e-03_rb, 2.54026e-03_rb/)
      selfrefo(:, 2) = (/ &
     & 9.76664e-03_rb, 8.47783e-03_rb, 7.35910e-03_rb, 6.38799e-03_rb, 5.54504e-03_rb, &
     & 4.81331e-03_rb, 4.17815e-03_rb, 3.62680e-03_rb, 3.14821e-03_rb, 2.73277e-03_rb/)
      selfrefo(:, 3) = (/ &
     & 9.53856e-03_rb, 8.23750e-03_rb, 7.11390e-03_rb, 6.14356e-03_rb, 5.30558e-03_rb, &
     & 4.58190e-03_rb, 3.95693e-03_rb, 3.41720e-03_rb, 2.95109e-03_rb, 2.54856e-03_rb/)
      selfrefo(:, 4) = (/ &
     & 8.47621e-03_rb, 7.29518e-03_rb, 6.27870e-03_rb, 5.40385e-03_rb, 4.65091e-03_rb, &
     & 4.00287e-03_rb, 3.44513e-03_rb, 2.96510e-03_rb, 2.55196e-03_rb, 2.19638e-03_rb/)
      selfrefo(:, 5) = (/ &
     & 6.71258e-03_rb, 5.95346e-03_rb, 5.28020e-03_rb, 4.68307e-03_rb, 4.15348e-03_rb, &
     & 3.68377e-03_rb, 3.26718e-03_rb, 2.89770e-03_rb, 2.57000e-03_rb, 2.27937e-03_rb/)
      selfrefo(:, 6) = (/ &
     & 6.29140e-03_rb, 5.55557e-03_rb, 4.90580e-03_rb, 4.33203e-03_rb, 3.82536e-03_rb, &
     & 3.37795e-03_rb, 2.98287e-03_rb, 2.63400e-03_rb, 2.32593e-03_rb, 2.05389e-03_rb/)
      selfrefo(:, 7) = (/ &
     & 6.00229e-03_rb, 5.28180e-03_rb, 4.64780e-03_rb, 4.08990e-03_rb, 3.59897e-03_rb, &
     & 3.16696e-03_rb, 2.78682e-03_rb, 2.45230e-03_rb, 2.15794e-03_rb, 1.89891e-03_rb/)
      selfrefo(:, 8) = (/ &
     & 5.78892e-03_rb, 5.07191e-03_rb, 4.44370e-03_rb, 3.89330e-03_rb, 3.41108e-03_rb, &
     & 2.98858e-03_rb, 2.61842e-03_rb, 2.29410e-03_rb, 2.00995e-03_rb, 1.76100e-03_rb/)
      selfrefo(:, 9) = (/ &
     & 4.96186e-03_rb, 4.56767e-03_rb, 4.20480e-03_rb, 3.87076e-03_rb, 3.56325e-03_rb, &
     & 3.28017e-03_rb, 3.01959e-03_rb, 2.77970e-03_rb, 2.55887e-03_rb, 2.35559e-03_rb/)
      selfrefo(:,10) = (/ &
     & 4.56849e-03_rb, 4.35527e-03_rb, 4.15200e-03_rb, 3.95822e-03_rb, 3.77348e-03_rb, &
     & 3.59736e-03_rb, 3.42946e-03_rb, 3.26940e-03_rb, 3.11681e-03_rb, 2.97134e-03_rb/)
      selfrefo(:,11) = (/ &
     & 4.47310e-03_rb, 4.32453e-03_rb, 4.18090e-03_rb, 4.04204e-03_rb, 3.90779e-03_rb, &
     & 3.77799e-03_rb, 3.65251e-03_rb, 3.53120e-03_rb, 3.41392e-03_rb, 3.30053e-03_rb/)
      selfrefo(:,12) = (/ &
     & 4.46459e-03_rb, 4.24031e-03_rb, 4.02730e-03_rb, 3.82499e-03_rb, 3.63284e-03_rb, &
     & 3.45035e-03_rb, 3.27702e-03_rb, 3.11240e-03_rb, 2.95605e-03_rb, 2.80755e-03_rb/)
      selfrefo(:,13) = (/ &
     & 4.43961e-03_rb, 4.35658e-03_rb, 4.27510e-03_rb, 4.19514e-03_rb, 4.11669e-03_rb, &
     & 4.03969e-03_rb, 3.96414e-03_rb, 3.89000e-03_rb, 3.81725e-03_rb, 3.74585e-03_rb/)
      selfrefo(:,14) = (/ &
     & 4.40512e-03_rb, 4.41515e-03_rb, 4.42520e-03_rb, 4.43527e-03_rb, 4.44537e-03_rb, &
     & 4.45549e-03_rb, 4.46563e-03_rb, 4.47580e-03_rb, 4.48599e-03_rb, 4.49620e-03_rb/)
      selfrefo(:,15) = (/ &
     & 3.21965e-03_rb, 3.42479e-03_rb, 3.64300e-03_rb, 3.87512e-03_rb, 4.12202e-03_rb, &
     & 4.38466e-03_rb, 4.66403e-03_rb, 4.96120e-03_rb, 5.27731e-03_rb, 5.61355e-03_rb/)
      selfrefo(:,16) = (/ &
     & 3.11402e-03_rb, 3.35870e-03_rb, 3.62260e-03_rb, 3.90724e-03_rb, 4.21424e-03_rb, &
     & 4.54536e-03_rb, 4.90250e-03_rb, 5.28770e-03_rb, 5.70317e-03_rb, 6.15128e-03_rb/)

      end subroutine lw_kgb138

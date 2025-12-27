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
!!prev      subroutine sw_kgb22
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine sw_kgb227
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg22, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg22, only : selfrefo, forrefo

      implicit none
      save

!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.351362e-07_rb, 0.341136e-07_rb, 0.181317e-06_rb /)
      forrefo(:, 2) = (/ 0.109648e-06_rb, 0.344240e-06_rb, 0.139709e-05_rb /)
      forrefo(:, 3) = (/ 0.374823e-06_rb, 0.103424e-05_rb, 0.188717e-05_rb /)
      forrefo(:, 4) = (/ 0.580041e-06_rb, 0.116876e-05_rb, 0.121183e-05_rb /)
      forrefo(:, 5) = (/ 0.115608e-05_rb, 0.148110e-05_rb, 0.836083e-06_rb /)
      forrefo(:, 6) = (/ 0.181460e-05_rb, 0.133313e-05_rb, 0.500167e-06_rb /)
      forrefo(:, 7) = (/ 0.199096e-05_rb, 0.115276e-05_rb, 0.432994e-06_rb /)
      forrefo(:, 8) = (/ 0.183730e-05_rb, 0.122260e-05_rb, 0.433248e-06_rb /)
      forrefo(:, 9) = (/ 0.198386e-05_rb, 0.100130e-05_rb, 0.269712e-06_rb /)
      forrefo(:,10) = (/ 0.276382e-05_rb, 0.749215e-06_rb, 0.236919e-06_rb /)
      forrefo(:,11) = (/ 0.298202e-05_rb, 0.629688e-06_rb, 0.228388e-06_rb /)
      forrefo(:,12) = (/ 0.364604e-05_rb, 0.455336e-06_rb, 0.206130e-06_rb /)
      forrefo(:,13) = (/ 0.373339e-05_rb, 0.245210e-06_rb, 0.201987e-06_rb /)
      forrefo(:,14) = (/ 0.480378e-05_rb, 0.177591e-06_rb, 0.171458e-06_rb /)
      forrefo(:,15) = (/ 0.521700e-05_rb, 0.203358e-06_rb, 0.189559e-06_rb /)
      forrefo(:,16) = (/ 0.542717e-05_rb, 0.219022e-06_rb, 0.218271e-06_rb /)

!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
        &  0.538526e-04_rb, 0.464603e-04_rb, 0.400828e-04_rb, 0.345807e-04_rb, 0.298339e-04_rb, &
        &  0.257386e-04_rb, 0.222055e-04_rb, 0.191574e-04_rb, 0.165277e-04_rb, 0.142590e-04_rb /)
      selfrefo(:, 2) = (/ &
        &  0.162409e-03_rb, 0.128347e-03_rb, 0.101430e-03_rb, 0.801571e-04_rb, 0.633460e-04_rb, &
        &  0.500607e-04_rb, 0.395616e-04_rb, 0.312645e-04_rb, 0.247075e-04_rb, 0.195257e-04_rb /)
      selfrefo(:, 3) = (/ &
        &  0.262882e-03_rb, 0.212793e-03_rb, 0.172247e-03_rb, 0.139427e-03_rb, 0.112860e-03_rb, &
        &  0.913557e-04_rb, 0.739487e-04_rb, 0.598584e-04_rb, 0.484529e-04_rb, 0.392206e-04_rb /)
      selfrefo(:, 4) = (/ &
        &  0.242873e-03_rb, 0.204225e-03_rb, 0.171726e-03_rb, 0.144399e-03_rb, 0.121421e-03_rb, &
        &  0.102099e-03_rb, 0.858516e-04_rb, 0.721899e-04_rb, 0.607022e-04_rb, 0.510426e-04_rb /)
      selfrefo(:, 5) = (/ &
        &  0.235614e-03_rb, 0.207814e-03_rb, 0.183293e-03_rb, 0.161666e-03_rb, 0.142591e-03_rb, &
        &  0.125766e-03_rb, 0.110927e-03_rb, 0.978381e-04_rb, 0.862939e-04_rb, 0.761119e-04_rb /)
      selfrefo(:, 6) = (/ &
        &  0.205508e-03_rb, 0.190174e-03_rb, 0.175985e-03_rb, 0.162854e-03_rb, 0.150702e-03_rb, &
        &  0.139458e-03_rb, 0.129052e-03_rb, 0.119423e-03_rb, 0.110513e-03_rb, 0.102267e-03_rb /)
      selfrefo(:, 7) = (/ &
        &  0.185027e-03_rb, 0.175148e-03_rb, 0.165796e-03_rb, 0.156944e-03_rb, 0.148565e-03_rb, &
        &  0.140633e-03_rb, 0.133124e-03_rb, 0.126016e-03_rb, 0.119288e-03_rb, 0.112919e-03_rb /)
      selfrefo(:, 8) = (/ &
        &  0.192634e-03_rb, 0.180192e-03_rb, 0.168554e-03_rb, 0.157668e-03_rb, 0.147484e-03_rb, &
        &  0.137959e-03_rb, 0.129048e-03_rb, 0.120713e-03_rb, 0.112917e-03_rb, 0.105624e-03_rb /)
      selfrefo(:, 9) = (/ &
        &  0.161632e-03_rb, 0.155919e-03_rb, 0.150408e-03_rb, 0.145092e-03_rb, 0.139963e-03_rb, &
        &  0.135016e-03_rb, 0.130244e-03_rb, 0.125640e-03_rb, 0.121199e-03_rb, 0.116915e-03_rb /)
      selfrefo(:,10) = (/ &
        &  0.120880e-03_rb, 0.125265e-03_rb, 0.129810e-03_rb, 0.134520e-03_rb, 0.139400e-03_rb, &
        &  0.144458e-03_rb, 0.149699e-03_rb, 0.155130e-03_rb, 0.160758e-03_rb, 0.166591e-03_rb /)
      selfrefo(:,11) = (/ &
        &  0.104705e-03_rb, 0.111761e-03_rb, 0.119291e-03_rb, 0.127330e-03_rb, 0.135910e-03_rb, &
        &  0.145068e-03_rb, 0.154843e-03_rb, 0.165277e-03_rb, 0.176414e-03_rb, 0.188302e-03_rb /)
      selfrefo(:,12) = (/ &
        &  0.846335e-04_rb, 0.951236e-04_rb, 0.106914e-03_rb, 0.120166e-03_rb, 0.135060e-03_rb, &
        &  0.151800e-03_rb, 0.170616e-03_rb, 0.191763e-03_rb, 0.215532e-03_rb, 0.242246e-03_rb /)
      selfrefo(:,13) = (/ &
        &  0.669754e-04_rb, 0.781902e-04_rb, 0.912829e-04_rb, 0.106568e-03_rb, 0.124413e-03_rb, &
        &  0.145245e-03_rb, 0.169566e-03_rb, 0.197959e-03_rb, 0.231107e-03_rb, 0.269805e-03_rb /)
      selfrefo(:,14) = (/ &
        &  0.597091e-04_rb, 0.722265e-04_rb, 0.873679e-04_rb, 0.105684e-03_rb, 0.127839e-03_rb, &
        &  0.154639e-03_rb, 0.187057e-03_rb, 0.226272e-03_rb, 0.273707e-03_rb, 0.331087e-03_rb /)
      selfrefo(:,15) = (/ &
        &  0.640410e-04_rb, 0.771879e-04_rb, 0.930338e-04_rb, 0.112133e-03_rb, 0.135152e-03_rb, &
        &  0.162897e-03_rb, 0.196338e-03_rb, 0.236644e-03_rb, 0.285225e-03_rb, 0.343778e-03_rb /)
      selfrefo(:,16) = (/ &
        &  0.666420e-04_rb, 0.801056e-04_rb, 0.962892e-04_rb, 0.115742e-03_rb, 0.139126e-03_rb, &
        &  0.167233e-03_rb, 0.201019e-03_rb, 0.241630e-03_rb, 0.290446e-03_rb, 0.349125e-03_rb /)
     
      end subroutine sw_kgb227

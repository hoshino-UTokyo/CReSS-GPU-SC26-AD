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
!!prev      subroutine sw_kgb23
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for  selfrefo(:,:), forrefo(:,:)
!!
      subroutine sw_kgb232
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
      use rrsw_kg23, only : selfrefo, forrefo

      implicit none
      save

!     -----------------------------------------------------------------

      forrefo(:, 1) = (/ 0.315770e-07_rb, 0.671978e-07_rb, 0.440649e-06_rb /)
      forrefo(:, 2) = (/ 0.313674e-06_rb, 0.285252e-06_rb, 0.421024e-05_rb /)
      forrefo(:, 3) = (/ 0.135818e-05_rb, 0.145071e-05_rb, 0.611285e-05_rb /)
      forrefo(:, 4) = (/ 0.534065e-05_rb, 0.586268e-05_rb, 0.933970e-05_rb /)
      forrefo(:, 5) = (/ 0.964007e-05_rb, 0.107110e-04_rb, 0.104486e-04_rb /)
      forrefo(:, 6) = (/ 0.302775e-04_rb, 0.357530e-04_rb, 0.340724e-04_rb /)
      forrefo(:, 7) = (/ 0.102437e-03_rb, 0.108475e-03_rb, 0.105245e-03_rb /)
      forrefo(:, 8) = (/ 0.146054e-03_rb, 0.141490e-03_rb, 0.133071e-03_rb /)
      forrefo(:, 9) = (/ 0.163978e-03_rb, 0.150208e-03_rb, 0.142864e-03_rb /)
      forrefo(:,10) = (/ 0.220412e-03_rb, 0.182943e-03_rb, 0.150941e-03_rb /)
      forrefo(:,11) = (/ 0.228877e-03_rb, 0.197679e-03_rb, 0.163220e-03_rb /)
      forrefo(:,12) = (/ 0.234177e-03_rb, 0.217734e-03_rb, 0.185038e-03_rb /)
      forrefo(:,13) = (/ 0.257187e-03_rb, 0.241570e-03_rb, 0.221178e-03_rb /)
      forrefo(:,14) = (/ 0.272455e-03_rb, 0.270637e-03_rb, 0.256269e-03_rb /)
      forrefo(:,15) = (/ 0.339445e-03_rb, 0.300268e-03_rb, 0.286574e-03_rb /)
      forrefo(:,16) = (/ 0.338841e-03_rb, 0.355428e-03_rb, 0.353794e-03_rb /)

!     -----------------------------------------------------------------
!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
        &  0.100945e-04_rb, 0.801113e-05_rb, 0.635771e-05_rb, 0.504554e-05_rb, 0.400419e-05_rb, &
        &  0.317777e-05_rb, 0.252191e-05_rb, 0.200141e-05_rb, 0.158834e-05_rb, 0.126052e-05_rb /)
      selfrefo(:, 2) = (/ &
        &  0.107573e-04_rb, 0.999809e-05_rb, 0.929245e-05_rb, 0.863661e-05_rb, 0.802706e-05_rb, &
        &  0.746053e-05_rb, 0.693399e-05_rb, 0.644460e-05_rb, 0.598976e-05_rb, 0.556702e-05_rb /)
      selfrefo(:, 3) = (/ &
        &  0.350389e-04_rb, 0.319234e-04_rb, 0.290850e-04_rb, 0.264989e-04_rb, 0.241428e-04_rb, &
        &  0.219962e-04_rb, 0.200404e-04_rb, 0.182586e-04_rb, 0.166351e-04_rb, 0.151560e-04_rb /)
      selfrefo(:, 4) = (/ &
        &  0.122993e-03_rb, 0.110885e-03_rb, 0.999691e-04_rb, 0.901277e-04_rb, 0.812551e-04_rb, &
        &  0.732559e-04_rb, 0.660443e-04_rb, 0.595426e-04_rb, 0.536809e-04_rb, 0.483963e-04_rb /)
      selfrefo(:, 5) = (/ &
        &  0.206434e-03_rb, 0.187435e-03_rb, 0.170185e-03_rb, 0.154522e-03_rb, 0.140301e-03_rb, &
        &  0.127388e-03_rb, 0.115664e-03_rb, 0.105019e-03_rb, 0.953540e-04_rb, 0.865783e-04_rb /)
      selfrefo(:, 6) = (/ &
        &  0.590645e-03_rb, 0.533109e-03_rb, 0.481177e-03_rb, 0.434305e-03_rb, 0.391998e-03_rb, &
        &  0.353812e-03_rb, 0.319346e-03_rb, 0.288238e-03_rb, 0.260160e-03_rb, 0.234817e-03_rb /)
      selfrefo(:, 7) = (/ &
        &  0.163029e-02_rb, 0.148773e-02_rb, 0.135763e-02_rb, 0.123891e-02_rb, 0.113057e-02_rb, &
        &  0.103170e-02_rb, 0.941483e-03_rb, 0.859153e-03_rb, 0.784023e-03_rb, 0.715462e-03_rb /)
      selfrefo(:, 8) = (/ &
        &  0.204528e-02_rb, 0.189258e-02_rb, 0.175128e-02_rb, 0.162053e-02_rb, 0.149954e-02_rb, &
        &  0.138758e-02_rb, 0.128398e-02_rb, 0.118812e-02_rb, 0.109941e-02_rb, 0.101733e-02_rb /)
      selfrefo(:, 9) = (/ &
        &  0.210589e-02_rb, 0.197078e-02_rb, 0.184434e-02_rb, 0.172601e-02_rb, 0.161528e-02_rb, &
        &  0.151164e-02_rb, 0.141466e-02_rb, 0.132390e-02_rb, 0.123896e-02_rb, 0.115947e-02_rb /)
      selfrefo(:,10) = (/ &
        &  0.245098e-02_rb, 0.233745e-02_rb, 0.222918e-02_rb, 0.212592e-02_rb, 0.202745e-02_rb, &
        &  0.193353e-02_rb, 0.184397e-02_rb, 0.175856e-02_rb, 0.167710e-02_rb, 0.159941e-02_rb /)
      selfrefo(:,11) = (/ &
        &  0.267460e-02_rb, 0.253325e-02_rb, 0.239936e-02_rb, 0.227255e-02_rb, 0.215244e-02_rb, &
        &  0.203868e-02_rb, 0.193093e-02_rb, 0.182888e-02_rb, 0.173222e-02_rb, 0.164067e-02_rb /)
      selfrefo(:,12) = (/ &
        &  0.304510e-02_rb, 0.283919e-02_rb, 0.264720e-02_rb, 0.246820e-02_rb, 0.230130e-02_rb, &
        &  0.214568e-02_rb, 0.200059e-02_rb, 0.186531e-02_rb, 0.173918e-02_rb, 0.162157e-02_rb /)
      selfrefo(:,13) = (/ &
        &  0.338445e-02_rb, 0.314719e-02_rb, 0.292655e-02_rb, 0.272139e-02_rb, 0.253060e-02_rb, &
        &  0.235319e-02_rb, 0.218822e-02_rb, 0.203482e-02_rb, 0.189217e-02_rb, 0.175952e-02_rb /)
      selfrefo(:,14) = (/ &
        &  0.388649e-02_rb, 0.357018e-02_rb, 0.327961e-02_rb, 0.301269e-02_rb, 0.276750e-02_rb, &
        &  0.254226e-02_rb, 0.233535e-02_rb, 0.214528e-02_rb, 0.197068e-02_rb, 0.181029e-02_rb /)
      selfrefo(:,15) = (/ &
        &  0.412547e-02_rb, 0.387413e-02_rb, 0.363810e-02_rb, 0.341646e-02_rb, 0.320831e-02_rb, &
        &  0.301285e-02_rb, 0.282930e-02_rb, 0.265693e-02_rb, 0.249506e-02_rb, 0.234305e-02_rb /)
      selfrefo(:,16) = (/ &
        &  0.534327e-02_rb, 0.482967e-02_rb, 0.436544e-02_rb, 0.394583e-02_rb, 0.356655e-02_rb, &
        &  0.322373e-02_rb, 0.291387e-02_rb, 0.263378e-02_rb, 0.238062e-02_rb, 0.215179e-02_rb /)
     
      end subroutine sw_kgb232

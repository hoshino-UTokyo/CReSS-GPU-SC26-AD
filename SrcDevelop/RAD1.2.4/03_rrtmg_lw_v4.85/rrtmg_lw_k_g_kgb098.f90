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
!!prev      subroutine lw_kgb09
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for kbo_mn2o(:,:,:)
!!
      subroutine lw_kgb098
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg09, only : fracrefao, fracrefbo, kao, kbo, kao_mn2o, &
!!prev                            kbo_mn2o, selfrefo, forrefo
      use rrlw_kg09, only : kbo_mn2o

      implicit none
      save

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

      kbo_mn2o(:, 1) = (/ &
     & 8.42688e-03_rb, 8.96787e-03_rb, 9.54358e-03_rb, 1.01563e-02_rb, 1.08083e-02_rb, &
     & 1.15021e-02_rb, 1.22405e-02_rb, 1.30263e-02_rb, 1.38626e-02_rb, 1.47525e-02_rb, &
     & 1.56996e-02_rb, 1.67075e-02_rb, 1.77800e-02_rb, 1.89215e-02_rb, 2.01362e-02_rb, &
     & 2.14289e-02_rb, 2.28045e-02_rb, 2.42685e-02_rb, 2.58265e-02_rb/)
      kbo_mn2o(:, 2) = (/ &
     & 2.24976e-02_rb, 2.38935e-02_rb, 2.53762e-02_rb, 2.69508e-02_rb, 2.86231e-02_rb, &
     & 3.03991e-02_rb, 3.22854e-02_rb, 3.42887e-02_rb, 3.64163e-02_rb, 3.86760e-02_rb, &
     & 4.10759e-02_rb, 4.36246e-02_rb, 4.63315e-02_rb, 4.92064e-02_rb, 5.22597e-02_rb, &
     & 5.55024e-02_rb, 5.89464e-02_rb, 6.26040e-02_rb, 6.64886e-02_rb/)
      kbo_mn2o(:, 3) = (/ &
     & 5.93542e-02_rb, 6.37312e-02_rb, 6.84310e-02_rb, 7.34774e-02_rb, 7.88960e-02_rb, &
     & 8.47141e-02_rb, 9.09613e-02_rb, 9.76692e-02_rb, 1.04872e-01_rb, 1.12605e-01_rb, &
     & 1.20910e-01_rb, 1.29826e-01_rb, 1.39400e-01_rb, 1.49680e-01_rb, 1.60718e-01_rb, &
     & 1.72570e-01_rb, 1.85296e-01_rb, 1.98961e-01_rb, 2.13633e-01_rb/)
      kbo_mn2o(:, 4) = (/ &
     & 1.98022e-01_rb, 2.05895e-01_rb, 2.14082e-01_rb, 2.22594e-01_rb, 2.31445e-01_rb, &
     & 2.40647e-01_rb, 2.50216e-01_rb, 2.60164e-01_rb, 2.70509e-01_rb, 2.81265e-01_rb, &
     & 2.92448e-01_rb, 3.04076e-01_rb, 3.16167e-01_rb, 3.28738e-01_rb, 3.41809e-01_rb, &
     & 3.55400e-01_rb, 3.69531e-01_rb, 3.84224e-01_rb, 3.99501e-01_rb/)
      kbo_mn2o(:, 5) = (/ &
     & 6.41413e-01_rb, 6.46239e-01_rb, 6.51101e-01_rb, 6.56000e-01_rb, 6.60936e-01_rb, &
     & 6.65910e-01_rb, 6.70920e-01_rb, 6.75968e-01_rb, 6.81054e-01_rb, 6.86179e-01_rb, &
     & 6.91342e-01_rb, 6.96544e-01_rb, 7.01785e-01_rb, 7.07065e-01_rb, 7.12385e-01_rb, &
     & 7.17746e-01_rb, 7.23146e-01_rb, 7.28587e-01_rb, 7.34070e-01_rb/)
      kbo_mn2o(:, 6) = (/ &
     & 1.47906e+00_rb, 1.48768e+00_rb, 1.49635e+00_rb, 1.50507e+00_rb, 1.51384e+00_rb, &
     & 1.52267e+00_rb, 1.53154e+00_rb, 1.54047e+00_rb, 1.54944e+00_rb, 1.55847e+00_rb, &
     & 1.56755e+00_rb, 1.57669e+00_rb, 1.58588e+00_rb, 1.59512e+00_rb, 1.60442e+00_rb, &
     & 1.61377e+00_rb, 1.62317e+00_rb, 1.63263e+00_rb, 1.64215e+00_rb/)
      kbo_mn2o(:, 7) = (/ &
     & 3.53152e+00_rb, 3.55492e+00_rb, 3.57848e+00_rb, 3.60219e+00_rb, 3.62606e+00_rb, &
     & 3.65008e+00_rb, 3.67427e+00_rb, 3.69862e+00_rb, 3.72313e+00_rb, 3.74780e+00_rb, &
     & 3.77263e+00_rb, 3.79763e+00_rb, 3.82279e+00_rb, 3.84812e+00_rb, 3.87362e+00_rb, &
     & 3.89929e+00_rb, 3.92513e+00_rb, 3.95114e+00_rb, 3.97732e+00_rb/)
      kbo_mn2o(:, 8) = (/ &
     & 9.06783e+00_rb, 9.04597e+00_rb, 9.02415e+00_rb, 9.00239e+00_rb, 8.98069e+00_rb, &
     & 8.95903e+00_rb, 8.93743e+00_rb, 8.91588e+00_rb, 8.89438e+00_rb, 8.87293e+00_rb, &
     & 8.85154e+00_rb, 8.83020e+00_rb, 8.80890e+00_rb, 8.78766e+00_rb, 8.76647e+00_rb, &
     & 8.74533e+00_rb, 8.72425e+00_rb, 8.70321e+00_rb, 8.68223e+00_rb/)
      kbo_mn2o(:, 9) = (/ &
     & 3.88220e+01_rb, 3.85805e+01_rb, 3.83405e+01_rb, 3.81019e+01_rb, 3.78649e+01_rb, &
     & 3.76293e+01_rb, 3.73952e+01_rb, 3.71625e+01_rb, 3.69313e+01_rb, 3.67016e+01_rb, &
     & 3.64732e+01_rb, 3.62463e+01_rb, 3.60208e+01_rb, 3.57967e+01_rb, 3.55740e+01_rb, &
     & 3.53527e+01_rb, 3.51327e+01_rb, 3.49142e+01_rb, 3.46970e+01_rb/)
      kbo_mn2o(:, 10) = (/ &
     & 1.14211e+02_rb, 1.13955e+02_rb, 1.13700e+02_rb, 1.13445e+02_rb, 1.13191e+02_rb, &
     & 1.12938e+02_rb, 1.12685e+02_rb, 1.12433e+02_rb, 1.12181e+02_rb, 1.11930e+02_rb, &
     & 1.11679e+02_rb, 1.11429e+02_rb, 1.11180e+02_rb, 1.10931e+02_rb, 1.10682e+02_rb, &
     & 1.10434e+02_rb, 1.10187e+02_rb, 1.09940e+02_rb, 1.09694e+02_rb/)
      kbo_mn2o(:, 11) = (/ &
     & 1.60513e+02_rb, 1.60857e+02_rb, 1.61201e+02_rb, 1.61547e+02_rb, 1.61893e+02_rb, &
     & 1.62240e+02_rb, 1.62587e+02_rb, 1.62936e+02_rb, 1.63285e+02_rb, 1.63635e+02_rb, &
     & 1.63985e+02_rb, 1.64337e+02_rb, 1.64689e+02_rb, 1.65041e+02_rb, 1.65395e+02_rb, &
     & 1.65749e+02_rb, 1.66105e+02_rb, 1.66460e+02_rb, 1.66817e+02_rb/)
      kbo_mn2o(:, 12) = (/ &
     & 1.71473e+02_rb, 1.72766e+02_rb, 1.74068e+02_rb, 1.75381e+02_rb, 1.76703e+02_rb, &
     & 1.78035e+02_rb, 1.79377e+02_rb, 1.80729e+02_rb, 1.82091e+02_rb, 1.83464e+02_rb, &
     & 1.84847e+02_rb, 1.86240e+02_rb, 1.87644e+02_rb, 1.89059e+02_rb, 1.90484e+02_rb, &
     & 1.91920e+02_rb, 1.93367e+02_rb, 1.94824e+02_rb, 1.96293e+02_rb/)
      kbo_mn2o(:, 13) = (/ &
     & 2.71287e+01_rb, 2.75538e+01_rb, 2.79856e+01_rb, 2.84241e+01_rb, 2.88695e+01_rb, &
     & 2.93219e+01_rb, 2.97814e+01_rb, 3.02480e+01_rb, 3.07220e+01_rb, 3.12035e+01_rb, &
     & 3.16924e+01_rb, 3.21890e+01_rb, 3.26934e+01_rb, 3.32058e+01_rb, 3.37261e+01_rb, &
     & 3.42546e+01_rb, 3.47914e+01_rb, 3.53365e+01_rb, 3.58903e+01_rb/)
      kbo_mn2o(:, 14) = (/ &
     & 1.70389e+01_rb, 1.70899e+01_rb, 1.71411e+01_rb, 1.71924e+01_rb, 1.72439e+01_rb, &
     & 1.72955e+01_rb, 1.73473e+01_rb, 1.73992e+01_rb, 1.74513e+01_rb, 1.75035e+01_rb, &
     & 1.75559e+01_rb, 1.76085e+01_rb, 1.76612e+01_rb, 1.77141e+01_rb, 1.77671e+01_rb, &
     & 1.78203e+01_rb, 1.78736e+01_rb, 1.79271e+01_rb, 1.79808e+01_rb/)
      kbo_mn2o(:, 15) = (/ &
     & 2.49725e+00_rb, 2.66861e+00_rb, 2.85174e+00_rb, 3.04743e+00_rb, 3.25655e+00_rb, &
     & 3.48003e+00_rb, 3.71883e+00_rb, 3.97403e+00_rb, 4.24673e+00_rb, 4.53815e+00_rb, &
     & 4.84957e+00_rb, 5.18236e+00_rb, 5.53798e+00_rb, 5.91801e+00_rb, 6.32412e+00_rb, &
     & 6.75809e+00_rb, 7.22185e+00_rb, 7.71742e+00_rb, 8.24701e+00_rb/)
      kbo_mn2o(:, 16) = (/ &
     & 1.82935e-03_rb, 2.58912e-03_rb, 3.66444e-03_rb, 5.18637e-03_rb, 7.34039e-03_rb, &
     & 1.03890e-02_rb, 1.47038e-02_rb, 2.08106e-02_rb, 2.94538e-02_rb, 4.16865e-02_rb, &
     & 5.89999e-02_rb, 8.35040e-02_rb, 1.18185e-01_rb, 1.67270e-01_rb, 2.36741e-01_rb, &
     & 3.35065e-01_rb, 4.74225e-01_rb, 6.71180e-01_rb, 9.49936e-01_rb/)

      end subroutine lw_kgb098

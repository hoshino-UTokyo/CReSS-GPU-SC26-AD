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
!!          set parameters for kbo_mo3(:,:)
!!
      subroutine lw_kgb137
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg13, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
!!prev                            kbo_mo3, selfrefo, forrefo
      use rrlw_kg13, only : kbo_mo3

      implicit none
      save

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

      kbo_mo3(:, 1) = (/ &
     & 1.07596e-02_rb, 1.12146e-02_rb, 1.16887e-02_rb, 1.21830e-02_rb, 1.26981e-02_rb, &
     & 1.32350e-02_rb, 1.37946e-02_rb, 1.43779e-02_rb, 1.49858e-02_rb, 1.56194e-02_rb, &
     & 1.62799e-02_rb, 1.69682e-02_rb, 1.76857e-02_rb, 1.84334e-02_rb, 1.92129e-02_rb, &
     & 2.00252e-02_rb, 2.08719e-02_rb, 2.17544e-02_rb, 2.26743e-02_rb/)
      kbo_mo3(:, 2) = (/ &
     & 9.48276e-02_rb, 9.66591e-02_rb, 9.85260e-02_rb, 1.00429e-01_rb, 1.02369e-01_rb, &
     & 1.04346e-01_rb, 1.06361e-01_rb, 1.08416e-01_rb, 1.10510e-01_rb, 1.12644e-01_rb, &
     & 1.14820e-01_rb, 1.17037e-01_rb, 1.19298e-01_rb, 1.21602e-01_rb, 1.23951e-01_rb, &
     & 1.26345e-01_rb, 1.28785e-01_rb, 1.31273e-01_rb, 1.33808e-01_rb/)
      kbo_mo3(:, 3) = (/ &
     & 3.54721e-01_rb, 3.55779e-01_rb, 3.56841e-01_rb, 3.57906e-01_rb, 3.58973e-01_rb, &
     & 3.60044e-01_rb, 3.61119e-01_rb, 3.62196e-01_rb, 3.63277e-01_rb, 3.64360e-01_rb, &
     & 3.65448e-01_rb, 3.66538e-01_rb, 3.67631e-01_rb, 3.68728e-01_rb, 3.69828e-01_rb, &
     & 3.70932e-01_rb, 3.72038e-01_rb, 3.73148e-01_rb, 3.74262e-01_rb/)
      kbo_mo3(:, 4) = (/ &
     & 6.46454e-01_rb, 6.43823e-01_rb, 6.41202e-01_rb, 6.38593e-01_rb, 6.35994e-01_rb, &
     & 6.33405e-01_rb, 6.30827e-01_rb, 6.28260e-01_rb, 6.25703e-01_rb, 6.23156e-01_rb, &
     & 6.20620e-01_rb, 6.18094e-01_rb, 6.15578e-01_rb, 6.13073e-01_rb, 6.10578e-01_rb, &
     & 6.08093e-01_rb, 6.05618e-01_rb, 6.03153e-01_rb, 6.00698e-01_rb/)
      kbo_mo3(:, 5) = (/ &
     & 9.29832e-01_rb, 9.22877e-01_rb, 9.15975e-01_rb, 9.09124e-01_rb, 9.02324e-01_rb, &
     & 8.95576e-01_rb, 8.88877e-01_rb, 8.82229e-01_rb, 8.75631e-01_rb, 8.69082e-01_rb, &
     & 8.62582e-01_rb, 8.56130e-01_rb, 8.49727e-01_rb, 8.43372e-01_rb, 8.37064e-01_rb, &
     & 8.30803e-01_rb, 8.24589e-01_rb, 8.18422e-01_rb, 8.12301e-01_rb/)
      kbo_mo3(:, 6) = (/ &
     & 1.43531e+00_rb, 1.42616e+00_rb, 1.41706e+00_rb, 1.40802e+00_rb, 1.39903e+00_rb, &
     & 1.39010e+00_rb, 1.38124e+00_rb, 1.37242e+00_rb, 1.36367e+00_rb, 1.35496e+00_rb, &
     & 1.34632e+00_rb, 1.33773e+00_rb, 1.32919e+00_rb, 1.32071e+00_rb, 1.31229e+00_rb, &
     & 1.30391e+00_rb, 1.29559e+00_rb, 1.28733e+00_rb, 1.27911e+00_rb/)
      kbo_mo3(:, 7) = (/ &
     & 2.68664e+00_rb, 2.67196e+00_rb, 2.65736e+00_rb, 2.64284e+00_rb, 2.62840e+00_rb, &
     & 2.61404e+00_rb, 2.59975e+00_rb, 2.58555e+00_rb, 2.57142e+00_rb, 2.55737e+00_rb, &
     & 2.54340e+00_rb, 2.52950e+00_rb, 2.51568e+00_rb, 2.50193e+00_rb, 2.48826e+00_rb, &
     & 2.47466e+00_rb, 2.46114e+00_rb, 2.44769e+00_rb, 2.43432e+00_rb/)
      kbo_mo3(:, 8) = (/ &
     & 2.45343e+00_rb, 2.43442e+00_rb, 2.41556e+00_rb, 2.39684e+00_rb, 2.37827e+00_rb, &
     & 2.35984e+00_rb, 2.34156e+00_rb, 2.32342e+00_rb, 2.30541e+00_rb, 2.28755e+00_rb, &
     & 2.26983e+00_rb, 2.25224e+00_rb, 2.23479e+00_rb, 2.21747e+00_rb, 2.20029e+00_rb, &
     & 2.18324e+00_rb, 2.16633e+00_rb, 2.14954e+00_rb, 2.13289e+00_rb/)
      kbo_mo3(:, 9) = (/ &
     & 1.55879e-01_rb, 1.55998e-01_rb, 1.56118e-01_rb, 1.56238e-01_rb, 1.56358e-01_rb, &
     & 1.56478e-01_rb, 1.56599e-01_rb, 1.56719e-01_rb, 1.56840e-01_rb, 1.56960e-01_rb, &
     & 1.57081e-01_rb, 1.57201e-01_rb, 1.57322e-01_rb, 1.57443e-01_rb, 1.57564e-01_rb, &
     & 1.57685e-01_rb, 1.57806e-01_rb, 1.57928e-01_rb, 1.58049e-01_rb/)
      kbo_mo3(:,10) = (/ &
     & 8.75149e-03_rb, 8.88794e-03_rb, 9.02651e-03_rb, 9.16725e-03_rb, 9.31018e-03_rb, &
     & 9.45534e-03_rb, 9.60276e-03_rb, 9.75248e-03_rb, 9.90454e-03_rb, 1.00590e-02_rb, &
     & 1.02158e-02_rb, 1.03751e-02_rb, 1.05368e-02_rb, 1.07011e-02_rb, 1.08680e-02_rb, &
     & 1.10374e-02_rb, 1.12095e-02_rb, 1.13843e-02_rb, 1.15618e-02_rb/)
      kbo_mo3(:,11) = (/ &
     & 8.83874e-03_rb, 8.97926e-03_rb, 9.12201e-03_rb, 9.26703e-03_rb, 9.41436e-03_rb, &
     & 9.56403e-03_rb, 9.71608e-03_rb, 9.87055e-03_rb, 1.00275e-02_rb, 1.01869e-02_rb, &
     & 1.03488e-02_rb, 1.05134e-02_rb, 1.06805e-02_rb, 1.08503e-02_rb, 1.10228e-02_rb, &
     & 1.11980e-02_rb, 1.13761e-02_rb, 1.15569e-02_rb, 1.17407e-02_rb/)
      kbo_mo3(:,12) = (/ &
     & 9.59461e-03_rb, 9.70417e-03_rb, 9.81498e-03_rb, 9.92705e-03_rb, 1.00404e-02_rb, &
     & 1.01550e-02_rb, 1.02710e-02_rb, 1.03883e-02_rb, 1.05069e-02_rb, 1.06269e-02_rb, &
     & 1.07482e-02_rb, 1.08709e-02_rb, 1.09951e-02_rb, 1.11206e-02_rb, 1.12476e-02_rb, &
     & 1.13760e-02_rb, 1.15059e-02_rb, 1.16373e-02_rb, 1.17702e-02_rb/)
      kbo_mo3(:,13) = (/ &
     & 1.13077e-02_rb, 1.14079e-02_rb, 1.15089e-02_rb, 1.16109e-02_rb, 1.17138e-02_rb, &
     & 1.18176e-02_rb, 1.19223e-02_rb, 1.20279e-02_rb, 1.21344e-02_rb, 1.22419e-02_rb, &
     & 1.23504e-02_rb, 1.24598e-02_rb, 1.25702e-02_rb, 1.26816e-02_rb, 1.27939e-02_rb, &
     & 1.29073e-02_rb, 1.30216e-02_rb, 1.31370e-02_rb, 1.32534e-02_rb/)
      kbo_mo3(:,14) = (/ &
     & 6.74844e-03_rb, 6.82637e-03_rb, 6.90519e-03_rb, 6.98493e-03_rb, 7.06558e-03_rb, &
     & 7.14717e-03_rb, 7.22970e-03_rb, 7.31318e-03_rb, 7.39762e-03_rb, 7.48304e-03_rb, &
     & 7.56945e-03_rb, 7.65686e-03_rb, 7.74527e-03_rb, 7.83470e-03_rb, 7.92517e-03_rb, &
     & 8.01668e-03_rb, 8.10925e-03_rb, 8.20289e-03_rb, 8.29761e-03_rb/)
      kbo_mo3(:,15) = (/ &
     & 7.94595e-03_rb, 8.00015e-03_rb, 8.05472e-03_rb, 8.10966e-03_rb, 8.16497e-03_rb, &
     & 8.22067e-03_rb, 8.27674e-03_rb, 8.33320e-03_rb, 8.39004e-03_rb, 8.44727e-03_rb, &
     & 8.50489e-03_rb, 8.56290e-03_rb, 8.62130e-03_rb, 8.68011e-03_rb, 8.73932e-03_rb, &
     & 8.79893e-03_rb, 8.85895e-03_rb, 8.91937e-03_rb, 8.98021e-03_rb/)
      kbo_mo3(:,16) = (/ &
     & 1.85967e-03_rb, 1.86082e-03_rb, 1.86197e-03_rb, 1.86312e-03_rb, 1.86428e-03_rb, &
     & 1.86543e-03_rb, 1.86658e-03_rb, 1.86774e-03_rb, 1.86889e-03_rb, 1.87005e-03_rb, &
     & 1.87121e-03_rb, 1.87236e-03_rb, 1.87352e-03_rb, 1.87468e-03_rb, 1.87584e-03_rb, &
     & 1.87700e-03_rb, 1.87816e-03_rb, 1.87932e-03_rb, 1.88049e-03_rb/)

      end subroutine lw_kgb137

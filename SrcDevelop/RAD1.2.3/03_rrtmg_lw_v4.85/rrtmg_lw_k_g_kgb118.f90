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
!!prev      subroutine lw_kgb11
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for kbo_mo2(:,:)
!!
      subroutine lw_kgb118
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg11, only : fracrefao, fracrefbo, kao, kbo, kao_mo2, &
!!prev                            kbo_mo2, selfrefo, forrefo
      use rrlw_kg11, only :  kbo_mo2

      implicit none
      save

!     The array KBO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

      kbo_mo2(:, 1) = (/ &
     & 4.97626e-07_rb, 5.05955e-07_rb, 5.14424e-07_rb, 5.23034e-07_rb, 5.31789e-07_rb, &
     & 5.40690e-07_rb, 5.49739e-07_rb, 5.58941e-07_rb, 5.68296e-07_rb, 5.77808e-07_rb, &
     & 5.87479e-07_rb, 5.97312e-07_rb, 6.07310e-07_rb, 6.17475e-07_rb, 6.27810e-07_rb, &
     & 6.38318e-07_rb, 6.49002e-07_rb, 6.59865e-07_rb, 6.70910e-07_rb/)
      kbo_mo2(:, 2) = (/ &
     & 3.10232e-06_rb, 3.06339e-06_rb, 3.02496e-06_rb, 2.98700e-06_rb, 2.94952e-06_rb, &
     & 2.91252e-06_rb, 2.87597e-06_rb, 2.83989e-06_rb, 2.80426e-06_rb, 2.76907e-06_rb, &
     & 2.73433e-06_rb, 2.70002e-06_rb, 2.66614e-06_rb, 2.63269e-06_rb, 2.59966e-06_rb, &
     & 2.56704e-06_rb, 2.53483e-06_rb, 2.50303e-06_rb, 2.47162e-06_rb/)
      kbo_mo2(:, 3) = (/ &
     & 2.91635e-06_rb, 2.88637e-06_rb, 2.85669e-06_rb, 2.82733e-06_rb, 2.79826e-06_rb, &
     & 2.76949e-06_rb, 2.74102e-06_rb, 2.71284e-06_rb, 2.68495e-06_rb, 2.65735e-06_rb, &
     & 2.63003e-06_rb, 2.60299e-06_rb, 2.57623e-06_rb, 2.54975e-06_rb, 2.52353e-06_rb, &
     & 2.49759e-06_rb, 2.47191e-06_rb, 2.44650e-06_rb, 2.42135e-06_rb/)
      kbo_mo2(:, 4) = (/ &
     & 3.15584e-06_rb, 3.11986e-06_rb, 3.08430e-06_rb, 3.04914e-06_rb, 3.01438e-06_rb, &
     & 2.98002e-06_rb, 2.94605e-06_rb, 2.91247e-06_rb, 2.87927e-06_rb, 2.84645e-06_rb, &
     & 2.81400e-06_rb, 2.78192e-06_rb, 2.75021e-06_rb, 2.71886e-06_rb, 2.68787e-06_rb, &
     & 2.65723e-06_rb, 2.62694e-06_rb, 2.59699e-06_rb, 2.56739e-06_rb/)
      kbo_mo2(:, 5) = (/ &
     & 2.52067e-06_rb, 2.50127e-06_rb, 2.48202e-06_rb, 2.46291e-06_rb, 2.44396e-06_rb, &
     & 2.42515e-06_rb, 2.40648e-06_rb, 2.38796e-06_rb, 2.36958e-06_rb, 2.35134e-06_rb, &
     & 2.33324e-06_rb, 2.31529e-06_rb, 2.29747e-06_rb, 2.27978e-06_rb, 2.26224e-06_rb, &
     & 2.24482e-06_rb, 2.22755e-06_rb, 2.21040e-06_rb, 2.19339e-06_rb/)
      kbo_mo2(:, 6) = (/ &
     & 2.37304e-06_rb, 2.36340e-06_rb, 2.35380e-06_rb, 2.34423e-06_rb, 2.33471e-06_rb, &
     & 2.32522e-06_rb, 2.31578e-06_rb, 2.30637e-06_rb, 2.29700e-06_rb, 2.28766e-06_rb, &
     & 2.27837e-06_rb, 2.26911e-06_rb, 2.25989e-06_rb, 2.25071e-06_rb, 2.24157e-06_rb, &
     & 2.23246e-06_rb, 2.22339e-06_rb, 2.21436e-06_rb, 2.20536e-06_rb/)
      kbo_mo2(:, 7) = (/ &
     & 2.56366e-06_rb, 2.56395e-06_rb, 2.56424e-06_rb, 2.56453e-06_rb, 2.56482e-06_rb, &
     & 2.56510e-06_rb, 2.56539e-06_rb, 2.56568e-06_rb, 2.56597e-06_rb, 2.56625e-06_rb, &
     & 2.56654e-06_rb, 2.56683e-06_rb, 2.56712e-06_rb, 2.56741e-06_rb, 2.56769e-06_rb, &
     & 2.56798e-06_rb, 2.56827e-06_rb, 2.56856e-06_rb, 2.56885e-06_rb/)
      kbo_mo2(:, 8) = (/ &
     & 2.54502e-06_rb, 2.55393e-06_rb, 2.56287e-06_rb, 2.57185e-06_rb, 2.58085e-06_rb, &
     & 2.58989e-06_rb, 2.59896e-06_rb, 2.60806e-06_rb, 2.61719e-06_rb, 2.62636e-06_rb, &
     & 2.63555e-06_rb, 2.64478e-06_rb, 2.65404e-06_rb, 2.66334e-06_rb, 2.67266e-06_rb, &
     & 2.68202e-06_rb, 2.69141e-06_rb, 2.70084e-06_rb, 2.71030e-06_rb/)
      kbo_mo2(:, 9) = (/ &
     & 1.84106e-06_rb, 1.83922e-06_rb, 1.83737e-06_rb, 1.83553e-06_rb, 1.83369e-06_rb, &
     & 1.83186e-06_rb, 1.83002e-06_rb, 1.82819e-06_rb, 1.82636e-06_rb, 1.82453e-06_rb, &
     & 1.82270e-06_rb, 1.82087e-06_rb, 1.81905e-06_rb, 1.81723e-06_rb, 1.81541e-06_rb, &
     & 1.81359e-06_rb, 1.81177e-06_rb, 1.80996e-06_rb, 1.80814e-06_rb/)
      kbo_mo2(:,10) = (/ &
     & 1.83886e-06_rb, 1.83632e-06_rb, 1.83379e-06_rb, 1.83126e-06_rb, 1.82874e-06_rb, &
     & 1.82622e-06_rb, 1.82370e-06_rb, 1.82119e-06_rb, 1.81868e-06_rb, 1.81617e-06_rb, &
     & 1.81367e-06_rb, 1.81117e-06_rb, 1.80867e-06_rb, 1.80618e-06_rb, 1.80369e-06_rb, &
     & 1.80120e-06_rb, 1.79872e-06_rb, 1.79624e-06_rb, 1.79377e-06_rb/)
      kbo_mo2(:,11) = (/ &
     & 2.30390e-06_rb, 2.30269e-06_rb, 2.30148e-06_rb, 2.30028e-06_rb, 2.29907e-06_rb, &
     & 2.29787e-06_rb, 2.29667e-06_rb, 2.29546e-06_rb, 2.29426e-06_rb, 2.29306e-06_rb, &
     & 2.29186e-06_rb, 2.29066e-06_rb, 2.28946e-06_rb, 2.28826e-06_rb, 2.28706e-06_rb, &
     & 2.28586e-06_rb, 2.28466e-06_rb, 2.28347e-06_rb, 2.28227e-06_rb/)
      kbo_mo2(:,12) = (/ &
     & 2.38201e-06_rb, 2.36536e-06_rb, 2.34882e-06_rb, 2.33240e-06_rb, 2.31609e-06_rb, &
     & 2.29990e-06_rb, 2.28382e-06_rb, 2.26785e-06_rb, 2.25199e-06_rb, 2.23625e-06_rb, &
     & 2.22061e-06_rb, 2.20508e-06_rb, 2.18967e-06_rb, 2.17436e-06_rb, 2.15915e-06_rb, &
     & 2.14406e-06_rb, 2.12907e-06_rb, 2.11418e-06_rb, 2.09940e-06_rb/)
      kbo_mo2(:,13) = (/ &
     & 2.33326e-06_rb, 2.32549e-06_rb, 2.31775e-06_rb, 2.31003e-06_rb, 2.30234e-06_rb, &
     & 2.29467e-06_rb, 2.28703e-06_rb, 2.27941e-06_rb, 2.27182e-06_rb, 2.26426e-06_rb, &
     & 2.25672e-06_rb, 2.24920e-06_rb, 2.24171e-06_rb, 2.23424e-06_rb, 2.22680e-06_rb, &
     & 2.21939e-06_rb, 2.21200e-06_rb, 2.20463e-06_rb, 2.19729e-06_rb/)
      kbo_mo2(:,14) = (/ &
     & 2.75292e-06_rb, 2.75210e-06_rb, 2.75129e-06_rb, 2.75047e-06_rb, 2.74965e-06_rb, &
     & 2.74883e-06_rb, 2.74801e-06_rb, 2.74720e-06_rb, 2.74638e-06_rb, 2.74556e-06_rb, &
     & 2.74475e-06_rb, 2.74393e-06_rb, 2.74311e-06_rb, 2.74230e-06_rb, 2.74148e-06_rb, &
     & 2.74067e-06_rb, 2.73985e-06_rb, 2.73904e-06_rb, 2.73822e-06_rb/)
      kbo_mo2(:,15) = (/ &
     & 2.55262e-06_rb, 2.53364e-06_rb, 2.51480e-06_rb, 2.49611e-06_rb, 2.47755e-06_rb, &
     & 2.45913e-06_rb, 2.44084e-06_rb, 2.42269e-06_rb, 2.40468e-06_rb, 2.38680e-06_rb, &
     & 2.36906e-06_rb, 2.35144e-06_rb, 2.33396e-06_rb, 2.31660e-06_rb, 2.29938e-06_rb, &
     & 2.28228e-06_rb, 2.26531e-06_rb, 2.24847e-06_rb, 2.23175e-06_rb/)
      kbo_mo2(:,16) = (/ &
     & 3.11382e-06_rb, 3.08751e-06_rb, 3.06141e-06_rb, 3.03554e-06_rb, 3.00989e-06_rb, &
     & 2.98445e-06_rb, 2.95923e-06_rb, 2.93422e-06_rb, 2.90942e-06_rb, 2.88483e-06_rb, &
     & 2.86045e-06_rb, 2.83628e-06_rb, 2.81231e-06_rb, 2.78854e-06_rb, 2.76498e-06_rb, &
     & 2.74161e-06_rb, 2.71844e-06_rb, 2.69547e-06_rb, 2.67269e-06_rb/)

      end subroutine lw_kgb118

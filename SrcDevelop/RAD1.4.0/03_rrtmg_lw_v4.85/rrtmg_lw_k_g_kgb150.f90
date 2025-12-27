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
!!prev      subroutine lw_kgb15
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao
!!
      subroutine lw_kgb150
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg15, only : fracrefao, kao, kao_mn2, selfrefo, forrefo
      use rrlw_kg15, only : fracrefao

      implicit none
      save

! Planck fraction mapping level : P = 1053. mb, T = 294.2 K
      fracrefao(:, 1) = (/ &
     &  1.0689e-01_rb,1.1563e-01_rb,1.2447e-01_rb,1.2921e-01_rb,1.2840e-01_rb,1.2113e-01_rb, &
     &  1.0643e-01_rb,8.4987e-02_rb,6.0142e-02_rb,6.6798e-03_rb,5.5293e-03_rb,4.3700e-03_rb, &
     &  3.2061e-03_rb,2.0476e-03_rb,7.7366e-04_rb,1.0897e-04_rb/)
      fracrefao(:, 2) = (/ &
     &  1.0782e-01_rb,1.1637e-01_rb,1.2290e-01_rb,1.2911e-01_rb,1.2841e-01_rb,1.2113e-01_rb, &
     &  1.0643e-01_rb,8.4987e-02_rb,6.0142e-02_rb,6.6798e-03_rb,5.5293e-03_rb,4.3700e-03_rb, &
     &  3.2061e-03_rb,2.0476e-03_rb,7.7366e-04_rb,1.0897e-04_rb/)
      fracrefao(:, 3) = (/ &
     &  1.0858e-01_rb,1.1860e-01_rb,1.2237e-01_rb,1.2665e-01_rb,1.2841e-01_rb,1.2111e-01_rb, &
     &  1.0642e-01_rb,8.4987e-02_rb,6.0142e-02_rb,6.6798e-03_rb,5.5293e-03_rb,4.3700e-03_rb, &
     &  3.2061e-03_rb,2.0476e-03_rb,7.7366e-04_rb,1.0897e-04_rb/)
      fracrefao(:, 4) = (/ &
     &  1.1022e-01_rb,1.1965e-01_rb,1.2334e-01_rb,1.2383e-01_rb,1.2761e-01_rb,1.2109e-01_rb, &
     &  1.0642e-01_rb,8.4987e-02_rb,6.0142e-02_rb,6.6798e-03_rb,5.5293e-03_rb,4.3700e-03_rb, &
     &  3.2061e-03_rb,2.0476e-03_rb,7.7366e-04_rb,1.0897e-04_rb/)
      fracrefao(:, 5) = (/ &
     &  1.1342e-01_rb,1.2069e-01_rb,1.2360e-01_rb,1.2447e-01_rb,1.2340e-01_rb,1.2020e-01_rb, &
     &  1.0639e-01_rb,8.4987e-02_rb,6.0142e-02_rb,6.6798e-03_rb,5.5293e-03_rb,4.3700e-03_rb, &
     &  3.2061e-03_rb,2.0476e-03_rb,7.7366e-04_rb,1.0897e-04_rb/)
      fracrefao(:, 6) = (/ &
     &  1.1771e-01_rb,1.2280e-01_rb,1.2177e-01_rb,1.2672e-01_rb,1.2398e-01_rb,1.1787e-01_rb, &
     &  1.0131e-01_rb,8.4987e-02_rb,6.0142e-02_rb,6.6798e-03_rb,5.5293e-03_rb,4.3700e-03_rb, &
     &  3.2061e-03_rb,2.0476e-03_rb,7.7366e-04_rb,1.0897e-04_rb/)
      fracrefao(:, 7) = (/ &
     &  1.2320e-01_rb,1.2491e-01_rb,1.2001e-01_rb,1.2936e-01_rb,1.2653e-01_rb,1.1929e-01_rb, &
     &  9.8955e-02_rb,7.4887e-02_rb,6.0142e-02_rb,6.6798e-03_rb,5.5293e-03_rb,4.3700e-03_rb, &
     &  3.2061e-03_rb,2.0476e-03_rb,7.7366e-04_rb,1.0897e-04_rb/)
      fracrefao(:, 8) = (/ &
     &  1.3105e-01_rb,1.2563e-01_rb,1.3055e-01_rb,1.2854e-01_rb,1.3402e-01_rb,1.1571e-01_rb, &
     &  9.4876e-02_rb,6.0459e-02_rb,5.6457e-02_rb,6.6798e-03_rb,5.5293e-03_rb,4.3700e-03_rb, &
     &  3.2061e-03_rb,2.0476e-03_rb,7.7366e-04_rb,1.0897e-04_rb/)
      fracrefao(:, 9) = (/ &
     &  1.1375e-01_rb,1.2090e-01_rb,1.2348e-01_rb,1.2458e-01_rb,1.2406e-01_rb,1.1921e-01_rb, &
     &  1.0802e-01_rb,8.6613e-02_rb,5.8125e-02_rb,6.2984e-03_rb,5.2359e-03_rb,4.0641e-03_rb, &
     &  2.9379e-03_rb,1.9001e-03_rb,7.2646e-04_rb,1.0553e-04_rb/)

      end subroutine lw_kgb150

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
!!          set parameters for fracrefao(:,:), fracrefbo(:,:)
!!
      subroutine lw_kgb110
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg11, only : fracrefao, fracrefbo, kao, kbo, kao_mo2, &
!!prev                            kbo_mo2, selfrefo, forrefo
      use rrlw_kg11, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level : P=1053.63 mb, T= 294.2 K
      fracrefao(:) = (/ &
     &  1.4601e-01_rb,1.3824e-01_rb,1.4240e-01_rb,1.3463e-01_rb,1.1948e-01_rb,1.0440e-01_rb, &
     &  8.8667e-02_rb,6.5792e-02_rb,4.3893e-02_rb,4.7941e-03_rb,4.0760e-03_rb,3.3207e-03_rb, &
     &  2.4087e-03_rb,1.3912e-03_rb,4.3482e-04_rb,6.0932e-05_rb/)

! Planck fraction mapping level : P=0.353 mb, T = 262.11 K
      fracrefbo(:) = (/ &
     &  7.2928e-02_rb,1.4900e-01_rb,1.6156e-01_rb,1.5603e-01_rb,1.3934e-01_rb,1.1394e-01_rb, &
     &  8.8783e-02_rb,6.2411e-02_rb,4.0191e-02_rb,4.4587e-03_rb,3.9533e-03_rb,3.0847e-03_rb, &
     &  2.2317e-03_rb,1.4410e-03_rb,5.6722e-04_rb,7.7933e-05_rb/)

      end subroutine lw_kgb110

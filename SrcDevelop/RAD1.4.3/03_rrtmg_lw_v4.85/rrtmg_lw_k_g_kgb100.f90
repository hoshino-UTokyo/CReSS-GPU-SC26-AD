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
!!prev      subroutine lw_kgb10
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:,:), fracrefbo(:,:)
!!
      subroutine lw_kgb100
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg10, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg10, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level : P = 212.7250, T = 223.06 K
      fracrefao(:) = (/ &
     &  1.6909e-01_rb, 1.5419e-01_rb, 1.3999e-01_rb, 1.2637e-01_rb, &
     &  1.1429e-01_rb, 9.9676e-02_rb, 8.0093e-02_rb, 6.0283e-02_rb, &
     &  4.1077e-02_rb, 4.4857e-03_rb, 3.6545e-03_rb, 2.9243e-03_rb, &
     &  2.0407e-03_rb, 1.2891e-03_rb, 4.8767e-04_rb, 6.7748e-05_rb/)

! Planck fraction mapping level : P = 95.58350 mb, T = 215.70 K
      fracrefbo(:) = (/ &
     &  1.7391e-01_rb, 1.5680e-01_rb, 1.4419e-01_rb, 1.2672e-01_rb, &
     &  1.0708e-01_rb, 9.7034e-02_rb, 7.8545e-02_rb, 5.9784e-02_rb, &
     &  4.0879e-02_rb, 4.4704e-03_rb, 3.7150e-03_rb, 2.9038e-03_rb, &
     &  2.1454e-03_rb, 1.2802e-03_rb, 4.8328e-04_rb, 6.7378e-05_rb/)

      end subroutine lw_kgb100

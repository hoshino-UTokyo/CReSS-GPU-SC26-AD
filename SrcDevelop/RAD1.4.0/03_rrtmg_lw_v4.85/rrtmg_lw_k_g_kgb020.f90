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
!!prev      subroutine lw_kgb02
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:), fracrefbo(:)
!!
      subroutine lw_kgb020
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg02, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg02, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level: P = 1053.630 mbar, T = 294.2 K
      fracrefao(:) = (/ &
        1.6388e-01_rb, 1.5241e-01_rb, 1.4290e-01_rb, 1.2864e-01_rb, &
        1.1615e-01_rb, 1.0047e-01_rb, 8.0013e-02_rb, 6.0445e-02_rb, &
        4.0530e-02_rb, 4.3879e-03_rb, 3.5726e-03_rb, 2.7669e-03_rb, &
        2.0078e-03_rb, 1.2864e-03_rb, 4.7630e-04_rb, 6.9109e-05_rb/)

! Planck fraction mapping level: P = 3.206e-2 mb, T = 197.92 K
      fracrefbo(:) = (/ &
        1.4697e-01_rb, 1.4826e-01_rb, 1.4278e-01_rb, 1.3320e-01_rb, &
        1.1965e-01_rb, 1.0297e-01_rb, 8.4170e-02_rb, 6.3282e-02_rb, &
        4.2868e-02_rb, 4.6644e-03_rb, 3.8619e-03_rb, 3.0533e-03_rb, &
        2.2359e-03_rb, 1.4226e-03_rb, 5.3642e-04_rb, 7.6316e-05_rb/)

      end subroutine lw_kgb020

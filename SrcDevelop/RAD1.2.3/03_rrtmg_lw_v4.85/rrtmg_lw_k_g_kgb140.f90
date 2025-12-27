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
!!prev      subroutine lw_kgb14
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:), fracrefbo(:)
!!
      subroutine lw_kgb140
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg14, only : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo
      use rrlw_kg14, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level : P = 142.5940 mb, T = 215.70 K
      fracrefao(:) = (/ &
     &  1.9360e-01_rb, 1.7276e-01_rb, 1.4811e-01_rb, 1.2238e-01_rb, &
     &  1.0242e-01_rb, 8.6830e-02_rb, 7.1890e-02_rb, 5.4030e-02_rb, &
     &  3.5075e-02_rb, 3.8052e-03_rb, 3.1458e-03_rb, 2.4873e-03_rb, &
     &  1.8182e-03_rb, 1.1563e-03_rb, 4.3251e-04_rb, 5.7744e-05_rb/)

! Planck fraction mapping level : P = 4.758820mb, T = 250.85 K
      fracrefbo(:) = (/ &
     &  1.8599e-01_rb, 1.6646e-01_rb, 1.4264e-01_rb, 1.2231e-01_rb, &
     &  1.0603e-01_rb, 9.2014e-02_rb, 7.5287e-02_rb, 5.6758e-02_rb, &
     &  3.8386e-02_rb, 4.2139e-03_rb, 3.5399e-03_rb, 2.7381e-03_rb, &
     &  1.9202e-03_rb, 1.2083e-03_rb, 4.5395e-04_rb, 6.2699e-05_rb/)

      end subroutine lw_kgb140

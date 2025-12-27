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
!!prev      subroutine lw_kgb06
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:), cfc11adjo(:), cfc12o(:)
!!
      subroutine lw_kgb060
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg06, only : fracrefao, kao, kao_mco2, selfrefo, forrefo, &
!!prev                            cfc11adjo, cfc12o
      use rrlw_kg06, only : fracrefao, &
                            cfc11adjo, cfc12o

      implicit none
      save

! Planck fraction mapping level : P = 473.4280 mb, T = 259.83 K
      fracrefao(:) = (/ &
     &   1.4353e-01_rb,1.4774e-01_rb,1.4467e-01_rb,1.3785e-01_rb,1.2376e-01_rb,1.0214e-01_rb, &
     &   8.1984e-02_rb,6.1152e-02_rb,4.0987e-02_rb,4.5067e-03_rb,4.0020e-03_rb,3.1772e-03_rb, &
     &   2.3458e-03_rb,1.5025e-03_rb,5.7415e-04_rb,8.2970e-05_rb/)

! Minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12

!      cfc11(:) = (/ &
!     &    0.,      0.,      26.5435, 108.850, &
!     &    58.7804, 54.0875, 41.1065, 35.6120, &
!     &    41.2328, 47.7402, 79.1026, 64.3005, &
!     &    108.206, 141.617, 186.565, 58.4782/)
! Original cfc11 is multiplied by 1.385 to account for the 1060-1107 cm-1 band.

      cfc11adjo(:) = (/ &
     &   0._rb,      0._rb,      36.7627_rb, 150.757_rb, &
     &   81.4109_rb, 74.9112_rb, 56.9325_rb, 49.3226_rb, &
     &   57.1074_rb, 66.1202_rb, 109.557_rb, 89.0562_rb, & 
     &   149.865_rb, 196.140_rb, 258.393_rb, 80.9923_rb/)
      cfc12o(:) = (/ &
     &   62.8368_rb, 43.2626_rb, 26.7549_rb, 22.2487_rb, &
     &   23.5029_rb, 34.8323_rb, 26.2335_rb, 23.2306_rb, &
     &   18.4062_rb, 13.9534_rb, 22.6268_rb, 24.2604_rb, &
     &   30.0088_rb, 26.3634_rb, 15.8237_rb, 57.5050_rb/)

      end subroutine lw_kgb060

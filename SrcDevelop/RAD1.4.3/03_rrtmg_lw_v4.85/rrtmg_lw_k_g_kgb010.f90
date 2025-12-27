!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 21:04:30 $
!
!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
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
!!prev      subroutine lw_kgb01
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao, fracrefbo
!!
      subroutine lw_kgb010
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg01, only : fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
!!prev                            selfrefo, forrefo
      use rrlw_kg01, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level: P = 212.7250 mbar, T = 223.06 K
      fracrefao(:) = (/ &
        2.1227E-01_rb,1.8897E-01_rb,1.3934E-01_rb,1.1557E-01_rb,9.5282E-02_rb,8.3359E-02_rb, &
        6.5333E-02_rb,5.2016E-02_rb,3.4272E-02_rb,4.0257E-03_rb,3.1857E-03_rb,2.6014E-03_rb, &
        1.9141E-03_rb,1.2612E-03_rb,5.3169E-04_rb,7.6476E-05_rb/)

! Planck fraction mapping level: P = 212.7250 mbar, T = 223.06 K
! These Planck fractions were calculated using lower atmosphere
! parameters.
      fracrefbo(:) = (/ &
        2.1227E-01_rb,1.8897E-01_rb,1.3934E-01_rb,1.1557E-01_rb,9.5282E-02_rb,8.3359E-02_rb, &
        6.5333E-02_rb,5.2016E-02_rb,3.4272E-02_rb,4.0257E-03_rb,3.1857E-03_rb,2.6014E-03_rb, &
        1.9141E-03_rb,1.2612E-03_rb,5.3169E-04_rb,7.6476E-05_rb/)

      end subroutine lw_kgb010

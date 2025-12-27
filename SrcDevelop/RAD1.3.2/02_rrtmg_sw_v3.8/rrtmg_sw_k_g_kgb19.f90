!     path:      $Source:
!     /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 22:22:22 $

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
!      subroutine sw_kgbnn
! **************************************************************************
!  RRTM Shortwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original by J.Delamere, Atmospheric & Environmental Research.
!  Reformatted for F90: JJMorcrette, ECMWF
!  Further F90 and GCM revisions:  MJIacono, AER, July 2002
!
!  This file contains 14 subroutines that include the 
!  absorption coefficients and other data for each of the 14 shortwave
!  spectral bands used in RRTM_SW.  Here, the data are defined for 16
!  g-points, or sub-intervals, per band.  These data are combined and
!  weighted using a mapping procedure in routine RRTMG_SW_INIT to reduce
!  the total number of g-points from 224 to 112 for use in the GCM.
! **************************************************************************
      subroutine sw_kgb19
! **************************************************************************

      implicit none
      save

      call sw_kgb190
      call sw_kgb191
      call sw_kgb192
      call sw_kgb193
      call sw_kgb194
      call sw_kgb195
      call sw_kgb196
      call sw_kgb197
  
      end subroutine sw_kgb19

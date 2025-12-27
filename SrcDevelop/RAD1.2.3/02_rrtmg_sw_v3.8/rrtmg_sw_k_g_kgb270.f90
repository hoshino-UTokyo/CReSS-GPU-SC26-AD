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
!!prev      subroutine sw_kgb27
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for  sfluxrefo(:), absch4o(:), raylo(:)
!!
      subroutine sw_kgb270
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg27, only : kao, kbo, sfluxrefo, raylo
      use rrsw_kg27, only : sfluxrefo, raylo

      implicit none
      save

! Kurucz solar source function
! The following values were obtained using the "low resolution"
! version of the Kurucz solar source function.  For unknown reasons,
! the total irradiance in this band differs from the corresponding
! total in the "high-resolution" version of the Kurucz function.
! Therefore, these values are scaled later by the factor SCALEKUR.
      sfluxrefo(:) = (/ &
        & 14.0526_rb    , 11.4794_rb    , 8.72590_rb    , 5.56966_rb    , &
        & 3.80927_rb    , 1.57690_rb    , 1.15099_rb    , 1.10012_rb    , &
        & 0.658212_rb   , 5.86859e-02_rb, 5.56186e-02_rb, 4.68040e-02_rb, &
        & 3.64897e-02_rb, 3.58053e-02_rb, 1.38130e-02_rb, 1.90193e-03_rb /)

! Rayleigh extinction coefficient at v = 2925 cm-1.
      raylo(:) = (/ &
        & 3.44534e-06_rb,4.14480e-06_rb,4.95069e-06_rb,5.81204e-06_rb, &
        & 6.69748e-06_rb,7.56488e-06_rb,8.36344e-06_rb,9.04135e-06_rb, &
        & 9.58324e-06_rb,9.81542e-06_rb,9.75119e-06_rb,9.74533e-06_rb, &
        & 9.74139e-06_rb,9.73525e-06_rb,9.73577e-06_rb,9.73618e-06_rb /)

      end subroutine sw_kgb270

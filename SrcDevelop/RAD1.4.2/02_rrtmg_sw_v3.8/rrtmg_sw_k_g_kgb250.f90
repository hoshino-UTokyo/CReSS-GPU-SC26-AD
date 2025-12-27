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
!!prev      subroutine sw_kgb25
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for sfluxrefo(:), raylo(:), abso3ao(:), abso3bo(:)
!!
      subroutine sw_kgb250
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg25, only : kao, sfluxrefo, &
!!prev                            raylo, abso3ao, abso3bo
      use rrsw_kg25, only : sfluxrefo, &
                            raylo, abso3ao, abso3bo

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:) = (/ &
        & 42.6858_rb , 45.7720_rb, 44.9872_rb, 45.9662_rb    , &
        & 46.5458_rb , 41.6926_rb, 32.2893_rb, 24.0928_rb    , &
        & 16.7686_rb , 1.86048_rb, 1.54057_rb, 1.23503_rb    , &
        & 0.915085_rb,0.590099_rb,0.218622_rb, 3.21287e-02_rb /)

! Rayleigh extinction coefficient at v = 2925 cm-1.
      raylo(:) = (/ &
        & 9.81132e-07_rb,8.25605e-07_rb,6.71302e-07_rb,5.53556e-07_rb, & 
        & 3.97383e-07_rb,3.68206e-07_rb,4.42379e-07_rb,4.57799e-07_rb, &
        & 4.22683e-07_rb,3.87113e-07_rb,3.79810e-07_rb,3.63192e-07_rb, &
        & 3.51921e-07_rb,3.34231e-07_rb,3.34294e-07_rb,3.32673e-07_rb /)
     
      abso3ao(:) = (/ &
        & 2.32664e-02_rb,5.76154e-02_rb,0.125389_rb,0.250158_rb, &
        & 0.378756_rb   ,0.402196_rb   ,0.352026_rb,0.352036_rb, &
        & 0.386253_rb   ,0.414598_rb   ,0.420079_rb,0.435471_rb, &
        & 0.445487_rb   ,0.459549_rb   ,0.452920_rb,0.456838_rb /)

      abso3bo(:) = (/ &     
        & 1.76917e-02_rb,4.64185e-02_rb,1.03640e-01_rb,0.189469_rb, &
        & 0.303858_rb   ,0.400248_rb   ,0.447357_rb   ,0.470009_rb, &
        & 0.498673_rb   ,0.515696_rb   ,0.517053_rb   ,0.517930_rb, &
        & 0.518345_rb   ,0.524952_rb   ,0.508244_rb   ,0.468981_rb /)

      end subroutine sw_kgb250

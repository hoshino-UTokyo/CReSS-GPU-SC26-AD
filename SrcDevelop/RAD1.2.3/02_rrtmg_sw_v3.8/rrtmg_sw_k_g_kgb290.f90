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
!!prev      subroutine sw_kgb29
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for sfluxrefo(:), absh2oo(:), absco2o(:), rayl
!!
      subroutine sw_kgb290
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg29, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
!!prev                            absh2oo, absco2o, rayl
      use rrsw_kg29, only : sfluxrefo, &
                            absh2oo, absco2o, rayl

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:) = (/ &
        & 1.32880_rb    , 2.14018_rb    , 1.97612_rb    , 1.79000_rb    , &
        & 1.51242_rb    , 1.22977_rb    , 1.06052_rb    , 0.800996_rb   , &
        & 0.748053_rb   , 8.64369e-02_rb, 7.10675e-02_rb, 5.62425e-02_rb, &
        & 4.46988e-02_rb, 3.07441e-02_rb, 1.16728e-02_rb, 1.65573e-03_rb /)

      absco2o(:) = (/ &
        & 2.90073e-06_rb, 2.12382e-05_rb, 1.03032e-04_rb, 1.86481e-04_rb, &
        & 4.31997e-04_rb, 6.08238e-04_rb, 2.17603e-03_rb, 4.64479e-02_rb, &
        & 2.96956_rb    , 14.9569_rb    , 28.4831_rb    , 61.3998_rb    , &
        & 164.129_rb    , 832.282_rb    , 4995.02_rb    , 12678.1_rb     /)

      absh2oo(:) = (/ &
        & 2.99508e-04_rb, 3.95012e-03_rb, 1.49316e-02_rb, 3.24384e-02_rb, &
        & 6.92879e-02_rb, 0.123523_rb   , 0.360985_rb   , 1.86434_rb    , &
        & 10.38157_rb   , 0.214129_rb   , 0.213914_rb   , 0.212781_rb   , &
        & 0.215562_rb   , 0.218087_rb   , 0.220918_rb   , 0.218546_rb    /)
           
! Rayleigh extinction coefficient at v = 2200 cm-1.
      rayl = 9.30e-11_rb

      end subroutine sw_kgb290

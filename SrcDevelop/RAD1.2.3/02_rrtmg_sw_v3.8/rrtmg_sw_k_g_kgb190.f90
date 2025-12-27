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
!!prev      subroutine sw_kgb19
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for sfluxrefo(:,:), rayl
!!
      subroutine sw_kgb190
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg19, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg19, only : sfluxrefo, rayl

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:,1) = (/ &
        & 3.25791_rb    , 3.29697_rb    , 3.16031_rb    , 2.96115_rb    , &
        & 2.69238_rb    , 2.33819_rb    , 1.92760_rb    , 1.44918_rb    , &
        & 0.979764_rb   , 0.107336_rb   , 8.94523e-02_rb, 6.98325e-02_rb, &
        & 5.12051e-02_rb, 3.23645e-02_rb, 1.23401e-02_rb, 1.71339e-03_rb /)
      sfluxrefo(:,2) = (/ &
        & 3.22769_rb    , 3.28817_rb    , 3.16687_rb    , 2.97662_rb    , &
        & 2.69495_rb    , 2.34392_rb    , 1.92900_rb    , 1.45391_rb    , &
        & 0.982522_rb   , 0.107638_rb   , 8.92458e-02_rb, 6.99885e-02_rb, &
        & 5.09679e-02_rb, 3.23789e-02_rb, 1.22673e-02_rb, 1.56040e-03_rb /)
      sfluxrefo(:,3) = (/ &
        & 3.22294_rb    , 3.27780_rb    , 3.17424_rb    , 2.97143_rb    , &
        & 2.69785_rb    , 2.34993_rb    , 1.93155_rb    , 1.45196_rb    , &
        & 0.985329_rb   , 0.108027_rb   , 8.93552e-02_rb, 6.99937e-02_rb, &
        & 5.11678e-02_rb, 3.24846e-02_rb, 1.20636e-02_rb, 1.56040e-03_rb /)
      sfluxrefo(:,4) = (/ &
        & 3.22445_rb    , 3.26113_rb    , 3.18438_rb    , 2.96921_rb    , &
        & 2.69579_rb    , 2.35586_rb    , 1.93454_rb    , 1.44949_rb    , &
        & 0.987347_rb   , 0.108611_rb   , 8.91643e-02_rb, 7.02236e-02_rb, &
        & 5.12980e-02_rb, 3.25282e-02_rb, 1.21189e-02_rb, 1.56040e-03_rb /)
      sfluxrefo(:,5) = (/ &
        & 3.22497_rb    , 3.25109_rb    , 3.18741_rb    , 2.96970_rb    , &
        & 2.69460_rb    , 2.36020_rb    , 1.93301_rb    , 1.45224_rb    , &
        & 0.988564_rb   , 0.108255_rb   , 8.93830e-02_rb, 7.03655e-02_rb, &
        & 5.13017e-02_rb, 3.29414e-02_rb, 1.21189e-02_rb, 1.56040e-03_rb /)
      sfluxrefo(:,6) = (/ &
        & 3.22632_rb    , 3.24174_rb    , 3.18524_rb    , 2.97402_rb    , &
        & 2.69807_rb    , 2.35742_rb    , 1.93377_rb    , 1.45621_rb    , &
        & 0.988132_rb   , 0.108344_rb   , 8.93188e-02_rb, 7.04907e-02_rb, &
        & 5.17938e-02_rb, 3.31465e-02_rb, 1.21155e-02_rb, 1.56040e-03_rb /)
      sfluxrefo(:,7) = (/ &
        & 3.22793_rb    , 3.23589_rb    , 3.17720_rb    , 2.97869_rb    , &
        & 2.70293_rb    , 2.35436_rb    , 1.93557_rb    , 1.45868_rb    , &
        & 0.988654_rb   , 0.108198_rb   , 8.93375e-02_rb, 7.09790e-02_rb, &
        & 5.24733e-02_rb, 3.31298e-02_rb, 1.21126e-02_rb, 1.56040e-03_rb /)
      sfluxrefo(:,8) = (/ &
        & 3.22966_rb    , 3.24087_rb    , 3.15676_rb    , 2.98171_rb    , &
        & 2.70894_rb    , 2.34975_rb    , 1.93855_rb    , 1.46354_rb    , &
        & 0.988544_rb   , 0.108574_rb   , 9.02522e-02_rb, 7.12908e-02_rb, &
        & 5.24844e-02_rb, 3.31084e-02_rb, 1.21060e-02_rb, 1.56040e-03_rb /)
      sfluxrefo(:,9) = (/ &
        & 3.27240_rb    , 3.24666_rb    , 3.13886_rb    , 2.95238_rb    , &
        & 2.70190_rb    , 2.34460_rb    , 1.93948_rb    , 1.47111_rb    , &
        & 0.990821_rb   , 0.108730_rb   , 9.01625e-02_rb, 7.13261e-02_rb, &
        & 5.24813e-02_rb, 3.31083e-02_rb, 1.21126e-02_rb, 1.56040e-03_rb /)

! Rayleigh extinction coefficient at v = 4900 cm-1.
      rayl = 2.29e-09_rb

      end subroutine sw_kgb190

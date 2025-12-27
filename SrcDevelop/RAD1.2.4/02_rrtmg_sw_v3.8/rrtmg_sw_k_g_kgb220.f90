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
!!prev      subroutine sw_kgb22
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for  sfluxrefo(:,:), rayl
!!
      subroutine sw_kgb220
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg22, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg22, only : sfluxrefo, rayl

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:, 1) = (/ &
        & 3.71641_rb    ,3.63190_rb    ,3.44795_rb    ,3.17936_rb    , &
        & 2.86071_rb    ,2.48490_rb    ,2.02471_rb    ,1.52475_rb    , &
        & 1.03811_rb    ,0.113272_rb   ,9.37115e-02_rb,7.38969e-02_rb, &
        & 5.44713e-02_rb,3.45905e-02_rb,1.30293e-02_rb,1.84198e-03_rb /)
      sfluxrefo(:, 2) = (/ &
        & 3.73933_rb    ,3.60360_rb    ,3.43370_rb    ,3.19749_rb    , & 
        & 2.87747_rb    ,2.47926_rb    ,2.02175_rb    ,1.52010_rb    , &
        & 1.03612_rb    ,0.113265_rb   ,9.37145e-02_rb,7.38951e-02_rb, &
        & 5.44714e-02_rb,3.45906e-02_rb,1.30293e-02_rb,1.84198e-03_rb /)
      sfluxrefo(:, 3) = (/ &
        & 3.73889_rb    ,3.60279_rb    ,3.43404_rb    ,3.20560_rb    , &
        & 2.87367_rb    ,2.47515_rb    ,2.02412_rb    ,1.52315_rb    , &
        & 1.03146_rb    ,0.113272_rb   ,9.36707e-02_rb,7.39080e-02_rb, &
        & 5.44598e-02_rb,3.45906e-02_rb,1.30293e-02_rb,1.84198e-03_rb /)
      sfluxrefo(:, 4) = (/ &
        & 3.73801_rb    ,3.60530_rb    ,3.43659_rb    ,3.20640_rb    , &
        & 2.87039_rb    ,2.47330_rb    ,2.02428_rb    ,1.52509_rb    , &
        & 1.03037_rb    ,0.112553_rb   ,9.35352e-02_rb,7.39675e-02_rb, &
        & 5.43951e-02_rb,3.45669e-02_rb,1.30292e-02_rb,1.84198e-03_rb /)
      sfluxrefo(:, 5) = (/ &
        & 3.73809_rb    ,3.60996_rb    ,3.43602_rb    ,3.20364_rb    , &
        & 2.87005_rb    ,2.47343_rb    ,2.02353_rb    ,1.52617_rb    , &
        & 1.03138_rb    ,0.111172_rb   ,9.29885e-02_rb,7.35034e-02_rb, &
        & 5.42427e-02_rb,3.45732e-02_rb,1.30169e-02_rb,1.84550e-03_rb /)
      sfluxrefo(:, 6) = (/ &
        & 3.73872_rb    ,3.62054_rb    ,3.42934_rb    ,3.20110_rb    , &
        & 2.86886_rb    ,2.47379_rb    ,2.02237_rb    ,1.52754_rb    , & 
        & 1.03228_rb    ,0.111597_rb   ,9.12252e-02_rb,7.33115e-02_rb, &
        & 5.35600e-02_rb,3.45187e-02_rb,1.30184e-02_rb,1.84551e-03_rb /)
      sfluxrefo(:, 7) = (/ &
        & 3.73969_rb    ,3.65461_rb    ,3.40646_rb    ,3.19082_rb    , &
        & 2.86919_rb    ,2.47289_rb    ,2.02312_rb    ,1.52629_rb    , &
        & 1.03329_rb    ,0.111611_rb   ,9.16275e-02_rb,7.14731e-02_rb, &
        & 5.31771e-02_rb,3.44980e-02_rb,1.30190e-02_rb,1.84551e-03_rb /)
      sfluxrefo(:, 8) = (/ &
        & 3.73995_rb    ,3.65348_rb    ,3.43707_rb    ,3.16351_rb    , &
        & 2.87003_rb    ,2.47392_rb    ,2.02114_rb    ,1.52548_rb    , & 
        & 1.03306_rb    ,0.111088_rb   ,9.12422e-02_rb,7.11146e-02_rb, &
        & 5.31333e-02_rb,3.45302e-02_rb,1.30209e-02_rb,1.84554e-03_rb /)
      sfluxrefo(:, 9) = (/ &
        & 3.73788_rb    ,3.65004_rb    ,3.46938_rb    ,3.15236_rb    , &
        & 2.86381_rb    ,2.47393_rb    ,2.01715_rb    ,1.52134_rb    , &
        & 1.03163_rb    ,0.111259_rb   ,9.12948e-02_rb,7.09999e-02_rb, &
        & 5.31792e-02_rb,3.44955e-02_rb,1.30189e-02_rb,1.84551e-03_rb /)

! Rayleigh extinction coefficient at v = 8000 cm-1.
      rayl = 1.54e-08_rb

      end subroutine sw_kgb220

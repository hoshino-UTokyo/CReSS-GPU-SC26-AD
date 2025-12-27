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
!!prev      subroutine sw_kgb18
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for sfluxrefo(:,:), rayl
!!
      subroutine sw_kgb180
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg18, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg18, only : sfluxrefo, rayl

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:,1) = (/ &
        & 3.65840_rb    , 3.54375_rb    , 3.34481_rb    , 3.10534_rb    , &
        & 2.79879_rb    , 2.42841_rb    , 1.98748_rb    , 1.49377_rb    , &
        & 1.00196_rb    , 0.108342_rb   , 8.95099e-02_rb, 7.05199e-02_rb, &
        & 5.16432e-02_rb, 3.27635e-02_rb, 1.25133e-02_rb, 1.73001e-03_rb /)  
      sfluxrefo(:,2) = (/ &
        & 3.86372_rb    , 3.48521_rb    , 3.30790_rb    , 3.08103_rb    , &
        & 2.77552_rb    , 2.40722_rb    , 1.97307_rb    , 1.48023_rb    , &
        & 0.993055_rb   , 0.107691_rb   , 8.84430e-02_rb, 6.99354e-02_rb, &
        & 5.07881e-02_rb, 3.24121e-02_rb, 1.19442e-02_rb, 1.57612e-03_rb /)
      sfluxrefo(:,3) = (/ &
        & 3.90370_rb    , 3.50657_rb    , 3.30629_rb    , 3.06046_rb    , &
        & 2.76982_rb    , 2.39907_rb    , 1.96358_rb    , 1.47458_rb    , &
        & 0.988475_rb   , 0.106698_rb   , 8.75242e-02_rb, 6.85898e-02_rb, &
        & 5.04798e-02_rb, 3.13718e-02_rb, 1.09533e-02_rb, 1.57612e-03_rb /)
      sfluxrefo(:,4) = (/ &
        & 3.93165_rb    , 3.52058_rb    , 3.31346_rb    , 3.04944_rb    , &
        & 2.76074_rb    , 2.39433_rb    , 1.95556_rb    , 1.46712_rb    , &
        & 0.984056_rb   , 0.105885_rb   , 8.73062e-02_rb, 6.84054e-02_rb, &
        & 4.87443e-02_rb, 2.99295e-02_rb, 1.09533e-02_rb, 1.57612e-03_rb /)
      sfluxrefo(:,5) = (/ &
        & 3.94082_rb    , 3.55221_rb    , 3.31863_rb    , 3.04730_rb    , &
        & 2.74918_rb    , 2.38328_rb    , 1.95212_rb    , 1.45889_rb    , &
        & 0.978888_rb   , 0.105102_rb   , 8.65732e-02_rb, 6.74563e-02_rb, &
        & 4.76592e-02_rb, 2.91017e-02_rb, 1.09533e-02_rb, 1.57612e-03_rb /)
      sfluxrefo(:,6) = (/ &
        & 3.94198_rb    , 3.58743_rb    , 3.32106_rb    , 3.05866_rb    , &
        & 2.74115_rb    , 2.36939_rb    , 1.94305_rb    , 1.45180_rb    , &
        & 0.971784_rb   , 1.04045e-01_rb, 8.53731e-02_rb, 6.60654e-02_rb, &
        & 4.63228e-02_rb, 2.91016e-02_rb, 1.09552e-02_rb, 1.57612e-03_rb /)
      sfluxrefo(:,7) = (/ &
        & 3.93596_rb    , 3.63366_rb    , 3.33144_rb    , 3.06252_rb    , &
        & 2.74054_rb    , 2.35492_rb    , 1.92769_rb    , 1.44300_rb    , &
        & 0.961809_rb   , 1.02867e-01_rb, 8.34164e-02_rb, 6.41005e-02_rb, &
        & 4.61826e-02_rb, 2.91006e-02_rb, 1.09553e-02_rb, 1.57612e-03_rb /)
      sfluxrefo(:,8) = (/ &
        & 3.92520_rb    , 3.69078_rb    , 3.35656_rb    , 3.07055_rb    , &
        & 2.73862_rb    , 2.34430_rb    , 1.90187_rb    , 1.42242_rb    , &
        & 0.946676_rb   , 9.96302e-02_rb, 8.14421e-02_rb, 6.38622e-02_rb, &
        & 4.61794e-02_rb, 2.91017e-02_rb, 1.09553e-02_rb, 1.57612e-03_rb /)
      sfluxrefo(:,9) = (/ &
        & 3.80721_rb    , 3.74437_rb    , 3.50205_rb    , 3.18009_rb    , &
        & 2.75757_rb    , 2.29188_rb    , 1.84382_rb    , 1.35694_rb    , &
        & 0.914040_rb   , 9.86811e-02_rb, 8.14321e-02_rb, 6.38541e-02_rb, &
        & 4.61795e-02_rb, 2.90960e-02_rb, 1.09613e-02_rb, 1.57612e-03_rb /)

! Rayleigh extinction coefficient at v = 4325 cm-1.
      rayl = 1.39e-09_rb

      end subroutine sw_kgb180

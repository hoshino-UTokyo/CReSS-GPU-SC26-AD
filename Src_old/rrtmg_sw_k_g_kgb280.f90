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
!!prev      subroutine sw_kgb28
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for  sfluxrefo(:), rayl(:)
!!
      subroutine sw_kgb280
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg28, only : kao, kbo, sfluxrefo, rayl
      use rrsw_kg28, only : sfluxrefo, rayl

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:,1) = (/ &
        & 1.06156_rb    , 0.599910_rb   , 0.422462_rb   , 0.400077_rb   , &
        & 0.282221_rb   , 0.187893_rb   , 6.77357e-02_rb, 3.04572e-02_rb, &
        & 2.00442e-02_rb, 2.30786e-03_rb, 2.08824e-03_rb, 1.42604e-03_rb, &
        & 9.67384e-04_rb, 6.35362e-04_rb, 1.47727e-04_rb, 6.87639e-06_rb /)
      sfluxrefo(:,2) = (/ &
        & 1.07598_rb    , 0.585099_rb   , 0.422852_rb   , 0.400077_rb   , &
        & 0.282221_rb   , 0.187893_rb   , 6.69686e-02_rb, 3.09070e-02_rb, &
        & 2.02400e-02_rb, 2.47760e-03_rb, 1.89411e-03_rb, 1.41122e-03_rb, &
        & 1.12449e-03_rb, 5.73505e-04_rb, 2.04160e-04_rb, 1.58371e-05_rb /)
      sfluxrefo(:,3) = (/ &
        & 0.461647_rb   , 0.406113_rb   , 0.332506_rb   , 0.307508_rb   , &
        & 0.211167_rb   , 0.235457_rb   , 0.495886_rb   , 0.363921_rb   , &
        & 0.192700_rb   , 2.04678e-02_rb, 1.55407e-02_rb, 1.03882e-02_rb, &
        & 1.10778e-02_rb, 1.00504e-02_rb, 4.93497e-03_rb, 5.73410e-04_rb /)
      sfluxrefo(:,4) = (/ &
        & 0.132669_rb   , 0.175058_rb   , 0.359263_rb   , 0.388142_rb   , &
        & 0.350359_rb   , 0.475892_rb   , 0.489593_rb   , 0.408437_rb   , &
        & 0.221049_rb   , 1.94514e-02_rb, 1.54848e-02_rb, 1.44999e-02_rb, &
        & 1.44568e-02_rb, 1.00527e-02_rb, 4.95897e-03_rb, 5.73327e-04_rb /)
      sfluxrefo(:,5) = (/ &
        & 7.54800e-02_rb, 0.232246_rb   , 0.359263_rb   , 0.388142_rb   , &
        & 0.350359_rb   , 0.426317_rb   , 0.493485_rb   , 0.432016_rb   , &
        & 0.239203_rb   , 1.74951e-02_rb, 1.74477e-02_rb, 1.83566e-02_rb, &
        & 1.44818e-02_rb, 1.01048e-02_rb, 4.97487e-03_rb, 5.66831e-04_rb /)

! Rayleigh extinction coefficient at v = ????? cm-1.
      rayl = 2.02e-05_rb

    end subroutine sw_kgb280

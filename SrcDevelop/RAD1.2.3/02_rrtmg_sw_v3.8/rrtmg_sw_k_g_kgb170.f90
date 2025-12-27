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
!!prev      subroutine sw_kgb17
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for sfluxrefo(:,:), rayl
!!
      subroutine sw_kgb170
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg17, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg17, only : sfluxrefo, rayl

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:,1) = (/ &
        & 3.15613_rb  ,  3.03449_rb  ,  2.92069_rb  ,  2.63874_rb   , &
        & 2.34581_rb  ,  2.06999_rb  ,  1.70906_rb  ,  1.29085_rb   , &
        & 0.874851_rb ,  0.0955392_rb,  0.0787813_rb,  0.0621951_rb , &
        & 0.0459076_rb,  0.0294129_rb,  0.0110387_rb,  0.00159668_rb /)
      sfluxrefo(:,2) = (/ &
        & 2.83147_rb  ,  2.95919_rb  ,  2.96674_rb  ,  2.77677_rb   , &
        & 2.46826_rb  ,  2.11481_rb  ,  1.73243_rb  ,  1.30279_rb   , &
        & 0.882714_rb ,  0.0962350_rb,  0.0802122_rb,  0.0636194_rb , &
        & 0.0472620_rb,  0.0299051_rb,  0.0110785_rb,  0.00159668_rb /)
      sfluxrefo(:,3) = (/ &
        & 2.82300_rb  ,  2.94845_rb  ,  2.95887_rb  ,  2.77593_rb   , &
        & 2.47096_rb  ,  2.12596_rb  ,  1.73847_rb  ,  1.30796_rb   , &
        & 0.884395_rb ,  0.0966936_rb,  0.0801996_rb,  0.0640199_rb , &
        & 0.0472803_rb,  0.0300515_rb,  0.0112366_rb,  0.00160814_rb /)
      sfluxrefo(:,4) = (/ &
        & 2.81715_rb  ,  2.93789_rb  ,  2.95091_rb  ,  2.77046_rb   , &
        & 2.47716_rb  ,  2.13591_rb  ,  1.74365_rb  ,  1.31277_rb   , &
        & 0.887443_rb ,  0.0967016_rb,  0.0803391_rb,  0.0642442_rb , &
        & 0.0472909_rb,  0.0300720_rb,  0.0114817_rb,  0.00161875_rb /)
      sfluxrefo(:,5) = (/ &
        & 2.82335_rb  ,  2.93168_rb  ,  2.91455_rb  ,  2.75213_rb   , &
        & 2.49168_rb  ,  2.14408_rb  ,  1.75726_rb  ,  1.32401_rb   , &
        & 0.893644_rb ,  0.0969523_rb,  0.0805197_rb,  0.0639936_rb , &
        & 0.0475099_rb,  0.0305667_rb,  0.0115372_rb,  0.00161875_rb /)

! Rayleigh extinction coefficient at v = 3625 cm-1.
      rayl = 6.86e-10_rb

      end subroutine sw_kgb170

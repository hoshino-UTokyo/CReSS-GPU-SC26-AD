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
!!prev      subroutine sw_kgb21
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for  sfluxrefo(:), rayl
!!
      subroutine sw_kgb210
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg21, only : kao, kbo, selfrefo, forrefo, sfluxrefo, rayl
      use rrsw_kg21, only : forrefo, sfluxrefo, rayl

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:, 1) = (/ &
        & 16.1643_rb , 15.5806_rb, 14.7254_rb    , 13.5541_rb    , &
        & 11.9519_rb ,10.44410_rb, 8.37884_rb    , 6.26384_rb    , &
        & 4.28435_rb ,0.465228_rb, 0.385095_rb   ,0.304226_rb    , &
        & 0.222479_rb,0.143286_rb, 5.58046e-02_rb, 7.84856e-03_rb /)
      sfluxrefo(:, 2) = (/ &
        & 15.6451_rb , 15.3170_rb, 14.6987_rb    , 13.7350_rb    , &
        & 12.2267_rb ,10.51646_rb, 8.47150_rb    , 6.38873_rb    , &
        & 4.33536_rb ,0.470610_rb,0.389426_rb    ,0.306461_rb    , &
        & 0.223537_rb,0.143273_rb, 5.58179e-02_rb, 7.84856e-03_rb /)
      sfluxrefo(:, 3) = (/ &
        & 15.6092_rb , 15.3293_rb, 14.6881_rb    , 13.6693_rb    , &
        & 12.2342_rb ,10.52010_rb, 8.49442_rb    , 6.42138_rb    , &
        & 4.35865_rb ,0.473349_rb,0.391349_rb    ,0.308861_rb    , &
        & 0.224666_rb,0.144799_rb, 5.58176e-02_rb, 7.84881e-03_rb /)
      sfluxrefo(:, 4) = (/ &
        & 15.5786_rb , 15.3422_rb, 14.6894_rb    , 13.6040_rb    , &
        & 12.2567_rb ,10.49400_rb, 8.53521_rb    , 6.44427_rb    , &
        & 4.37208_rb ,0.475709_rb,0.392956_rb    ,0.309737_rb    , &
        & 0.226274_rb,0.146483_rb, 5.59325e-02_rb, 7.84881e-03_rb /)
      sfluxrefo(:, 5) = (/ &
        & 15.5380_rb , 15.3826_rb, 14.6575_rb    , 13.5722_rb    , &
        & 12.2646_rb ,10.47672_rb, 8.57158_rb    , 6.46343_rb    , &
        & 4.38259_rb ,0.477647_rb,0.393982_rb    ,0.310686_rb    , &
        & 0.227620_rb,0.148376_rb, 5.60398e-02_rb, 7.83925e-03_rb /)
      sfluxrefo(:, 6) = (/ &
        & 15.5124_rb , 15.3986_rb, 14.6240_rb    , 13.5535_rb    , &
        & 12.2468_rb ,10.48891_rb, 8.60434_rb    , 6.47985_rb    , &
        & 4.39448_rb ,0.478267_rb,0.395618_rb    ,0.311043_rb    , &
        & 0.230927_rb,0.148774_rb, 5.61189e-02_rb, 7.83925e-03_rb /)
      sfluxrefo(:, 7) = (/ &
        & 15.4910_rb , 15.4028_rb, 14.5772_rb    , 13.5507_rb    , &
        & 12.2122_rb ,10.52735_rb, 8.62650_rb    , 6.49644_rb    , &
        & 4.41173_rb ,0.478627_rb,0.396433_rb    ,0.314199_rb    , & 
        & 0.233125_rb,0.149052_rb, 5.62309e-02_rb, 7.83925e-03_rb /)
      sfluxrefo(:, 8) = (/ &
        & 15.4562_rb , 15.3928_rb, 14.5510_rb    , 13.5122_rb    , &
        & 12.1890_rb , 10.5826_rb, 8.65842_rb    , 6.51558_rb    , &
        & 4.42747_rb ,0.480669_rb,0.400143_rb    ,0.318144_rb    , &
        & 0.233937_rb,0.149119_rb, 5.62309e-02_rb, 7.83925e-03_rb /)
      sfluxrefo(:, 9) = (/ &
        & 15.0069_rb , 15.1479_rb, 14.7802_rb    , 13.6085_rb    , &
        & 12.2793_rb , 10.6929_rb, 8.72723_rb    , 6.57114_rb    , &
        & 4.46330_rb ,0.486724_rb,0.401446_rb    ,0.318879_rb    , &
        & 0.233959_rb,0.149119_rb, 5.62309e-02_rb, 7.83925e-03_rb /)

! Rayleigh extinction coefficient at v = 6925 cm-1.
      rayl = 9.41e-09_rb

      end subroutine sw_kgb210

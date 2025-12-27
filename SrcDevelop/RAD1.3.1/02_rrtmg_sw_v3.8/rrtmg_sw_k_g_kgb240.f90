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
!!prev      subroutine sw_kgb24
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for sfluxrefo(:,:), raylao(:,:), raylbo(:), abso3ao(:), abso3bo(:)
!!
      subroutine sw_kgb240
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrsw_kg24, only : kao, kbo, selfrefo, forrefo, sfluxrefo, &
!!prev                            raylao, raylbo, abso3ao, abso3bo
      use rrsw_kg24, only : sfluxrefo, &
                            raylao, raylbo, abso3ao, abso3bo

      implicit none
      save

! Kurucz solar source function
      sfluxrefo(:,1) = (/ &
        & 34.3610_rb , 33.1240_rb, 31.3948_rb, 28.7248_rb, &
        & 24.7884_rb , 21.4892_rb, 17.3972_rb, 13.7928_rb, &
        & 9.54462_rb , 1.05002_rb,0.867332_rb,0.685753_rb, &
        & 0.504718_rb,0.323112_rb,0.122183_rb, 1.70288e-02_rb /)
      sfluxrefo(:,2) = (/ &
        & 34.2367_rb , 32.4327_rb, 30.0863_rb, 28.2085_rb, & 
        & 25.6533_rb , 22.3412_rb, 18.3112_rb, 13.8521_rb, &
        & 9.51035_rb , 1.04138_rb,0.863493_rb,0.682790_rb, &
        & 0.504721_rb,0.323102_rb,0.122193_rb, 1.70288e-02_rb /)
      sfluxrefo(:,3) = (/ &
        & 34.1883_rb , 32.2479_rb, 30.2650_rb, 28.2914_rb, &
        & 25.6626_rb , 22.3163_rb, 18.3327_rb, 13.8508_rb, &
        & 9.49190_rb , 1.03672_rb,0.858272_rb,0.681485_rb, &
        & 0.501363_rb,0.323110_rb,0.122183_rb, 1.70288e-02_rb /)
      sfluxrefo(:,4) = (/ &
        & 34.1365_rb , 32.2316_rb, 30.3325_rb, 28.3305_rb, &
        & 25.6420_rb , 22.3223_rb, 18.3411_rb, 13.8471_rb, &
        & 9.47492_rb , 1.03376_rb,0.855380_rb,0.679085_rb, &
        & 0.497998_rb,0.323053_rb,0.122183_rb, 1.70288e-02_rb /)
      sfluxrefo(:,5) = (/ &
        & 34.0460_rb , 32.2795_rb, 30.4147_rb, 28.3123_rb, &
        & 25.6438_rb , 22.3238_rb, 18.3441_rb, 13.8528_rb, &
        & 9.45222_rb , 1.03058_rb,0.854037_rb,0.675554_rb, &
        & 0.498344_rb,0.320072_rb,0.122193_rb, 1.70288e-02_rb /)
      sfluxrefo(:,6) = (/ &
        & 33.9909_rb , 32.3127_rb, 30.4854_rb, 28.3005_rb, &
        & 25.6310_rb , 22.3294_rb, 18.3459_rb, 13.8488_rb, &
        & 9.43336_rb , 1.02901_rb,0.852728_rb,0.672322_rb, &
        & 0.498056_rb,0.317753_rb,0.122183_rb, 1.70288e-02_rb /)
      sfluxrefo(:,7) = (/ &
        & 33.9225_rb , 32.4097_rb, 30.5125_rb, 28.2810_rb, &
        & 25.6387_rb , 22.3080_rb, 18.3715_rb, 13.8248_rb, &
        & 9.41834_rb , 1.02735_rb,0.850807_rb,0.671379_rb, &
        & 0.496975_rb,0.317158_rb,0.119297_rb, 1.70207e-02_rb /)
      sfluxrefo(:,8) = (/ &
        & 33.8940_rb , 32.4951_rb, 30.5494_rb, 28.2788_rb, &
        & 25.5975_rb , 22.3225_rb, 18.3358_rb, 13.8199_rb, &
        & 9.40283_rb , 1.02751_rb,0.850729_rb,0.670152_rb, &
        & 0.494294_rb,0.315829_rb,0.116195_rb, 1.64138e-02_rb /)
      sfluxrefo(:,9) = (/ &
        & 34.6501_rb , 32.6690_rb, 30.2872_rb, 28.0955_rb, &
        & 25.4662_rb , 22.1446_rb, 18.2754_rb, 13.7573_rb, &
        & 9.36645_rb , 1.02356_rb,0.847154_rb,0.668519_rb, &
        & 0.489186_rb,0.313790_rb,0.117074_rb, 1.60943e-02_rb /)

! Rayleigh extinction coefficient at all v
      raylao(:,1) = (/ &
        & 1.28405e-07_rb,1.45501e-07_rb,1.67272e-07_rb,1.94856e-07_rb, &
        & 2.15248e-07_rb,2.34920e-07_rb,2.48558e-07_rb,1.80004e-07_rb, &
        & 1.46504e-07_rb,1.31355e-07_rb,1.33562e-07_rb,1.35618e-07_rb, &
        & 1.22412e-07_rb,1.19842e-07_rb,1.19924e-07_rb,1.20264e-07_rb /)
      raylao(:,2) = (/ &
        & 1.41622e-07_rb,1.93436e-07_rb,2.25057e-07_rb,2.01025e-07_rb, &
        & 1.85138e-07_rb,1.72672e-07_rb,1.64771e-07_rb,1.59312e-07_rb, &
        & 1.44961e-07_rb,1.37448e-07_rb,1.37506e-07_rb,1.38081e-07_rb, &
        & 1.22432e-07_rb,1.19844e-07_rb,1.19921e-07_rb,1.20287e-07_rb /)
      raylao(:,3) = (/ &
        & 1.45382e-07_rb,1.97020e-07_rb,2.22781e-07_rb,1.96062e-07_rb, &
        & 1.83495e-07_rb,1.72495e-07_rb,1.64910e-07_rb,1.58797e-07_rb, &
        & 1.46208e-07_rb,1.42274e-07_rb,1.40445e-07_rb,1.39496e-07_rb, &
        & 1.26940e-07_rb,1.19844e-07_rb,1.19921e-07_rb,1.20287e-07_rb /)
      raylao(:,4) = (/ &
        & 1.48247e-07_rb,1.99958e-07_rb,2.18048e-07_rb,1.93896e-07_rb, &
        & 1.83125e-07_rb,1.73244e-07_rb,1.64320e-07_rb,1.58298e-07_rb, &
        & 1.48428e-07_rb,1.44769e-07_rb,1.43704e-07_rb,1.38498e-07_rb, &
        & 1.31732e-07_rb,1.22299e-07_rb,1.19921e-07_rb,1.20287e-07_rb /)
      raylao(:,5) = (/ &
        & 1.51343e-07_rb,1.99621e-07_rb,2.14563e-07_rb,1.93824e-07_rb, &
        & 1.82992e-07_rb,1.73143e-07_rb,1.64587e-07_rb,1.57355e-07_rb, &
        & 1.51198e-07_rb,1.46373e-07_rb,1.45438e-07_rb,1.38095e-07_rb, &
        & 1.35026e-07_rb,1.27504e-07_rb,1.19921e-07_rb,1.20287e-07_rb /)
      raylao(:,6) = (/ &
        & 1.54462e-07_rb,1.97610e-07_rb,2.11992e-07_rb,1.93831e-07_rb, &
        & 1.83900e-07_rb,1.73125e-07_rb,1.64093e-07_rb,1.57651e-07_rb, &
        & 1.53158e-07_rb,1.46843e-07_rb,1.44733e-07_rb,1.40611e-07_rb, &
        & 1.37320e-07_rb,1.33932e-07_rb,1.20423e-07_rb,1.20287e-07_rb /)
      raylao(:,7) = (/ &
        & 1.59068e-07_rb,1.92757e-07_rb,2.09865e-07_rb,1.95132e-07_rb, &
        & 1.83641e-07_rb,1.73778e-07_rb,1.63215e-07_rb,1.59462e-07_rb, &
        & 1.54331e-07_rb,1.46177e-07_rb,1.45819e-07_rb,1.43177e-07_rb, &
        & 1.39797e-07_rb,1.36780e-07_rb,1.33385e-07_rb,1.20287e-07_rb /)
      raylao(:,8) = (/ &
        & 1.62066e-07_rb,1.87529e-07_rb,2.07191e-07_rb,1.97788e-07_rb, &
        & 1.84920e-07_rb,1.72951e-07_rb,1.65450e-07_rb,1.60344e-07_rb, &
        & 1.54403e-07_rb,1.47679e-07_rb,1.47287e-07_rb,1.44951e-07_rb, &
        & 1.42517e-07_rb,1.41107e-07_rb,1.48688e-07_rb,1.51127e-07_rb /)
      raylao(:,9) = (/ &
        & 1.19177e-07_rb,1.86522e-07_rb,2.20324e-07_rb,2.13543e-07_rb, &
        & 1.92198e-07_rb,1.81641e-07_rb,1.70092e-07_rb,1.65072e-07_rb, &
        & 1.59804e-07_rb,1.56745e-07_rb,1.51235e-07_rb,1.51400e-07_rb, &
        & 1.49635e-07_rb,1.48056e-07_rb,1.49046e-07_rb,1.51010e-07_rb /)

      raylbo(:) = (/ &
        & 1.23766e-07_rb,1.40524e-07_rb,1.61610e-07_rb,1.83232e-07_rb, &
        & 2.02951e-07_rb,2.21367e-07_rb,2.38367e-07_rb,2.53019e-07_rb, &
        & 2.12202e-07_rb,1.36977e-07_rb,1.39118e-07_rb,1.37097e-07_rb, &
        & 1.33223e-07_rb,1.38695e-07_rb,1.19868e-07_rb,1.20062e-07_rb /)

      abso3ao(:) = (/ &
        & 8.03067e-02_rb,0.180926_rb   ,0.227484_rb   ,0.168015_rb   , &
        & 0.138284_rb   ,0.114537_rb   ,9.50114e-02_rb,8.06816e-02_rb, &
        & 6.76406e-02_rb,5.69802e-02_rb,5.63283e-02_rb,4.57592e-02_rb, &
        & 4.21862e-02_rb,3.47949e-02_rb,2.65731e-02_rb,2.67628e-02_rb /)

      abso3bo(:) = (/ &
        & 2.94848e-02_rb,4.33642e-02_rb,6.70197e-02_rb,0.104990_rb   , &
        & 0.156180_rb   ,0.214638_rb   ,0.266281_rb   ,0.317941_rb   , &
        & 0.355327_rb   ,0.371241_rb   ,0.374396_rb   ,0.326847_rb   , &
        & 0.126497_rb   ,6.95264e-02_rb,2.58175e-02_rb,2.52862e-02_rb /)

      end subroutine sw_kgb240

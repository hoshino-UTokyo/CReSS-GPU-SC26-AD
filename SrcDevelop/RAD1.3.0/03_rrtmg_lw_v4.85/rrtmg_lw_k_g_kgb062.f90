!     path:      $Source:
!     /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_k_g.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 21:04:30 $
!
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
!      subroutine lw_kgbnn
! **************************************************************************
!  RRTM Longwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original version:   E. J. Mlawer, et al.
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90 formatting:  Michael J. Iacono; June 2006
!
!  This file contains 16 subroutines that include the 
!  absorption coefficients and other data for each of the 16 longwave
!  spectral bands used in RRTM.  Here, the data are defined for 16
!  g-points, or sub-intervals, per band.  These data are combined and
!  weighted using a mapping procedure in routine RRTMG_LW_INIT to reduce
!  the total number of g-points from 256 to 140 for use in the GCM.
! **************************************************************************
!!prev      subroutine lw_kgb06
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for kao_mco2(:,:)
!!
      subroutine lw_kgb062
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg06, only : fracrefao, kao, kao_mco2, selfrefo, forrefo, &
!!prev                            cfc11adjo, cfc12o
      use rrlw_kg06, only : kao_mco2

      implicit none
      save

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

      kao_mco2(:, 1) = (/ &
     & 1.45661e-05_rb, 1.73337e-05_rb, 2.06273e-05_rb, 2.45466e-05_rb, 2.92105e-05_rb, &
     & 3.47607e-05_rb, 4.13654e-05_rb, 4.92251e-05_rb, 5.85781e-05_rb, 6.97083e-05_rb, &
     & 8.29533e-05_rb, 9.87149e-05_rb, 1.17471e-04_rb, 1.39792e-04_rb, 1.66353e-04_rb, &
     & 1.97961e-04_rb, 2.35574e-04_rb, 2.80335e-04_rb, 3.33600e-04_rb/)
      kao_mco2(:, 2) = (/ &
     & 9.96332e-06_rb, 1.21229e-05_rb, 1.47506e-05_rb, 1.79478e-05_rb, 2.18381e-05_rb, &
     & 2.65716e-05_rb, 3.23310e-05_rb, 3.93389e-05_rb, 4.78658e-05_rb, 5.82408e-05_rb, &
     & 7.08647e-05_rb, 8.62250e-05_rb, 1.04914e-04_rb, 1.27655e-04_rb, 1.55325e-04_rb, &
     & 1.88992e-04_rb, 2.29957e-04_rb, 2.79801e-04_rb, 3.40448e-04_rb/)
      kao_mco2(:, 3) = (/ &
     & 1.14968e-05_rb, 1.39890e-05_rb, 1.70215e-05_rb, 2.07115e-05_rb, 2.52013e-05_rb, &
     & 3.06644e-05_rb, 3.73118e-05_rb, 4.54002e-05_rb, 5.52420e-05_rb, 6.72173e-05_rb, &
     & 8.17887e-05_rb, 9.95188e-05_rb, 1.21092e-04_rb, 1.47343e-04_rb, 1.79283e-04_rb, &
     & 2.18148e-04_rb, 2.65438e-04_rb, 3.22980e-04_rb, 3.92995e-04_rb/)
      kao_mco2(:, 4) = (/ &
     & 1.02186e-05_rb, 1.23232e-05_rb, 1.48613e-05_rb, 1.79222e-05_rb, 2.16134e-05_rb, &
     & 2.60649e-05_rb, 3.14332e-05_rb, 3.79071e-05_rb, 4.57145e-05_rb, 5.51297e-05_rb, &
     & 6.64843e-05_rb, 8.01773e-05_rb, 9.66905e-05_rb, 1.16605e-04_rb, 1.40621e-04_rb, &
     & 1.69583e-04_rb, 2.04510e-04_rb, 2.46631e-04_rb, 2.97426e-04_rb/)
      kao_mco2(:, 5) = (/ &
     & 1.03469e-05_rb, 1.24680e-05_rb, 1.50239e-05_rb, 1.81037e-05_rb, 2.18149e-05_rb, &
     & 2.62869e-05_rb, 3.16756e-05_rb, 3.81690e-05_rb, 4.59935e-05_rb, 5.54220e-05_rb, &
     & 6.67833e-05_rb, 8.04737e-05_rb, 9.69704e-05_rb, 1.16849e-04_rb, 1.40803e-04_rb, &
     & 1.69667e-04_rb, 2.04448e-04_rb, 2.46359e-04_rb, 2.96861e-04_rb/)
      kao_mco2(:, 6) = (/ &
     & 1.71660e-05_rb, 2.07334e-05_rb, 2.50420e-05_rb, 3.02461e-05_rb, 3.65317e-05_rb, &
     & 4.41235e-05_rb, 5.32930e-05_rb, 6.43680e-05_rb, 7.77446e-05_rb, 9.39010e-05_rb, &
     & 1.13415e-04_rb, 1.36984e-04_rb, 1.65451e-04_rb, 1.99835e-04_rb, 2.41363e-04_rb, &
     & 2.91522e-04_rb, 3.52104e-04_rb, 4.25276e-04_rb, 5.13654e-04_rb/)
      kao_mco2(:, 7) = (/ &
     & 4.78803e-05_rb, 5.79395e-05_rb, 7.01119e-05_rb, 8.48418e-05_rb, 1.02666e-04_rb, &
     & 1.24235e-04_rb, 1.50336e-04_rb, 1.81920e-04_rb, 2.20139e-04_rb, 2.66388e-04_rb, &
     & 3.22354e-04_rb, 3.90077e-04_rb, 4.72028e-04_rb, 5.71197e-04_rb, 6.91199e-04_rb, &
     & 8.36413e-04_rb, 1.01214e-03_rb, 1.22477e-03_rb, 1.48209e-03_rb/)
      kao_mco2(:, 8) = (/ &
     & 1.27954e-04_rb, 1.55281e-04_rb, 1.88445e-04_rb, 2.28692e-04_rb, 2.77534e-04_rb, &
     & 3.36808e-04_rb, 4.08741e-04_rb, 4.96037e-04_rb, 6.01977e-04_rb, 7.30542e-04_rb, &
     & 8.86566e-04_rb, 1.07591e-03_rb, 1.30570e-03_rb, 1.58456e-03_rb, 1.92298e-03_rb, &
     & 2.33367e-03_rb, 2.83208e-03_rb, 3.43694e-03_rb, 4.17097e-03_rb/)
      kao_mco2(:, 9) = (/ &
     & 2.93792e-05_rb, 3.55109e-05_rb, 4.29223e-05_rb, 5.18805e-05_rb, 6.27083e-05_rb, &
     & 7.57960e-05_rb, 9.16151e-05_rb, 1.10736e-04_rb, 1.33847e-04_rb, 1.61782e-04_rb, &
     & 1.95547e-04_rb, 2.36359e-04_rb, 2.85689e-04_rb, 3.45315e-04_rb, 4.17384e-04_rb, &
     & 5.04495e-04_rb, 6.09787e-04_rb, 7.37054e-04_rb, 8.90882e-04_rb/)
      kao_mco2(:,10) = (/ &
     & 5.08569e-05_rb, 6.24700e-05_rb, 7.67350e-05_rb, 9.42574e-05_rb, 1.15781e-04_rb, &
     & 1.42220e-04_rb, 1.74695e-04_rb, 2.14587e-04_rb, 2.63588e-04_rb, 3.23778e-04_rb, &
     & 3.97712e-04_rb, 4.88530e-04_rb, 6.00085e-04_rb, 7.37114e-04_rb, 9.05433e-04_rb, &
     & 1.11219e-03_rb, 1.36616e-03_rb, 1.67812e-03_rb, 2.06131e-03_rb/)
      kao_mco2(:,11) = (/ &
     & 4.82546e-06_rb, 6.21462e-06_rb, 8.00369e-06_rb, 1.03078e-05_rb, 1.32752e-05_rb, &
     & 1.70969e-05_rb, 2.20188e-05_rb, 2.83575e-05_rb, 3.65211e-05_rb, 4.70348e-05_rb, &
     & 6.05753e-05_rb, 7.80138e-05_rb, 1.00472e-04_rb, 1.29397e-04_rb, 1.66647e-04_rb, &
     & 2.14622e-04_rb, 2.76407e-04_rb, 3.55980e-04_rb, 4.58459e-04_rb/)
      kao_mco2(:,12) = (/ &
     & 2.41346e-06_rb, 2.96282e-06_rb, 3.63723e-06_rb, 4.46516e-06_rb, 5.48153e-06_rb, &
     & 6.72926e-06_rb, 8.26100e-06_rb, 1.01414e-05_rb, 1.24498e-05_rb, 1.52837e-05_rb, &
     & 1.87627e-05_rb, 2.30335e-05_rb, 2.82765e-05_rb, 3.47129e-05_rb, 4.26144e-05_rb, &
     & 5.23144e-05_rb, 6.42225e-05_rb, 7.88410e-05_rb, 9.67871e-05_rb/)
      kao_mco2(:,13) = (/ &
     & 2.76412e-06_rb, 3.46195e-06_rb, 4.33596e-06_rb, 5.43062e-06_rb, 6.80164e-06_rb, &
     & 8.51879e-06_rb, 1.06695e-05_rb, 1.33631e-05_rb, 1.67367e-05_rb, 2.09621e-05_rb, &
     & 2.62542e-05_rb, 3.28824e-05_rb, 4.11839e-05_rb, 5.15813e-05_rb, 6.46035e-05_rb, &
     & 8.09134e-05_rb, 1.01341e-04_rb, 1.26925e-04_rb, 1.58969e-04_rb/)
      kao_mco2(:,14) = (/ &
     & 1.25126e-06_rb, 1.54971e-06_rb, 1.91935e-06_rb, 2.37715e-06_rb, 2.94416e-06_rb, &
     & 3.64640e-06_rb, 4.51615e-06_rb, 5.59335e-06_rb, 6.92749e-06_rb, 8.57985e-06_rb, &
     & 1.06263e-05_rb, 1.31610e-05_rb, 1.63001e-05_rb, 2.01881e-05_rb, 2.50034e-05_rb, &
     & 3.09672e-05_rb, 3.83536e-05_rb, 4.75018e-05_rb, 5.88319e-05_rb/)
      kao_mco2(:,15) = (/ &
     & 1.59748e-06_rb, 2.08378e-06_rb, 2.71812e-06_rb, 3.54557e-06_rb, 4.62491e-06_rb, &
     & 6.03282e-06_rb, 7.86932e-06_rb, 1.02649e-05_rb, 1.33897e-05_rb, 1.74658e-05_rb, &
     & 2.27827e-05_rb, 2.97182e-05_rb, 3.87649e-05_rb, 5.05657e-05_rb, 6.59589e-05_rb, &
     & 8.60380e-05_rb, 1.12230e-04_rb, 1.46394e-04_rb, 1.90959e-04_rb/)
      kao_mco2(:,16) = (/ &
     & 1.68148e-06_rb, 2.17133e-06_rb, 2.80388e-06_rb, 3.62071e-06_rb, 4.67549e-06_rb, &
     & 6.03756e-06_rb, 7.79642e-06_rb, 1.00677e-05_rb, 1.30006e-05_rb, 1.67879e-05_rb, &
     & 2.16786e-05_rb, 2.79941e-05_rb, 3.61493e-05_rb, 4.66803e-05_rb, 6.02792e-05_rb, &
     & 7.78398e-05_rb, 1.00516e-04_rb, 1.29799e-04_rb, 1.67612e-04_rb/)

      end subroutine lw_kgb062

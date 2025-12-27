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
!!prev      subroutine lw_kgb07
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:,:), fracrefbo(:,:)
!!
      subroutine lw_kgb070
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg07, only : fracrefao, fracrefbo, kao, kbo, kao_mco2, &
!!prev                            kbo_mco2, selfrefo, forrefo
      use rrlw_kg07, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level : P = 706.27 mb, T = 278.94 K
      fracrefao(:, 1) = (/ &
        1.6312e-01_rb,1.4949e-01_rb,1.4305e-01_rb,1.3161e-01_rb,1.1684e-01_rb,9.9900e-02_rb, &
        8.0912e-02_rb,6.0203e-02_rb,4.0149e-02_rb,4.3365e-03_rb,3.5844e-03_rb,2.8019e-03_rb, &
        2.0756e-03_rb,1.3449e-03_rb,5.0492e-04_rb,7.1194e-05_rb/)
      fracrefao(:, 2) = (/ &
        1.6329e-01_rb,1.4989e-01_rb,1.4328e-01_rb,1.3101e-01_rb,1.1691e-01_rb,9.9754e-02_rb, &
        8.0956e-02_rb,5.9912e-02_rb,4.0271e-02_rb,4.3298e-03_rb,3.5626e-03_rb,2.8421e-03_rb, &
        2.1031e-03_rb,1.3360e-03_rb,4.8965e-04_rb,6.8900e-05_rb/)
      fracrefao(:, 3) = (/ &
        1.6236e-01_rb,1.5081e-01_rb,1.4341e-01_rb,1.3083e-01_rb,1.1684e-01_rb,9.9701e-02_rb, &
        8.0956e-02_rb,5.9884e-02_rb,4.0245e-02_rb,4.3837e-03_rb,3.6683e-03_rb,2.9250e-03_rb, &
        2.0969e-03_rb,1.3320e-03_rb,4.8965e-04_rb,6.8900e-05_rb/)
      fracrefao(:, 4) = (/ &
        1.6096e-01_rb,1.5183e-01_rb,1.4354e-01_rb,1.3081e-01_rb,1.1687e-01_rb,9.9619e-02_rb, &
        8.0947e-02_rb,5.9899e-02_rb,4.0416e-02_rb,4.4389e-03_rb,3.7280e-03_rb,2.9548e-03_rb, &
        2.0977e-03_rb,1.3305e-03_rb,4.8965e-04_rb,6.8900e-05_rb/)
      fracrefao(:, 5) = (/ &
        1.5661e-01_rb,1.5478e-01_rb,1.4414e-01_rb,1.3097e-01_rb,1.1695e-01_rb,9.9823e-02_rb, &
        8.0750e-02_rb,6.0100e-02_rb,4.0741e-02_rb,4.4598e-03_rb,3.7366e-03_rb,2.9521e-03_rb, &
        2.0980e-03_rb,1.3297e-03_rb,4.8965e-04_rb,6.8900e-05_rb/)
      fracrefao(:, 6) = (/ &
        1.4879e-01_rb,1.5853e-01_rb,1.4586e-01_rb,1.3162e-01_rb,1.1729e-01_rb,1.0031e-01_rb, &
        8.0908e-02_rb,6.0460e-02_rb,4.1100e-02_rb,4.4578e-03_rb,3.7388e-03_rb,2.9508e-03_rb, &
        2.0986e-03_rb,1.3288e-03_rb,4.8965e-04_rb,6.8900e-05_rb/)
      fracrefao(:, 7) = (/ &
        1.4117e-01_rb,1.4838e-01_rb,1.4807e-01_rb,1.3759e-01_rb,1.2218e-01_rb,1.0228e-01_rb, &
        8.2130e-02_rb,6.1546e-02_rb,4.1522e-02_rb,4.4577e-03_rb,3.7428e-03_rb,2.9475e-03_rb, &
        2.0997e-03_rb,1.3277e-03_rb,4.8965e-04_rb,6.8900e-05_rb/)
      fracrefao(:, 8) = (/ &
        1.4018e-01_rb,1.4207e-01_rb,1.3919e-01_rb,1.3332e-01_rb,1.2325e-01_rb,1.0915e-01_rb, &
        9.0280e-02_rb,6.5554e-02_rb,4.1852e-02_rb,4.4707e-03_rb,3.7572e-03_rb,2.9364e-03_rb, &
        2.1023e-03_rb,1.3249e-03_rb,4.8965e-04_rb,6.8900e-05_rb/)
      fracrefao(:, 9) = (/ &
        1.4863e-01_rb,1.4926e-01_rb,1.4740e-01_rb,1.3558e-01_rb,1.1999e-01_rb,1.0044e-01_rb, &
        8.1927e-02_rb,6.0989e-02_rb,4.0665e-02_rb,4.4481e-03_rb,3.7369e-03_rb,2.9482e-03_rb, &
        2.0976e-03_rb,1.3281e-03_rb,4.8965e-04_rb,6.8900e-05_rb/)

! Planck fraction mapping level : P=95.58 mbar, T= 215.70 K
      fracrefbo(:) = (/ &
        1.5872e-01_rb,1.5443e-01_rb,1.4413e-01_rb,1.3147e-01_rb,1.1634e-01_rb,9.8914e-02_rb, &
        8.0236e-02_rb,6.0197e-02_rb,4.0624e-02_rb,4.4225e-03_rb,3.6688e-03_rb,2.9074e-03_rb, &
        2.0862e-03_rb,1.3039e-03_rb,4.8561e-04_rb,6.8854e-05_rb/)

      end subroutine lw_kgb070

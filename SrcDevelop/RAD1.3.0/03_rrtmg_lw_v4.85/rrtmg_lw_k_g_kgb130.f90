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
!!prev      subroutine lw_kgb13
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for fracrefao(:,:), fracrefbo(:)
!!
      subroutine lw_kgb130
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg13, only : fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
!!prev                            kbo_mo3, selfrefo, forrefo
      use rrlw_kg13, only : fracrefao, fracrefbo

      implicit none
      save

! Planck fraction mapping level : P=473.4280 mb, T = 259.83 K      
      fracrefao(:, 1) = (/ &
     &  1.7534e-01_rb,1.7394e-01_rb,1.6089e-01_rb,1.3782e-01_rb,1.0696e-01_rb,8.5853e-02_rb, &
     &  6.6548e-02_rb,4.9053e-02_rb,3.2064e-02_rb,3.4820e-03_rb,2.8763e-03_rb,2.2204e-03_rb, &
     &  1.5612e-03_rb,9.8572e-04_rb,3.6853e-04_rb,5.1612e-05_rb/)
      fracrefao(:, 2) = (/ &
     &  1.7489e-01_rb,1.7309e-01_rb,1.5981e-01_rb,1.3782e-01_rb,1.0797e-01_rb,8.6367e-02_rb, &
     &  6.7042e-02_rb,4.9257e-02_rb,3.2207e-02_rb,3.4820e-03_rb,2.8767e-03_rb,2.2203e-03_rb, &
     &  1.5613e-03_rb,9.8571e-04_rb,3.6853e-04_rb,5.1612e-05_rb/)
      fracrefao(:, 3) = (/ &
     &  1.7459e-01_rb,1.7259e-01_rb,1.5948e-01_rb,1.3694e-01_rb,1.0815e-01_rb,8.7376e-02_rb, &
     &  6.7339e-02_rb,4.9541e-02_rb,3.2333e-02_rb,3.5019e-03_rb,2.8958e-03_rb,2.2527e-03_rb, &
     &  1.6099e-03_rb,9.8574e-04_rb,3.6853e-04_rb,5.1612e-05_rb/)
      fracrefao(:, 4) = (/ &
     &  1.7391e-01_rb,1.7244e-01_rb,1.5921e-01_rb,1.3644e-01_rb,1.0787e-01_rb,8.7776e-02_rb, &
     &  6.8361e-02_rb,4.9628e-02_rb,3.2578e-02_rb,3.5117e-03_rb,2.9064e-03_rb,2.2571e-03_rb, &
     &  1.6887e-03_rb,1.0045e-03_rb,3.6853e-04_rb,5.1612e-05_rb/)
      fracrefao(:, 5) = (/ &
     &  1.7338e-01_rb,1.7157e-01_rb,1.5957e-01_rb,1.3571e-01_rb,1.0773e-01_rb,8.7966e-02_rb, &
     &  6.9000e-02_rb,5.0300e-02_rb,3.2813e-02_rb,3.5470e-03_rb,2.9425e-03_rb,2.2552e-03_rb, &
     &  1.7038e-03_rb,1.1025e-03_rb,3.6853e-04_rb,5.1612e-05_rb/)
      fracrefao(:, 6) = (/ &
     &  1.7230e-01_rb,1.7082e-01_rb,1.5917e-01_rb,1.3562e-01_rb,1.0806e-01_rb,8.7635e-02_rb, &
     &  6.9815e-02_rb,5.1155e-02_rb,3.3139e-02_rb,3.6264e-03_rb,2.9436e-03_rb,2.3417e-03_rb, &
     &  1.7731e-03_rb,1.1156e-03_rb,4.4533e-04_rb,5.1612e-05_rb/)
      fracrefao(:, 7) = (/ &
     &  1.7073e-01_rb,1.6961e-01_rb,1.5844e-01_rb,1.3594e-01_rb,1.0821e-01_rb,8.7791e-02_rb, &
     &  7.0502e-02_rb,5.1904e-02_rb,3.4107e-02_rb,3.5888e-03_rb,2.9574e-03_rb,2.5851e-03_rb, &
     &  1.9127e-03_rb,1.1537e-03_rb,4.7789e-04_rb,1.0016e-04_rb/)
      fracrefao(:, 8) = (/ &
     &  1.6700e-01_rb,1.6848e-01_rb,1.5628e-01_rb,1.3448e-01_rb,1.1011e-01_rb,8.9016e-02_rb, &
     &  7.1973e-02_rb,5.2798e-02_rb,3.5650e-02_rb,3.8534e-03_rb,3.4142e-03_rb,2.7799e-03_rb, &
     &  2.1288e-03_rb,1.3043e-03_rb,6.2858e-04_rb,1.0016e-04_rb/)
      fracrefao(:, 9) = (/ &
     &  1.6338e-01_rb,1.5565e-01_rb,1.4470e-01_rb,1.3500e-01_rb,1.1909e-01_rb,9.8312e-02_rb, &
     &  7.9023e-02_rb,5.5728e-02_rb,3.6831e-02_rb,3.6569e-03_rb,3.0552e-03_rb,2.3431e-03_rb, &
     &  1.7088e-03_rb,1.1082e-03_rb,3.6829e-04_rb,5.1612e-05_rb/)

! Planck fraction mapping level : P=4.758820 mb, T = 250.85 K
      fracrefbo(:) = (/ &
     &  1.5411e-01_rb,1.3573e-01_rb,1.2527e-01_rb,1.2698e-01_rb,1.2394e-01_rb,1.0876e-01_rb, &
     &  8.9906e-02_rb,6.9551e-02_rb,4.8240e-02_rb,5.2434e-03_rb,4.3630e-03_rb,3.4262e-03_rb, &
     &  2.5124e-03_rb,1.5479e-03_rb,3.7294e-04_rb,5.1050e-05_rb/)

      end subroutine lw_kgb130

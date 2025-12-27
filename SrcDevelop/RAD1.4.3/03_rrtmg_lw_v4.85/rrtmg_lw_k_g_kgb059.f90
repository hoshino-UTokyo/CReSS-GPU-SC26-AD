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
!!prev      subroutine lw_kgb05
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for selfrefo(:,:), forrefo(:,:)
!!
      subroutine lw_kgb059
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg05, only : fracrefao, fracrefbo, kao, kbo, kao_mo3, &
!!prev                            selfrefo, forrefo, ccl4o
      use rrlw_kg05, only : selfrefo, forrefo

      implicit none
      save

!     The array FORREFO contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  The first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  The second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.0689e-05_rb,1.6987e-05_rb,1.8993e-05_rb,3.4470e-05_rb,4.0873e-05_rb,4.8275e-05_rb, &
     &6.1178e-05_rb,6.4035e-05_rb,6.6253e-05_rb,7.8914e-05_rb,8.1640e-05_rb,7.9738e-05_rb, &
     &7.8492e-05_rb,9.1565e-05_rb,1.0262e-04_rb,1.0368e-04_rb/)
      forrefo(2,:) = (/ &
     &1.1194e-05_rb,1.6128e-05_rb,1.7213e-05_rb,2.6845e-05_rb,4.1361e-05_rb,5.1508e-05_rb, &
     &6.8245e-05_rb,7.4063e-05_rb,7.6273e-05_rb,8.4061e-05_rb,8.2492e-05_rb,8.1720e-05_rb, &
     &7.7626e-05_rb,1.0096e-04_rb,1.0519e-04_rb,1.0631e-04_rb/)
      forrefo(3,:) = (/ &
     &1.0891e-05_rb,1.4933e-05_rb,1.7964e-05_rb,2.2577e-05_rb,4.4290e-05_rb,5.4675e-05_rb, &
     &7.2494e-05_rb,7.8410e-05_rb,7.6948e-05_rb,7.5742e-05_rb,7.7654e-05_rb,8.2760e-05_rb, &
     &7.8443e-05_rb,9.8384e-05_rb,1.0634e-04_rb,1.0838e-04_rb/)
      forrefo(4,:) = (/ &
     &1.1316e-05_rb,1.5470e-05_rb,2.1246e-05_rb,3.3349e-05_rb,4.8704e-05_rb,5.6424e-05_rb, &
     &5.8569e-05_rb,5.8780e-05_rb,6.0358e-05_rb,6.1586e-05_rb,6.4281e-05_rb,6.9333e-05_rb, &
     &7.2763e-05_rb,7.2675e-05_rb,7.3754e-05_rb,1.0131e-04_rb/)

!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 1.27686e-01_rb, 1.09347e-01_rb, 9.36410e-02_rb, 8.01912e-02_rb, 6.86732e-02_rb, &
     & 5.88096e-02_rb, 5.03627e-02_rb, 4.31290e-02_rb, 3.69343e-02_rb, 3.16294e-02_rb/)
      selfrefo(:, 2) = (/ &
     & 1.40051e-01_rb, 1.20785e-01_rb, 1.04170e-01_rb, 8.98402e-02_rb, 7.74816e-02_rb, &
     & 6.68231e-02_rb, 5.76308e-02_rb, 4.97030e-02_rb, 4.28658e-02_rb, 3.69691e-02_rb/)
      selfrefo(:, 3) = (/ &
     & 1.42322e-01_rb, 1.22872e-01_rb, 1.06080e-01_rb, 9.15829e-02_rb, 7.90671e-02_rb, &
     & 6.82616e-02_rb, 5.89329e-02_rb, 5.08790e-02_rb, 4.39258e-02_rb, 3.79228e-02_rb/)
      selfrefo(:, 4) = (/ &
     & 1.53244e-01_rb, 1.33057e-01_rb, 1.15530e-01_rb, 1.00311e-01_rb, 8.70977e-02_rb, &
     & 7.56244e-02_rb, 6.56626e-02_rb, 5.70130e-02_rb, 4.95028e-02_rb, 4.29819e-02_rb/)
      selfrefo(:, 5) = (/ &
     & 1.71011e-01_rb, 1.46680e-01_rb, 1.25810e-01_rb, 1.07910e-01_rb, 9.25563e-02_rb, &
     & 7.93874e-02_rb, 6.80922e-02_rb, 5.84040e-02_rb, 5.00943e-02_rb, 4.29669e-02_rb/)
      selfrefo(:, 6) = (/ &
     & 1.76012e-01_rb, 1.51010e-01_rb, 1.29560e-01_rb, 1.11157e-01_rb, 9.53672e-02_rb, &
     & 8.18207e-02_rb, 7.01984e-02_rb, 6.02270e-02_rb, 5.16720e-02_rb, 4.43322e-02_rb/)
      selfrefo(:, 7) = (/ &
     & 1.85600e-01_rb, 1.59051e-01_rb, 1.36300e-01_rb, 1.16803e-01_rb, 1.00095e-01_rb, &
     & 8.57776e-02_rb, 7.35077e-02_rb, 6.29930e-02_rb, 5.39823e-02_rb, 4.62606e-02_rb/)
      selfrefo(:, 8) = (/ &
     & 1.88931e-01_rb, 1.61727e-01_rb, 1.38440e-01_rb, 1.18506e-01_rb, 1.01442e-01_rb, &
     & 8.68356e-02_rb, 7.43321e-02_rb, 6.36290e-02_rb, 5.44670e-02_rb, 4.66243e-02_rb/)
      selfrefo(:, 9) = (/ &
     & 1.91122e-01_rb, 1.63407e-01_rb, 1.39710e-01_rb, 1.19450e-01_rb, 1.02128e-01_rb, &
     & 8.73176e-02_rb, 7.46552e-02_rb, 6.38290e-02_rb, 5.45728e-02_rb, 4.66589e-02_rb/)
      selfrefo(:,10) = (/ &
     & 1.91334e-01_rb, 1.64872e-01_rb, 1.42070e-01_rb, 1.22421e-01_rb, 1.05490e-01_rb, &
     & 9.09008e-02_rb, 7.83291e-02_rb, 6.74960e-02_rb, 5.81612e-02_rb, 5.01174e-02_rb/)
      selfrefo(:,11) = (/ &
     & 1.89858e-01_rb, 1.63934e-01_rb, 1.41550e-01_rb, 1.22222e-01_rb, 1.05534e-01_rb, &
     & 9.11237e-02_rb, 7.86814e-02_rb, 6.79380e-02_rb, 5.86615e-02_rb, 5.06517e-02_rb/)
      selfrefo(:,12) = (/ &
     & 1.89783e-01_rb, 1.63757e-01_rb, 1.41300e-01_rb, 1.21923e-01_rb, 1.05203e-01_rb, &
     & 9.07760e-02_rb, 7.83274e-02_rb, 6.75860e-02_rb, 5.83176e-02_rb, 5.03202e-02_rb/)
      selfrefo(:,13) = (/ &
     & 1.87534e-01_rb, 1.62016e-01_rb, 1.39970e-01_rb, 1.20924e-01_rb, 1.04470e-01_rb, &
     & 9.02541e-02_rb, 7.79730e-02_rb, 6.73630e-02_rb, 5.81967e-02_rb, 5.02778e-02_rb/)
      selfrefo(:,14) = (/ &
     & 1.99128e-01_rb, 1.71410e-01_rb, 1.47550e-01_rb, 1.27011e-01_rb, 1.09332e-01_rb, &
     & 9.41131e-02_rb, 8.10128e-02_rb, 6.97360e-02_rb, 6.00289e-02_rb, 5.16731e-02_rb/)
      selfrefo(:,15) = (/ &
     & 1.99460e-01_rb, 1.72342e-01_rb, 1.48910e-01_rb, 1.28664e-01_rb, 1.11171e-01_rb, &
     & 9.60560e-02_rb, 8.29962e-02_rb, 7.17120e-02_rb, 6.19620e-02_rb, 5.35376e-02_rb/)
      selfrefo(:,16) = (/ &
     & 1.99906e-01_rb, 1.72737e-01_rb, 1.49260e-01_rb, 1.28974e-01_rb, 1.11445e-01_rb, &
     & 9.62982e-02_rb, 8.32102e-02_rb, 7.19010e-02_rb, 6.21288e-02_rb, 5.36848e-02_rb/)

      end subroutine lw_kgb059

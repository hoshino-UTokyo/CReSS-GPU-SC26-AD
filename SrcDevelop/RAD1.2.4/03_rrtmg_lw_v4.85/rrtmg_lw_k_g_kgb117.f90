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
!!prev      subroutine lw_kgb11
!!
!!          divided by M.Yoshioka 2015/02/12
!!          set parameters for kao_mo2(:,:)
!!
      subroutine lw_kgb117
! **************************************************************************

      use parkind, only : im => kind_im, rb => kind_rb 
!!prev      use rrlw_kg11, only : fracrefao, fracrefbo, kao, kbo, kao_mo2, &
!!prev                            kbo_mo2, selfrefo, forrefo
      use rrlw_kg11, only :  kao_mo2

      implicit none
      save

!     The array KAO_Mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   The first index refers to temperature 
!     in 7.2 degree increments.  For instance, JT = 1 refers to a 
!     temperature of 188.0, JT = 2 refers to 195.2, etc. The second index 
!     runs over the g-channel (1 to 16).

      kao_mo2(:, 1) = (/ &
     & 2.31723e-06_rb, 2.28697e-06_rb, 2.25710e-06_rb, 2.22762e-06_rb, 2.19852e-06_rb, &
     & 2.16981e-06_rb, 2.14147e-06_rb, 2.11350e-06_rb, 2.08590e-06_rb, 2.05865e-06_rb, &
     & 2.03176e-06_rb, 2.00523e-06_rb, 1.97904e-06_rb, 1.95319e-06_rb, 1.92768e-06_rb, &
     & 1.90250e-06_rb, 1.87765e-06_rb, 1.85313e-06_rb, 1.82893e-06_rb/)
      kao_mo2(:, 2) = (/ &
     & 1.81980e-06_rb, 1.81352e-06_rb, 1.80726e-06_rb, 1.80101e-06_rb, 1.79479e-06_rb, &
     & 1.78860e-06_rb, 1.78242e-06_rb, 1.77626e-06_rb, 1.77013e-06_rb, 1.76402e-06_rb, &
     & 1.75792e-06_rb, 1.75185e-06_rb, 1.74580e-06_rb, 1.73977e-06_rb, 1.73377e-06_rb, &
     & 1.72778e-06_rb, 1.72181e-06_rb, 1.71587e-06_rb, 1.70994e-06_rb/)
      kao_mo2(:, 3) = (/ &
     & 2.26922e-06_rb, 2.25413e-06_rb, 2.23914e-06_rb, 2.22425e-06_rb, 2.20945e-06_rb, &
     & 2.19476e-06_rb, 2.18016e-06_rb, 2.16566e-06_rb, 2.15126e-06_rb, 2.13695e-06_rb, &
     & 2.12274e-06_rb, 2.10862e-06_rb, 2.09459e-06_rb, 2.08066e-06_rb, 2.06683e-06_rb, &
     & 2.05308e-06_rb, 2.03942e-06_rb, 2.02586e-06_rb, 2.01239e-06_rb/)
      kao_mo2(:, 4) = (/ &
     & 2.15555e-06_rb, 2.14539e-06_rb, 2.13527e-06_rb, 2.12520e-06_rb, 2.11517e-06_rb, &
     & 2.10520e-06_rb, 2.09527e-06_rb, 2.08538e-06_rb, 2.07555e-06_rb, 2.06576e-06_rb, &
     & 2.05601e-06_rb, 2.04631e-06_rb, 2.03666e-06_rb, 2.02706e-06_rb, 2.01749e-06_rb, &
     & 2.00798e-06_rb, 1.99851e-06_rb, 1.98908e-06_rb, 1.97970e-06_rb/)
      kao_mo2(:, 5) = (/ &
     & 2.05821e-06_rb, 2.04914e-06_rb, 2.04011e-06_rb, 2.03111e-06_rb, 2.02216e-06_rb, &
     & 2.01324e-06_rb, 2.00437e-06_rb, 1.99553e-06_rb, 1.98673e-06_rb, 1.97798e-06_rb, &
     & 1.96926e-06_rb, 1.96057e-06_rb, 1.95193e-06_rb, 1.94333e-06_rb, 1.93476e-06_rb, &
     & 1.92623e-06_rb, 1.91774e-06_rb, 1.90928e-06_rb, 1.90087e-06_rb/)
      kao_mo2(:, 6) = (/ &
     & 2.20148e-06_rb, 2.18998e-06_rb, 2.17854e-06_rb, 2.16717e-06_rb, 2.15585e-06_rb, &
     & 2.14459e-06_rb, 2.13339e-06_rb, 2.12225e-06_rb, 2.11117e-06_rb, 2.10014e-06_rb, &
     & 2.08918e-06_rb, 2.07827e-06_rb, 2.06741e-06_rb, 2.05662e-06_rb, 2.04588e-06_rb, &
     & 2.03519e-06_rb, 2.02457e-06_rb, 2.01399e-06_rb, 2.00348e-06_rb/)
      kao_mo2(:, 7) = (/ &
     & 2.28960e-06_rb, 2.27651e-06_rb, 2.26349e-06_rb, 2.25054e-06_rb, 2.23767e-06_rb, &
     & 2.22487e-06_rb, 2.21215e-06_rb, 2.19950e-06_rb, 2.18692e-06_rb, 2.17441e-06_rb, &
     & 2.16198e-06_rb, 2.14961e-06_rb, 2.13732e-06_rb, 2.12509e-06_rb, 2.11294e-06_rb, &
     & 2.10085e-06_rb, 2.08884e-06_rb, 2.07689e-06_rb, 2.06501e-06_rb/)
      kao_mo2(:, 8) = (/ &
     & 2.28564e-06_rb, 2.27363e-06_rb, 2.26168e-06_rb, 2.24980e-06_rb, 2.23798e-06_rb, &
     & 2.22622e-06_rb, 2.21452e-06_rb, 2.20288e-06_rb, 2.19131e-06_rb, 2.17980e-06_rb, &
     & 2.16834e-06_rb, 2.15695e-06_rb, 2.14562e-06_rb, 2.13434e-06_rb, 2.12313e-06_rb, &
     & 2.11197e-06_rb, 2.10087e-06_rb, 2.08984e-06_rb, 2.07886e-06_rb/)
      kao_mo2(:, 9) = (/ &
     & 2.28505e-06_rb, 2.27395e-06_rb, 2.26291e-06_rb, 2.25192e-06_rb, 2.24099e-06_rb, &
     & 2.23011e-06_rb, 2.21928e-06_rb, 2.20850e-06_rb, 2.19778e-06_rb, 2.18711e-06_rb, &
     & 2.17649e-06_rb, 2.16592e-06_rb, 2.15540e-06_rb, 2.14494e-06_rb, 2.13452e-06_rb, &
     & 2.12416e-06_rb, 2.11385e-06_rb, 2.10358e-06_rb, 2.09337e-06_rb/)
      kao_mo2(:,10) = (/ &
     & 2.25915e-06_rb, 2.24938e-06_rb, 2.23965e-06_rb, 2.22997e-06_rb, 2.22032e-06_rb, &
     & 2.21072e-06_rb, 2.20116e-06_rb, 2.19164e-06_rb, 2.18216e-06_rb, 2.17272e-06_rb, &
     & 2.16333e-06_rb, 2.15397e-06_rb, 2.14465e-06_rb, 2.13538e-06_rb, 2.12614e-06_rb, &
     & 2.11695e-06_rb, 2.10779e-06_rb, 2.09868e-06_rb, 2.08960e-06_rb/)
      kao_mo2(:,11) = (/ &
     & 2.52025e-06_rb, 2.50423e-06_rb, 2.48831e-06_rb, 2.47249e-06_rb, 2.45677e-06_rb, &
     & 2.44115e-06_rb, 2.42563e-06_rb, 2.41021e-06_rb, 2.39489e-06_rb, 2.37967e-06_rb, &
     & 2.36454e-06_rb, 2.34951e-06_rb, 2.33457e-06_rb, 2.31973e-06_rb, 2.30498e-06_rb, &
     & 2.29033e-06_rb, 2.27577e-06_rb, 2.26130e-06_rb, 2.24692e-06_rb/)
      kao_mo2(:,12) = (/ &
     & 2.52634e-06_rb, 2.51180e-06_rb, 2.49735e-06_rb, 2.48299e-06_rb, 2.46871e-06_rb, &
     & 2.45451e-06_rb, 2.44039e-06_rb, 2.42635e-06_rb, 2.41239e-06_rb, 2.39851e-06_rb, &
     & 2.38472e-06_rb, 2.37100e-06_rb, 2.35736e-06_rb, 2.34380e-06_rb, 2.33032e-06_rb, &
     & 2.31691e-06_rb, 2.30358e-06_rb, 2.29033e-06_rb, 2.27716e-06_rb/)
      kao_mo2(:,13) = (/ &
     & 2.66614e-06_rb, 2.64897e-06_rb, 2.63191e-06_rb, 2.61496e-06_rb, 2.59812e-06_rb, &
     & 2.58138e-06_rb, 2.56476e-06_rb, 2.54824e-06_rb, 2.53183e-06_rb, 2.51552e-06_rb, &
     & 2.49932e-06_rb, 2.48322e-06_rb, 2.46723e-06_rb, 2.45134e-06_rb, 2.43555e-06_rb, &
     & 2.41987e-06_rb, 2.40428e-06_rb, 2.38880e-06_rb, 2.37341e-06_rb/)
      kao_mo2(:,14) = (/ &
     & 2.96755e-06_rb, 2.94803e-06_rb, 2.92864e-06_rb, 2.90937e-06_rb, 2.89023e-06_rb, &
     & 2.87122e-06_rb, 2.85233e-06_rb, 2.83357e-06_rb, 2.81493e-06_rb, 2.79641e-06_rb, &
     & 2.77802e-06_rb, 2.75974e-06_rb, 2.74159e-06_rb, 2.72355e-06_rb, 2.70563e-06_rb, &
     & 2.68784e-06_rb, 2.67015e-06_rb, 2.65259e-06_rb, 2.63514e-06_rb/)
      kao_mo2(:,15) = (/ &
     & 1.30668e-06_rb, 1.31378e-06_rb, 1.32091e-06_rb, 1.32808e-06_rb, 1.33530e-06_rb, &
     & 1.34255e-06_rb, 1.34984e-06_rb, 1.35717e-06_rb, 1.36454e-06_rb, 1.37195e-06_rb, &
     & 1.37941e-06_rb, 1.38690e-06_rb, 1.39443e-06_rb, 1.40200e-06_rb, 1.40962e-06_rb, &
     & 1.41727e-06_rb, 1.42497e-06_rb, 1.43271e-06_rb, 1.44049e-06_rb/)
      kao_mo2(:,16) = (/ &
     & 5.99001e-07_rb, 6.16844e-07_rb, 6.35219e-07_rb, 6.54141e-07_rb, 6.73626e-07_rb, &
     & 6.93692e-07_rb, 7.14356e-07_rb, 7.35635e-07_rb, 7.57548e-07_rb, 7.80114e-07_rb, &
     & 8.03352e-07_rb, 8.27282e-07_rb, 8.51925e-07_rb, 8.77302e-07_rb, 9.03435e-07_rb, &
     & 9.30347e-07_rb, 9.58060e-07_rb, 9.86599e-07_rb, 1.01599e-06_rb/)

      end subroutine lw_kgb117

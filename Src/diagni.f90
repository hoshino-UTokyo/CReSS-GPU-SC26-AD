!***********************************************************************
      module m_diagni
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2001/10/18, 2001/11/20, 2002/01/07,
!                   2002/01/15, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2003/12/12, 2004/04/01, 2004/05/31,
!                   2004/06/10, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2005/04/04, 2005/09/30, 2005/10/05,
!                   2006/02/13, 2007/10/19, 2007/11/26, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2009/11/05, 2011/03/18,
!                   2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the diagnostic concentrations of all categories of the ice
!     hydrometeor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_comphy
      use m_getiname
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diagni, s_diagni

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diagni

        module procedure s_diagni

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max
      intrinsic min
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_diagni(fphaiopt,ni,nj,nk,nqi,nni,rbr,qice,nidia)
!***********************************************************************

! Input variables

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

! Output variable

      real, intent(out) :: nidia(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Diagnostic concentrations of ice hydrometeor

! Internal shared variables

      integer haiopt   ! Option for additional hail processes

      real miiv        ! 4.0 / (3.0 x mimax)

      real msmiv2      ! 0.01 / msmax
      real ms0iv2      ! 100.0 / ms0

      real mgmiv2      ! 0.01 / mgmax
      real mg0iv2      ! 100.0 / mg0

      real mhmiv2      ! 0.01 / mhmax
      real mh0iv2      ! 100.0 / mh0

      real cdiaqs      ! Coefficient of mean diameter of snow
      real cdiaqg      ! Coefficient of mean diameter of graupel
      real cdiaqh      ! Coefficient of mean diameter of hail

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real rbv         ! Inverse of base state density


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_diagni = 0
      integer, parameter :: DUMP_TARGET_diagni = 1
      logical, save :: dump_done_diagni = .false.

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fphaiopt,haiopt)

! -----

! Set the common used variables.

      miiv=4.e0/(3.e0*mimax)

      msmiv2=1.e-2/msmax
      ms0iv2=1.e2/ms0

      mgmiv2=1.e-2/mgmax
      mg0iv2=1.e2/mg0

      mhmiv2=1.e-2/mhmax
      mh0iv2=1.e2/mh0

      cdiaqs=ns0*ns0*ns0/(cc*rhos)
      cdiaqg=ng0*ng0*ng0/(cc*rhog)
      cdiaqh=nh0*nh0*nh0/(cc*rhoh)

! -----

!!! Get the diagnostic concentrations of all categories of the ice
!!! hydrometeor.

!@llm start meta_info ----------------------------------------------------
! Location: diagni.f90 :: s_diagni
! Summary : Calculate diagnostic concentrations for all ice hydrometeor
!           categories (cloud ice, snow, graupel, hail) from mixing ratios.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic max, min, sqrt functions (GPU-compatible)
!   - Private variable k for outer loop; rbv for local scalar
!   - Writes to nidia output array for multiple ice categories
!   - Conditional branch based on haiopt (3 vs 4 categories)
!   - Independent point-wise operations per grid cell
! Next:
!   - Direct conversion to OpenACC with collapsed loops
!   - Handle haiopt conditional outside kernel or use single kernel with masking
!   - Data managed automatically via Unified Memory
! Runtime:
!   - Calls: 1
!   - AvgLoops: 102.4M
!   - TotalTime: 0.012s (0.00%)
!   - AvgTime: 11.720ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('diagni.f90', 's_diagni', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)

! Dump input data at target call
dump_call_count_diagni = dump_call_count_diagni + 1
if (dump_call_count_diagni == DUMP_TARGET_diagni .and. .not. dump_done_diagni) then
  call dump_init('diagni')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('nqi', nqi)
  call dump_scalar_i('nni', nni)
  call dump_scalar_i('haiopt', haiopt)
  call dump_array_3d('rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_4d('qice.bin', qice, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)
  call dump_scalar_r('cdiaqg', cdiaqg)
  call dump_scalar_r('cdiaqh', cdiaqh)
  call dump_scalar_r('cdiaqs', cdiaqs)
  call dump_scalar_r('mg0iv2', mg0iv2)
  call dump_scalar_r('mgmiv2', mgmiv2)
  call dump_scalar_r('mh0iv2', mh0iv2)
  call dump_scalar_r('mhmiv2', mhmiv2)
  call dump_scalar_r('miiv', miiv)
  call dump_scalar_r('ms0iv2', ms0iv2)
  call dump_scalar_r('msmiv2', msmiv2)
end if

call profile_start(prof_id1)

#if defined(USE_GPU) && !defined(DISABLE_GPU_083)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    ! Local variables

    ! Derived constants

    ! Set derived constants
    miiv = 4.e0 / (3.e0 * mimax)
    msmiv2 = 1.e-2 / msmax
    ms0iv2 = 1.e2 / ms0
    mgmiv2 = 1.e-2 / mgmax
    mg0iv2 = 1.e2 / mg0
    mhmiv2 = 1.e-2 / mhmax
    mh0iv2 = 1.e2 / mh0
    cdiaqs = ns0*ns0*ns0 / (cc*rhos)
    cdiaqg = ng0*ng0*ng0 / (cc*rhog)
    cdiaqh = nh0*nh0*nh0 / (cc*rhoh)

    if (haiopt == 0) then

      !$acc kernels
      !$acc loop independent collapse(3) private(rbv)
      do k = 1, nk-1
        do j = 1, nj-1
          do i = 1, ni-1
            rbv = 1.e0 / rbr(i,j,k)

            ! Cloud ice
            nidia(i,j,k,1) = miiv * qice(i,j,k,1)

            ! Snow
            nidia(i,j,k,2) = sqrt(sqrt(cdiaqs*rbr(i,j,k)*qice(i,j,k,2))) * rbv
            nidia(i,j,k,2) = min(max(nidia(i,j,k,2), msmiv2*qice(i,j,k,2)), ms0iv2*qice(i,j,k,2))

            ! Graupel
            nidia(i,j,k,3) = sqrt(sqrt(cdiaqg*rbr(i,j,k)*qice(i,j,k,3))) * rbv
            nidia(i,j,k,3) = min(max(nidia(i,j,k,3), mgmiv2*qice(i,j,k,3)), mg0iv2*qice(i,j,k,3))
          end do
        end do
      end do
      !$acc end kernels

    else

      !$acc kernels
      !$acc loop independent collapse(3) private(rbv)
      do k = 1, nk-1
        do j = 1, nj-1
          do i = 1, ni-1
            rbv = 1.e0 / rbr(i,j,k)

            ! Cloud ice
            nidia(i,j,k,1) = miiv * qice(i,j,k,1)

            ! Snow
            nidia(i,j,k,2) = sqrt(sqrt(cdiaqs*rbr(i,j,k)*qice(i,j,k,2))) * rbv
            nidia(i,j,k,2) = min(max(nidia(i,j,k,2), msmiv2*qice(i,j,k,2)), ms0iv2*qice(i,j,k,2))

            ! Graupel
            nidia(i,j,k,3) = sqrt(sqrt(cdiaqg*rbr(i,j,k)*qice(i,j,k,3))) * rbv
            nidia(i,j,k,3) = min(max(nidia(i,j,k,3), mgmiv2*qice(i,j,k,3)), mg0iv2*qice(i,j,k,3))

            ! Hail
            nidia(i,j,k,4) = sqrt(sqrt(cdiaqh*rbr(i,j,k)*qice(i,j,k,4))) * rbv
            nidia(i,j,k,4) = min(max(nidia(i,j,k,4), mhmiv2*qice(i,j,k,4)), mh0iv2*qice(i,j,k,4))
          end do
        end do
      end do
      !$acc end kernels

    end if

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared) private(k)

!! Get the diagnostic concentrations of the cloud ice, snow and graupel.

      if(haiopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,rbv)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the inverse of base state density.

            rbv=1.e0/rbr(i,j,k)

! -----

! Get the diagnostic concentrations of cloud ice.

            nidia(i,j,k,1)=miiv*qice(i,j,k,1)

! -----

! Get the diagnostic concentrations of snow.

            nidia(i,j,k,2)=sqrt(                                        &
     &        sqrt(cdiaqs*rbr(i,j,k)*qice(i,j,k,2)))*rbv

            nidia(i,j,k,2)=min(max(                                     &
     &        nidia(i,j,k,2),msmiv2*qice(i,j,k,2)),ms0iv2*qice(i,j,k,2))

! -----

! Get the diagnostic concentrations of graupel.

            nidia(i,j,k,3)=sqrt(                                        &
     &        sqrt(cdiaqg*rbr(i,j,k)*qice(i,j,k,3)))*rbv

            nidia(i,j,k,3)=min(max(                                     &
     &        nidia(i,j,k,3),mgmiv2*qice(i,j,k,3)),mg0iv2*qice(i,j,k,3))

! -----

          end do
          end do

!$omp end do

        end do

!! -----

!! Get the diagnostic concentrations of the cloud ice, snow, graupel and
!! hail.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,rbv)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the inverse of base state density.

            rbv=1.e0/rbr(i,j,k)

! -----

! Get the diagnostic concentrations of cloud ice.

            nidia(i,j,k,1)=miiv*qice(i,j,k,1)

! -----

! Get the diagnostic concentrations of snow.

            nidia(i,j,k,2)=sqrt(                                        &
     &        sqrt(cdiaqs*rbr(i,j,k)*qice(i,j,k,2)))*rbv

            nidia(i,j,k,2)=min(max(                                     &
     &        nidia(i,j,k,2),msmiv2*qice(i,j,k,2)),ms0iv2*qice(i,j,k,2))

! -----

! Get the diagnostic concentrations of graupel.

            nidia(i,j,k,3)=sqrt(                                        &
     &        sqrt(cdiaqg*rbr(i,j,k)*qice(i,j,k,3)))*rbv

            nidia(i,j,k,3)=min(max(                                     &
     &        nidia(i,j,k,3),mgmiv2*qice(i,j,k,3)),mg0iv2*qice(i,j,k,3))

! -----

! Get the diagnostic concentrations of hail.

            nidia(i,j,k,4)=sqrt(                                        &
     &        sqrt(cdiaqh*rbr(i,j,k)*qice(i,j,k,4)))*rbv

            nidia(i,j,k,4)=min(max(                                     &
     &        nidia(i,j,k,4),mhmiv2*qice(i,j,k,4)),mh0iv2*qice(i,j,k,4))

! -----

          end do
          end do

!$omp end do

        end do

      end if

!! -----

!$omp end parallel
#endif

call profile_stop(prof_id1, loop_len)

! Dump output data at target call
if (dump_call_count_diagni == DUMP_TARGET_diagni .and. .not. dump_done_diagni) then
  call dump_array_4d('nidia_ref.bin', nidia, 0, ni+1, 0, nj+1, 1, nk, 1, nni)
  call dump_finalize()
  dump_done_diagni = .true.
end if

!!! -----

      end subroutine s_diagni

!-----7--------------------------------------------------------------7--

      end module m_diagni

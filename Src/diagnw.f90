!***********************************************************************
      module m_diagnw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/11/26
!     Modification: 2008/05/02, 2008/08/25, 2009/01/30, 2009/02/27,
!                   2009/11/05, 2011/03/18, 2011/09/22, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the diagnostic concentrations of the water hydrometeor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_dump_kernel
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diagnw, s_diagnw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diagnw

        module procedure s_diagnw

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
      subroutine s_diagnw(ni,nj,nk,nqw,nnw,rbr,qwtr,nwdia)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

! Output variable

      real, intent(out) :: nwdia(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Diagnostic concentrations of water hydrometeor

! Internal shared variables

      real mrmiv2      ! 0.01 / mrmax
      real mr0iv2      ! 100.0 / mr0

      real cdiaqr      ! Coefficient of mean diameter of rain water

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real rbv         ! Inverse of base state density


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_diagnw = 0
      integer, parameter :: DUMP_TARGET_diagnw = 360
      logical, save :: dump_done_diagnw = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variables.

      mrmiv2=1.e-2/mrmax
      mr0iv2=1.e2/mr0

      cdiaqr=nr0*nr0*nr0/(cc*rhow)

! -----

!! Get the diagnostic concentrations of the water hydrometeor.

!@llm start meta_info ----------------------------------------------------
! Location: diagnw.f90 :: s_diagnw
! Summary : Compute diagnostic concentrations of cloud water and rain water
!           based on base state density and water hydrometeor mixing ratios.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region (only intrinsics: sqrt, min, max)
!   - No global/module variable writes
!   - No synchronization constructs (barriers, critical, atomic)
!   - Simple 3D loop with k-loop outside, j-i loops inside with schedule(runtime)
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - Consider collapse(3) after loop restructuring for better GPU utilization
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 2.073s (0.07%)
!   - AvgTime: 5.759ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('diagnw.f90', 's_diagnw', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_diagnw = dump_call_count_diagnw + 1
if (dump_call_count_diagnw == DUMP_TARGET_diagnw .and. .not. dump_done_diagnw) then
  call dump_init('diagnw')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('nqw', nqw)
  call dump_scalar_i('nnw', nnw)
  call dump_array_3d('rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_4d('qwtr.bin', qwtr, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
  call dump_scalar_r('cdiaqr', cdiaqr)
  call dump_scalar_r('mr0iv2', mr0iv2)
  call dump_scalar_r('mrmiv2', mrmiv2)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_085)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------

    ! Set the common used variables
    mrmiv2 = 1.0e-2 / mrmax
    mr0iv2 = 1.0e2 / mr0
    cdiaqr = nr0 * nr0 * nr0 / (cc * rhow)

    !$acc kernels
    !$acc loop independent collapse(3) private(rbv)
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          ! Calculate the inverse of base state density
          rbv = 1.0 / rbr(i,j,k)

          ! Get the diagnostic concentrations of cloud water
          nwdia(i,j,k,1) = nclcst * rbv

          ! Get the diagnostic concentrations of rain water
          nwdia(i,j,k,2) = sqrt(sqrt(cdiaqr * rbr(i,j,k) * qwtr(i,j,k,2))) * rbv
          nwdia(i,j,k,2) = min(max(nwdia(i,j,k,2), mrmiv2*qwtr(i,j,k,2)), &
                               mr0iv2*qwtr(i,j,k,2))
        end do
      end do
    end do
    !$acc end kernels

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,rbv)

        do j=1,nj-1
        do i=1,ni-1

! Calculate the inverse of base state density.

          rbv=1.e0/rbr(i,j,k)

! -----

! Get the diagnostic concentrations of cloud water.

          nwdia(i,j,k,1)=nclcst*rbv

! -----

! Get the diagnostic concentrations of rain water.

          nwdia(i,j,k,2)=sqrt(                                          &
     &      sqrt(cdiaqr*rbr(i,j,k)*qwtr(i,j,k,2)))*rbv

          nwdia(i,j,k,2)=min(max(                                       &
     &      nwdia(i,j,k,2),mrmiv2*qwtr(i,j,k,2)),mr0iv2*qwtr(i,j,k,2))

! -----

        end do
        end do

!$omp end do

      end do

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_diagnw == DUMP_TARGET_diagnw .and. .not. dump_done_diagnw) then
  call dump_array_4d('nwdia_ref.bin', nwdia, 0, ni+1, 0, nj+1, 1, nk, 1, nnw)
  call dump_finalize()
  dump_done_diagnw = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_diagnw

!-----7--------------------------------------------------------------7--

      end module m_diagnw

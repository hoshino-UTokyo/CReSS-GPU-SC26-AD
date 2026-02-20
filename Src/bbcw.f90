!***********************************************************************
      module m_bbcw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/07/01
!     Modification: 2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the bottom boundary conditions for the z components of
!     velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bbcw, s_bbcw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bbcw

        module procedure s_bbcw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_bbcw(fpmpopt,fpmfcopt,ni,nj,nk,j31,j32,mf,aa,        &
     &                  uf,vf,wf,j31u2,j32v2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: aa(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(in) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(in) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

! Input and output variable

      real, intent(inout) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      real, intent(inout) :: j31u2(0:ni+1,0:nj+1)
                       ! 2.0 x j31 x u

      real, intent(inout) :: j32v2(0:ni+1,0:nj+1)
                       ! 2.0 x j32 x v

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_bbcw = 0
      integer, parameter :: DUMP_TARGET_bbcw = 14400
      logical, save :: dump_done_bbcw = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)

! -----

! Set the bottom boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: bbcw.f90 :: s_bbcw
! Summary : Set bottom boundary conditions for vertical velocity (wf)
!           using terrain-following coordinate transformations
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to intent(inout) arrays wf, j31u2, j32v2
!   - Multiple conditional branches based on mfcopt and mpopt options
!   - Uses work arrays j31u2, j32v2 for intermediate calculations
! Next:
!   - GPU port requires handling conditional branches
!   - Consider separating compute kernels by mfcopt/mpopt case
!   - Work arrays j31u2, j32v2 should be device-resident
!   - Branch divergence from mfcopt/mpopt may impact performance
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 807.3K
!   - TotalTime: 1.749s (0.06%)
!   - AvgTime: 0.121ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('bbcw.f90', 's_bbcw', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_bbcw = dump_call_count_bbcw + 1
if (dump_call_count_bbcw == DUMP_TARGET_bbcw .and. .not. dump_done_bbcw) then
  call dump_init('bbcw')
  call dump_scalar_i('mpopt', mpopt)
  call dump_scalar_i('mfcopt', mfcopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('mf.bin', mf, 0, ni+1, 0, nj+1)
  call dump_array_3d('aa.bin', aa, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('uf.bin', uf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vf.bin', vf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wf_in.bin', wf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('j31u2_in.bin', j31u2, 0, ni+1, 0, nj+1)
  call dump_array_2d('j32v2_in.bin', j32v2, 0, ni+1, 0, nj+1)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_025)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni
        j31u2(i,j) = (uf(i,j,1) + uf(i,j,2)) * j31(i,j,2)
      end do
    end do
    !$acc end kernels

    ! Compute j32v2
    !$acc kernels
    !$acc loop independent
    do j = 1, nj
      !$acc loop independent
      do i = 1, ni-1
        j32v2(i,j) = (vf(i,j,1) + vf(i,j,2)) * j32(i,j,2)
      end do
    end do
    !$acc end kernels

    ! Compute wf based on mfcopt and mpopt
    if (mfcopt == 0) then
      !$acc kernels
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          wf(i,j,3) = wf(i,j,3) + 0.25e0 * aa(i,j,3) &
               * ((j31u2(i,j) + j31u2(i+1,j)) + (j32v2(i,j) + j32v2(i,j+1)))
        end do
      end do
      !$acc end kernels
    else
      if (mpopt == 0 .or. mpopt == 10) then
        !$acc kernels
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            wf(i,j,3) = wf(i,j,3) &
                 + 0.25e0 * aa(i,j,3) * (mf(i,j) * (j31u2(i,j) + j31u2(i+1,j)) &
                 + (j32v2(i,j) + j32v2(i,j+1)))
          end do
        end do
        !$acc end kernels
      else if (mpopt == 5) then
        !$acc kernels
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            wf(i,j,3) = wf(i,j,3) &
                 + 0.25e0 * aa(i,j,3) * ((j31u2(i,j) + j31u2(i+1,j)) &
                 + mf(i,j) * (j32v2(i,j) + j32v2(i,j+1)))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 2, ni-2
            wf(i,j,3) = wf(i,j,3) + 0.25e0 * mf(i,j) * aa(i,j,3) &
                 * ((j31u2(i,j) + j31u2(i+1,j)) + (j32v2(i,j) + j32v2(i,j+1)))
          end do
        end do
        !$acc end kernels
      end if
    end if

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni
        j31u2(i,j)=(uf(i,j,1)+uf(i,j,2))*j31(i,j,2)
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj
      do i=1,ni-1
        j32v2(i,j)=(vf(i,j,1)+vf(i,j,2))*j32(i,j,2)
      end do
      end do

!$omp end do

      if(mfcopt.eq.0) then

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wf(i,j,3)=wf(i,j,3)+.25e0*aa(i,j,3)                           &
     &      *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
        end do
        end do

!$omp end do

      else

        if(mpopt.eq.0.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wf(i,j,3)=wf(i,j,3)                                         &
     &        +.25e0*aa(i,j,3)*(mf(i,j)*(j31u2(i,j)+j31u2(i+1,j))       &
     &        +(j32v2(i,j)+j32v2(i,j+1)))
          end do
          end do

!$omp end do

        else if(mpopt.eq.5) then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wf(i,j,3)=wf(i,j,3)                                         &
     &        +.25e0*aa(i,j,3)*((j31u2(i,j)+j31u2(i+1,j))               &
     &        +mf(i,j)*(j32v2(i,j)+j32v2(i,j+1)))
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wf(i,j,3)=wf(i,j,3)+.25e0*mf(i,j)*aa(i,j,3)                 &
     &        *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
          end do
          end do

!$omp end do

        end if

      end if

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_bbcw == DUMP_TARGET_bbcw .and. .not. dump_done_bbcw) then
  call dump_array_3d('wf_ref.bin', wf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('j31u2_ref.bin', j31u2, 0, ni+1, 0, nj+1)
  call dump_array_2d('j32v2_ref.bin', j32v2, 0, ni+1, 0, nj+1)
  call dump_finalize()
  dump_done_bbcw = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_bbcw

!-----7--------------------------------------------------------------7--

      end module m_bbcw

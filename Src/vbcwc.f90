!***********************************************************************
      module m_vbcwc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/08/18,
!                   1999/08/23, 1999/09/06, 1999/10/12, 2000/01/17,
!                   2001/06/06, 2001/12/11, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the vertical boundary conditions for the zeta components of
!     contravariant velocity.

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

      public :: vbcwc, s_vbcwc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vbcwc

        module procedure s_vbcwc

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
      subroutine s_vbcwc(fpbbc,fptbc,ni,nj,nk,wc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpbbc
                       ! Formal parameter of unique index of bbc

      integer, intent(in) :: fptbc
                       ! Formal parameter of unique index of tbc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Internal shared variables

      integer bbc      ! Option for bottom boundary conditions
      integer tbc      ! Option for top boundary conditions

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_vbcwc = 0
      integer, parameter :: DUMP_TARGET_vbcwc = 15121
      logical, save :: dump_done_vbcwc = .false.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpbbc,bbc)
      call getiname(fptbc,tbc)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

! -----

!! Set the bottom and top boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: vbcwc.f90 :: s_vbcwc
! Summary : Set vertical boundary conditions (bottom/top) for zeta contravariant velocity
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple array writes to wc (inout) with conditional branches
!   - Multiple !$omp do regions with private(i,j) and schedule(runtime)
!   - No synchronization constructs beyond implicit barriers at end do
! Next:
!   - Direct OpenACC with Unified Memory (no explicit data transfer needed)
!   - Consider collapsing i,j loops and using teams distribute
! Runtime:
!   - Calls: 15121
!   - AvgLoops: 806.4K
!   - TotalTime: 0.833s (0.03%)
!   - AvgTime: 0.055ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('vbcwc.f90', 's_vbcwc', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)

! Dump input data at target call
dump_call_count_vbcwc = dump_call_count_vbcwc + 1
if (dump_call_count_vbcwc == DUMP_TARGET_vbcwc .and. .not. dump_done_vbcwc) then
  call dump_init('vbcwc')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('bbc', bbc)
  call dump_scalar_i('tbc', tbc)
  call dump_array_3d('wc_in.bin', wc, 0, ni+1, 0, nj+1, 1, nk)
end if

call profile_start(prof_id1)

!$omp parallel default(shared)

! Set the bottom boundary conditions.

      if(bbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wc(i,j,1)=-wc(i,j,3)
          wc(i,j,2)=0.e0
        end do
        end do

!$omp end do

      else if(bbc.eq.3) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wc(i,j,1)=wc(i,j,2)
        end do
        end do

!$omp end do

      end if

! -----

! Set the top boundary conditions.

      if(tbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wc(i,j,nk)=-wc(i,j,nkm2)
          wc(i,j,nkm1)=0.e0
        end do
        end do

!$omp end do

      else if(tbc.ge.3) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wc(i,j,nk)=wc(i,j,nkm1)
        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! Dump output data at target call
if (dump_call_count_vbcwc == DUMP_TARGET_vbcwc .and. .not. dump_done_vbcwc) then
  call dump_array_3d('wc_ref.bin', wc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_vbcwc = .true.
end if

!! -----

      end subroutine s_vbcwc

!-----7--------------------------------------------------------------7--

      end module m_vbcwc

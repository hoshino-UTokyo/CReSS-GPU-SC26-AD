!***********************************************************************
      module m_vbcp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/08/18,
!                   1999/08/23, 1999/09/06, 1999/10/12, 1999/12/06,
!                   2000/01/17, 2000/03/17, 2001/04/15, 2001/06/29,
!                   2001/07/13, 2001/08/07, 2001/12/11, 2002/04/02,
!                   2002/07/23, 2002/08/15, 2003/04/30, 2003/05/19,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the vertical boundary conditions for the pressure.

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

      public :: vbcp, s_vbcp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vbcp

        module procedure s_vbcp

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
      subroutine s_vbcp(fpbbc,ni,nj,nk,ppf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpbbc
                       ! Formal parameter of unique index of bbc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

! Internal shared variables

      integer bbc      ! Option for bottom boundary conditions

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_vbcp = 0
      integer, parameter :: DUMP_TARGET_vbcp = 14401
      logical, save :: dump_done_vbcp = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpbbc,bbc)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

! -----

!! Set the bottom and top boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: vbcp.f90 :: s_vbcp
! Summary : Sets vertical boundary conditions for pressure perturbation at
!           bottom and top boundaries with extrapolation or copying.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global/module variable writes, only local array writes
!   - No synchronization constructs (barrier, critical, atomic)
!   - Conditional branching on bbc value (executed by all threads)
!   - Multiple separate omp do regions within single parallel region
! Next:
!   - Direct conversion to OpenACC parallel loop or OpenACC
!   - Consider using OpenACC kernels directive for multiple loops
!   - Conditionals can remain as they are data-independent
! Runtime:
!   - Calls: 14401
!   - AvgLoops: 806.4K
!   - TotalTime: 0.658s (0.02%)
!   - AvgTime: 0.046ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('vbcp.f90', 's_vbcp', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_vbcp = dump_call_count_vbcp + 1
if (dump_call_count_vbcp == DUMP_TARGET_vbcp .and. .not. dump_done_vbcp) then
  call dump_init('vbcp')
  call dump_scalar_i('bbc', bbc)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('ppf_in.bin', ppf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_i('nkm1', nkm1)
  call dump_scalar_i('nkm2', nkm2)
end if

!$omp parallel default(shared)

! Set the bottom boundary conditions.

      if(bbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ppf(i,j,1)=2.e0*ppf(i,j,2)-ppf(i,j,3)
        end do
        end do

!$omp end do

      else if(bbc.eq.3) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ppf(i,j,1)=ppf(i,j,2)
        end do
        end do

!$omp end do

      end if

! -----

! Set the top boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        ppf(i,j,nkm1)=ppf(i,j,nkm2)
      end do
      end do

!$omp end do

! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_vbcp == DUMP_TARGET_vbcp .and. .not. dump_done_vbcp) then
  call dump_array_3d('ppf_ref.bin', ppf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_vbcp = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_vbcp

!-----7--------------------------------------------------------------7--

      end module m_vbcp

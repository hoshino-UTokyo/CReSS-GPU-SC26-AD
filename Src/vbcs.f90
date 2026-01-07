!***********************************************************************
      module m_vbcs
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/20,
!                   1999/08/18, 1999/09/06, 1999/10/12, 1999/11/01,
!                   1999/12/06, 2000/01/17, 2000/03/17, 2001/09/13,
!                   2001/11/20, 2001/12/11, 2002/04/02, 2002/07/23,
!                   2002/08/15, 2003/04/30, 2003/05/19, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the vertical boundary condition for optional scalar variable.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing
      use m_comprofile
      use m_dump_kernel

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vbcs, s_vbcs

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vbcs

        module procedure s_vbcs

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
      subroutine s_vbcs(ni,nj,nk,sf)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: sf(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at future

! Internal shared variables

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_vbcs = 0
      integer, parameter :: DUMP_TARGET_vbcs = 3962
      logical, save :: dump_done_vbcs = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

! -----

!! Set the bottom and top boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: vbcs.f90 :: s_vbcs
! Summary : Sets vertical boundary conditions for scalar variable by copying
!           values from adjacent levels at bottom and top boundaries.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global/module variable writes, only local array writes
!   - No synchronization constructs (barrier, critical, atomic)
!   - Two separate omp do regions for bottom and top boundaries
! Next:
!   - Direct conversion to OpenACC parallel loop or OpenACC
!   - Both loops are independent and can run concurrently on GPU
! Runtime:
!   - Calls: 3962
!   - AvgLoops: 806.4K
!   - TotalTime: 0.164s (0.01%)
!   - AvgTime: 0.042ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('vbcs.f90', 's_vbcs', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_vbcs = dump_call_count_vbcs + 1
if (dump_call_count_vbcs == DUMP_TARGET_vbcs .and. .not. dump_done_vbcs) then
  call dump_init('vbcs')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('sf_in.bin', sf, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared)

! Set the bottom boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        sf(i,j,1)=sf(i,j,2)
      end do
      end do

!$omp end do

! -----

! Set the top boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        sf(i,j,nkm1)=sf(i,j,nkm2)
      end do
      end do

!$omp end do

! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_vbcs == DUMP_TARGET_vbcs .and. .not. dump_done_vbcs) then
  call dump_array_3d('sf_ref.bin', sf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_vbcs = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_vbcs

!-----7--------------------------------------------------------------7--

      end module m_vbcs

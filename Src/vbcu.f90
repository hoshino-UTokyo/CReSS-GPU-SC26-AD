!***********************************************************************
      module m_vbcu
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/05/20, 1999/08/18, 1999/09/06, 1999/10/12,
!                   1999/11/01, 1999/12/06, 2000/01/17, 2000/03/17,
!                   2000/03/23, 2001/01/15, 2001/04/15, 2001/07/13,
!                   2001/08/07, 2001/09/13, 2001/11/20, 2001/12/11,
!                   2002/04/02, 2002/07/23, 2002/08/15, 2003/04/30,
!                   2003/05/19, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the vertical boundary conditions for the x components of
!     velocity.

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

      public :: vbcu, s_vbcu

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vbcu

        module procedure s_vbcu

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
      subroutine s_vbcu(ni,nj,nk,uf)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

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
      integer, save :: dump_call_count_vbcu = 0
      integer, parameter :: DUMP_TARGET_vbcu = 14401
      logical, save :: dump_done_vbcu = .false.

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

! -----

!! Set the bottom and top boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: vbcu.f90 :: s_vbcu
! Summary : Sets vertical boundary conditions for x-velocity component by
!           copying values from adjacent levels at bottom and top boundaries.
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
!   - Calls: 14401
!   - AvgLoops: 807.3K
!   - TotalTime: 0.590s (0.02%)
!   - AvgTime: 0.041ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('vbcu.f90', 's_vbcu', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni)-(1)+1,8)

! Dump input data at target call
dump_call_count_vbcu = dump_call_count_vbcu + 1
if (dump_call_count_vbcu == DUMP_TARGET_vbcu .and. .not. dump_done_vbcu) then
  call dump_init('vbcu')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('uf_in.bin', uf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_i('nkm1', nkm1)
  call dump_scalar_i('nkm2', nkm2)
end if

call profile_start(prof_id1)

#if defined(USE_GPU) && !defined(DISABLE_GPU_361)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    ! Bottom boundary
    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni
        uf(i,j,1) = uf(i,j,2)
      end do
    end do
    !$acc end kernels

    ! Top boundary
    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni
        uf(i,j,nkm1) = uf(i,j,nkm2)
      end do
    end do
    !$acc end kernels

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared)

! Set the bottom boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni
        uf(i,j,1)=uf(i,j,2)
      end do
      end do

!$omp end do

! -----

! Set the top boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni
        uf(i,j,nkm1)=uf(i,j,nkm2)
      end do
      end do

!$omp end do

! -----

!$omp end parallel

#endif

call profile_stop(prof_id1, loop_len)

! Dump output data at target call
if (dump_call_count_vbcu == DUMP_TARGET_vbcu .and. .not. dump_done_vbcu) then
  call dump_array_3d('uf_ref.bin', uf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_vbcu = .true.
end if

!! -----

      end subroutine s_vbcu

!-----7--------------------------------------------------------------7--

      end module m_vbcu

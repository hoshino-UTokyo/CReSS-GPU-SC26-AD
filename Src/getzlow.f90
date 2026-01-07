!***********************************************************************
      module m_getzlow
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/04/01
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the z physical coordinates at lowest plane.

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

      public :: getzlow, s_getzlow

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getzlow

        module procedure s_getzlow

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
      subroutine s_getzlow(ni,nj,nk,zph,za)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Output variable

      real, intent(out) :: za(0:ni+1,0:nj+1)
                       ! z physical coordinates at lowest plane

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_getzlow = 0
      integer, parameter :: DUMP_TARGET_getzlow = 361
      logical, save :: dump_done_getzlow = .false.


!-----7--------------------------------------------------------------7--

! Calculate the z physical coordinates at lowest plane.

!@llm start meta_info ----------------------------------------------------
! Location: getzlow.f90 :: s_getzlow
! Summary : Calculate height of lowest model level above terrain by
!           averaging vertical spacing between levels 2 and 3.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Simple 2D loop with element-wise arithmetic
!   - Reads from 3D zph array at fixed k indices (2,3)
!   - No global writes, only output array za is modified
!   - No synchronization constructs
! Next:
!   - Direct translation to OpenACC with teams distribute
!   - Consider loop collapse for j,i dimensions
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.012s (0.00%)
!   - AvgTime: 0.032ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('getzlow.f90', 's_getzlow', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_getzlow = dump_call_count_getzlow + 1
if (dump_call_count_getzlow == DUMP_TARGET_getzlow .and. .not. dump_done_getzlow) then
  call dump_init('getzlow')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('zph.bin', zph, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        za(i,j)=.5e0*(zph(i,j,3)-zph(i,j,2))
      end do
      end do

!$omp end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_getzlow == DUMP_TARGET_getzlow .and. .not. dump_done_getzlow) then
  call dump_array_2d('za_ref.bin', za, 0, ni+1, 0, nj+1)
  call dump_finalize()
  dump_done_getzlow = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_getzlow

!-----7--------------------------------------------------------------7--

      end module m_getzlow

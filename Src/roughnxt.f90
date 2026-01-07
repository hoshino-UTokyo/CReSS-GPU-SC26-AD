!***********************************************************************
      module m_roughnxt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/15
!     Modification: 2001/10/18, 2002/04/02, 2003/02/05, 2003/04/30,
!                   2003/05/19, 2003/07/15, 2003/12/12, 2004/05/07,
!                   2004/08/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/06/16, 2009/11/13, 2011/06/01,
!                   2011/11/10, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     reset the roughness parameter on the sea surface to the next time
!     step.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: roughnxt, s_roughnxt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface roughnxt

        module procedure s_roughnxt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_roughnxt(ni,nj,land,va,cm,z0m,z0h)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: va(0:ni+1,0:nj+1)
                       ! Magnitude of velocity at lowest plane

      real, intent(in) :: cm(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity

! Input and output variables

      real, intent(inout) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(inout) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real ust         ! Friction velocity


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_roughnxt = 0
      integer, parameter :: DUMP_TARGET_roughnxt = 361
      logical, save :: dump_done_roughnxt = .false.


!-----7--------------------------------------------------------------7--

! Calculate the roughness parameter on the sea surface to the next time
! step.

!@llm start meta_info ----------------------------------------------------
! Location: roughnxt.f90 :: s_roughnxt
! Summary : Update sea surface roughness length for next time step based on friction velocity
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses intrinsic max function
!   - Simple 2D loop with conditional on land use (land < 3)
!   - Conditional on ust threshold for different roughness formulas
!   - Writes to z0m and z0h arrays
!   - No synchronization constructs
! Next:
!   - Convert to OpenACC with teams distribute parallel for
!   - Straightforward GPU port with collapse(2) clause
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.012s (0.00%)
!   - AvgTime: 0.033ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('roughnxt.f90', 's_roughnxt', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_roughnxt = dump_call_count_roughnxt + 1
if (dump_call_count_roughnxt == DUMP_TARGET_roughnxt .and. .not. dump_done_roughnxt) then
  call dump_init('roughnxt')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_array_2d_int('land.bin', land, 0, ni+1, 0, nj+1)
  call dump_array_2d('va.bin', va, 0, ni+1, 0, nj+1)
  call dump_array_2d('cm.bin', cm, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0m_in.bin', z0m, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0h_in.bin', z0h, 0, ni+1, 0, nj+1)
end if

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,ust)

      do j=1,nj-1
      do i=1,ni-1

        if(land(i,j).lt.3) then

          ust=cm(i,j)*va(i,j)

          if(ust.lt.1.08e0) then

            z0m(i,j)=max(-34.7e-6+8.28e-4*ust,z0min)

          else

            z0m(i,j)=max(-.277e-2+3.39e-3*ust,z0min)

          end if

          z0h(i,j)=z0m(i,j)

        end if

      end do
      end do

!$omp end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_roughnxt == DUMP_TARGET_roughnxt .and. .not. dump_done_roughnxt) then
  call dump_array_2d('z0m_ref.bin', z0m, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0h_ref.bin', z0h, 0, ni+1, 0, nj+1)
  call dump_finalize()
  dump_done_roughnxt = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_roughnxt

!-----7--------------------------------------------------------------7--

      end module m_roughnxt

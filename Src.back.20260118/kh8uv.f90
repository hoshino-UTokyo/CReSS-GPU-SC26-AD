!***********************************************************************
      module m_kh8uv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/06/10
!     Modification: 2006/02/13, 2006/11/06, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/09, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the horizontal eddy diffusivity at the u and v points.

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

      public :: kh8uv, s_kh8uv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface kh8uv

        module procedure s_kh8uv

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
      subroutine s_kh8uv(fpmpopt,fpmfcopt,ni,nj,nk,rmf,rkh,rkh8u,rkh8v)
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

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

! Input and output variable

      real, intent(inout) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy diffusivity / Jacobian

! Output variables

      real, intent(out) :: rkh8u(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity
                       ! / Jacobian at u points

      real, intent(out) :: rkh8v(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity
                       ! / Jacobian at v points

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

! Remark

!     rkh: This variable is also temporary, because it is not used
!          again.


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_kh8uv = 0
      integer, parameter :: DUMP_TARGET_kh8uv = 360
      logical, save :: dump_done_kh8uv = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)

! -----

! Set the horizontal eddy diffusivity at the u and v points.

!@llm start meta_info ----------------------------------------------------
! Location: kh8uv.f90 :: s_kh8uv
! Summary : Set horizontal eddy diffusivity at u and v points (rkh8u, rkh8v)
!           by averaging rkh values, with map scale factor corrections
! GPU diff: Medium
! Findings:
!   - Multiple conditional branches (mfcopt, mpopt)
!   - Simple arithmetic operations (addition, multiplication)
!   - Private variables: k, i, j
!   - Reads from rkh, rmf arrays
!   - Writes to rkh, rkh8u, rkh8v arrays
!   - rkh is modified in-place then used (potential ordering concern)
!   - Multiple omp do regions within single parallel block
!   - No sync constructs between threads
! Next:
!   - Can be ported to GPU with OpenACC parallel loop
!   - Careful attention needed for rkh modification ordering
!   - Consider separating different mpopt/mfcopt cases into different kernels
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.3M
!   - TotalTime: 3.068s (0.10%)
!   - AvgTime: 8.522ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('kh8uv.f90', 's_kh8uv', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(2)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_kh8uv = dump_call_count_kh8uv + 1
if (dump_call_count_kh8uv == DUMP_TARGET_kh8uv .and. .not. dump_done_kh8uv) then
  call dump_init('kh8uv')
  call dump_scalar_i('mpopt', mpopt)
  call dump_scalar_i('mfcopt', mfcopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call dump_array_3d('rkh_in.bin', rkh, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

      if(mfcopt.eq.1.and.(mpopt.eq.0.or.mpopt.eq.10)) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-1
            rkh8u(i,j,k)=rkh(i-1,j,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            rkh(i,j,k)=rmf(i,j,2)*rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=1,ni-1
            rkh8v(i,j,k)=rkh(i,j-1,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

      else if(mfcopt.eq.1.and.mpopt.eq.5) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=1,ni-1
            rkh8v(i,j,k)=rkh(i,j-1,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            rkh(i,j,k)=rmf(i,j,2)*rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-1
            rkh8u(i,j,k)=rkh(i-1,j,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-1
            rkh8u(i,j,k)=rkh(i-1,j,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=1,ni-1
            rkh8v(i,j,k)=rkh(i,j-1,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! Dump output data at target call
if (dump_call_count_kh8uv == DUMP_TARGET_kh8uv .and. .not. dump_done_kh8uv) then
  call dump_array_3d('rkh8u_ref.bin', rkh8u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rkh8v_ref.bin', rkh8v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rkh_ref.bin', rkh, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_kh8uv = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_kh8uv

!-----7--------------------------------------------------------------7--

      end module m_kh8uv

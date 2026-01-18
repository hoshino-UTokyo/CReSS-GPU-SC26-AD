!***********************************************************************
      module m_buoytke
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/10/12
!     Modification: 1999/11/01, 1999/11/24, 1999/11/30, 2000/01/17,
!                   2000/03/08, 2000/04/18, 2000/07/05, 2000/10/18,
!                   2000/11/18, 2000/12/19, 2001/04/15, 2001/05/29,
!                   2001/06/29, 2001/12/11, 2002/01/21, 2002/04/02,
!                   2002/06/06, 2002/12/02, 2003/01/04, 2003/02/13,
!                   2003/03/13, 2003/04/30, 2003/05/19, 2003/11/28,
!                   2003/12/12, 2004/02/01, 2004/03/05, 2004/04/15,
!                   2004/06/10, 2004/07/01, 2004/09/01, 2004/09/10,
!                   2005/04/04, 2005/06/10, 2006/01/10, 2006/02/13,
!                   2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/11/10, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the buoyancy production in the turbulent kinetic energy
!     equation.

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

      public :: buoytke, s_buoytke

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface buoytke

        module procedure s_buoytke

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
      subroutine s_buoytke(ni,nj,nk,jcb8w,nsq8w,rkv8s,tkefrc,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: nsq8w(0:ni+1,0:nj+1,1:nk)
                       ! Half value of Brunt Vaisala frequency squared
                       ! at w points

      real, intent(in) :: rkv8s(0:ni+1,0:nj+1,1:nk)
                       ! Half value of rbr x vertical eddy diffusivity
                       ! at scalar points

! Input and output variable

      real, intent(inout) :: tkefrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in turbulent kinetic energy equation

! Internal shared variable

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_buoytke = 0
      integer, parameter :: DUMP_TARGET_buoytke = 360
      logical, save :: dump_done_buoytke = .false.


!-----7--------------------------------------------------------------7--

! Calculate the buoyancy production.

!@llm start meta_info ----------------------------------------------------
! Location: buoytke.f90 :: s_buoytke
! Summary : Calculates buoyancy production term for turbulent kinetic
!           energy equation using Brunt-Vaisala frequency and diffusivity.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple arithmetic operations only
!   - Two sequential loop nests: first computes tmp1, second updates tkefrc
!   - Data dependency: tmp1 must be computed before tkefrc update
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Use collapse(2) for nested i,j loops within each k-loop
!   - Keep two separate target regions or use explicit barrier between phases
!   - Consider fusing loops if tmp1 dependency can be restructured
! Runtime:
!   - Calls: 360
!   - AvgLoops: 101.2M
!   - TotalTime: 3.027s (0.10%)
!   - AvgTime: 8.409ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('buoytke.f90', 's_buoytke', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(2)+1,8) &
     & * int((nj-2)-(2)+1,8) &
     & * int((ni-2)-(2)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_buoytke = dump_call_count_buoytke + 1
if (dump_call_count_buoytke == DUMP_TARGET_buoytke .and. .not. dump_done_buoytke) then
  call dump_init('buoytke')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('jcb8w.bin', jcb8w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nsq8w.bin', nsq8w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rkv8s.bin', rkv8s, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tkefrc_in.bin', tkefrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_in.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            tmp1(i,j,k)=jcb8w(i,j,k)*nsq8w(i,j,k)                       &
     &        *(rkv8s(i,j,k-1)+rkv8s(i,j,k))
          end do
          end do

!$omp end do

        end do

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            tkefrc(i,j,k)=tkefrc(i,j,k)-(tmp1(i,j,k)+tmp1(i,j,k+1))
          end do
          end do

!$omp end do

        end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_buoytke == DUMP_TARGET_buoytke .and. .not. dump_done_buoytke) then
  call dump_array_3d('tkefrc_ref.bin', tkefrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_ref.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_buoytke = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_buoytke

!-----7--------------------------------------------------------------7--

      end module m_buoytke

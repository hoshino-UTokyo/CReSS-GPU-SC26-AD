!***********************************************************************
      module m_getexner
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/03/08
!     Modification: 2000/04/18, 2000/07/05, 2001/05/29, 2002/04/02,
!                   2003/04/30, 2003/05/19, 2004/04/15, 2007/10/19,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the total pressure variable and Exner function.

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

      public :: getexner, s_getexner

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getexner

        module procedure s_getexner

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getexner(ni,nj,nk,pbr,pp,pi,p)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

! Output variables

      real, intent(out) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(out) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exner function

! Internal shared variables

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_getexner = 0
      integer, parameter :: DUMP_TARGET_getexner = 1080
      logical, save :: dump_done_getexner = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Calculate the total pressure variable and Exner function.

!@llm start meta_info ----------------------------------------------------
! Location: getexner.f90 :: s_getexner
! Summary : Calculates total pressure (p = pbr + pp) and Exner function
!           (pi = (p/p0)^(rd/cp)) at each grid point
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - Uses intrinsic exp() and log() functions - GPU compatible
!   - Simple element-wise computation, no dependencies between iterations
!   - Module constants rd, cp, p0 used from m_comphy
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Direct port to OpenACC parallel loop
!   - Ensure module constants are accessible on device
!   - Consider collapsing k,j,i loops for better GPU occupancy
! Runtime:
!   - Calls: 1080
!   - AvgLoops: 102.4M
!   - TotalTime: 6.873s (0.23%)
!   - AvgTime: 6.364ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('getexner.f90', 's_getexner', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_getexner = dump_call_count_getexner + 1
if (dump_call_count_getexner == DUMP_TARGET_getexner .and. .not. dump_done_getexner) then
  call dump_init('getexner')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('p.bin', p, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pbr.bin', pbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pp.bin', pp, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          p(i,j,k)=pbr(i,j,k)+pp(i,j,k)

          pi(i,j,k)=exp(rddvcp*log(p0iv*p(i,j,k)))

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_getexner == DUMP_TARGET_getexner .and. .not. dump_done_getexner) then
  call dump_array_3d('pi_ref.bin', pi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_getexner = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_getexner

!-----7--------------------------------------------------------------7--

      end module m_getexner

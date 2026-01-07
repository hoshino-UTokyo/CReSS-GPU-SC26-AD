!***********************************************************************
      module m_copy1d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/08/23, 1999/10/12,
!                   2000/01/17, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2006/01/10, 2007/01/31, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     copy 1 dimensional invar to outvar.

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

      public :: copy1d, s_copy1d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface copy1d

        module procedure s_copy1d

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
      subroutine s_copy1d(kmin,kmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: invar(kmin:kmax)
                       ! Coping variable

! Output variable

      real, intent(out) :: outvar(kmin:kmax)
                       ! Copied variable

! Internal private variable

      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_copy1d = 0
      integer, parameter :: DUMP_TARGET_copy1d = 1
      logical, save :: dump_done_copy1d = .false.


!-----7--------------------------------------------------------------7--

! Copy the invar to the outvar.

!@llm start meta_info ----------------------------------------------------
! Location: copy1d.f90 :: s_copy1d
! Summary : Simple 1D array copy from invar to outvar.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls
!   - Trivial memory copy operation
!   - Independent element-wise operations
! Next:
!   - Straightforward GPU port; consider using cudaMemcpy or similar
!   - For small arrays, overhead may exceed benefit of GPU execution
!   - May be better to keep data resident on GPU and avoid copy
! Runtime:
!   - Calls: 1
!   - AvgLoops: 128
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.012ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('copy1d.f90', 's_copy1d', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(kmin)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_copy1d = dump_call_count_copy1d + 1
if (dump_call_count_copy1d == DUMP_TARGET_copy1d .and. .not. dump_done_copy1d) then
  call dump_init('copy1d')
  call dump_scalar_i('kmin', kmin)
  call dump_scalar_i('kmax', kmax)
end if

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=kmin,kmax
        outvar(k)=invar(k)
      end do

!$omp end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_copy1d == DUMP_TARGET_copy1d .and. .not. dump_done_copy1d) then
  call dump_finalize()
  dump_done_copy1d = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_copy1d

!-----7--------------------------------------------------------------7--

      end module m_copy1d

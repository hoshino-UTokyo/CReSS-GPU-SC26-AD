!***********************************************************************
      module m_setcst2d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/08/23, 1999/09/30,
!                   1999/10/12, 2000/01/17, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2006/01/10, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2010/12/13, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     fill in the 2 dimensional array with the constant value.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setcst2d, s_setcst2d, s_setcst2d_r8

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setcst2d

        module procedure s_setcst2d, s_setcst2d_r8

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
      subroutine s_setcst2d(imin,imax,jmin,jmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      real, intent(in) :: invar
                       ! Constant variable

! Output variable

      real, intent(out) :: outvar(imin:imax,jmin:jmax)
                       ! Copied variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_setcst2d = 0
      integer, parameter :: DUMP_TARGET_setcst2d = 41
      logical, save :: dump_done_setcst2d = .false.


!-----7--------------------------------------------------------------7--

! Fill in the array with the constant value.

!@llm start meta_info ----------------------------------------------------
! Location: setcst2d.f90 :: s_setcst2d
! Summary : Fill 2D array with a constant value (real type)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls
!   - Simple 2D loop with direct assignment
!   - No synchronization constructs
!   - Trivially parallelizable
! Next:
!   - Convert to OpenACC with parallel loop collapse(2)
!   - Consider using memset or array assignment for better performance
! Runtime:
!   - Calls: 41
!   - AvgLoops: 811.8K
!   - TotalTime: 0.002s (0.00%)
!   - AvgTime: 0.054ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('setcst2d.f90', 's_setcst2d', &
   & 'OMP section 1')
end if
loop_len = int((jmax)-(jmin)+1,8) * int((imax)-(imin)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_setcst2d = dump_call_count_setcst2d + 1
if (dump_call_count_setcst2d == DUMP_TARGET_setcst2d .and. .not. dump_done_setcst2d) then
  call dump_init('setcst2d')
  call dump_scalar_i('imin', imin)
  call dump_scalar_i('imax', imax)
  call dump_scalar_i('jmin', jmin)
  call dump_scalar_i('jmax', jmax)
  ! FIXME: invar is array - call dump_scalar_r('invar', invar)
end if

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

      do j=jmin,jmax
      do i=imin,imax
        outvar(i,j)=invar
      end do
      end do

!$omp end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_setcst2d == DUMP_TARGET_setcst2d .and. .not. dump_done_setcst2d) then
  call dump_array_2d('outvar_ref.bin', outvar, imin, imax, jmin, jmax)
  call dump_finalize()
  dump_done_setcst2d = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_setcst2d

!***********************************************************************
      subroutine s_setcst2d_r8(imin,imax,jmin,jmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      real, intent(in) :: invar
                       ! Constant variable

! Output variable

      real(kind=r8), intent(out) :: outvar(imin:imax,jmin:jmax)
                       ! Copied variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Fill in the array with the constant value.

!@llm start meta_info ----------------------------------------------------
! Location: setcst2d.f90 :: s_setcst2d_r8
! Summary : Fill 2D array with a constant value (real*8 type)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls
!   - Simple 2D loop with direct assignment
!   - No synchronization constructs
!   - Trivially parallelizable
! Next:
!   - Convert to OpenACC with parallel loop collapse(2)
!   - Consider using memset or array assignment for better performance
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('setcst2d.f90', 's_setcst2d_r8', &
   & 'OMP section 1')
end if
loop_len = int((jmax)-(jmin)+1,8) * int((imax)-(imin)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

      do j=jmin,jmax
      do i=imin,imax
        outvar(i,j)=invar
      end do
      end do

!$omp end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_setcst2d_r8

!-----7--------------------------------------------------------------7--

      end module m_setcst2d

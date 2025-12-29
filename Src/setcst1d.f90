!***********************************************************************
      module m_setcst1d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/08/23, 1999/09/30,
!                   1999/10/12, 2000/01/17, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2006/01/10, 2007/01/31, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2010/12/13,
!                   2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     fill in the 1 dimensional array with the constant value.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_comprofile

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setcst1d, s_setcst1d, s_setcst1d_r8

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setcst1d

        module procedure s_setcst1d, s_setcst1d_r8

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
      subroutine s_setcst1d(kmin,kmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: invar
                       ! Constant variable

! Output variable

      real, intent(out) :: outvar(kmin:kmax)
                       ! Copied variable

! Internal private variable

      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Fill in the array with the constant value.

!@llm start meta_info ----------------------------------------------------
! Location: setcst1d.f90 :: s_setcst1d
! Summary : Fill 1D array with a constant value (real type)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls
!   - Simple 1D loop with direct assignment
!   - No synchronization constructs
!   - Trivially parallelizable
! Next:
!   - Convert to OpenACC with parallel loop
!   - Consider using memset or array assignment for better performance
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('setcst1d.f90', 's_setcst1d', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(kmin)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=kmin,kmax
        outvar(k)=invar
      end do

!$omp end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_setcst1d

!***********************************************************************
      subroutine s_setcst1d_r8(kmin,kmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: invar
                       ! Constant variable

! Output variable

      real(kind=r8), intent(out) :: outvar(kmin:kmax)
                       ! Copied variable

! Internal private variable

      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Fill in the array with the constant value.

!@llm start meta_info ----------------------------------------------------
! Location: setcst1d.f90 :: s_setcst1d_r8
! Summary : Fill 1D array with a constant value (real*8 type)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls
!   - Simple 1D loop with direct assignment
!   - No synchronization constructs
!   - Trivially parallelizable
! Next:
!   - Convert to OpenACC with parallel loop
!   - Consider using memset or array assignment for better performance
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('setcst1d.f90', 's_setcst1d_r8', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(kmin)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=kmin,kmax
        outvar(k)=invar
      end do

!$omp end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_setcst1d_r8

!-----7--------------------------------------------------------------7--

      end module m_setcst1d

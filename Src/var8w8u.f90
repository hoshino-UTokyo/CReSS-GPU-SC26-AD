!***********************************************************************
      module m_var8w8u
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/05/07
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     optional variable at w points be averaged to u points.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing
      use m_comprofile

      implicit none

! Default access control

      private

! Exceptional access control

      public :: var8w8u, s_var8w8u

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface var8w8u

        module procedure s_var8w8u

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
      subroutine s_var8w8u(ni,nj,nk,var8w,var8u)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: var8w(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at w points

! Output variable

      real, intent(out) :: var8u(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at u points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Be averaged to u points.

!@llm start meta_info ----------------------------------------------------
! Location: var8w8u.f90 :: s_var8w8u
! Summary : Averages variable at w points to u points using 4-point averaging
!           in x and z directions.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - No global/module variable writes, only local array writes
!   - No synchronization constructs (barrier, critical, atomic)
!   - Simple loop structure with private loop indices
! Next:
!   - Direct conversion to OpenACC parallel loop or OpenACC
!   - Consider collapsing nested loops for better GPU utilization
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('var8w8u.f90', 's_var8w8u', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) * int((nj)-(0)+1,8) * int((ni)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=1,ni
          var8u(i,j,k)=.25e0*((var8w(i-1,j,k)+var8w(i-1,j,k+1))         &
     &      +(var8w(i,j,k)+var8w(i,j,k+1)))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_var8w8u

!-----7--------------------------------------------------------------7--

      end module m_var8w8u

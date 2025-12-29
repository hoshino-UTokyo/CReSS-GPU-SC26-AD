!***********************************************************************
      module m_var8u8s
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/08/18
!     Modification: 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     optional variable at u points be averaged to scalar points.

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

      public :: var8u8s, s_var8u8s

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface var8u8s

        module procedure s_var8u8s

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
      subroutine s_var8u8s(ni,nj,nk,var8u,var8s)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: var8u(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at u points

! Output variable

      real, intent(out) :: var8s(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at scalar points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Be averaged to scalar points.

!@llm start meta_info ----------------------------------------------------
! Location: var8u8s.f90 :: s_var8u8s
! Summary : Average optional variable from u-points to scalar points
!           using simple 2-point stencil in x-direction
! GPU diff: Easy
! Findings:
!   - Serial k-loop wrapping parallel i,j loops (private(k))
!   - Simple 2-point averaging: (var8u(i,j,k)+var8u(i+1,j,k))*0.5
!   - No function calls or complex operations
!   - Output array is independent of input (no race condition)
! Next:
!   - Convert to OpenACC with collapse clause
!   - Straightforward GPU port with good memory access pattern
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('var8u8s.f90', 's_var8u8s', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) * int((nj)-(0)+1,8) * int((ni)-(0)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          var8s(i,j,k)=.5e0*(var8u(i,j,k)+var8u(i+1,j,k))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_var8u8s

!-----7--------------------------------------------------------------7--

      end module m_var8u8s

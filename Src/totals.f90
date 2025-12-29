!***********************************************************************
      module m_totals
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/08/18
!     Modification: 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the total value of optional scalar variable from base state
!     and perturbation value.

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

      public :: totals, s_totals

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface totals

        module procedure s_totals

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
      subroutine s_totals(ni,nj,nk,sbr,sp,s)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: sbr(0:ni+1,0:nj+1,1:nk)
                       ! Optional base state scalar variable

      real, intent(in) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar perturbation variable

! Output variable

      real, intent(out) :: s(0:ni+1,0:nj+1,1:nk)
                       ! Total value of optional scalar variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the total value of optional scalar variable from base state and
! perturbation value.

!@llm start meta_info ----------------------------------------------------
! Location: totals.f90 :: s_totals
! Summary : Add base state and perturbation values to get total scalar
!           variable (s = sbr + sp)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls
!   - Single output array (s)
!   - Simple element-wise addition
!   - No synchronization constructs
! Next:
!   - Very simple GPU port - ideal candidate
!   - Consider fusing with other scalar operations
!   - Memory bandwidth bound operation
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('totals.f90', 's_totals', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) * int((nj)-(0)+1,8) * int((ni)-(0)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          s(i,j,k)=sbr(i,j,k)+sp(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_totals

!-----7--------------------------------------------------------------7--

      end module m_totals

!***********************************************************************
      module m_adjstqa
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the aerosol mixing ratio more than user specified value.

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

      public :: adjstqa, s_adjstqa

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstqa

        module procedure s_adjstqa

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
      subroutine s_adjstqa(ni,nj,nk,nqa,qasl)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

! Input and output variable

      real, intent(inout) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Force the aerosol mixing ratio more than user specified value.

!@llm start meta_info ----------------------------------------------------
! Location: adjstqa.f90 :: subroutine s_adjstqa
! Summary : Forces aerosol mixing ratio to be non-negative using max(0).
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Pure max() operation, GPU compatible.
!   - 4D array with aerosol categories in outer loop.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse for (n,k,j,i).
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('adjstqa.f90', 's_adjstqa', &
   & 'OMP section 1')
end if
loop_len = int((nqa(0))-(1)+1,8) &
     & * int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k,n)

      do n=1,nqa(0)

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            qasl(i,j,k,n)=max(qasl(i,j,k,n),0.e0)
          end do
          end do

!$omp end do

        end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_adjstqa

!-----7--------------------------------------------------------------7--

      end module m_adjstqa

!***********************************************************************
      module m_diagnci
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2001/10/18, 2001/11/20, 2002/01/07,
!                   2002/01/15, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2003/12/12, 2004/04/01, 2004/05/31,
!                   2004/06/10, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2005/04/04, 2005/09/30, 2005/10/05,
!                   2006/01/10, 2006/02/13, 2006/04/03, 2007/10/19,
!                   2007/11/26, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the diagnostic concentrations of the cloud ice.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diagnci, s_diagnci

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diagnci

        module procedure s_diagnci

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
      subroutine s_diagnci(ni,nj,nk,nqi,nni,qice,nidia)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

! Output variable

      real, intent(out) :: nidia(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Diagnostic concentrations of ice hydrometeor

! Internal shared variable

      real miiv        ! 4.0 / (3.0 x mimax)

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      miiv=4.e0/(3.e0*mimax)

! -----

! Get the diagnostic concentrations of the cloud ice.

!@llm start meta_info ----------------------------------------------------
! Location: diagnci.f90 :: s_diagnci
! Summary : Calculate diagnostic cloud ice concentrations from ice mixing
!           ratio using a simple linear scaling with inverse max mass.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Private variable k for outer loop
!   - Writes to nidia output array
!   - Simple point-wise multiplication operation
!   - Completely independent iterations
! Next:
!   - Direct conversion to OpenACC with collapsed loops
!   - Minimal data transfer: input qice, output nidia
!   - Excellent GPU candidate due to simplicity
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('diagnci.f90', 's_diagnci', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          nidia(i,j,k,1)=miiv*qice(i,j,k,1)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_diagnci

!-----7--------------------------------------------------------------7--

      end module m_diagnci

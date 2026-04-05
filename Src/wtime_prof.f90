!***********************************************************************
      module m_wtime_prof
!***********************************************************************

! In this module,
!     declare the wall-clock profiling accumulators.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_defmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      public

!-----7--------------------------------------------------------------7--

! Module variables

      double precision wt_comm
                       ! Accumulated wall time for MPI communication

!-----7--------------------------------------------------------------7--

      contains

!-----7--------------------------------------------------------------7--

      subroutine wtime_prof_init()

! Initialize profiling accumulators.

      wt_comm = 0.0d0

      end subroutine wtime_prof_init

!-----7--------------------------------------------------------------7--

      subroutine wtime_prof_report(mype,root)

! Report profiling results from rank 0.

      integer, intent(in) :: mype
      integer, intent(in) :: root

      if(mype.eq.root) then
        write(6,'(a)')       ' ================================================'
        write(6,'(a)')       '  Wall-clock profiling summary'
        write(6,'(a)')       ' ================================================'
        write(6,'(a,f12.3,a)') '  MPI communication : ',wt_comm,' sec'
        write(6,'(a)')       ' ================================================'
      end if

      end subroutine wtime_prof_report

!-----7--------------------------------------------------------------7--

      end module m_wtime_prof

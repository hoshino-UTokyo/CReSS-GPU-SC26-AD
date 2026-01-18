!***********************************************************************
      module m_chkitr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/06/01
!     Modification: 2000/07/05, 2001/04/15, 2001/05/29, 2001/08/07,
!                   2002/04/02, 2003/04/30, 2003/05/19, 2003/11/05,
!                   2004/08/20, 2004/09/10, 2006/12/04, 2007/01/20,
!                   2007/05/14, 2007/10/19, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2009/02/27, 2009/11/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the convergence of the iteration.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comprofile
      use m_dump_kernel
      use m_defmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkitr, s_chkitr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkitr

        module procedure s_chkitr

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chkitr(fproc,istr,iend,jstr,jend,kstr,kend,itcon,    &
     &                    ni,nj,nk,dvar)
!***********************************************************************

! Input variables

      character(len=6), intent(in) :: fproc
                       ! Control flag of parallel processing

      integer, intent(in) :: istr
                       ! Minimum do loops index in x direction

      integer, intent(in) :: iend
                       ! Maximum do loops index in x direction

      integer, intent(in) :: jstr
                       ! Minimum do loops index in y direction

      integer, intent(in) :: jend
                       ! Maximum do loops index in y direction

      integer, intent(in) :: kstr
                       ! Minimum do loops index in z direction

      integer, intent(in) :: kend
                       ! Maximum do loops index in z direction

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dvar(0:ni+1,0:nj+1,1:nk)
                       ! Variations of iteration

! Output variable

      real, intent(out) :: itcon
                       ! Control flag of continuation of interation

! Internal shared variables

      integer ierr     ! Error descriptor

      real intitc      ! Internal processed control flag of
                       ! continuation of interation

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_chkitr = 0
      integer, parameter :: DUMP_TARGET_chkitr = 16
      logical, save :: dump_done_chkitr = .false.


!-----7--------------------------------------------------------------7--

! Initialize the processed variable, intitc.

      intitc=0.e0

! -----

! Check the convergence of the iteration in each processor element.

!@llm start meta_info ----------------------------------------------------
! Location: chkitr.f90 :: s_chkitr
! Summary : Find maximum absolute value of iteration variations (dvar)
!           for convergence checking with max reduction
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function/subroutine calls inside parallel region
!   - Uses reduction(max:) on single variable intitc
!   - Uses intrinsic functions (abs, max)
!   - Simple 3D loop with element-wise max computation
!   - MPI_allreduce called after parallel region (not inside)
! Next:
!   - Direct OpenACC with collapse(3) and reduction(max:)
!   - GPU reduction primitives well-suited for this pattern
! Runtime:
!   - Calls: 16
!   - AvgLoops: 806.4K
!   - TotalTime: 0.002s (0.00%)
!   - AvgTime: 0.133ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('chkitr.f90', 's_chkitr', &
   & 'OMP section 1')
end if
loop_len = int((kend)-(kstr)+1,8) &
     & * int((jend)-(jstr)+1,8) &
     & * int((iend)-(istr)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_chkitr = dump_call_count_chkitr + 1
if (dump_call_count_chkitr == DUMP_TARGET_chkitr .and. .not. dump_done_chkitr) then
  call dump_init('chkitr')
  call dump_scalar_i('istr', istr)
  call dump_scalar_i('iend', iend)
  call dump_scalar_i('jstr', jstr)
  call dump_scalar_i('jend', jend)
  call dump_scalar_i('kstr', kstr)
  call dump_scalar_i('kend', kend)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('itcon', itcon)
  call dump_array_3d('dvar.bin', dvar, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,k) reduction(max: intitc)

      do k=kstr,kend
      do j=jstr,jend
      do i=istr,iend
        intitc=max(abs(dvar(i,j,k)),intitc)
      end do
      end do
      end do

!$omp end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_chkitr == DUMP_TARGET_chkitr .and. .not. dump_done_chkitr) then
  call dump_finalize()
  dump_done_chkitr = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

! Hold the common control flag itcon.

      if(fproc(1:6).eq.'common') then

        call mpi_allreduce(intitc,itcon,1,mpi_real,mpi_max,             &
     &                     mpi_comm_cress,ierr)

      end if

! -----

      end subroutine s_chkitr

!-----7--------------------------------------------------------------7--

      end module m_chkitr

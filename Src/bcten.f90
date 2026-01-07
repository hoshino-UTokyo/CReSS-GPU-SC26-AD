!***********************************************************************
      module m_bcten
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/07/05
!     Modification: 1999/08/23, 1999/09/30, 1999/10/12, 1999/11/01,
!                   1999/11/19, 2000/01/17, 2001/06/06, 2001/12/11,
!                   2002/04/02, 2002/07/23, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for optional tensor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bcten, s_bcten

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcten

        module procedure s_bcten

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
      subroutine s_bcten(fpbbc,fptbc,ni,nj,nk,ten)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpbbc
                       ! Formal parameter of unique index of bbc

      integer, intent(in) :: fptbc
                       ! Formal parameter of unique index of tbc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: ten(0:ni+1,0:nj+1,1:nk)
                       ! Optional tensor

! Internal shared variables

      integer bbc      ! Option for bottom boundary conditions
      integer tbc      ! Option for top boundary conditions

      integer nkm2     ! nk - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_bcten = 0
      integer, parameter :: DUMP_TARGET_bcten = 1440
      logical, save :: dump_done_bcten = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpbbc,bbc)
      call getiname(fptbc,tbc)

! -----

! Set the common used variable.

      nkm2=nk-2

! -----

!! Set the bottom and top boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: bcten.f90 :: s_bcten
! Summary : Sets bottom and top boundary conditions for optional tensor array
!           by copying or negating values at boundary layers.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple array assignments to boundary planes (k=1, k=nk)
!   - No synchronization constructs beyond implicit barriers at omp end do
!   - Conditional branches based on bbc/tbc values (control flow divergence)
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Use collapse(2) for nested i,j loops to increase parallelism
! Runtime:
!   - Calls: 1440
!   - AvgLoops: 811.8K
!   - TotalTime: 0.068s (0.00%)
!   - AvgTime: 0.047ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('bcten.f90', 's_bcten', &
   & 'OMP section 1')
end if
loop_len = int((nj+1)-(0)+1,8) * int((ni+1)-(0)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_bcten = dump_call_count_bcten + 1
if (dump_call_count_bcten == DUMP_TARGET_bcten .and. .not. dump_done_bcten) then
  call dump_init('bcten')
  call dump_scalar_i('fpbbc', fpbbc)
  call dump_scalar_i('fptbc', fptbc)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('ten_in.bin', ten, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared)

! Set the bottom boundary conditions.

      if(bbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj+1
        do i=0,ni+1
          ten(i,j,1)=-ten(i,j,3)
        end do
        end do

!$omp end do

      else if(bbc.ge.3) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj+1
        do i=0,ni+1
          ten(i,j,1)=ten(i,j,3)
        end do
        end do

!$omp end do

      end if

! -----

! Set the top boundary conditions.

      if(tbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj+1
        do i=0,ni+1
          ten(i,j,nk)=-ten(i,j,nkm2)
        end do
        end do

!$omp end do

      else if(tbc.ge.3) then

!$omp do schedule(runtime) private(i,j)

        do j=0,nj+1
        do i=0,ni+1
          ten(i,j,nk)=ten(i,j,nkm2)
        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_bcten == DUMP_TARGET_bcten .and. .not. dump_done_bcten) then
  call dump_array_3d('ten_ref.bin', ten, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_bcten = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_bcten

!-----7--------------------------------------------------------------7--

      end module m_bcten

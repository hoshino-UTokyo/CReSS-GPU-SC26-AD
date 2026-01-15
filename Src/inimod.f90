!***********************************************************************
      module m_inimod
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/11/05
!     Modification: 2006/12/04, 2007/01/05, 2007/01/20, 2007/01/31,
!                   2007/07/30, 2007/10/19, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2009/01/30, 2009/02/27, 2011/09/22, 2011/11/10,
!                   2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     initialize the common used module variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comdmp
      use m_comprofile
      use m_dump_kernel
      use m_comerr
      use m_comname
      use m_comsave
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: inimod, s_inimod

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface inimod

        module procedure s_inimod

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
      subroutine s_inimod
!***********************************************************************

! Internal shared variables

      integer in       ! Namelist table index

      integer ierr     ! Index of table to check namelist error

! Internal private variable

      integer in_sub   ! Substitute for in


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_inimod = 0
      integer, parameter :: DUMP_TARGET_inimod = 1
      logical, save :: dump_done_inimod = .false.


!-----7--------------------------------------------------------------7--

!!! Initialize the common used module variables.

! For the module file, m_comdmp.

      call inichar(fl3c)
      call inichar(fl3d)

      call inichar(fl2c)
      call inichar(fl2d)

      write(fdmp(1:3),'(a3)') 'off'
      write(fmon(1:3),'(a3)') 'off'

      write(border(1:7),'(a7)') 'unknown'

      nc3d=0
      nc3c=0

      nc2d=0
      nc2c=0

      io3c=0
      io3d=0

      io2c=0
      io2d=0

      rec3d=0
      rec2d=0

      cnt3d=0
      cnt2d=0

! -----

! For the module file, m_comerr.

      do ierr=1,nerr

        write(errlst(ierr)(1:14),'(a14)') '              '

      end do

! -----

!! For the module file, m_comname.

! For the extra defined variables.

      iname(-1)=0
      iname(0)=0

! -----

! For the character variables.

      do in=1,ncn

        call inichar(cname(in))
        call inichar(rcname(in))

      end do

! -----

! For the integer and real variables.

!@llm start meta_info ----------------------------------------------------
! Location: inimod.f90 :: s_inimod
! Summary : Initialize integer (iname, riname) and real (rname, rrname) namelist
!           table arrays to zero
! GPU diff: Easy
! Findings:
!   - Two simple loops initializing namelist arrays to zero
!   - No function calls within the parallel region
!   - No global writes beyond array initialization
!   - No sync constructs or thread-dependent logic
! Next:
!   - Can be directly ported to GPU with OpenACC parallel loop
!   - Consider using array syntax for simpler GPU offload
! Runtime:
!   - Calls: 1
!   - AvgLoops: 219
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.021ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('inimod.f90', 's_inimod', &
   & 'OMP section 1')
end if
loop_len = int((nin)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_inimod = dump_call_count_inimod + 1
if (dump_call_count_inimod == DUMP_TARGET_inimod .and. .not. dump_done_inimod) then
  call dump_init('inimod')
  ! FIXME: iname is array - call dump_scalar_r('iname', iname)
  ! FIXME: riname is array - call dump_scalar_r('riname', riname)
  ! FIXME: rname is array - call dump_scalar_r('rname', rname)
  ! FIXME: rrname is array - call dump_scalar_r('rrname', rrname)
end if

!$omp parallel default(shared)

!$omp do schedule(runtime) private(in_sub)

      do in_sub=1,nin
        iname(in_sub)=0
        riname(in_sub)=0
      end do

!$omp end do

!$omp do schedule(runtime) private(in_sub)

      do in_sub=1,nrn
        rname(in_sub)=0.e0
        rrname(in_sub)=0.e0
      end do

!$omp end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_inimod == DUMP_TARGET_inimod .and. .not. dump_done_inimod) then
  call dump_finalize()
  dump_done_inimod = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

!! -----

! For the module file, m_comsave.

      nxtgpv=0

      nxtasl=0

      nxtrdr(1)=0
      nxtrdr(2)=0

      nxtsst=0

      extcom=0

! -----

!!! -----

      end subroutine s_inimod

!-----7--------------------------------------------------------------7--

      end module m_inimod

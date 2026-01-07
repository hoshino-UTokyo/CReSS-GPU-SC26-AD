!***********************************************************************
      module m_bcyclex
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/02/10
!     Modification: 2006/12/04, 2007/01/05, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the periodic boundary conditions in x direction.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comprofile
      use m_dump_kernel
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bcyclex, s_bcyclex

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcyclex

        module procedure s_bcyclex

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
      subroutine s_bcyclex(fpwbc,fpebc,iwsnd,iwrcv,iesnd,iercv,         &
     &                     ni,nj,kmax,var)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: iwsnd
                       ! Array index of west sended value

      integer, intent(in) :: iwrcv
                       ! Array index of west received value

      integer, intent(in) :: iesnd
                       ! Array index of east sended value

      integer, intent(in) :: iercv
                       ! Array index of east received value

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

! Input and output variable

      real, intent(inout) :: var(0:ni+1,0:nj+1,1:kmax)
                       ! Optional variable

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

! Internal private variables

      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_bcyclex = 0
      integer, parameter :: DUMP_TARGET_bcyclex = 4
      logical, save :: dump_done_bcyclex = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)

! -----

! Set the periodic boundary conditions in x direction.

!@llm start meta_info ----------------------------------------------------
! Location: bcyclex.f90 :: s_bcyclex
! Summary : Sets periodic boundary conditions in x direction by copying
!           values between west and east boundaries for cyclic domains.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Module variable nisub from m_commpi used (read-only)
!   - Simple 1D array copy operations along j-dimension
!   - Sequential k-loop with parallel j loops inside
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Consider collapsing k-loop with j-loop for better GPU utilization
!   - Ensure nisub is mapped or use firstprivate
! Runtime:
!   - Calls: 4
!   - AvgLoops: 86.7K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.004ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('bcyclex.f90', 's_bcyclex', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(1)+1,8) * int((nj+1)-(0)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_bcyclex = dump_call_count_bcyclex + 1
if (dump_call_count_bcyclex == DUMP_TARGET_bcyclex .and. .not. dump_done_bcyclex) then
  call dump_init('bcyclex')
  call dump_scalar_i('fpwbc', fpwbc)
  call dump_scalar_i('fpebc', fpebc)
  call dump_scalar_i('iwsnd', iwsnd)
  call dump_scalar_i('iwrcv', iwrcv)
  call dump_scalar_i('iesnd', iesnd)
  call dump_scalar_i('iercv', iercv)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('kmax', kmax)
  call dump_array_3d('var_in.bin', var, 0, ni+1, 0, nj+1, 1, kmax)
end if

!$omp parallel default(shared) private(k)

      if(nisub.eq.1) then

        if(wbc.eq.-1.and.ebc.eq.-1) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var(iwrcv,j,k)=var(iesnd,j,k)
            end do

!$omp end do

          end do

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var(iercv,j,k)=var(iwsnd,j,k)
            end do

!$omp end do

          end do

        end if

      end if

!$omp end parallel

! Dump output data at target call
if (dump_call_count_bcyclex == DUMP_TARGET_bcyclex .and. .not. dump_done_bcyclex) then
  call dump_array_3d('var_ref.bin', var, 0, ni+1, 0, nj+1, 1, kmax)
  call dump_finalize()
  dump_done_bcyclex = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_bcyclex

!-----7--------------------------------------------------------------7--

      end module m_bcyclex

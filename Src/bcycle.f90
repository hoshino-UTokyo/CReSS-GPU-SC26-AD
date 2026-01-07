!***********************************************************************
      module m_bcycle
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 1999/03/25, 1999/04/06, 1999/07/05, 1999/08/18,
!                   1999/08/23, 1999/10/07, 1999/10/12, 1999/11/01,
!                   1999/11/19, 2000/01/17, 2000/04/18, 2000/07/05,
!                   2001/05/29, 2002/04/02, 2002/06/06, 2003/04/30,
!                   2003/05/19, 2004/08/20, 2005/02/10, 2006/12/04,
!                   2007/01/05, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the periodic boundary conditions.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comprofile
      use m_getiname
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bcycle, s_bcycle

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcycle

        module procedure s_bcycle

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
      subroutine s_bcycle(fpwbc,fpebc,fpsbc,fpnbc,                      &
     &                    iwsnd,iwrcv,iesnd,iercv,                      &
     &                    jssnd,jsrcv,jnsnd,jnrcv,ni,nj,kmax,var)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: iwsnd
                       ! Array index of west sended value

      integer, intent(in) :: iwrcv
                       ! Array index of west received value

      integer, intent(in) :: iesnd
                       ! Array index of east sended value

      integer, intent(in) :: iercv
                       ! Array index of east received value

      integer, intent(in) :: jssnd
                       ! Array index of south sended value

      integer, intent(in) :: jsrcv
                       ! Array index of south received value

      integer, intent(in) :: jnsnd
                       ! Array index of north sended value

      integer, intent(in) :: jnrcv
                       ! Array index of north received value

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
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_bcycle = 0
      integer, parameter :: DUMP_TARGET_bcycle = 72374
      logical, save :: dump_done_bcycle = .false.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

!! Set the periodic boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: bcycle.f90 :: s_bcycle
! Summary : Sets periodic boundary conditions by copying values between
!           west/east and south/north boundaries for cyclic domains.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Module variables nisub, njsub from m_commpi used (read-only)
!   - Simple array copy operations for boundary planes
!   - Sequential k-loop with parallel j or i loops inside
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Convert to OpenACC with Unified Memory (no explicit data transfer needed)
!   - Consider collapsing k-loop with inner loop for better GPU utilization
!   - Ensure nisub, njsub are mapped or use firstprivate
! Runtime:
!   - Calls: 72374
!   - AvgLoops: 115.3K
!   - TotalTime: 0.301s (0.01%)
!   - AvgTime: 0.004ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('bcycle.f90', 's_bcycle', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(1)+1,8) * int((nj+1)-(0)+1,8)

! Dump input data at target call
dump_call_count_bcycle = dump_call_count_bcycle + 1
if (dump_call_count_bcycle == DUMP_TARGET_bcycle .and. .not. dump_done_bcycle) then
  call dump_init('bcycle')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('kmax', kmax)
  call dump_scalar_i('wbc', wbc)
  call dump_scalar_i('ebc', ebc)
  call dump_scalar_i('sbc', sbc)
  call dump_scalar_i('nbc', nbc)
  call dump_scalar_i('iwsnd', iwsnd)
  call dump_scalar_i('iwrcv', iwrcv)
  call dump_scalar_i('iesnd', iesnd)
  call dump_scalar_i('iercv', iercv)
  call dump_scalar_i('jssnd', jssnd)
  call dump_scalar_i('jsrcv', jsrcv)
  call dump_scalar_i('jnsnd', jnsnd)
  call dump_scalar_i('jnrcv', jnrcv)
  call dump_array_3d('var_in.bin', var, 0, ni+1, 0, nj+1, 1, kmax)
end if

call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

! Set the west and east boundary conditions.

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

! -----

! Set the south and north boundary conditions.

      if(njsub.eq.1) then

        if(sbc.eq.-1.and.nbc.eq.-1) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var(i,jsrcv,k)=var(i,jnsnd,k)
            end do

!$omp end do

          end do

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var(i,jnrcv,k)=var(i,jssnd,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! Dump output data at target call
if (dump_call_count_bcycle == DUMP_TARGET_bcycle .and. .not. dump_done_bcycle) then
  call dump_array_3d('var_ref.bin', var, 0, ni+1, 0, nj+1, 1, kmax)
  call dump_finalize()
  dump_done_bcycle = .true.
end if

!! -----

      end subroutine s_bcycle

!-----7--------------------------------------------------------------7--

      end module m_bcycle

!***********************************************************************
      module m_bc8u
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/07/28, 1999/08/03, 1999/08/18, 1999/08/23,
!                   1999/09/30, 1999/10/07, 1999/11/01, 2000/01/17,
!                   2001/12/11, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2004/08/20, 2006/12/04, 2007/01/05, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for optional variable at the u points.

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

      public :: bc8u, s_bc8u

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bc8u

        module procedure s_bc8u

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
      subroutine s_bc8u(fpwbc,fpebc,ni,nj,kmax,var8u)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

! Input and output variable

      real, intent(inout) :: var8u(0:ni+1,0:nj+1,1:kmax)
                       ! Optional variable at u points

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2

! Internal private variables

      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_bc8u = 0
      integer, parameter :: DUMP_TARGET_bc8u = 720
      logical, save :: dump_done_bc8u = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2

! -----

!! Set the west and east boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: bc8u.f90 :: s_bc8u
! Summary : Sets west and east boundary conditions for optional variable at u points
!           by copying from adjacent interior points based on BC type.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Nested loops with outer k-loop serial, inner j-loop parallelized
!   - Uses module variables from m_commpi (ebw, ebe, isub, nisub)
!   - Conditional execution based on BC type (wbc, ebc) and subdomain position
! Next:
!   - Convert to OpenACC with collapsed j,k loops
!   - Restructure loops to have k as inner loop for better GPU coalescing
! Runtime:
!   - Calls: 720
!   - AvgLoops: 115.3K
!   - TotalTime: 1.712s (0.06%)
!   - AvgTime: 2.378ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('bc8u.f90', 's_bc8u', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(1)+1,8) * int((nj+1)-(0)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_bc8u = dump_call_count_bc8u + 1
if (dump_call_count_bc8u == DUMP_TARGET_bc8u .and. .not. dump_done_bc8u) then
  call dump_init('bc8u')
  call dump_scalar_i('wbc', wbc)
  call dump_scalar_i('ebc', ebc)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('kmax', kmax)
  call dump_array_3d('var8u_in.bin', var8u, 0, ni+1, 0, nj+1, 1, kmax)
  call dump_scalar_i('ebe', ebe)
  call dump_scalar_i('ebw', ebw)
  call dump_scalar_i('isub', isub)
  call dump_scalar_i('nim1', nim1)
  call dump_scalar_i('nim2', nim2)
  call dump_scalar_i('nisub', nisub)
end if

!$omp parallel default(shared) private(k)

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0) then

        if(wbc.eq.2) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var8u(1,j,k)=var8u(3,j,k)
            end do

!$omp end do

          end do

        else if(wbc.ge.3) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var8u(1,j,k)=var8u(2,j,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(ebc.eq.2) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var8u(ni,j,k)=var8u(nim2,j,k)
            end do

!$omp end do

          end do

        else if(ebc.ge.3) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var8u(ni,j,k)=var8u(nim1,j,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_bc8u == DUMP_TARGET_bc8u .and. .not. dump_done_bc8u) then
  call dump_array_3d('var8u_ref.bin', var8u, 0, ni+1, 0, nj+1, 1, kmax)
  call dump_finalize()
  dump_done_bc8u = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_bc8u

!-----7--------------------------------------------------------------7--

      end module m_bc8u

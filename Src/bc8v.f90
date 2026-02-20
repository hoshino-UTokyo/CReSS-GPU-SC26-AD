!***********************************************************************
      module m_bc8v
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
!     set the boundary conditions for optional variable at the v points.

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

      public :: bc8v, s_bc8v

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bc8v

        module procedure s_bc8v

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
      subroutine s_bc8v(fpsbc,fpnbc,ni,nj,kmax,var8v)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

! Input and output variable

      real, intent(inout) :: var8v(0:ni+1,0:nj+1,1:kmax)
                       ! Optional variable at v points

! Internal shared variables

      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer njm1     ! nj - 1
      integer njm2     ! nj - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_bc8v = 0
      integer, parameter :: DUMP_TARGET_bc8v = 720
      logical, save :: dump_done_bc8v = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

! Set the common used variables.

      njm1=nj-1
      njm2=nj-2

! -----

!! Set the south and north boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: bc8v.f90 :: s_bc8v
! Summary : Sets south and north boundary conditions for optional variable at v points
!           by copying from adjacent interior points based on BC type.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Nested loops with outer k-loop serial, inner i-loop parallelized
!   - Uses module variables from m_commpi (ebs, ebn, jsub, njsub)
!   - Conditional execution based on BC type (sbc, nbc) and subdomain position
! Next:
!   - Convert to OpenACC with collapsed i,k loops
!   - Restructure loops to have k as inner loop for better GPU coalescing
! Runtime:
!   - Calls: 720
!   - AvgLoops: 115.3K
!   - TotalTime: 1.731s (0.06%)
!   - AvgTime: 2.404ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('bc8v.f90', 's_bc8v', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(1)+1,8) * int((ni+1)-(0)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_bc8v = dump_call_count_bc8v + 1
if (dump_call_count_bc8v == DUMP_TARGET_bc8v .and. .not. dump_done_bc8v) then
  call dump_init('bc8v')
  call dump_scalar_i('sbc', sbc)
  call dump_scalar_i('nbc', nbc)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('kmax', kmax)
  call dump_array_3d('var8v_in.bin', var8v, 0, ni+1, 0, nj+1, 1, kmax)
  call dump_scalar_i('ebn', ebn)
  call dump_scalar_i('ebs', ebs)
  call dump_scalar_i('jsub', jsub)
  call dump_scalar_i('njm1', njm1)
  call dump_scalar_i('njm2', njm2)
  call dump_scalar_i('njsub', njsub)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_029)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    if (ebs == 1 .and. jsub == 0) then

      if (sbc == 2) then
        !$acc kernels
        !$acc loop independent
        do k = 1, kmax
          !$acc loop independent
          do i = 0, ni+1
            var8v(i,1,k) = var8v(i,3,k)
          end do
        end do
        !$acc end kernels

      else if (sbc >= 3) then
        !$acc kernels
        !$acc loop independent
        do k = 1, kmax
          !$acc loop independent
          do i = 0, ni+1
            var8v(i,1,k) = var8v(i,2,k)
          end do
        end do
        !$acc end kernels

      end if

    end if

    ! Set the north boundary conditions
    if (ebn == 1 .and. jsub == njsub-1) then

      if (nbc == 2) then
        !$acc kernels
        !$acc loop independent
        do k = 1, kmax
          !$acc loop independent
          do i = 0, ni+1
            var8v(i,nj,k) = var8v(i,njm2,k)
          end do
        end do
        !$acc end kernels

      else if (nbc >= 3) then
        !$acc kernels
        !$acc loop independent
        do k = 1, kmax
          !$acc loop independent
          do i = 0, ni+1
            var8v(i,nj,k) = var8v(i,njm1,k)
          end do
        end do
        !$acc end kernels

      end if

    end if

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared) private(k)

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.eq.2) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var8v(i,1,k)=var8v(i,3,k)
            end do

!$omp end do

          end do

        else if(sbc.ge.3) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var8v(i,1,k)=var8v(i,2,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(nbc.eq.2) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var8v(i,nj,k)=var8v(i,njm2,k)
            end do

!$omp end do

          end do

        else if(nbc.ge.3) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var8v(i,nj,k)=var8v(i,njm1,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_bc8v == DUMP_TARGET_bc8v .and. .not. dump_done_bc8v) then
  call dump_array_3d('var8v_ref.bin', var8v, 0, ni+1, 0, nj+1, 1, kmax)
  call dump_finalize()
  dump_done_bc8v = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_bc8v

!-----7--------------------------------------------------------------7--

      end module m_bc8v

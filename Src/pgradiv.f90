!***********************************************************************
      module m_pgradiv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/12/20
!     Modification: 2000/01/17, 2001/06/29, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2003/12/12, 2003/12/26, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the pressure gradient force vertically with the
!     horizontally explicit and vertically implicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: pgradiv, s_pgradiv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface pgradiv

        module procedure s_pgradiv

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
      subroutine s_pgradiv(fpdziv,fpweicoe,dts,ni,nj,nk,jcb,wfrc,fp,fw, &
     &                     fpdvj)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: fpweicoe
                       ! Formal parameter of unique index of weicoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation in large time steps

      real, intent(in) :: fp(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

! Input and output variable

      real, intent(inout) :: fw(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      real dziv        ! Inverce of dz

      real weicoe      ! Weighting coefficient for implicit method

      real dtdzw       ! dts x dziv x weicoe

      real, intent(inout) :: fpdvj(0:ni+1,0:nj+1,1:nk)
                       ! fp / jcb

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_pgradiv = 0
      integer, parameter :: DUMP_TARGET_pgradiv = 14400
      logical, save :: dump_done_pgradiv = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpdziv,dziv)
      call getrname(fpweicoe,weicoe)

! -----

! Set the common used variable.

      dtdzw=dts*dziv*weicoe

! -----

! Calculate the pressure gradient force vertically.

!@llm start meta_info ----------------------------------------------------
! Location: pgradiv.f90 :: s_pgradiv
! Summary : Compute vertical pressure gradient force using implicit method,
!           updating forcing term for w equation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region (pure arithmetic)
!   - No global/module variable writes (only intent(inout) fw, fpdvj arrays)
!   - No sync constructs
!   - Two separate k-loops: first computes fpdvj, second updates fw
!   - Second loop has k-1 dependency on fpdvj (read from previous level)
! Next:
!   - Convert to OpenACC with collapse for k,j,i loops
!   - First loop is independent; second loop needs fpdvj from k-1 level
!   - Can fuse loops or ensure proper synchronization between them
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 100.4M
!   - TotalTime: 117.847s (3.95%)
!   - AvgTime: 8.184ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('pgradiv.f90', 's_pgradiv', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) &
     & * int((nj-2)-(2)+1,8) &
     & * int((ni-2)-(2)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_pgradiv = dump_call_count_pgradiv + 1
if (dump_call_count_pgradiv == DUMP_TARGET_pgradiv .and. .not. dump_done_pgradiv) then
  call dump_init('pgradiv')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('dts', dts)
  call dump_scalar_r('dziv', dziv)
  call dump_scalar_r('weicoe', weicoe)
  call dump_array_3d('jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wfrc.bin', wfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('fp.bin', fp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('fw_in.bin', fw, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('fpdvj_in.bin', fpdvj, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('dtdzw', dtdzw)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_231)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    ! First loop: compute fpdvj
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          fpdvj(i,j,k) = fp(i,j,k) / jcb(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Second loop: update fw (depends on fpdvj, so separate kernels region)
    !$acc kernels
    !$acc loop independent
    do k = 3, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          fw(i,j,k) = wfrc(i,j,k) + fw(i,j,k) &
               + (fpdvj(i,j,k-1) - fpdvj(i,j,k)) * dtdzw
        end do
      end do
    end do
    !$acc end kernels

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          fpdvj(i,j,k)=fp(i,j,k)/jcb(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=3,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          fw(i,j,k)                                                     &
     &      =wfrc(i,j,k)+fw(i,j,k)+(fpdvj(i,j,k-1)-fpdvj(i,j,k))*dtdzw
        end do
        end do

!$omp end do

      end do

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_pgradiv == DUMP_TARGET_pgradiv .and. .not. dump_done_pgradiv) then
  call dump_array_3d('fw_ref.bin', fw, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('fpdvj_ref.bin', fpdvj, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_pgradiv = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_pgradiv

!-----7--------------------------------------------------------------7--

      end module m_pgradiv

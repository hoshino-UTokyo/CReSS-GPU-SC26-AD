!***********************************************************************
      module m_coriuv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/04/18, 2002/04/02, 2002/09/02,
!                   2003/01/04, 2003/04/30, 2003/05/19, 2003/11/28,
!                   2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the Coriolis force in the x and the y components of
!     velocity equation.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing
      use m_comprofile
      use m_dump_kernel

      implicit none

! Default access control

      private

! Exceptional access control

      public :: coriuv, s_coriuv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface coriuv

        module procedure s_coriuv

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
      subroutine s_coriuv(ni,nj,nk,fc,rst,u,v,ufrc,vfrc,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

! Internal shared variable

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_coriuv = 0
      integer, parameter :: DUMP_TARGET_coriuv = 360
      logical, save :: dump_done_coriuv = .false.


!-----7--------------------------------------------------------------7--

!! Calculate the Coriolis force in the x and the y components of
!! velocity equation.

!@llm start meta_info ----------------------------------------------------
! Location: coriuv.f90 :: subroutine s_coriuv
! Summary : Calculates Coriolis force terms for u and v velocity equations
!           using f-plane or beta-plane approximation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Two-stage calculation: (1) compute tmp1, (2) add to forcing.
!   - Stencil averaging of velocity for Coriolis term.
!   - All grid points are independent within each loop nest.
! Next:
!   - Direct OpenACC kernels with collapse(2) on i,j loops.
!   - Consider fusing the two stages into single kernel per component.
!   - Temporary array tmp1 needed for staggered grid averaging.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.6M
!   - TotalTime: 5.402s (0.18%)
!   - AvgTime: 15.006ms
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('coriuv.f90', 's_coriuv', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) &
     & * int((nj-2)-(2)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_coriuv = dump_call_count_coriuv + 1
if (dump_call_count_coriuv == DUMP_TARGET_coriuv .and. .not. dump_done_coriuv) then
  call dump_init('coriuv')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('fc.bin', fc, 0, ni+1, 0, nj+1, 1, 2)
  call dump_array_3d('rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ufrc_in.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc_in.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_in.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_069)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------

! Calculate the Coriolis force in the x components of velocity equation.

    !$acc kernels
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 1, ni-1
          tmp1(i,j,k) = fc(i,j,1) * rst(i,j,k) * (v(i,j,k) + v(i,j+1,k))
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-1
          ufrc(i,j,k) = ufrc(i,j,k) + (tmp1(i-1,j,k) + tmp1(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

! Calculate the Coriolis force in the y components of velocity equation.

    !$acc kernels
    do k = 2, nk-2
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 2, ni-2
          tmp1(i,j,k) = -fc(i,j,1) * rst(i,j,k) * (u(i,j,k) + u(i+1,j,k))
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-1
        !$acc loop independent
        do i = 2, ni-2
          vfrc(i,j,k) = vfrc(i,j,k) + (tmp1(i,j-1,k) + tmp1(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared) private(k)

! Calculate the Coriolis force in the x components of velocity equation.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=1,ni-1
          tmp1(i,j,k)=fc(i,j,1)*rst(i,j,k)*(v(i,j,k)+v(i,j+1,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          ufrc(i,j,k)=ufrc(i,j,k)+(tmp1(i-1,j,k)+tmp1(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the Coriolis force in the y components of velocity equation.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=2,ni-2
          tmp1(i,j,k)=-fc(i,j,1)*rst(i,j,k)*(u(i,j,k)+u(i+1,j,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          vfrc(i,j,k)=vfrc(i,j,k)+(tmp1(i,j-1,k)+tmp1(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_coriuv == DUMP_TARGET_coriuv .and. .not. dump_done_coriuv) then
  call dump_array_3d('ufrc_ref.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc_ref.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_ref.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_coriuv = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_coriuv

!-----7--------------------------------------------------------------7--

      end module m_coriuv

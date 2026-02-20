!***********************************************************************
      module m_diver2d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/12/17
!     Modification: 2000/01/17, 2001/02/24, 2001/06/29, 2002/04/02,
!                   2003/01/04, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2004/06/10, 2006/11/06, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/09, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the negative divergence horizontally.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname
      use m_comprofile
      use m_getrname
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diver2d, s_diver2d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diver2d

        module procedure s_diver2d

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
      subroutine s_diver2d(fpmpopt,fpmfcopt,fpdxiv,fpdyiv,ni,nj,nk,     &
     &                     mf,rmf,rmf8u,rmf8v,var8u,var8v,u,v,          &
     &                     div2d,tmp1,tmp2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: var8u(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at u points

      real, intent(in) :: var8v(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at v points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

! Output variable

      real, intent(out) :: div2d(0:ni+1,0:nj+1,1:nk)
                       ! 2 dimensional divergence

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_diver2d = 0
      integer, parameter :: DUMP_TARGET_diver2d = 14400
      logical, save :: dump_done_diver2d = .false.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)

! -----

!! Calculate the negative divergence horizontally.

!@llm start meta_info ----------------------------------------------------
! Location: diver2d.f90 :: s_diver2d
! Summary : Calculate 2D horizontal negative divergence using velocity
!           components u,v weighted by Jacobian and map scale factors.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - No global/module variable writes
!   - No synchronization constructs
!   - Multiple branches (mfcopt, mpopt) but all are simple data-parallel loops
!   - Two-phase computation: multiply then difference
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - tmp1, tmp2 are temporary arrays that can be fused or kept on GPU
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 102.5M
!   - TotalTime: 162.333s (5.45%)
!   - AvgTime: 11.273ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('diver2d.f90', 's_diver2d', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni)-(1)+1,8)

! Dump input data at target call
dump_call_count_diver2d = dump_call_count_diver2d + 1
if (dump_call_count_diver2d == DUMP_TARGET_diver2d .and. .not. dump_done_diver2d) then
  call dump_init('diver2d')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('mpopt', mpopt)
  call dump_scalar_i('mfcopt', mfcopt)
  call dump_scalar_r('dxiv', dxiv)
  call dump_scalar_r('dyiv', dyiv)
  call dump_array_2d('mf.bin', mf, 0, ni+1, 0, nj+1)
  call dump_array_3d('rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call dump_array_3d('rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('var8u.bin', var8u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('var8v.bin', var8v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  ! FIXME: tmp1 is an array, not scalar
  ! ! FIXME: tmp1 is array - call dump_scalar_i('tmp1', tmp1)
  ! FIXME: tmp2 is an array, not scalar
  ! ! FIXME: tmp2 is array - call dump_scalar_r('tmp2', tmp2)
end if

call profile_start(prof_id1)

#if defined(USE_GPU) && !defined(DISABLE_GPU_090)
! GPU version (OpenACC)

! Optional variables at u, v and w points are multiplyed by u, v and wc.

      if(mfcopt.eq.0) then

        !$acc kernels
        do k=1,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1,ni
              tmp1(i,j,k)=var8u(i,j,k)*u(i,j,k)
            end do
          end do
          !$acc loop independent
          do j=1,nj
            !$acc loop independent
            do i=1,ni-1
              tmp2(i,j,k)=var8v(i,j,k)*v(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else

        if(mpopt.eq.0.or.mpopt.eq.10) then

          !$acc kernels
          do k=1,nk-1
            !$acc loop independent
            do j=1,nj-1
              !$acc loop independent
              do i=1,ni
                tmp1(i,j,k)=var8u(i,j,k)*u(i,j,k)
              end do
            end do
            !$acc loop independent
            do j=1,nj
              !$acc loop independent
              do i=1,ni-1
                tmp2(i,j,k)=rmf8v(i,j,2)*var8v(i,j,k)*v(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

        else if(mpopt.eq.5) then

          !$acc kernels
          do k=1,nk-1
            !$acc loop independent
            do j=1,nj-1
              !$acc loop independent
              do i=1,ni
                tmp1(i,j,k)=rmf8u(i,j,2)*var8u(i,j,k)*u(i,j,k)
              end do
            end do
            !$acc loop independent
            do j=1,nj
              !$acc loop independent
              do i=1,ni-1
                tmp2(i,j,k)=var8v(i,j,k)*v(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

        else

          !$acc kernels
          do k=1,nk-1
            !$acc loop independent
            do j=1,nj-1
              !$acc loop independent
              do i=1,ni
                tmp1(i,j,k)=rmf8u(i,j,2)*var8u(i,j,k)*u(i,j,k)
              end do
            end do
            !$acc loop independent
            do j=1,nj
              !$acc loop independent
              do i=1,ni-1
                tmp2(i,j,k)=rmf8v(i,j,2)*var8v(i,j,k)*v(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

        end if

      end if

! -----

! Calculate the divergence horizontally.

      if(mfcopt.eq.0) then

        !$acc kernels
        do k=1,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1,ni-1
              div2d(i,j,k)=(tmp1(i,j,k)-tmp1(i+1,j,k))*dxiv               &
     &          +(tmp2(i,j,k)-tmp2(i,j+1,k))*dyiv
            end do
          end do
        end do
        !$acc end kernels

      else

        if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

          !$acc kernels
          do k=1,nk-1
            !$acc loop independent
            do j=1,nj-1
              !$acc loop independent
              do i=1,ni-1
                div2d(i,j,k)=mf(i,j)*((tmp1(i,j,k)-tmp1(i+1,j,k))*dxiv    &
     &            +(tmp2(i,j,k)-tmp2(i,j+1,k))*dyiv)
              end do
            end do
          end do
          !$acc end kernels

        else

          !$acc kernels
          do k=1,nk-1
            !$acc loop independent
            do j=1,nj-1
              !$acc loop independent
              do i=1,ni-1
                div2d(i,j,k)=rmf(i,j,1)*((tmp1(i,j,k)-tmp1(i+1,j,k))*dxiv &
     &            +(tmp2(i,j,k)-tmp2(i,j+1,k))*dyiv)
              end do
            end do
          end do
          !$acc end kernels

        end if

      end if

#else
! CPU version (OpenMP) - Original code preserved

!$omp parallel default(shared) private(k)

! Optional variables at u, v and w points are multiplyed by u, v and wc.

      if(mfcopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni
            tmp1(i,j,k)=var8u(i,j,k)*u(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=1,nj
          do i=1,ni-1
            tmp2(i,j,k)=var8v(i,j,k)*v(i,j,k)
          end do
          end do

!$omp end do

        end do

      else

        if(mpopt.eq.0.or.mpopt.eq.10) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni
              tmp1(i,j,k)=var8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni-1
              tmp2(i,j,k)=rmf8v(i,j,2)*var8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

        else if(mpopt.eq.5) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni
              tmp1(i,j,k)=rmf8u(i,j,2)*var8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni-1
              tmp2(i,j,k)=var8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni
              tmp1(i,j,k)=rmf8u(i,j,2)*var8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni-1
              tmp2(i,j,k)=rmf8v(i,j,2)*var8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Calculate the divergence horizontally.

      if(mfcopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            div2d(i,j,k)=(tmp1(i,j,k)-tmp1(i+1,j,k))*dxiv               &
     &        +(tmp2(i,j,k)-tmp2(i,j+1,k))*dyiv
          end do
          end do

!$omp end do

        end do

      else

        if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              div2d(i,j,k)=mf(i,j)*((tmp1(i,j,k)-tmp1(i+1,j,k))*dxiv    &
     &          +(tmp2(i,j,k)-tmp2(i,j+1,k))*dyiv)
            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              div2d(i,j,k)=rmf(i,j,1)*((tmp1(i,j,k)-tmp1(i+1,j,k))*dxiv &
     &          +(tmp2(i,j,k)-tmp2(i,j+1,k))*dyiv)
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel
#endif

call profile_stop(prof_id1, loop_len)

! Dump output data at target call
if (dump_call_count_diver2d == DUMP_TARGET_diver2d .and. .not. dump_done_diver2d) then
  call dump_array_3d('div2d_ref.bin', div2d, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_diver2d = .true.
end if

!! -----

      end subroutine s_diver2d

!-----7--------------------------------------------------------------7--

      end module m_diver2d

!***********************************************************************
      module m_phy2cnt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/18, 1999/08/23, 1999/09/14, 1999/09/30,
!                   1999/10/12, 1999/11/01, 1999/11/19, 1999/12/20,
!                   2000/01/17, 2000/12/18, 2001/06/06, 2001/11/20,
!                   2002/04/02, 2003/01/04, 2003/03/21, 2003/04/30,
!                   2003/05/19, 2003/11/05, 2003/12/12, 2004/06/10,
!                   2006/04/03, 2006/06/21, 2006/11/06, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/10/19, 2008/05/02,
!                   2008/06/09, 2008/08/25, 2009/02/27, 2011/07/15,
!                   2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the zeta components of contravariant velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_comprofile
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_lbcwc
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
      use m_vbcwc
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: phy2cnt, s_phy2cnt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface phy2cnt

        module procedure s_phy2cnt

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
      subroutine s_phy2cnt(fpsthopt,fptrnopt,fpmpopt,fpmfcopt,fpadvopt, &
     &                     ni,nj,nk,j31,j32,jcb8w,mf,u,v,w,wc,          &
     &                     mf25,j31u2,j32v2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsthopt
                       ! Formal parameter of unique index of sthopt

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Output variable

      real, intent(out) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Internal shared variables

      integer sthopt   ! Option for vertical grid stretching

      integer trnopt   ! Option for terrain height setting

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      integer advopt   ! Option for advection scheme

      real, intent(inout) :: mf25(0:ni+1,0:nj+1)
                       ! 0.25 x mf

      real, intent(inout) :: j31u2(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x j31 x u

      real, intent(inout) :: j32v2(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x j32 x v

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_phy2cnt = 0
      integer, parameter :: DUMP_TARGET_phy2cnt = 15121
      logical, save :: dump_done_phy2cnt = .false.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsthopt,sthopt)
      call getiname(fptrnopt,trnopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpadvopt,advopt)

! -----

!! Calculate the zeta components of contravariant velocity.

! Get the zeta components of contravariant velocity.

!@llm start meta_info ----------------------------------------------------
! Location: phy2cnt.f90 :: s_phy2cnt
! Summary : Calculate zeta components of contravariant velocity from
!           physical velocity components with terrain-following coordinates.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to output array wc, work arrays j31u2, j32v2, mf25
!   - Simple 3D stencil computations with straightforward data access
!   - Multiple conditional paths based on trnopt, sthopt, mfcopt, mpopt
!   - No explicit barriers but implicit at !$omp end do
! Next:
!   - Use OpenACC parallel loop with collapse(3) for 3D loops
!   - Straightforward GPU port with data region for arrays
!   - Consider kernel fusion for consecutive loops
! Runtime:
!   - Calls: 15121
!   - AvgLoops: 101.6M
!   - TotalTime: 42.305s (1.42%)
!   - AvgTime: 2.798ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('phy2cnt.f90', 's_phy2cnt', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(2)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)

! Dump input data at target call
dump_call_count_phy2cnt = dump_call_count_phy2cnt + 1
if (dump_call_count_phy2cnt == DUMP_TARGET_phy2cnt .and. .not. dump_done_phy2cnt) then
  call dump_init('phy2cnt')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('sthopt', sthopt)
  call dump_scalar_i('trnopt', trnopt)
  call dump_scalar_i('mpopt', mpopt)
  call dump_scalar_i('mfcopt', mfcopt)
  call dump_array_3d('j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('jcb8w.bin', jcb8w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('mf.bin', mf, 0, ni+1, 0, nj+1)
  call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  ! FIXME: j31u2 is an array, not scalar
  ! ! FIXME: j31u2 is array - call dump_scalar_r('j31u2', j31u2)
  ! FIXME: j32v2 is an array, not scalar
  ! ! FIXME: j32v2 is array - call dump_scalar_r('j32v2', j32v2)
  ! FIXME: mf25 is an array, not scalar
  ! ! FIXME: mf25 is array - call dump_scalar_i('mf25', mf25)
end if

call profile_start(prof_id1)

#if defined(USE_GPU) && !defined(DISABLE_GPU_236)
! GPU version (OpenACC)
      if(trnopt.eq.0) then

        if(sthopt.eq.0) then

!$acc kernels
!$acc loop independent
          do k=2,nk-1
!$acc loop independent
            do j=1,nj-1
!$acc loop independent
            do i=1,ni-1
              wc(i,j,k)=w(i,j,k)
            end do
            end do
          end do
!$acc end kernels

        else if(sthopt.ge.1) then

!$acc kernels
!$acc loop independent
          do k=2,nk-1
!$acc loop independent
            do j=1,nj-1
!$acc loop independent
            do i=1,ni-1
              wc(i,j,k)=w(i,j,k)/jcb8w(i,j,k)
            end do
            end do
          end do
!$acc end kernels

        end if

      else

!$acc kernels
!$acc loop independent
        do k=2,nk-1
!$acc loop independent
          do j=1,nj-1
!$acc loop independent
          do i=1,ni
            j31u2(i,j,k)=(u(i,j,k-1)+u(i,j,k))*j31(i,j,k)
          end do
          end do
        end do
!$acc end kernels

!$acc kernels
!$acc loop independent
        do k=2,nk-1
!$acc loop independent
          do j=1,nj
!$acc loop independent
          do i=1,ni-1
            j32v2(i,j,k)=(v(i,j,k-1)+v(i,j,k))*j32(i,j,k)
          end do
          end do
        end do
!$acc end kernels

        if(mfcopt.eq.0) then

!$acc kernels
!$acc loop independent
          do k=2,nk-1
!$acc loop independent
            do j=1,nj-1
!$acc loop independent
            do i=1,ni-1
              wc(i,j,k)=(.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k))           &
     &          +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
            end do
            end do
          end do
!$acc end kernels

        else

          if(mpopt.eq.0.or.mpopt.eq.10) then

!$acc kernels
!$acc loop independent
            do k=2,nk-1
!$acc loop independent
              do j=1,nj-1
!$acc loop independent
              do i=1,ni-1
                wc(i,j,k)=(.25e0*(mf(i,j)*(j31u2(i,j,k)+j31u2(i+1,j,k)) &
     &            +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
              end do
              end do
            end do
!$acc end kernels

          else if(mpopt.eq.5) then

!$acc kernels
!$acc loop independent
            do k=2,nk-1
!$acc loop independent
              do j=1,nj-1
!$acc loop independent
              do i=1,ni-1
                wc(i,j,k)=(.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k))         &
     &            +mf(i,j)*(j32v2(i,j,k)+j32v2(i,j+1,k)))               &
     &            +w(i,j,k))/jcb8w(i,j,k)
              end do
              end do
            end do
!$acc end kernels

          else

!$acc kernels
!$acc loop independent
            do j=1,nj-1
!$acc loop independent
            do i=1,ni-1
              mf25(i,j)=.25e0*mf(i,j)
            end do
            end do
!$acc end kernels

!$acc kernels
!$acc loop independent
            do k=2,nk-1
!$acc loop independent
              do j=1,nj-1
!$acc loop independent
              do i=1,ni-1
                wc(i,j,k)=(mf25(i,j)*((j31u2(i,j,k)+j31u2(i+1,j,k))     &
     &            +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
              end do
              end do
            end do
!$acc end kernels

          end if

        end if

      end if

#else
! CPU version (OpenMP) - Original code preserved
!$omp parallel default(shared) private(k)

      if(trnopt.eq.0) then

        if(sthopt.eq.0) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wc(i,j,k)=w(i,j,k)
            end do
            end do

!$omp end do

          end do

        else if(sthopt.ge.1) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wc(i,j,k)=w(i,j,k)/jcb8w(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      else

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni
            j31u2(i,j,k)=(u(i,j,k-1)+u(i,j,k))*j31(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=1,nj
          do i=1,ni-1
            j32v2(i,j,k)=(v(i,j,k-1)+v(i,j,k))*j32(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.0) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wc(i,j,k)=(.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k))           &
     &          +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                wc(i,j,k)=(.25e0*(mf(i,j)*(j31u2(i,j,k)+j31u2(i+1,j,k)) &
     &            +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
              end do
              end do

!$omp end do

            end do

          else if(mpopt.eq.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                wc(i,j,k)=(.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k))         &
     &            +mf(i,j)*(j32v2(i,j,k)+j32v2(i,j+1,k)))               &
     &            +w(i,j,k))/jcb8w(i,j,k)
              end do
              end do

!$omp end do

            end do

          else

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              mf25(i,j)=.25e0*mf(i,j)
            end do
            end do

!$omp end do

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                wc(i,j,k)=(mf25(i,j)*((j31u2(i,j,k)+j31u2(i+1,j,k))     &
     &            +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
              end do
              end do

!$omp end do

            end do

          end if

        end if

      end if

!$omp end parallel
#endif

call profile_stop(prof_id1, loop_len)

! Dump output data at target call
if (dump_call_count_phy2cnt == DUMP_TARGET_phy2cnt .and. .not. dump_done_phy2cnt) then
  call dump_array_3d('wc_ref.bin', wc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_phy2cnt = .true.
end if

! -----

!! -----

!! Exchange the value horizontally.

      if(advopt.ge.4) then

! Exchange the value horizontally between sub domain.

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,wc,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,wc,1,1,rbuf)

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,wc,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,wc,1,1,rbuf)

! -----

! Exchange the value horizontally betweein group domain.

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,wc,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,wc,1,1,rbuf)

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,wc,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,wc,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,wc,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,wc,1,1,rbuf)

! -----

      end if

!! -----

! Set the lateral boundary conditions.

      if(advopt.ge.4) then

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,wc)

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,ni,nj,nk,wc)

        call lbcwc(idwbc,idebc,idsbc,idnbc,ni,nj,nk,wc)

      end if

! -----

! Set the bottom and the top boundary conditions.

      call vbcwc(idbbc,idtbc,ni,nj,nk,wc)

! -----

      end subroutine s_phy2cnt

!-----7--------------------------------------------------------------7--

      end module m_phy2cnt

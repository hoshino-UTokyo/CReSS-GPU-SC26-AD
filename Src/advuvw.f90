!***********************************************************************
      module m_advuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/18, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/04/18, 2000/12/18, 2001/04/15,
!                   2001/05/29, 2001/06/06, 2001/08/07, 2001/11/20,
!                   2002/04/02, 2003/01/04, 2003/04/30, 2003/05/19,
!                   2003/11/28, 2004/02/01, 2004/03/05, 2004/04/15,
!                   2004/06/10, 2004/08/01, 2004/09/10, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the velocity advection.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_dump_kernel
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: advuvw, s_advuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface advuvw

        module procedure s_advuvw

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
      subroutine s_advuvw(fpadvopt,fpiwest,fpieast,fpjsouth,fpjnorth,   &
     &                    fpdxiv,fpdyiv,fpdziv,ni,nj,nk,rstxu,rstxv,    &
     &                    rstxwc,u,v,w,ufrc,vfrc,wfrc,hadv,vadv,        &
     &                    tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpiwest
                       ! Formal parameter of unique index of iwest

      integer, intent(in) :: fpieast
                       ! Formal parameter of unique index of ieast

      integer, intent(in) :: fpjsouth
                       ! Formal parameter of unique index of jsouth

      integer, intent(in) :: fpjnorth
                       ! Formal parameter of unique index of jnorth

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rstxu(0:ni+1,0:nj+1,1:nk)
                       ! u x base state density x Jacobian

      real, intent(in) :: rstxv(0:ni+1,0:nj+1,1:nk)
                       ! v x base state density x Jacobian

      real, intent(in) :: rstxwc(0:ni+1,0:nj+1,1:nk)
                       ! wc x base state density x Jacobian

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      integer advopt   ! Option for advection scheme

      integer iwest    ! Added index on west boundary
      integer jsouth   ! Added index on south boundary

      integer ieast    ! Subtracted index on east boundary
      integer jnorth   ! Subtracted index on north boundary

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real dxv25n      ! - 0.25 x dxiv
      real dyv25n      ! - 0.25 x dyiv
      real dzv25n      ! - 0.25 x dziv

      real dxv24       ! dxiv / 24.0
      real dyv24       ! dyiv / 24.0
      real dzv24       ! dziv / 24.0

      real dxv48       ! dxiv / 48.0
      real dyv48       ! dyiv / 48.0
      real dzv48       ! dziv / 48.0

      real, intent(inout) :: hadv(0:ni+1,0:nj+1,1:nk)
                       ! Advection value in x and y direction

      real, intent(inout) :: vadv(0:ni+1,0:nj+1,1:nk)
                       ! Advection value in z direction

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_advuvw = 0
      integer, parameter :: DUMP_TARGET_advuvw = 360
      logical, save :: dump_done_advuvw = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpadvopt,advopt)
      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      dxv25n=-.25e0*dxiv
      dyv25n=-.25e0*dyiv
      dzv25n=-.25e0*dziv

      dxv24=oned24*dxiv
      dyv24=oned24*dyiv
      dzv24=oned24*dziv

      dxv48=.5e0*dxv24
      dyv48=.5e0*dyv24
      dzv48=.5e0*dzv24

! -----

!!! Calculate the velocity advection.

!@llm start meta_info ----------------------------------------------------
! Location: advuvw.f90 :: subroutine s_advuvw
! Summary : Calculates velocity (u, v, w) advection using 2nd or 4th order
!           centered finite difference schemes with mass-weighted fluxes.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constants (oned24, fourd3) from commath.
!   - No synchronization constructs.
!   - Multiple code paths based on advopt (1=2nd order, 2=2nd+4th, 3=separated).
!   - Multi-stage stencil with temporary arrays (tmp1, tmp2, tmp3, hadv, vadv).
!   - Three separate velocity components (u, v, w) computed sequentially.
!   - All loops are embarrassingly parallel within each stage.
! Next:
!   - Can use multiple OpenACC kernels matching the loop structure.
!   - Consider computing u, v, w advection in parallel if temporary arrays
!     are independent.
!   - Temporary arrays already allocated - good for GPU data management.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.6M
!   - TotalTime: 43.681s (1.47%)
!   - AvgTime: 121.336ms
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('advuvw.f90', 's_advuvw', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) &
     & * int((nj-2)-(2)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_advuvw = dump_call_count_advuvw + 1
if (dump_call_count_advuvw == DUMP_TARGET_advuvw .and. .not. dump_done_advuvw) then
  call dump_init('advuvw')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('advopt', advopt)
  call dump_scalar_i('iwest', iwest)
  call dump_scalar_i('ieast', ieast)
  call dump_scalar_i('jsouth', jsouth)
  call dump_scalar_i('jnorth', jnorth)
  call dump_scalar_r('dxv24', dxv24)
  call dump_scalar_r('dxv25n', dxv25n)
  call dump_scalar_r('dxv48', dxv48)
  call dump_scalar_r('dyv24', dyv24)
  call dump_scalar_r('dyv25n', dyv25n)
  call dump_scalar_r('dyv48', dyv48)
  call dump_scalar_r('dzv24', dzv24)
  call dump_scalar_r('dzv25n', dzv25n)
  call dump_scalar_r('dzv48', dzv48)
  call dump_array_3d('rstxu.bin', rstxu, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rstxv.bin', rstxv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rstxwc.bin', rstxwc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ufrc_in.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc_in.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wfrc_in.bin', wfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('hadv_in.bin', hadv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vadv_in.bin', vadv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_in.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp2_in.bin', tmp2, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp3_in.bin', tmp3, 0, ni+1, 0, nj+1, 1, nk)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_016)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 1, ni-1
          tmp1(i,j,k) = (rstxu(i,j,k) + rstxu(i+1,j,k)) * (u(i+1,j,k) - u(i,j,k)) * dxv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-1
        do i = 2, ni-1
          tmp2(i,j,k) = (rstxv(i-1,j,k) + rstxv(i,j,k)) * (u(i,j,k) - u(i,j-1,k)) * dyv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-2
        do i = 2, ni-1
          tmp3(i,j,k) = (rstxwc(i-1,j,k) + rstxwc(i,j,k)) * (u(i,j,k) - u(i,j,k-1)) * dzv25n
        end do
      end do
    end do
    !$acc end kernels

    if (advopt == 1) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2, nj-2
          do i = 2, ni-1
            ufrc(i,j,k) = ufrc(i,j,k) + ((tmp3(i,j,k) + tmp3(i,j,k+1)) &
                 + ((tmp1(i-1,j,k) + tmp1(i,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))))
          end do
        end do
      end do
      !$acc end kernels
    else
      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              hadv(i,j,k) = (tmp1(i-1,j,k) + tmp1(i,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              hadv(i,j,k) = (tmp1(i-1,j,k) + tmp1(i,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))
              vadv(i,j,k) = tmp3(i,j,k) + tmp3(i,j,k+1)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! 4th order u advection
    if (advopt == 2 .or. advopt == 3) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-2-jnorth
          do i = 1+iwest, ni-ieast
            tmp1(i,j,k) = (rstxu(i-1,j,k) + rstxu(i+1,j,k)) * (u(i+1,j,k) - u(i-1,j,k)) * dxv24
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 1+jsouth, nj-1-jnorth
          do i = 2+iwest, ni-1-ieast
            tmp2(i,j,k) = ((rstxv(i-1,j,k) + rstxv(i,j,k)) + (rstxv(i-1,j+1,k) + rstxv(i,j+1,k))) &
                 * (u(i,j+1,k) - u(i,j-1,k)) * dyv48
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-2-jnorth
          do i = 2+iwest, ni-1-ieast
            hadv(i,j,k) = fourd3 * hadv(i,j,k) &
                 + ((tmp1(i-1,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j+1,k)))
          end do
        end do
      end do
      !$acc end kernels

      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) + hadv(i,j,k) + (tmp3(i,j,k) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-1-ieast
              tmp3(i,j,k) = ((rstxwc(i-1,j,k) + rstxwc(i,j,k)) + (rstxwc(i-1,j,k+1) + rstxwc(i,j,k+1))) &
                   * (u(i,j,k+1) - u(i,j,k-1)) * dzv48
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 3, nk-3
          do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-1-ieast
              vadv(i,j,k) = fourd3 * vadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) + hadv(i,j,k) + vadv(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    !! Calculate the v advection (2nd order)
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-1
        do i = 2, ni-1
          tmp1(i,j,k) = (rstxu(i,j-1,k) + rstxu(i,j,k)) * (v(i,j,k) - v(i-1,j,k)) * dxv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 1, nj-1
        do i = 2, ni-2
          tmp2(i,j,k) = (rstxv(i,j,k) + rstxv(i,j+1,k)) * (v(i,j+1,k) - v(i,j,k)) * dyv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-1
        do i = 2, ni-2
          tmp3(i,j,k) = (rstxwc(i,j-1,k) + rstxwc(i,j,k)) * (v(i,j,k) - v(i,j,k-1)) * dzv25n
        end do
      end do
    end do
    !$acc end kernels

    if (advopt == 1) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2, nj-1
          do i = 2, ni-2
            vfrc(i,j,k) = vfrc(i,j,k) + ((tmp3(i,j,k) + tmp3(i,j,k+1)) &
                 + ((tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j,k))))
          end do
        end do
      end do
      !$acc end kernels
    else
      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j,k))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j,k))
              vadv(i,j,k) = tmp3(i,j,k) + tmp3(i,j,k+1)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! 4th order v advection
    if (advopt == 2 .or. advopt == 3) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-1-jnorth
          do i = 1+iwest, ni-1-ieast
            tmp1(i,j,k) = ((rstxu(i,j-1,k) + rstxu(i,j,k)) + (rstxu(i+1,j-1,k) + rstxu(i+1,j,k))) &
                 * (v(i+1,j,k) - v(i-1,j,k)) * dxv48
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 1+jsouth, nj-jnorth
          do i = 2+iwest, ni-2-ieast
            tmp2(i,j,k) = (rstxv(i,j-1,k) + rstxv(i,j+1,k)) * (v(i,j+1,k) - v(i,j-1,k)) * dyv24
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-1-jnorth
          do i = 2+iwest, ni-2-ieast
            hadv(i,j,k) = fourd3 * hadv(i,j,k) &
                 + ((tmp1(i-1,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j+1,k)))
          end do
        end do
      end do
      !$acc end kernels

      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) + hadv(i,j,k) + (tmp3(i,j,k) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2+jsouth, nj-1-jnorth
            do i = 2+iwest, ni-2-ieast
              tmp3(i,j,k) = ((rstxwc(i,j-1,k) + rstxwc(i,j,k)) + (rstxwc(i,j-1,k+1) + rstxwc(i,j,k+1))) &
                   * (v(i,j,k+1) - v(i,j,k-1)) * dzv48
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 3, nk-3
          do j = 2+jsouth, nj-1-jnorth
            do i = 2+iwest, ni-2-ieast
              vadv(i,j,k) = fourd3 * vadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) + hadv(i,j,k) + vadv(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    !! Calculate the w advection (2nd order)
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-2
        do i = 2, ni-1
          tmp1(i,j,k) = (rstxu(i,j,k-1) + rstxu(i,j,k)) * (w(i,j,k) - w(i-1,j,k)) * dxv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-1
        do i = 2, ni-2
          tmp2(i,j,k) = (rstxv(i,j,k-1) + rstxv(i,j,k)) * (w(i,j,k) - w(i,j-1,k)) * dyv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-1
      do j = 2, nj-2
        do i = 2, ni-2
          tmp3(i,j,k) = (rstxwc(i,j,k) + rstxwc(i,j,k+1)) * (w(i,j,k+1) - w(i,j,k)) * dzv25n
        end do
      end do
    end do
    !$acc end kernels

    if (advopt == 1) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-1
        do j = 2, nj-2
          do i = 2, ni-2
            wfrc(i,j,k) = wfrc(i,j,k) + ((tmp3(i,j,k-1) + tmp3(i,j,k)) &
                 + ((tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))))
          end do
        end do
      end do
      !$acc end kernels
    else
      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2, nj-2
            do i = 2, ni-2
              hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2, nj-2
            do i = 2, ni-2
              hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))
              vadv(i,j,k) = tmp3(i,j,k-1) + tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! 4th order w advection
    if (advopt == 2 .or. advopt == 3) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-1
        do j = 2+jsouth, nj-2-jnorth
          do i = 1+iwest, ni-1-ieast
            tmp1(i,j,k) = ((rstxu(i,j,k-1) + rstxu(i,j,k)) + (rstxu(i+1,j,k-1) + rstxu(i+1,j,k))) &
                 * (w(i+1,j,k) - w(i-1,j,k)) * dxv48
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-1
        do j = 1+jsouth, nj-1-jnorth
          do i = 2+iwest, ni-2-ieast
            tmp2(i,j,k) = ((rstxv(i,j,k-1) + rstxv(i,j,k)) + (rstxv(i,j+1,k-1) + rstxv(i,j+1,k))) &
                 * (w(i,j+1,k) - w(i,j-1,k)) * dyv48
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-1
        do j = 2+jsouth, nj-2-jnorth
          do i = 2+iwest, ni-2-ieast
            hadv(i,j,k) = fourd3 * hadv(i,j,k) &
                 + ((tmp1(i-1,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j+1,k)))
          end do
        end do
      end do
      !$acc end kernels

      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2, nj-2
            do i = 2, ni-2
              wfrc(i,j,k) = wfrc(i,j,k) + hadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-2-ieast
              tmp3(i,j,k) = (rstxwc(i,j,k-1) + rstxwc(i,j,k+1)) * (w(i,j,k+1) - w(i,j,k-1)) * dzv24
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 3, nk-2
          do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-2-ieast
              vadv(i,j,k) = fourd3 * vadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2, nj-2
            do i = 2, ni-2
              wfrc(i,j,k) = wfrc(i,j,k) + hadv(i,j,k) + vadv(i,j,k)
            end do
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

!! Calculate the u advection.

! Calculate the 2nd order u advection.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=1,ni-1
          tmp1(i,j,k)                                                   &
     &      =(rstxu(i,j,k)+rstxu(i+1,j,k))*(u(i+1,j,k)-u(i,j,k))*dxv25n
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-1
          tmp2(i,j,k)                                                   &
     &      =(rstxv(i-1,j,k)+rstxv(i,j,k))*(u(i,j,k)-u(i,j-1,k))*dyv25n
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          tmp3(i,j,k)=(rstxwc(i-1,j,k)+rstxwc(i,j,k))                   &
     &      *(u(i,j,k)-u(i,j,k-1))*dzv25n
        end do
        end do

!$omp end do

      end do

      if(advopt.eq.1) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            ufrc(i,j,k)=ufrc(i,j,k)+((tmp3(i,j,k)+tmp3(i,j,k+1))        &
     &        +((tmp1(i-1,j,k)+tmp1(i,j,k))                             &
     &        +(tmp2(i,j,k)+tmp2(i,j+1,k))))
          end do
          end do

!$omp end do

        end do

      else

        if(advopt.eq.2) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              hadv(i,j,k)                                               &
     &          =(tmp1(i-1,j,k)+tmp1(i,j,k))+(tmp2(i,j,k)+tmp2(i,j+1,k))
            end do
            end do

!$omp end do

          end do

        else if(advopt.eq.3) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              hadv(i,j,k)                                               &
     &          =(tmp1(i-1,j,k)+tmp1(i,j,k))+(tmp2(i,j,k)+tmp2(i,j+1,k))

              vadv(i,j,k)=tmp3(i,j,k)+tmp3(i,j,k+1)

            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Calculate the 4th order u advection.

      if(advopt.eq.2.or.advopt.eq.3) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=1+iwest,ni-ieast
            tmp1(i,j,k)=(rstxu(i-1,j,k)+rstxu(i+1,j,k))                 &
     &        *(u(i+1,j,k)-u(i-1,j,k))*dxv24
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=1+jsouth,nj-1-jnorth
          do i=2+iwest,ni-1-ieast
            tmp2(i,j,k)=((rstxv(i-1,j,k)+rstxv(i,j,k))                  &
     &        +(rstxv(i-1,j+1,k)+rstxv(i,j+1,k)))                       &
     &        *(u(i,j+1,k)-u(i,j-1,k))*dyv48
          end do
          end do

!$omp end do

        end do

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=2+iwest,ni-1-ieast
            hadv(i,j,k)=fourd3*hadv(i,j,k)                              &
     &        +((tmp1(i-1,j,k)+tmp1(i+1,j,k))                           &
     &        +(tmp2(i,j-1,k)+tmp2(i,j+1,k)))
          end do
          end do

!$omp end do

        end do

        if(advopt.eq.2) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)                                   &
     &          +hadv(i,j,k)+(tmp3(i,j,k)+tmp3(i,j,k+1))
            end do
            end do

!$omp end do

          end do

        else if(advopt.eq.3) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2+jsouth,nj-2-jnorth
            do i=2+iwest,ni-1-ieast
              tmp3(i,j,k)=((rstxwc(i-1,j,k)+rstxwc(i,j,k))              &
     &           +(rstxwc(i-1,j,k+1)+rstxwc(i,j,k+1)))                  &
     &           *(u(i,j,k+1)-u(i,j,k-1))*dzv48
            end do
            end do

!$omp end do

          end do

          do k=3,nk-3

!$omp do schedule(runtime) private(i,j)

            do j=2+jsouth,nj-2-jnorth
            do i=2+iwest,ni-1-ieast
              vadv(i,j,k)                                               &
     &          =fourd3*vadv(i,j,k)+(tmp3(i,j,k-1)+tmp3(i,j,k+1))
            end do
            end do

!$omp end do

          end do

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)+hadv(i,j,k)+vadv(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!! -----

!! Calculate the v advection.

! Calculate the 2nd order v advection.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-1
          tmp1(i,j,k)                                                   &
     &      =(rstxu(i,j-1,k)+rstxu(i,j,k))*(v(i,j,k)-v(i-1,j,k))*dxv25n
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=2,ni-2
          tmp2(i,j,k)                                                   &
     &      =(rstxv(i,j,k)+rstxv(i,j+1,k))*(v(i,j+1,k)-v(i,j,k))*dyv25n
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          tmp3(i,j,k)=(rstxwc(i,j-1,k)+rstxwc(i,j,k))                   &
     &      *(v(i,j,k)-v(i,j,k-1))*dzv25n
        end do
        end do

!$omp end do

      end do

      if(advopt.eq.1) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            vfrc(i,j,k)=vfrc(i,j,k)+((tmp3(i,j,k)+tmp3(i,j,k+1))        &
     &        +((tmp1(i,j,k)+tmp1(i+1,j,k))                             &
     &        +(tmp2(i,j-1,k)+tmp2(i,j,k))))
          end do
          end do

!$omp end do

        end do

      else

        if(advopt.eq.2) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              hadv(i,j,k)                                               &
     &          =(tmp1(i,j,k)+tmp1(i+1,j,k))+(tmp2(i,j-1,k)+tmp2(i,j,k))
            end do
            end do

!$omp end do

          end do

        else if(advopt.eq.3) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              hadv(i,j,k)                                               &
     &          =(tmp1(i,j,k)+tmp1(i+1,j,k))+(tmp2(i,j-1,k)+tmp2(i,j,k))

              vadv(i,j,k)=tmp3(i,j,k)+tmp3(i,j,k+1)

            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Calculate the 4th order v advection.

      if(advopt.eq.2.or.advopt.eq.3) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-1-jnorth
          do i=1+iwest,ni-1-ieast
            tmp1(i,j,k)=((rstxu(i,j-1,k)+rstxu(i,j,k))                  &
     &        +(rstxu(i+1,j-1,k)+rstxu(i+1,j,k)))                       &
     &        *(v(i+1,j,k)-v(i-1,j,k))*dxv48
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=1+jsouth,nj-jnorth
          do i=2+iwest,ni-2-ieast
            tmp2(i,j,k)=(rstxv(i,j-1,k)+rstxv(i,j+1,k))                 &
     &        *(v(i,j+1,k)-v(i,j-1,k))*dyv24
          end do
          end do

!$omp end do

        end do

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-1-jnorth
          do i=2+iwest,ni-2-ieast
            hadv(i,j,k)=fourd3*hadv(i,j,k)                              &
     &        +((tmp1(i-1,j,k)+tmp1(i+1,j,k))                           &
     &        +(tmp2(i,j-1,k)+tmp2(i,j+1,k)))
          end do
          end do

!$omp end do

        end do

        if(advopt.eq.2) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)                                   &
     &          +hadv(i,j,k)+(tmp3(i,j,k)+tmp3(i,j,k+1))
            end do
            end do

!$omp end do

          end do

        else if(advopt.eq.3) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2+jsouth,nj-1-jnorth
            do i=2+iwest,ni-2-ieast
              tmp3(i,j,k)=((rstxwc(i,j-1,k)+rstxwc(i,j,k))              &
     &          +(rstxwc(i,j-1,k+1)+rstxwc(i,j,k+1)))                   &
     &          *(v(i,j,k+1)-v(i,j,k-1))*dzv48
            end do
            end do

!$omp end do

          end do

          do k=3,nk-3

!$omp do schedule(runtime) private(i,j)

            do j=2+jsouth,nj-1-jnorth
            do i=2+iwest,ni-2-ieast
              vadv(i,j,k)                                               &
     &          =fourd3*vadv(i,j,k)+(tmp3(i,j,k-1)+tmp3(i,j,k+1))
            end do
            end do

!$omp end do

          end do

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)+hadv(i,j,k)+vadv(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!! -----

!! Calculate the w advection.

! Calculate the 2nd order w advection.

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          tmp1(i,j,k)                                                   &
     &      =(rstxu(i,j,k-1)+rstxu(i,j,k))*(w(i,j,k)-w(i-1,j,k))*dxv25n
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          tmp2(i,j,k)                                                   &
     &      =(rstxv(i,j,k-1)+rstxv(i,j,k))*(w(i,j,k)-w(i,j-1,k))*dyv25n
        end do
        end do

!$omp end do

      end do

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tmp3(i,j,k)=(rstxwc(i,j,k)+rstxwc(i,j,k+1))                   &
     &      *(w(i,j,k+1)-w(i,j,k))*dzv25n
        end do
        end do

!$omp end do

      end do

      if(advopt.eq.1) then

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wfrc(i,j,k)=wfrc(i,j,k)+((tmp3(i,j,k-1)+tmp3(i,j,k))        &
     &        +((tmp1(i,j,k)+tmp1(i+1,j,k))                             &
     &        +(tmp2(i,j,k)+tmp2(i,j+1,k))))
          end do
          end do

!$omp end do

        end do

      else

        if(advopt.eq.2) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              hadv(i,j,k)                                               &
     &          =(tmp1(i,j,k)+tmp1(i+1,j,k))+(tmp2(i,j,k)+tmp2(i,j+1,k))
            end do
            end do

!$omp end do

          end do

        else if(advopt.eq.3) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              hadv(i,j,k)                                               &
     &          =(tmp1(i,j,k)+tmp1(i+1,j,k))+(tmp2(i,j,k)+tmp2(i,j+1,k))

              vadv(i,j,k)=tmp3(i,j,k-1)+tmp3(i,j,k)

            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Calculate the 4th order w advection.

      if(advopt.eq.2.or.advopt.eq.3) then

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=1+iwest,ni-1-ieast
            tmp1(i,j,k)=((rstxu(i,j,k-1)+rstxu(i,j,k))                  &
     &        +(rstxu(i+1,j,k-1)+rstxu(i+1,j,k)))                       &
     &        *(w(i+1,j,k)-w(i-1,j,k))*dxv48
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=1+jsouth,nj-1-jnorth
          do i=2+iwest,ni-2-ieast
            tmp2(i,j,k)=((rstxv(i,j,k-1)+rstxv(i,j,k))                  &
     &        +(rstxv(i,j+1,k-1)+rstxv(i,j+1,k)))                       &
     &        *(w(i,j+1,k)-w(i,j-1,k))*dyv48
          end do
          end do

!$omp end do

        end do

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=2+iwest,ni-2-ieast
            hadv(i,j,k)=fourd3*hadv(i,j,k)                              &
     &        +((tmp1(i-1,j,k)+tmp1(i+1,j,k))                           &
     &        +(tmp2(i,j-1,k)+tmp2(i,j+1,k)))
          end do
          end do

!$omp end do

        end do

        if(advopt.eq.2) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              wfrc(i,j,k)=wfrc(i,j,k)                                   &
     &          +hadv(i,j,k)+(tmp3(i,j,k-1)+tmp3(i,j,k))
            end do
            end do

!$omp end do

          end do

        else if(advopt.eq.3) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2+jsouth,nj-2-jnorth
            do i=2+iwest,ni-2-ieast
              tmp3(i,j,k)=(rstxwc(i,j,k-1)+rstxwc(i,j,k+1))             &
     &          *(w(i,j,k+1)-w(i,j,k-1))*dzv24
            end do
            end do

!$omp end do

          end do

          do k=3,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2+jsouth,nj-2-jnorth
            do i=2+iwest,ni-2-ieast
              vadv(i,j,k)                                               &
     &          =fourd3*vadv(i,j,k)+(tmp3(i,j,k-1)+tmp3(i,j,k+1))
            end do
            end do

!$omp end do

          end do

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              wfrc(i,j,k)=wfrc(i,j,k)+hadv(i,j,k)+vadv(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!! -----

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_advuvw == DUMP_TARGET_advuvw .and. .not. dump_done_advuvw) then
  call dump_array_3d('ufrc_ref.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc_ref.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wfrc_ref.bin', wfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('hadv_ref.bin', hadv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vadv_ref.bin', vadv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_ref.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp2_ref.bin', tmp2, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp3_ref.bin', tmp3, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_advuvw = .true.
end if


call profile_stop(prof_id1, loop_len)

!!! ----

      end subroutine s_advuvw

!-----7--------------------------------------------------------------7--

      end module m_advuvw

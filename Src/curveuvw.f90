!***********************************************************************
      module m_curveuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/11/11
!     Modification: 2003/01/04, 2003/04/30, 2003/05/19, 2003/10/31,
!                   2003/11/28, 2003/12/12, 2004/04/15, 2004/06/10,
!                   2006/02/13, 2006/11/06, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/09, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the curvature of earth.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_dump_kernel
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: curveuvw, s_curveuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface curveuvw

        module procedure s_curveuvw

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
      subroutine s_curveuvw(fpmpopt,fpmfcopt,ni,nj,nk,rmf8u,rmf8v,      &
     &                      rst,u,v,w,ufrc,vfrc,wfrc,tmp1,tmp2,         &
     &                      tmp3,tmp4,tmp5)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

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

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      real rev125      ! 0.125 / rearth

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_curveuvw = 0
      integer, parameter :: DUMP_TARGET_curveuvw = 360
      logical, save :: dump_done_curveuvw = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)

! -----

! Set the common used variable.

      rev125=.125e0/rearth

! -----

!! Calculate the curvature of earth in the x and y components of
!! velocity equation.

!@llm start meta_info ----------------------------------------------------
! Location: curveuvw.f90 :: s_curveuvw
! Summary : Calculate earth curvature forcing terms for u, v, w velocity
!           equations using temporary arrays and map scale factors.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Private variable k used for outer loop
!   - Writes to ufrc, vfrc, wfrc (forcing terms) and tmp1-tmp5 (temporaries)
!   - Multiple sequential k-loops with data dependencies between them
!   - Conditional branches based on mpopt and mfcopt options
! Next:
!   - Convert to OpenACC with data region for all arrays
!   - May need to fuse some k-loops or restructure for better parallelism
!   - Handle conditional logic for map projection options on GPU
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 17.817s (0.60%)
!   - AvgTime: 49.492ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('curveuvw.f90', 's_curveuvw', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_curveuvw = dump_call_count_curveuvw + 1
if (dump_call_count_curveuvw == DUMP_TARGET_curveuvw .and. .not. dump_done_curveuvw) then
  call dump_init('curveuvw')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('mpopt', mpopt)
  call dump_scalar_i('mfcopt', mfcopt)
  call dump_scalar_r('rearth', rearth)
  call dump_array_3d('rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ufrc_in.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc_in.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wfrc_in.bin', wfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_in.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp2_in.bin', tmp2, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp3_in.bin', tmp3, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp4_in.bin', tmp4, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp5_in.bin', tmp5, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('rev125', rev125)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_073)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------

! Set the common used array - tmp1, tmp2
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          tmp1(i,j,k) = u(i,j,k) + u(i+1,j,k)
          tmp2(i,j,k) = v(i,j,k) + v(i,j+1,k)
        end do
      end do
    end do
    !$acc end kernels

! Set tmp3
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 1, nj-1
        do i = 1, ni-1
          tmp3(i,j,k) = rst(i,j,k) * (w(i,j,k) + w(i,j,k+1))
        end do
      end do
    end do
    !$acc end kernels

! For x and y components of velocity - tmp4, tmp5
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 1, ni-1
          tmp4(i,j,k) = tmp1(i,j,k) * tmp3(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 1, nj-1
        do i = 2, ni-2
          tmp5(i,j,k) = tmp2(i,j,k) * tmp3(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

! Update ufrc, vfrc
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 2, ni-1
          ufrc(i,j,k) = ufrc(i,j,k) - rev125 * (tmp4(i-1,j,k) + tmp4(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-1
        do i = 2, ni-2
          vfrc(i,j,k) = vfrc(i,j,k) - rev125 * (tmp5(i,j-1,k) + tmp5(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

! Calculate the curvature of earth in the z components of velocity equation
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-1
      do j = 2, nj-2
        do i = 2, ni-2
          tmp3(i,j,k) = rst(i,j,k) &
            * (tmp1(i,j,k) * tmp1(i,j,k) + tmp2(i,j,k) * tmp2(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-2
        do i = 2, ni-2
          wfrc(i,j,k) = wfrc(i,j,k) + rev125 * (tmp3(i,j,k-1) + tmp3(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

! Add the terms in the case of turning on the option for map scale factor
    if (mfcopt == 1) then

      if (mpopt == 0 .or. mpopt == 10) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 1, nj-1
            do i = 1, ni-1
              tmp3(i,j,k) = rmf8v(i,j,3) * rst(i,j,k) * tmp1(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else if (mpopt == 5) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 1, nj-1
            do i = 1, ni-1
              tmp3(i,j,k) = rmf8u(i,j,3) * rst(i,j,k) * tmp2(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 1, nj-1
            do i = 1, ni-1
              tmp3(i,j,k) = rst(i,j,k) &
                * (rmf8v(i,j,3) * tmp1(i,j,k) - rmf8u(i,j,3) * tmp2(i,j,k))
            end do
          end do
        end do
        !$acc end kernels
      end if

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 1, nj-1
          do i = 2, ni-2
            tmp1(i,j,k) = tmp1(i,j,k) * tmp3(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2, nj-2
          do i = 1, ni-1
            tmp2(i,j,k) = tmp2(i,j,k) * tmp3(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      if (mpopt == 5) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) - (tmp2(i-1,j,k) + tmp2(i,j,k))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) + (tmp1(i,j-1,k) + tmp1(i,j,k))
            end do
          end do
        end do
        !$acc end kernels

      else
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) + (tmp2(i-1,j,k) + tmp2(i,j,k))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) - (tmp1(i,j-1,k) + tmp1(i,j,k))
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

! Set the common used array.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          tmp1(i,j,k)=u(i,j,k)+u(i+1,j,k)
          tmp2(i,j,k)=v(i,j,k)+v(i,j+1,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          tmp3(i,j,k)=rst(i,j,k)*(w(i,j,k)+w(i,j,k+1))
        end do
        end do

!$omp end do

      end do

! -----

! For x and y components of velocity.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=1,ni-1
          tmp4(i,j,k)=tmp1(i,j,k)*tmp3(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=2,ni-2
          tmp5(i,j,k)=tmp2(i,j,k)*tmp3(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          ufrc(i,j,k)=ufrc(i,j,k)-rev125*(tmp4(i-1,j,k)+tmp4(i,j,k))
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          vfrc(i,j,k)=vfrc(i,j,k)-rev125*(tmp5(i,j-1,k)+tmp5(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the curvature of earth in the z components of velocity
! equation.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tmp3(i,j,k)=rst(i,j,k)                                        &
     &      *(tmp1(i,j,k)*tmp1(i,j,k)+tmp2(i,j,k)*tmp2(i,j,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wfrc(i,j,k)=wfrc(i,j,k)+rev125*(tmp3(i,j,k-1)+tmp3(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Add the terms in the case of turning on the option for map scale
! factor.

      if(mfcopt.eq.1) then

        if(mpopt.eq.0.or.mpopt.eq.10) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp3(i,j,k)=rmf8v(i,j,3)*rst(i,j,k)*tmp1(i,j,k)
            end do
            end do

!$omp end do

          end do

        else if(mpopt.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp3(i,j,k)=rmf8u(i,j,3)*rst(i,j,k)*tmp2(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp3(i,j,k)=rst(i,j,k)                                    &
     &          *(rmf8v(i,j,3)*tmp1(i,j,k)-rmf8u(i,j,3)*tmp2(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-2
            tmp1(i,j,k)=tmp1(i,j,k)*tmp3(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=1,ni-1
            tmp2(i,j,k)=tmp2(i,j,k)*tmp3(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(mpopt.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)-(tmp2(i-1,j,k)+tmp2(i,j,k))
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)+(tmp1(i,j-1,k)+tmp1(i,j,k))
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)+(tmp2(i-1,j,k)+tmp2(i,j,k))
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)-(tmp1(i,j-1,k)+tmp1(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

#endif

! Dump output data at target call
if (dump_call_count_curveuvw == DUMP_TARGET_curveuvw .and. .not. dump_done_curveuvw) then
  call dump_array_3d('ufrc_ref.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc_ref.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wfrc_ref.bin', wfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_ref.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp2_ref.bin', tmp2, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp3_ref.bin', tmp3, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp4_ref.bin', tmp4, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp5_ref.bin', tmp5, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_curveuvw = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_curveuvw

!-----7--------------------------------------------------------------7--

      end module m_curveuvw

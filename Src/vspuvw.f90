!***********************************************************************
      module m_vspuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/07/05, 1999/08/03, 1999/09/14, 1999/09/30,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2001/06/29,
!                   2001/11/20, 2002/04/02, 2002/08/15, 2003/04/30,
!                   2003/05/19, 2003/12/12, 2004/04/15, 2006/01/10,
!                   2006/09/21, 2006/11/27, 2007/05/07, 2007/07/30,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the vertical sponge damping for the velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_comprofile
      use m_dump_kernel
      use m_getiname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vspuvw, s_vspuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vspuvw

        module procedure s_vspuvw

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
      subroutine s_vspuvw(fpgpvvar,fpvspvar,fpvspopt,ksp0,gtinc,        &
     &                    ni,nj,nk,ubr,vbr,rst8u,rst8v,rst8w,up,vp,wp,  &
     &                    rbct,ugpv,utd,vgpv,vtd,wgpv,wtd,              &
     &                    ufrc,vfrc,wfrc,rbct8s)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpvspvar
                       ! Formal parameter of unique index of vspvar

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) vspvar
                       ! Control flag of
                       ! vertical sponge damped variables

      integer vspopt   ! Option for vertical sponge damping

      real, intent(inout) :: rbct8s(1:ni,1:nj,1:nk)
                       ! 0.5 x relaxed top sponge damping coefficients
                       ! at scalar points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_vspuvw = 0
      integer, parameter :: DUMP_TARGET_vspuvw = 360
      logical, save :: dump_done_vspuvw = .false.


!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(vspvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fpvspvar,vspvar)
      call getiname(fpvspopt,vspopt)

! -----

!!! Calculate the vertical sponge damping for the velocity.

!@llm start meta_info ----------------------------------------------------
! Location: vspuvw.f90 :: subroutine s_vspuvw
! Summary : Calculates vertical sponge damping for u, v, w velocity
!           components near model top to absorb outgoing waves.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - No writes to module/global variables.
!   - No synchronization constructs.
!   - Multiple code paths based on vspvar and vspopt flags.
!   - Loops start from ksp0 (sponge layer index) to top.
!   - All grid points are independent within each loop.
!   - Damping to GPV data or base state based on vspopt.
! Next:
!   - Direct OpenACC kernels for each component.
!   - Evaluate conditions outside kernel to select code path.
!   - Consider fusing loops for same component if beneficial.
! Runtime:
!   - Calls: 360
!   - AvgLoops: 101.6M
!   - TotalTime: 8.704s (0.29%)
!   - AvgTime: 24.177ms
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('vspuvw.f90', 's_vspuvw', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(ksp0(1)-1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_vspuvw = dump_call_count_vspuvw + 1
if (dump_call_count_vspuvw == DUMP_TARGET_vspuvw .and. .not. dump_done_vspuvw) then
  call dump_init('vspuvw')
  call dump_scalar_i('fpgpvvar', fpgpvvar)
  call dump_scalar_i('fpvspvar', fpvspvar)
  call dump_scalar_i('fpvspopt', fpvspopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('gtinc', gtinc)
  call dump_array_3d('ubr.bin', ubr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vbr.bin', vbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rst8u.bin', rst8u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rst8v.bin', rst8v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rst8w.bin', rst8w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('up.bin', up, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vp.bin', vp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wp.bin', wp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_4d('rbct.bin', rbct, 1, ni, 1, nj, 1, nk, 1, 2)
  call dump_array_3d('ugpv.bin', ugpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('utd.bin', utd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vgpv.bin', vgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vtd.bin', vtd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wgpv.bin', wgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wtd.bin', wtd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ufrc_in.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc_in.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wfrc_in.bin', wfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rbct8s_in.bin', rbct8s, 1, ni, 1, nj, 1, nk)
end if

!$omp parallel default(shared) private(k)

! Set the common used variable.

      if(vspvar(1:1).eq.'o'.or.vspvar(2:2).eq.'o') then

        if(vspopt.eq.1) then

          do k=ksp0(1)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              rbct8s(i,j,k)=.25e0*(rbct(i,j,k,1)+rbct(i,j,k+1,1))
            end do
            end do

!$omp end do

          end do

        else

          do k=ksp0(2)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              rbct8s(i,j,k)=.25e0*(rbct(i,j,k,2)+rbct(i,j,k+1,2))
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!! For the x components of velocity.

      if(vspvar(1:1).eq.'o') then

! Damp to the GPV data.

        if(vspopt.eq.1) then

          do k=ksp0(1)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)                                   &
     &          -rst8u(i,j,k)*(rbct8s(i-1,j,k)+rbct8s(i,j,k))           &
     &          *(up(i,j,k)-(ugpv(i,j,k)+utd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the base state value.

        else

          do k=ksp0(2)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)-rst8u(i,j,k)                      &
     &          *(rbct8s(i-1,j,k)+rbct8s(i,j,k))*(up(i,j,k)-ubr(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!! For the y components of velocity.

      if(vspvar(2:2).eq.'o') then

! Damp to the GPV data.

        if(vspopt.eq.1) then

          do k=ksp0(1)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)                                   &
     &          -rst8v(i,j,k)*(rbct8s(i,j-1,k)+rbct8s(i,j,k))           &
     &          *(vp(i,j,k)-(vgpv(i,j,k)+vtd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the base state value.

        else

          do k=ksp0(2)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)-rst8v(i,j,k)                      &
     &          *(rbct8s(i,j-1,k)+rbct8s(i,j,k))*(vp(i,j,k)-vbr(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!! For the z components of velocity.

      if(vspvar(3:3).eq.'o') then

! Damp to the GPV data.

        if(vspopt.eq.1.and.gpvvar(1:1).eq.'o') then

          do k=ksp0(1),nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              wfrc(i,j,k)=wfrc(i,j,k)-rbct(i,j,k,1)*rst8w(i,j,k)        &
     &          *(wp(i,j,k)-(wgpv(i,j,k)+wtd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the 0.

        else

          do k=ksp0(2),nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              wfrc(i,j,k)=wfrc(i,j,k)                                   &
     &          -rbct(i,j,k,2)*rst8w(i,j,k)*wp(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

!! -----

      end if

!$omp end parallel

! Dump output data at target call
if (dump_call_count_vspuvw == DUMP_TARGET_vspuvw .and. .not. dump_done_vspuvw) then
  call dump_array_3d('ufrc_ref.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc_ref.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wfrc_ref.bin', wfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rbct8s_ref.bin', rbct8s, 1, ni, 1, nj, 1, nk)
  call dump_finalize()
  dump_done_vspuvw = .true.
end if


call profile_stop(prof_id1, loop_len)

!!! -----

      end subroutine s_vspuvw

!-----7--------------------------------------------------------------7--

      end module m_vspuvw

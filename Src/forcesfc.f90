!***********************************************************************
      module m_forcesfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/03/13
!     Modification: 2003/04/30, 2003/05/19, 2004/02/01, 2004/03/05,
!                   2004/09/10, 2005/06/10, 2006/01/10, 2007/01/20,
!                   2007/05/21, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/09/22, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the surface flux to bottom boundary.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forcesfc, s_forcesfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forcesfc

        module procedure s_forcesfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_forcesfc(fmois,ni,nj,nk,j31,j32,ptbr,u,v,w,ptp,qv,   &
     &                      ptv,qvsfc,ce,ct,cq,ufrc,vfrc,ptfrc,qvfrc)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

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

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(in) :: qvsfc(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio on surface

      real, intent(in) :: ce(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface momentum flux

      real, intent(in) :: ct(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface heat flux

      real, intent(in) :: cq(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface moisture flux

! Output variables

      real, intent(out) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(out) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(out) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in potential temperature equation

      real, intent(out) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real j318u       ! j31 at u point
      real j328v       ! j32 at v point

      real xcomp       ! Temporary variable
      real ycomp       ! Temporary variable
      real zcomp       ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_forcesfc = 0
      integer, parameter :: DUMP_TARGET_forcesfc = 361
      logical, save :: dump_done_forcesfc = .false.


!-----7--------------------------------------------------------------7--

!! Get the surface flux to bottom boundary.

!@llm start meta_info ----------------------------------------------------
! Location: forcesfc.f90 :: s_forcesfc
! Summary : Compute surface flux forcing terms for potential temperature,
!           water vapor, and velocity components at bottom boundary.
! GPU diff: Easy
! Findings:
!   - Multiple omp do regions for ptfrc, qvfrc, ufrc, vfrc calculations
!   - Conditional branches based on fmois (dry/moist) flag
!   - Uses sqrt intrinsic function for velocity calculations
!   - All operations are on 2D surface layer (k=1 or k=2)
!   - No reductions or synchronization
! Next:
!   - Port as 2D GPU kernels for surface layer
!   - Handle dry/moist branching with separate kernels or compile-time flag
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.050s (0.00%)
!   - AvgTime: 0.139ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('forcesfc.f90', 's_forcesfc', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_forcesfc = dump_call_count_forcesfc + 1
if (dump_call_count_forcesfc == DUMP_TARGET_forcesfc .and. .not. dump_done_forcesfc) then
  call dump_init('forcesfc')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('epsav', epsav)
  call dump_array_3d('j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptv.bin', ptv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('qvsfc.bin', qvsfc, 0, ni+1, 0, nj+1)
  call dump_array_2d('ce.bin', ce, 0, ni+1, 0, nj+1)
  call dump_array_2d('ct.bin', ct, 0, ni+1, 0, nj+1)
  call dump_array_2d('cq.bin', cq, 0, ni+1, 0, nj+1)
end if

!$omp parallel default(shared)

! Get the surface flux for the potential tempeture.

      if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ptfrc(i,j,1)=ct(i,j)*(ptv(i,j,2)-ptv(i,j,1))
        end do
        end do

!$omp end do

      else if(fmois(1:5).eq.'moist') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ptfrc(i,j,1)=ct(i,j)*((ptbr(i,j,2)+ptp(i,j,2))                &
     &      -ptv(i,j,1)*(1.e0+qvsfc(i,j))/(1.e0+epsav*qvsfc(i,j)))
        end do
        end do

!$omp end do

      end if

! -----

! Get the surface flux for water vapor mixing ratio.

      if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qvfrc(i,j,1)=0.e0
        end do
        end do

!$omp end do

      else if(fmois(1:5).eq.'moist') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qvfrc(i,j,1)=cq(i,j)*(qv(i,j,2)-qvsfc(i,j))
        end do
        end do

!$omp end do

      end if

! -----

! Get the surface flux for the x components of velocity.

!$omp do schedule(runtime) private(i,j,j318u,xcomp,zcomp)

      do j=1,nj-1
      do i=2,ni-1
        j318u=j31(i,j,2)+j31(i,j,3)

        xcomp=1.e0/sqrt(4.e0+j318u*j318u)
        zcomp=.125e0*j318u*xcomp

        ufrc(i,j,1)=(ce(i-1,j)+ce(i,j))*(u(i,j,2)*xcomp                 &
     &    +((w(i-1,j,2)+w(i,j,3))+(w(i-1,j,3)+w(i,j,2)))*zcomp)

      end do
      end do

!$omp end do

! -----

! Get the surface flux for the x components of velocity.

!$omp do schedule(runtime) private(i,j,j328v,ycomp,zcomp)

      do j=2,nj-1
      do i=1,ni-1
        j328v=j32(i,j,2)+j32(i,j,3)

        ycomp=1.e0/sqrt(4.e0+j328v*j328v)
        zcomp=.125e0*j328v*ycomp

        vfrc(i,j,1)=(ce(i,j-1)+ce(i,j))*(v(i,j,2)*ycomp                 &
     &    +((w(i,j-1,2)+w(i,j,3))+(w(i,j-1,3)+w(i,j,2)))*zcomp)

      end do
      end do

!$omp end do

! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_forcesfc == DUMP_TARGET_forcesfc .and. .not. dump_done_forcesfc) then
  call dump_array_3d('ufrc_ref.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc_ref.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptfrc_ref.bin', ptfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvfrc_ref.bin', qvfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_forcesfc = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_forcesfc

!-----7--------------------------------------------------------------7--

      end module m_forcesfc

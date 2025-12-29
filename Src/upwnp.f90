!***********************************************************************
      module m_upwnp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/09/30, 2007/10/19, 2007/11/26, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2009/11/05, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the sedimentation for optional precipitation
!     concentrations.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname
      use m_comprofile

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: upwnp, s_upwnp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface upwnp

        module procedure s_upwnp

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_upwnp(fpdziv,dtp,ni,nj,nk,rbr,rst,un,ncf,ncflx)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtp
                       ! Time steps interval of fall out integration

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: un(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity
                       ! of optional precipitation concentrations

! Input and output variable

      real, intent(inout) :: ncf(0:ni+1,0:nj+1,1:nk)
                       ! Optional precipitation concentrations at future

! Internal shared variables

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      real dziv        ! Inverse of dz

      real dzvdt       ! dziv x dtp

      real, intent(inout) :: ncflx(0:ni+1,0:nj+1,1:nk)
                       ! Fallout flux of
                       ! optional precipitation concentrations

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

      dzvdt=dziv*dtp

! -----

! Calculate the sedimentation.

!@llm start meta_info ----------------------------------------------------
! Location: upwnp.f90 :: s_upwnp
! Summary : Calculate sedimentation for optional precipitation
!           concentrations using upwind flux divergence scheme
! GPU diff: Easy
! Findings:
!   - Serial k-loop wrapping parallel i,j loops (private(k))
!   - Three sequential loop nests: flux calc, update, boundary copy
!   - Contains max() intrinsic for non-negative concentration
!   - Simple arithmetic operations, no function calls
! Next:
!   - Collapse loops or use OpenACC kernels with loop directive
!   - Can potentially fuse kernels for better performance
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('upwnp.f90', 's_upwnp', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ncflx(i,j,k)=rbr(i,j,k)*un(i,j,k)*ncf(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=1,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ncf(i,j,k)=max(ncf(i,j,k)                                     &
     &      +(ncflx(i,j,k+1)-ncflx(i,j,k))/rst(i,j,k)*dzvdt,0.e0)
        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        ncf(i,j,nkm1)=ncf(i,j,nkm2)
      end do
      end do

!$omp end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_upwnp

!-----7--------------------------------------------------------------7--

      end module m_upwnp

!***********************************************************************
      module m_diverpiv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/12/17
!     Modification: 1999/12/20, 2000/01/17, 2001/06/29, 2002/04/02,
!                   2003/04/30, 2003/05/19, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the divergence vertically in the pressure equation
!     with the horizontally explicit and vertically implicit method.

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

      public :: diverpiv, s_diverpiv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diverpiv

        module procedure s_diverpiv

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
      subroutine s_diverpiv(fpdziv,ni,nj,nk,rcsq,w,pdiv)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: ni
                       ! Model dimemsion in x direction

      integer, intent(in) :: nj
                       ! Model dimemsion in y direction

      integer, intent(in) :: nk
                       ! Model dimemsion in z direction

      real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Output variable

      real, intent(out) :: pdiv(0:ni+1,0:nj+1,1:nk)
                       ! Vertical divergence value in pressure equation

! Internal shared variable

      real dziv        ! Inverse of dz

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

! Calculate the divergence vertically in the pressure equation.

!@llm start meta_info ----------------------------------------------------
! Location: diverpiv.f90 :: s_diverpiv
! Summary : Calculate vertical divergence for pressure equation (HEVI method),
!           computing rcsq * (w(k) - w(k+1)) * dziv at each grid point.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - No global/module variable writes
!   - No synchronization constructs
!   - Simple 3D loop with straightforward vertical differencing
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - Very simple kernel, good candidate for early GPU porting
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('diverpiv.f90', 's_diverpiv', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) &
     & * int((nj-2)-(2)+1,8) &
     & * int((ni-2)-(2)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          pdiv(i,j,k)=rcsq(i,j,k)*(w(i,j,k)-w(i,j,k+1))*dziv
        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_diverpiv

!-----7--------------------------------------------------------------7--

      end module m_diverpiv

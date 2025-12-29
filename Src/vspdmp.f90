!***********************************************************************
      module m_vspdmp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/07/05, 1999/08/03, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2000/04/18, 2000/07/05,
!                   2001/04/15, 2001/05/29, 2001/11/20, 2001/12/11,
!                   2002/04/02, 2002/06/18, 2002/09/09, 2003/04/30,
!                   2003/05/19, 2003/10/10, 2003/12/12, 2004/04/15,
!                   2004/08/31, 2004/09/10, 2006/11/27, 2007/01/31,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the relaxed vertical sponge damping coefficients.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vspdmp, s_vspdmp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vspdmp

        module procedure s_vspdmp

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_vspdmp(fpvspopt,fpvspgpv,fpvspbar,fpbotgpv,fpbotbar, &
     &                    ksp0,ni,nj,nk,zph,rbct,z1dmax)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpvspgpv
                       ! Formal parameter of unique index of vspgpv

      integer, intent(in) :: fpvspbar
                       ! Formal parameter of unique index of vspbar

      integer, intent(in) :: fpbotgpv
                       ! Formal parameter of unique index of botgpv

      integer, intent(in) :: fpbotbar
                       ! Formal parameter of unique index of botbar

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Output variables

      integer, intent(out) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      real, intent(out) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

! Internal shared variables

      integer vspopt   ! Option for vertical sponge damping

      integer nkm1     ! nk - 1

      real vspgpv      ! Vertical sponge damping coefficient
                       ! for external data
      real vspbar      ! Vertical sponge damping coefficient
                       ! for base state or 0

      real botgpv      ! Lowest height of vertical sponge damping
                       ! for external data
      real botbar      ! Lowest height of vertical sponge damping
                       ! for base state or 0

      real cgpv05      ! 0.5 x vspgpv
      real cbar05      ! 0.5 x vspbar

      real, intent(inout) :: z1dmax(1:nk)
                       ! Maximum z physical coordinates at each plane

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpvspopt,vspopt)
      call getrname(fpvspgpv,vspgpv)
      call getrname(fpvspbar,vspbar)
      call getrname(fpbotgpv,botgpv)
      call getrname(fpbotbar,botbar)

! -----

! Set the common used variables.

      nkm1=nk-1

      cgpv05=.5e0*vspgpv
      cbar05=.5e0*vspbar

! -----

! Initialize the processed variables.

      ksp0(1)=2
      ksp0(2)=2

! -----

!! Calculate the relaxed vertical sponge damping coefficients.

!@llm start meta_info ----------------------------------------------------
! Location: vspdmp.f90 :: s_vspdmp
! Summary : Calculates relaxed vertical sponge damping coefficients with
!           maximum z-coordinate search and cosine-based damping profiles.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Uses max() intrinsic for reduction-like operation on z1dmax
!   - !$omp single block for sequential ksp0 index search (do_k_1, do_k_2)
!   - Multiple k loops with different purposes (init, max-find, coef-calc)
!   - Conditional vspopt branches inside parallel region
!   - Potential race condition in z1dmax(k)=max(...) without proper reduction
! Next:
!   - z1dmax computation needs reduction or atomic operations for GPU
!   - Sequential ksp0 search should remain on CPU or use parallel reduction
!   - Split into separate kernels: max-find, ksp0-search, coef-calculation
!   - Use cosine from device math library
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('vspdmp.f90', 's_vspdmp', &
   & 'OMP section 1')
end if
loop_len = int((nk)-(3)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

! Set the maximum z physical coordinates at each plane.

!$omp do schedule(runtime)

      do k=3,nk
        z1dmax(k)=lim36n
      end do

!$omp end do

      do k=3,nk

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          z1dmax(k)=max(zph(i,j,k),z1dmax(k))
        end do
        end do

!$omp end do

      end do

! -----

! Get the lowest damping level.

!$omp single private(k)

      if(vspopt.eq.1) then

        do_k_1: do k=3,nk

          if(z1dmax(k).gt.botgpv) then

            ksp0(1)=k-1

            exit do_k_1

          end if

        end do do_k_1

      end if

      do_k_2: do k=3,nk

        if(z1dmax(k).gt.botbar) then

          ksp0(2)=k-1

          exit do_k_2

        end if

      end do do_k_2

!$omp end single

! -----

! Finally get the relaxed vertical sponge damping coefficients.

      if(vspopt.eq.1) then

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(zph(i,j,k).gt.botgpv) then

              rbct(i,j,k,1)=cgpv05*(1.e0-cos(cc*(zph(i,j,k)-botgpv)     &
     &          /(zph(i,j,nkm1)-botgpv)))

            else

              rbct(i,j,k,1)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

      end if

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(zph(i,j,k).gt.botbar) then

            rbct(i,j,k,2)=cbar05                                        &
     &        *(1.e0-cos(cc*(zph(i,j,k)-botbar)/(zph(i,j,nkm1)-botbar)))

          else

            rbct(i,j,k,2)=0.e0

          end if

        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_vspdmp

!-----7--------------------------------------------------------------7--

      end module m_vspdmp

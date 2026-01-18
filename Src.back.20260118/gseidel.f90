!***********************************************************************
      module m_gseidel
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/06/01
!     Modification: 2000/07/05, 2001/04/15, 2001/05/29, 2001/10/18,
!                   2001/11/20, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2003/12/12, 2004/03/05, 2004/09/01,
!                   2006/01/10, 2007/01/20, 2007/01/31, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the tridiagonal equation with the Gauss-Seidel method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkitr
      use m_comprofile
      use m_commath
      use m_copy2d
      use m_copy3d
      use m_getrname
      use m_setcst2d
      use m_setcst3d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: gseidel, s_gseidel

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface gseidel

        module procedure s_gseidel

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_gseidel(fpgsdeps,ni,nj,nk,rr,ss,tt,bb,ff,nr,nrp,dnr)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgsdeps
                       ! Formal parameter of unique index of gsdeps

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rr(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(in) :: ss(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(in) :: tt(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

! Input and output variable

      real, intent(inout) :: bb(0:ni+1,0:nj+1,1:nk)
                       ! Right hand and solved vector

! Internal shared variables

      integer nkm2     ! nk - 2
      integer nkm3     ! nk - 3

      real gsdeps      ! Value of convergence of iteration

      real itcon       ! Control flag of continuation of iteration

      real, intent(inout) :: ff(0:ni+1,0:nj+1,1:nk)
                       ! Solved vector

      real, intent(inout) :: nr(0:ni+1,0:nj+1)
                       ! Norm of ff at current

      real, intent(inout) :: nrp(0:ni+1,0:nj+1)
                       ! Norm of ff at past

      real, intent(inout) :: dnr(0:ni+1,0:nj+1)
                       ! Differencial of norm

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpgsdeps,gsdeps)

! -----

! Fill in the array with constant value.

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ff)

      call setcst2d(0,ni+1,0,nj+1,eps,nr)
      call setcst2d(0,ni+1,0,nj+1,eps,nrp)

      call setcst2d(0,ni+1,0,nj+1,1.e0,dnr)

! -----

! Set the common used variables.

      nkm2=nk-2
      nkm3=nk-3

! -----

!! Solve the tridiagonal equation with the Gauss-Seidel method.

      iterate: do

! Perform the Gauss-Seidel method.

!@llm start meta_info ----------------------------------------------------
! Location: gseidel.f90 :: s_gseidel
! Summary : Solve tridiagonal system using Gauss-Seidel iterative method
!           with convergence checking at each column.
! GPU diff: Hard
! Findings:
!   - No omp_get_thread usage
!   - No external function calls within parallel region
!   - Uses intrinsic abs function (GPU compatible)
!   - Sequential k-dependency in tridiagonal solve (ff(k) depends on ff(k-1))
!   - Conditional execution based on dnr(i,j) > gsdeps (divergent branches)
!   - Multiple omp do regions with implicit barriers between them
!   - Accumulation into nr array (potential race if k loop were parallelized)
!   - Iterative algorithm with external convergence check (chkitr)
! Next:
!   - Gauss-Seidel has inherent sequential dependency in k-direction
!   - Consider switching to Thomas algorithm (direct solve) for GPU
!   - Or use parallel cyclic reduction / PCR algorithm
!   - Column-wise parallelism (i,j) is safe but k must remain sequential
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('gseidel.f90', 's_gseidel', &
   & 'OMP section 1')
end if
loop_len = int((nj-2)-(2)+1,8) * int((ni-2)-(2)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2

          if(dnr(i,j).gt.gsdeps) then

            ff(i,j,3)=(bb(i,j,3)-tt(i,j,3)*ff(i,j,4))/ss(i,j,3)

          end if

        end do
        end do

!$omp end do

        do k=4,nk-3

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2

            if(dnr(i,j).gt.gsdeps) then

              ff(i,j,k)=(bb(i,j,k)                                      &
     &          -rr(i,j,k)*ff(i,j,k-1)-tt(i,j,k)*ff(i,j,k+1))/ss(i,j,k)

            end if

          end do
          end do

!$omp end do

        end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2

          if(dnr(i,j).gt.gsdeps) then

            ff(i,j,nkm2)                                                &
     &        =(bb(i,j,nkm2)-rr(i,j,nkm2)*ff(i,j,nkm3))/ss(i,j,nkm2)

          end if

        end do
        end do

!$omp end do

        do k=3,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2

            if(dnr(i,j).gt.gsdeps) then

              nr(i,j)=nr(i,j)+ff(i,j,k)*ff(i,j,k)

            end if

          end do
          end do

!$omp end do

        end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2

          if(dnr(i,j).gt.gsdeps) then

            dnr(i,j)=abs(nr(i,j)/nrp(i,j)-1.e0)

          end if

        end do
        end do

!$omp end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

! Check the continuity of the iteration.

        call s_chkitr('single',2,ni-2,2,nj-2,1,1,itcon,ni,nj,1,dnr)

        if(itcon.gt.gsdeps) then

          call copy2d(0,ni+1,0,nj+1,nr,nrp)

          call setcst2d(0,ni+1,0,nj+1,eps,nr)

        else

          exit iterate

        end if

! -----

      end do iterate

!! -----

! Copy the array ff to the bb.

      call copy3d(0,ni+1,0,nj+1,1,nk,ff,bb)

! -----

      end subroutine s_gseidel

!-----7--------------------------------------------------------------7--

      end module m_gseidel

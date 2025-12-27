!***********************************************************************
      module m_getxy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/04/06
!     Modification: 1999/05/10, 1999/07/05, 1999/07/23, 1999/09/30,
!                   2000/01/17, 2001/05/29, 2002/04/02, 2002/07/03,
!                   2003/04/30, 2003/05/19, 2003/09/01, 2006/09/21,
!                   2006/12/04, 2007/01/05, 2007/01/31, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the x and the y coordinates.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getxy, s_getxy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getxy

        module procedure s_getxy

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getxy(fpdx,fpdy,xo,imin,imax,jmin,jmax,x,y)
!***********************************************************************

! Input variables

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

! Output variables

      real, intent(out) :: x(imin:imax)
                       ! x coordinates

      real, intent(out) :: y(jmin:jmax)
                       ! y coordinates

! Internal shared variables

      integer ies2     ! Start index in entire domain in x d. - 2
      integer jes2     ! Start index in entire domain in y d. - 2

      integer ies23    ! 2 x start index in entire domain in x d. - 3
      integer jes23    ! 2 x start index in entire domain in y d. - 3

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpdx,dx)
      call getrname(fpdy,dy)

! -----

! Set the common used variables.

      if(imin.eq.0) then

        ies2=(imax-4)*(nisub*igrp+isub)-2
        ies23=2*(imax-4)*(nisub*igrp+isub)-3

      else

        ies2=-2
        ies23=-3

      end if

      if(jmin.eq.0) then

        jes2=(jmax-4)*(njsub*jgrp+jsub)-2
        jes23=2*(jmax-4)*(njsub*jgrp+jsub)-3

      else

        jes2=-2
        jes23=-3

      end if

! -----

!! Calculate the x and the y coordinates.

!@llm start meta_info ----------------------------------------------------
! Location: getxy.f90 :: s_getxy
! Summary : Calculate x and y coordinates for different grid staggering
!           (scalar, u, v, w points) based on grid spacing and MPI domain.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Uses intrinsic real function (GPU compatible)
!   - Multiple conditional branches for different stagger types
!   - Separate 1D loops for x and y arrays
!   - No global writes, only output arrays x and y are modified
! Next:
!   - Direct translation to OpenMP target with teams distribute
!   - 1D arrays are small, consider keeping on CPU or async transfer
!   - Separate kernels for x and y may be more efficient
!@llm end meta_info ------------------------------------------------------
!$omp parallel default(shared)

! Calculate the x and the y coordinates at the data grid points.

      if(xo(1:2).eq.'oo') then

!$omp do schedule(runtime) private(i)

        do i=imin,imax
          x(i)=real(ies2+i)*dx
        end do

!$omp end do

!$omp do schedule(runtime) private(j)

        do j=jmin,jmax
          y(j)=real(jes2+j)*dy
        end do

!$omp end do

! -----

! Calculate the x and the y coordinates at the u points.

      else if(xo(1:2).eq.'ox') then

!$omp do schedule(runtime) private(i)

        do i=imin,imax
          x(i)=real(ies2+i)*dx
        end do

!$omp end do

!$omp do schedule(runtime) private(j)

        do j=jmin,jmax-1
          y(j)=.5e0*real(jes23+2*j)*dy
        end do

!$omp end do

! -----

! Calculate the x and the y coordinates at the v points.

      else if(xo(1:2).eq.'xo') then

!$omp do schedule(runtime) private(i)

        do i=imin,imax-1
          x(i)=.5e0*real(ies23+2*i)*dx
        end do

!$omp end do

!$omp do schedule(runtime) private(j)

        do j=jmin,jmax
          y(j)=real(jes2+j)*dy
        end do

!$omp end do

! -----

! Calculate the x and the y coordinates at the w and the scalar points.

      else

!$omp do schedule(runtime) private(i)

        do i=imin,imax-1
          x(i)=.5e0*real(ies23+2*i)*dx
        end do

!$omp end do

!$omp do schedule(runtime) private(j)

        do j=jmin,jmax-1
          y(j)=.5e0*real(jes23+2*j)*dy
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_getxy

!-----7--------------------------------------------------------------7--

      end module m_getxy

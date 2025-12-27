!***********************************************************************
      module m_cg2mkbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/03/29, 2007/10/19, 2008/01/11, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/18, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     change measurement from [c] [g] to [m] [k] for warm bin cloud
!     physics.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: cg2mkbw, s_cg2mkbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface cg2mkbw

        module procedure s_cg2mkbw

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
      subroutine s_cg2mkbw(ni,nj,nk,nqw,nnw,rbv,mwbin,nwbin,prr)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density [cm^3/g]

! Input and output variables

      real, intent(inout) :: mwbin(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Total water mass [g/cm^3]

      real, intent(inout) :: nwbin(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations [1/cm^3]

      real, intent(inout) :: prr(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for rain
                       ! [cm/s], [cm]

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in bin categories

!-----7--------------------------------------------------------------7--

! Change measurement from [c] [g] to [m] [k].

!@llm start meta_info ----------------------------------------------------
! Location: cg2mkbw.f90 :: s_cg2mkbw
! Summary : Converts units from CGS [c][g] to MKS [m][k] for warm bin
!           cloud microphysics arrays (mass, concentration, precipitation).
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple scalar multiplication operations
!   - Nested loops over bin categories (n), vertical levels (k), and grid (i,j)
!   - Independent updates to mwbin, nwbin, prr arrays
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Convert to OpenMP target with data mapping for mwbin, nwbin, prr, rbv
!   - Use collapse for nested loops (consider collapse(3) or collapse(4))
!   - May need to handle n-loop separately if nqw varies at runtime
!   - Simple structure well-suited for GPU offload
!@llm end meta_info ------------------------------------------------------
!$omp parallel default(shared) private(k,n)

      do n=1,nqw

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            mwbin(i,j,k,n)=rbv(i,j,k)*mwbin(i,j,k,n)
            nwbin(i,j,k,n)=1.e3*rbv(i,j,k)*nwbin(i,j,k,n)
          end do
          end do

!$omp end do

        end do

      end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        prr(i,j,1)=1.e-2*prr(i,j,1)
        prr(i,j,2)=1.e-2*prr(i,j,2)
      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_cg2mkbw

!-----7--------------------------------------------------------------7--

      end module m_cg2mkbw

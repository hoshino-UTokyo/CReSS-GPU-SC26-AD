!***********************************************************************
      module m_specnudm
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2016/04/07
!     Modification: 2016/04/08, 2016/11/11, 2017/01/07, 2017/01/29,
!                   2017/02/01, 2017/06/09

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform spectral nudging.
!     1. Gather data from each process to root process.
!     2. Call splitwv.
!     3. Scatter the data in the root process to each process.
! This module is supported for multiple processing.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_cpondpe
      use m_defmpi
      use m_splitwv

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: specnudm, s_specnudm

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface specnudm

        module procedure s_specnudm

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
      subroutine s_specnudm( xdim, ydim, xsub, ysub, spnx, spny, ni, nj,&
     &                       nk, nsl, numpe, spnk2pe, spnpe2k, pval )
!***********************************************************************

! Input variables

      integer, intent(in) :: xdim
                       ! Model dimension in x direction

      integer, intent(in) :: ydim
                       ! Model dimension in y direction

      integer, intent(in) :: xsub
                       ! Number of sub domain in x direction

      integer, intent(in) :: ysub
                       ! Number of sub domain in y direction

      integer, intent(in) :: spnx
                       ! Truncation wave number in x direction

      integer, intent(in) :: spny
                       ! Truncation wave number in y direction

      integer, intent(in) :: ni
                       ! Model dimension of each sub domain in x direction

      integer, intent(in) :: nj
                       ! Model dimension of each sub domain in y direction

      integer, intent(in) :: nk
                       ! Model dimension of each sub domain in z direction

      integer, intent(in) :: nsl
                       ! ngglev+1 (actual lowest level)
                       ! Model dimension of each sub domain in z direction

      integer, intent(in) :: numpe
                       ! Total number of processor elements

      integer, intent(in) :: spnk2pe(nsl:nk-2)
                       ! Array to classify SPN processor number from k

      logical, intent(in) :: spnpe2k(1:numpe)
                       ! Array to classify k from SPN processor number

! Input and Output variables

      real, intent(inout) :: pval(0:ni+1,0:nj+1,nsl:nk-2)
                       ! Array which FFT will be performed

! Internal private variables

      integer IERROR
      integer i            ! Array index in x direction
      integer j            ! Array index in y direction
      integer k            ! Array index in z direction

      integer :: icounter  ! Array index in one-dimension
      integer :: ip        ! Array index in sub x direction
      integer :: jp        ! Array index in sub y direction

!-- each processor (level) array
      real :: tot3sp(1:(xdim-3)*(ydim-3))
      real :: tot2sp(1:(xdim-3+xsub)*(ydim-3+ysub))

!-- each processor (domain) array
      real :: pat3sp(1:(ni-3)*(nj-3),nsl:nk-2)
      real :: pat2sp(1:(ni-2)*(nj-2),nsl:nk-2)

      real :: ival(1:xdim-3,1:ydim-3)
      real :: tmpval(1:ni-2,1:nj-2,nsl:nk-2)
      real :: o3val(1:xdim-3,1:ydim-3)
      real :: o2val(1:xdim-3+xsub,1:ydim-3+ysub)

!-- Rearrange pval to tmpval (2:ni-1, 2:nj-1 -> 1:ni-2, 1:nj-2)

      do k=nsl,nk-2
      do j=1,nj-2
      do i=1,ni-2
        tmpval(i,j,k)=pval(i+1,j+1,k)
      end do
      end do
      end do

      pval=0.0e0

!-- Communicate patsp (ni-2,nj-2 -> ni-3,nj-3)
!-- for each processor (pe = nsl -- nk-2)

      do k=nsl,nk-2

        icounter=0

        do j=1,nj-3
        do i=1,ni-3
          icounter=icounter+1
          pat3sp(icounter,k)=tmpval(i,j,k)
        end do
        end do

        call cpondpe

!-- gathering the partial data to the whole data

        call MPI_GATHER( pat3sp(1,k), (ni-3)*(nj-3), MPI_REAL,  &
  &                      tot3sp(1), (ni-3)*(nj-3), MPI_REAL,  &
  &                      spnk2pe(k), mpi_comm_cress, IERROR )
                         ! pe == spnk2pe(k) (not root)

        call cpondpe

      end do

      if(spnpe2k(mype-root+1).eqv..true.)then
        write(*,*) "This processor is ", mype
!-- rearranging from comm array to data array

        icounter=0

        do jp=1,ysub
        do ip=1,xsub
        do j=1,nj-3
        do i=1,ni-3
          icounter=icounter+1
          ival((ip-1)*(ni-3)+i,(jp-1)*(nj-3)+j)=tot3sp(icounter)
        end do
        end do
        end do
        end do

        call splitwv( xdim, ydim, spnx, spny, ival, o3val )

!-- rearranging from data array to comm array (ni-3,nj-3 -> ni-2,nj-2)

! process 1 (document p.1)

!$omp parallel default(shared)
!$omp do schedule(runtime) private(i,j,ip,jp)

        do jp=1,ysub
        do ip=1,xsub
        do j=1,nj-3
        do i=1,ni-3
           o2val((ip-1)*(ni-2)+i,(jp-1)*(nj-2)+j)=o3val((ip-1)*(ni-3)+i,(jp-1)*(nj-3)+j)
        end do
        end do
        end do
        end do

!$omp end do
!$omp end parallel

! process 2-1 (document p.2)

        do jp=1,ysub-1
        do ip=1,xsub-1
          do j=1,nj-3
            o2val(ip*(ni-2),(jp-1)*(nj-2)+j)=o3val(ip*(ni-3)+1,(jp-1)*(nj-3)+j)
          end do
          do i=1,ni-3
            o2val((ip-1)*(ni-2)+i,jp*(nj-2))=o3val((ip-1)*(ni-3)+i,jp*(nj-3)+1)
          end do

          o2val(ip*(ni-2),jp*(nj-2))=o3val(ip*(ni-3)+1,jp*(nj-3)+1)

        end do
        end do

! process 2-2 (document p.2)

        do jp=1,ysub-1   ! == north-south lateral in the most westside domains
          do i=1,ni-3
            o2val((xsub-1)*(ni-2)+i,jp*(nj-2))=o3val((xsub-1)*(ni-3)+i,jp*(nj-3)+1)
          end do
        end do

        do ip=1,xsub-1   ! == east-west lateral in the most northside domains
          do j=1,nj-3
            o2val(ip*(ni-2),(ysub-1)*(nj-2)+j)=o3val(ip*(ni-3)+1,(ysub-1)*(nj-3)+j)
          end do
        end do

! process 3 (document p.3)

        do i=1,xdim-3+xsub-1   ! == xsub*(ni-2)-1
          o2val(i,ydim-3+ysub)=o2val(i,ydim-3+ysub-1)
        end do
        do j=1,ydim-3+ysub-1   ! == ysub*(nj-2)-1
          o2val(xdim-3+xsub,j)=o2val(xdim-3+xsub-1,j)
        end do

        o2val(xdim-3+xsub,ydim-3+ysub)=o2val(xdim-3+xsub-1,ydim-3+ysub-1)

        icounter=0

        do jp=1,ysub
        do ip=1,xsub
        do j=1,nj-2
        do i=1,ni-2
          icounter=icounter+1
          tot2sp(icounter)=o2val((ip-1)*(ni-2)+i,(jp-1)*(nj-2)+j)
        end do
        end do
        end do
        end do

      end if

      call cpondpe

!-- splitting the whole data to partial data

      do k=nsl,nk-2

        call cpondpe

        call MPI_SCATTER( tot2sp(1), (ni-2)*(nj-2), MPI_REAL,  &
  &                       pat2sp(1,k), (ni-2)*(nj-2), MPI_REAL,  &
  &                       spnk2pe(k), mpi_comm_cress, IERROR )
                          ! pe == spnk2pe(k) (not root)

        call cpondpe

!-- rearranging from comm array to data array

        icounter=0

        do j=1,nj-2
        do i=1,ni-2
          icounter=icounter+1
          tmpval(i,j,k)=pat2sp(icounter,k)
        end do
        end do

      end do

!-- Rearrange tmpval to pval

      do k=nsl,nk-2
      do j=1,nj-2
      do i=1,ni-2
        pval(i+1,j+1,k)=tmpval(i,j,k)
      end do
      end do
      end do

      end subroutine s_specnudm

!-----7--------------------------------------------------------------7--

      end module m_specnudm

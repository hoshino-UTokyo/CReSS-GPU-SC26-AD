!***********************************************************************
      module m_splitwv
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2016/12/23, 2017/06/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform FFT and isolate the long wave components for ISPACK library.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comspc

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: splitwv, s_splitwv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface splitwv

        module procedure s_splitwv

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
      subroutine s_splitwv( xdim, ydim, spnx, spny, ival, oval )
!***********************************************************************

      use m_comspc

      implicit none

! Input variables

      integer, intent(in) :: xdim
                       ! Model dimension in x direction

      integer, intent(in) :: ydim
                       ! Model dimension in y direction

      integer, intent(in) :: spnx
                       ! Truncation wave number in x direction

      integer, intent(in) :: spny
                       ! Truncation wave number in y direction

      real, intent(in) :: ival(1:xdim-3,1:ydim-3)
                       ! Array which FFT will be performed

! Output variables

      real, intent(out) :: oval(1:xdim-3,1:ydim-3)
                       ! Array which FFT was performed

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer wnx      ! Working array number of x direction
      integer wny      ! Working array number of y direction

      double precision :: cpival(1:ydim-3,1:xdim-3)
                       ! variable in real space before FFT
                       ! == cpival(1:ngy,1:ngx)
      double precision :: cpoval(-((ydim-3)/2-1):(ydim-3)/2-1,  &
  &                              -((xdim-3)/2-1):(xdim-3)/2-1)
                       ! variable in real space before FFT
                       ! == cpoval(-(ngy/2-1):ngy/2-1,-(ngx/2-1):ngx/2-1)

! In the future, the following arrays will be not allocate attribute.

      double precision :: sptmp(1:ydim-3,1:xdim-3)
                       ! Temporary array after FFT
      double precision ::   &
  &     sp(-((ydim-3)/2-1):(ydim-3)/2-1,-((xdim-3)/2-1):((xdim-3))/2-1)
                       ! Isolated wave number in spectral space

! Initialize and allocate internal arrays

      wnx=xdim-3
      wny=ydim-3

      cpival=0.0d0
      cpoval=0.0d0
      sptmp=0.0d0
      sp=0.0d0

! Substitute the input variable into internal private array

      do j=1,wny
         do i=1,wnx
            cpival(j,i)=dble(ival(i,j))
         end do
      end do

! Calculate RFFT (Real Space -> Spectral Space)

      call P2G2SA( wny/2-1, wnx/2-1, wny, wnx, cpival, cpoval, WGL,  &
  &                LITJ, LTJ, LITI, LTI )

! Isolate the long wave components from all waves

!-- wave number 0 (horizontal mean field)

      sp(0,0)=cpoval(0,0)

!-- boundary for each direction of wave

      if(spny>0)then

        do j=1,spny

          sp(j,0)=cpoval(j,0)
          sp(-j,0)=cpoval(-j,0)

        end do

      end if

      if(spnx>0)then

        do i=1,spnx

          sp(0,i)=cpoval(0,i)
          sp(0,-i)=cpoval(0,-i)

        end do

      end if

!-- internal region for each direction of wave

      if(spnx>0.and.spny>0)then

        do i=1,spnx
        do j=1,spny
          sp(j,i)=cpoval(j,i)
          sp(j,-i)=cpoval(j,-i)
          sp(-j,i)=cpoval(-j,i)
          sp(-j,-i)=cpoval(-j,-i)
        end do
        end do

      end if

! Calculate IFFT (Spectral Space -> Real Space)

      call P2S2GA( wny/2-1, wnx/2-1, wny, wnx, sp, sptmp, WGL,  &
  &                LITJ, LTJ, LITI, LTI )

! Substitute the internal private variable into output array

      do j=1,ydim-3
         do i=1,xdim-3
            oval(i,j)=real(sptmp(j,i))
         end do
      end do

  end subroutine s_splitwv

end module m_splitwv

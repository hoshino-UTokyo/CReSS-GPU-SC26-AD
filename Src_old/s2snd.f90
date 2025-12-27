!***********************************************************************
      module m_s2snd
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2016/06/25
!     Modification: 2016/08/18, 2017/01/29, 2017/02/01, 2017/06/09

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the spectral nudging to GPV data of optional scalar
!     variable for multiple processor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comspc
      use m_getcname
      use m_getiname
      use m_inichar
      use m_specnudm
      use m_specnuds

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: s2snd, s_s2snd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface s2snd

        module procedure s_s2snd

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
      subroutine s_s2snd(fpnggvar,fpxdim,fpydim,fpxsub,fpysub,fpspnx,   &
     &                   fpspny,fpngglev,fpnumpe,nggdmp,gtinc,ni,nj,    &
     &                   nk,rst,sp,sgpv,std,sfrc,fft_flag)
!***********************************************************************

      use m_comspc

! Input variables

      integer, intent(in) :: fpnggvar
                       ! Formal parameter of unique index of nggvar

      integer, intent(in) :: fpxdim
                       ! Formal parameter of unique index of xdim

      integer, intent(in) :: fpydim
                       ! Formal parameter of unique index of ydim

      integer, intent(in) :: fpxsub
                       ! Formal parameter of unique index of xsub

      integer, intent(in) :: fpysub
                       ! Formal parameter of unique index of ysub

      integer, intent(in) :: fpspnx
                       ! Formal parameter of unique index of spnx

      integer, intent(in) :: fpspny
                       ! Formal parameter of unique index of spny

      integer, intent(in) :: fpngglev
                       ! Formal parameter of unique index of ngglev

      integer, intent(in) :: fpnumpe
                       ! Formal parameter of unique index of numpe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian

      real, intent(in) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at past

      real, intent(in) :: sgpv(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable of GPV data
                       ! at marked time

      real, intent(in) :: std(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! optional scalar variable of GPV data

      logical, intent(in) :: fft_flag
                       ! Whether forcings for spectral nudging is updated or not

! Input and output variable

      real, intent(inout) :: sfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term of optional scalar

! Internal shared variable

      character(len=108) nggvar
                       ! Control flag of
                       ! analysis nudged variables to GPV

      integer xdim
                       ! Model dimension in x direction

      integer ydim
                       ! Model dimension in y direction

      integer xsub
                       ! Number of sub domain in x direction

      integer ysub
                       ! Number of sub domain in y direction

      integer spnx
                       ! Truncation wave number in x direction

      integer spny
                       ! Truncation wave number in y direction

      integer ngglev
                       ! The lowest level of spectral nudging to GPV

      integer numpe
                       ! Total number of processor elements

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction
      integer nsl      ! ngglev+1 (actual lowest level)

      real :: wcoe(1:nk)
                       ! weighted coefficient with varying height
                       ! If ngglev > 1,
                       ! wcoe = nggcoe * {k-(ngglev+1)}/{nk-2-(ngglev+1)}
                       ! If ngglev = 1, wcoe = 1

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(nggvar)

! -----

! Get the required namelist variable.

      call getcname(fpnggvar,nggvar)
      call getiname(fpxdim,xdim)
      call getiname(fpydim,ydim)
      call getiname(fpxsub,xsub)
      call getiname(fpysub,ysub)
      call getiname(fpspnx,spnx)
      call getiname(fpspny,spny)
      call getiname(fpngglev,ngglev)
      call getiname(fpnumpe,numpe)

! -----

      nsl=ngglev+1
      wcoe=0.0e0

      if(ngglev>1)then
        do k=nsl,nk-2
          wcoe(k)=real(k-nsl)/real(nk-2-nsl)
        end do
      else
        wcoe=1.0e0
      end if

! Calculate the analysis nudging terms for scalar variables.

! Update the forcing value.

      if(fft_flag.eqv..true.)then

!$omp parallel default(shared)
!$omp do schedule(runtime) private(i,j,k)

        do k=nsl,nk-2
        do j=2,nj-2
        do i=2,ni-2
          spnsfrc(i,j,k)=wcoe(k)*nggdmp*rst(i,j,k)                      &
     &                  *((sgpv(i,j,k)+std(i,j,k)*gtinc)-sp(i,j,k))
        end do
        end do
        end do

!$omp end do
!$omp end parallel

! Removing high frequency component in spnsfrc.

        if(spn_mul_flag.eqv..true.)then  ! SPN for multiple processing

          call specnudm( xdim, ydim, xsub, ysub, spnx, spny,            &
     &                   ni, nj, nk, nsl, numpe,                        &
     &                   spnk2pe(nsl:nk-2), spnpe2k(1:numpe),           &
     &                   spnsfrc(0:ni+1,0:nj+1,nsl:nk-2) )

        else  ! SPN for single processing

! NOTE : This loop must not be parallelizing by OpenMP.

          do k=nsl,nk-2
            call specnuds( xdim, ydim, xsub, ysub, spnx, spny, ni, nj,  &
     &                     spnsfrc(0:ni+1,0:nj+1,k) )
          end do

        end if

      end if

!$omp parallel default(shared)
!$omp do schedule(runtime) private(i,j,k)

      do k=nsl,nk-2
      do j=2,nj-2
      do i=2,ni-2
        sfrc(i,j,k)=sfrc(i,j,k)+spnsfrc(i,j,k)
      end do
      end do
      end do

!$omp end do
!$omp end parallel

! -----

      end subroutine s_s2snd

!-----7--------------------------------------------------------------7--

      end module m_s2snd

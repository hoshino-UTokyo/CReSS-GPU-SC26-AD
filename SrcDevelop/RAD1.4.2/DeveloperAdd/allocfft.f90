!***********************************************************************
      module m_allocfft
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2016/12/23, 2017/01/29, 2017/02/01, 2017/06/09

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for FFT for ISPACK library
!     for multiple processor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_comspc
      use m_cpondpe
      use m_destroy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocfft, s_allocfft

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocfft

        module procedure s_allocfft

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
      subroutine s_allocfft(fpnumpe,fpnggopt,fpngglev,fpxdim,fpydim,    &
     &                      fpzdim)
!***********************************************************************

      use m_comspc

! Input Valriables
      integer, intent(in) :: fpnumpe
                       ! Formal parameter of unique index of numpe

      integer, intent(in) :: fpnggopt
                       ! Formal parameter of unique index of nggopt

      integer, intent(in) :: fpngglev
                       ! Formal parameter of unique index of ngglev

      integer, intent(in) :: fpxdim
                       ! Formal parameter of unique index of xdim

      integer, intent(in) :: fpydim
                       ! Formal parameter of unique index of ydim

      integer, intent(in) :: fpzdim
                       ! Formal parameter of unique index of zdim

! Internal shared variables
      integer numpe    ! Total number of processor elements

      integer nggopt   ! Option for analysis nudging to GPV

      integer ngglev   ! The lowest level of spectral nudging to GPV

      integer xdim     ! Model number of x direction

      integer ydim     ! Model number of y direction

      integer zdim     ! Model number of y direction

! Internal Valriables
      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

      integer k, m     ! Loop index

      integer wnx      ! Working grid number of x direction

      integer wny      ! Working grid number of y direction

!-----7--------------------------------------------------------------7--

      call getiname(fpnumpe,numpe)
      call getiname(fpnggopt,nggopt)
      call getiname(fpngglev,ngglev)
      call getiname(fpzdim,zdim)

      if(numpe>zdim)then
        spn_mul_flag=.true.   ! SPN for multiple processing
      else
        spn_mul_flag=.false.  ! SPN for single processing
      end if

      if(nggopt.eq.2)then

        stat=0

        if(spn_mul_flag.eqv..true.)then   ! SPN for multiple processing

          allocate(spnk2pe(1:zdim),stat=cstat)

          stat=stat+abs(cstat)

          allocate(spnpe2k(1:numpe),stat=cstat)

          stat=stat+abs(cstat)

! Set processor number corresponding to k, and k corresponding to
! processor number.

          spnk2pe=0
          spnpe2k=.false.

          do k=ngglev+1,zdim-2   ! == nsl -- nk-2
            m=1+(numpe/zdim)*(k-1)  ! SPN process number
            spnk2pe(k)=m       ! SPN process number for k
            spnpe2k(m-root+1)=.true.  ! whether SPN process number for k or not
          end do

          if(spnpe2k(mype-root+1).eqv..true.)then

            call getiname(fpxdim,xdim)
            call getiname(fpydim,ydim)

            wnx=xdim-3
            wny=ydim-3

            allocate(LTI(wnx*2),stat=cstat)

            stat=stat+abs(cstat)

            allocate(LTJ(wny*2),stat=cstat)

            stat=stat+abs(cstat)

            allocate(WGL(wny,wnx),stat=cstat)

            stat=stat+abs(cstat)

! Calculate basic arrays for prime factors (2,3,5,7)

            call P2INIT( wny, wnx, LITJ, LTJ, LITI, LTI )

! Initialize allocated arrays

            WGL=0.0d0

          end if

        else   ! SPN for single processing

          if(mype==root)then

            call getiname(fpxdim,xdim)
            call getiname(fpydim,ydim)

            wnx=xdim-3
            wny=ydim-3

            allocate(LTI(wnx*2),stat=cstat)

            stat=stat+abs(cstat)

            allocate(LTJ(wny*2),stat=cstat)

            stat=stat+abs(cstat)

            allocate(WGL(wny,wnx),stat=cstat)

            stat=stat+abs(cstat)

! Calculate basic arrays for prime factors (2,3,5,7)

            call P2INIT( wny, wnx, LITJ, LTJ, LITI, LTI )

! Initialize allocated arrays

            WGL=0.0d0

          end if

        end if

! If error occured, call the procedure destroy.

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('allocfft',8,'cont',5,'              ',14,101, &
     &                    stat)

          end if

          call cpondpe

          call destroy('allocfft',8,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

      else

        write(*,*) "### Message (allocfft) ### : Not allocate arrays "  &
     &           //"for spectral nudging."

      end if

! -----

!! -----

      end subroutine s_allocfft

!-----7--------------------------------------------------------------7--

      end module m_allocfft

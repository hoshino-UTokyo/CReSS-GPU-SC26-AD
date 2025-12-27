!***********************************************************************
      module m_allocspc
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2016/04/08
!     Modification: 2016/08/18, 2017/06/09

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for spectral nudging.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_comspc
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_getcname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocspc, s_allocspc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocspc

        module procedure s_allocspc

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
      subroutine s_allocspc(idi1,idc1,ni,nj,nk)
!***********************************************************************

! Input Valriables
      integer, intent(in) :: idi1
                       ! Unique index in integer namelist table

      integer, intent(in) :: idc1
                       ! Unique index in character namelist table
                       ! Not active

      integer, intent(in) :: ni
                       ! Model dimension in sub x direction

      integer, intent(in) :: nj
                       ! Model dimension in sub y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Internal Valriables
      integer nggopt   ! Option for analysis nudging to GPV

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

      character(len=108) nggvar
                       ! Control flag of
                       ! analysis nudged variables to GPV

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(idi1,nggopt)
      call getcname(idc1,nggvar)

      if(nggopt.eq.2)then

        stat=0

        allocate(spnufrc(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(spnvfrc(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(spnsfrc(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

! Initialize allocated arrays

        spnufrc=0.0e0
        spnvfrc=0.0e0
        spnsfrc=0.0e0

! If error occured, call the procedure destroy.

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('allocspc',8,'cont',5,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('allocspc',8,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

      else

        write(*,*) "### Message (allocspc) ### : Not allocate arrays "  &
     &           //"for spectral nudging."

      end if

! -----

!! -----

      end subroutine s_allocspc

!-----7--------------------------------------------------------------7--

      end module m_allocspc

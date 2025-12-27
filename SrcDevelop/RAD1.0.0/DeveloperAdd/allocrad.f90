!***********************************************************************
      module m_allocrad
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/08
!     Modification: 2010/12/17, 2010/12/21, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for mstranx radiation scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_commstrn
      use m_comrad
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_setcst1d
      use m_setcst2d
      use m_setcst3d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocrad, s_allocrad

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocrad

        module procedure s_allocrad

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
      subroutine s_allocrad(fpcphopt,fpradopt,ni,nj,nk)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpradopt
                       ! Formal parameter of unique index of radopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer radopt   ! Option for turning on mstranx radiation scheme

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fpradopt,radopt)

! ------

!! Allocate the array for mstranx radiation scheme.

! Perform allocate.

      stat=0

      if(abs(cphopt).le.4.and.radopt.eq.1) then

        allocate(zl(1:nk-3),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pl(1:nk-3),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tl(1:nk-3),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pb(1:nk-2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tb(1:nk-2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cpcl(1:nk-3,1:kpcl,1:kpclc),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cpcl_data(1:kln,1:kpcl,1:kpclc),stat=cstat)

        stat=stat+abs(cstat)

        allocate(gdcfrc(1:nk-3),stat=cstat)

        stat=stat+abs(cstat)

        allocate(gdcfrc_data(1:kln),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cgas(1:nk-3,1:kmol),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cgas_data(1:kln,1:kmol),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ccfc(1:kcfc),stat=cstat)

        stat=stat+abs(cstat)

        allocate(prg(1:kprg),stat=cstat)

        stat=stat+abs(cstat)

        allocate(fd(1:nk-2,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(fu(1:nk-2,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(htrsd(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(htrsu(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(htrld(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(htrlu(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocrad',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocrad',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

! Fill in all array for the program mstranx radiation scheme with 0.

      if(abs(cphopt).le.4.and.radopt.eq.1) then

        call setcst1d(1,nk-3,0.e0,zl)
        call setcst1d(1,nk-3,0.e0,pl)
        call setcst1d(1,nk-3,0.e0,tl)

        call setcst1d(1,nk-2,0.e0,pb)
        call setcst1d(1,nk-2,0.e0,tb)

        call setcst3d(1,nk-3,1,kpcl,1,kpclc,0.e0,cpcl)
        call setcst3d(1,kln,1,kpcl,1,kpclc,0.e0,cpcl_data)

        call setcst1d(1,nk-3,0.e0,gdcfrc)
        call setcst1d(1,kln,0.e0,gdcfrc_data)

        call setcst2d(1,nk-3,1,kmol,0.e0,cgas)
        call setcst2d(1,kln,1,kmol,0.e0,cgas_data)

        call setcst1d(1,kcfc,0.e0,ccfc)

        call setcst1d(1,kprg,0.e0,prg)

        call setcst2d(1,nk-2,1,2,0.e0,fd)
        call setcst2d(1,nk-2,1,2,0.e0,fu)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,htrsd)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,htrsu)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,htrld)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,htrlu)

      end if

! -----

      end subroutine s_allocrad

!-----7--------------------------------------------------------------7--

      end module m_allocrad

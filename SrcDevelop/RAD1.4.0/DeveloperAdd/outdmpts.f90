!***********************************************************************
      module m_outdmpts
!***********************************************************************

!     Author      : Hasegawa Koichi
!     Date        : 2014/06/10
!     Modification:

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the 3 dimensional variables (soil temperature)
!     to the dumped file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comdmp
      use m_comindx
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_outcap
      use m_outstd09

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outdmpts, s_outdmpts

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outdmpts

        module procedure s_outdmpts

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
      subroutine s_outdmpts(vname,ncvn,vcap,ncvc,xo,ni,nj,nund,var3d)
!***********************************************************************

! Input variables

      character(len=6), intent(in) :: vname
                       ! Optional variable name

      character(len=60), intent(in) :: vcap
                       ! Caption for dumped variable

      character(len=3), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: ncvn
                       ! Number of character of vname

      integer, intent(in) :: ncvc
                       ! Number of character of vcap

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nund
                       ! Model dimension in z direction

      real, intent(in) :: var3d(0:ni+1,0:nj+1,1:nund)
                       ! 3 dimensional optional variable

! Internal shared variables

      integer dmpfmt   ! Option for dumped file format
      integer dmpmon   ! Option for monitor variables output

      integer stat     ! Runtime status

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer cntdmp   ! Counter of dumped variables

!-----7--------------------------------------------------------------7--

!! Read in the 3 dimensional variables to the dumped file.

      if(fdmp(1:3).eq.'act') then

! Get the required namelist variables.

        call getiname(iddmpfmt,dmpfmt)
        call getiname(iddmpmon,dmpmon)

! -----

! Initialize the runtime status.

        stat=1

! -----

! Increse the cnt3d.

        cnt3d=cnt3d+1

! -----

! Read in the variable at the u points to the dumped file.

        if(xo(1:3).eq.'oxx') then

          if(dmpfmt.eq.1) then

            do k=1,nund

              write(io3d,*,iostat=stat,err=100)                         &
     &          ((.5e0*(var3d(i,j,k)+var3d(i+1,j,k)),i=2,ni-2),j=2,nj-2)

            end do

          else if(dmpfmt.eq.2) then

            do k=1,nund

              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &          ((.5e0*(var3d(i,j,k)+var3d(i+1,j,k)),i=2,ni-2),j=2,nj-2)

            end do

          end if

! -----

! Read in the variable at the v points to the dumped file.

        else if(xo(1:3).eq.'xox') then

          if(dmpfmt.eq.1) then

            do k=1,nund

              write(io3d,*,iostat=stat,err=100)                         &
     &          ((.5e0*(var3d(i,j,k)+var3d(i,j+1,k)),i=2,ni-2),j=2,nj-2)

            end do

          else if(dmpfmt.eq.2) then

            do k=1,nund

              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &          ((.5e0*(var3d(i,j,k)+var3d(i,j+1,k)),i=2,ni-2),j=2,nj-2)

            end do

          end if

! -----

! Read in the variable at the w points to the dumped file.

        else if(xo(1:3).eq.'xxo') then

          if(dmpfmt.eq.1) then

            do k=1,nund

              write(io3d,*,iostat=stat,err=100)                         &
     &          ((.5e0*(var3d(i,j,k)+var3d(i,j,k+1)),i=2,ni-2),j=2,nj-2)

            end do

          else if(dmpfmt.eq.2) then

            do k=1,nund

              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &          ((.5e0*(var3d(i,j,k)+var3d(i,j,k+1)),i=2,ni-2),j=2,nj-2)

            end do

          end if

! -----

! Read in the variable at the scalar points to the dumped file.

        else if(xo(1:3).eq.'xxx') then

          if(dmpfmt.eq.1) then

            do k=1,nund

              write(io3d,*,iostat=stat,err=100)                         &
     &             ((var3d(i,j,k),i=2,ni-2),j=2,nj-2)

            end do

          else if(dmpfmt.eq.2) then

            do k=1,nund

              rec3d=rec3d+1

              write(io3d,rec=rec3d,iostat=stat,err=100)                 &
     &             ((var3d(i,j,k),i=2,ni-2),j=2,nj-2)

            end do

          end if

        end if

! -----

! If error occured, call the procedure destroy.

  100   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outdmpts',8,'cont',4,'              ',14,     &
     &                   io3d,stat)

          end if

          call cpondpe

          call destroy('outdmpts',8,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

! -----

! Read in the messages to the dumped data checking file.

        call outcap('dmp',vname,vcap,ncvc,io3c,nund+3)

! -----

! Read in the messages to the standard i/o.

        if(mype.eq.root) then

          if(dmpmon.eq.0) then

            call outstd09(vname,ncvn,vcap,ncvc,3,cnt3d)

          else

            cntdmp=cnt3d+cnt2d

            call outstd09(vname,ncvn,vcap,ncvc,3,cntdmp)

          end if

        end if

        call cpondpe

! -----

      end if

!! -----

      end subroutine s_outdmpts

!-----7--------------------------------------------------------------7--

      end module m_outdmpts

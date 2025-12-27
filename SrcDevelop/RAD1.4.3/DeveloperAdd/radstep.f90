!***********************************************************************
      module m_radstep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/20
!     Modification: 2010/12/17, 2013/02/13

!     Author      : Hasegawa Koichi
!     Modification: 2013/11/28, 2014/06/10, 2015/07/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the control flag of radiation scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: radstep, s_radstep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface radstep

        module procedure s_radstep

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_radstep(fpradopt,fpraddlt,ctime,radon)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpradopt
                       ! Formal parameter of unique index of radopt

      integer, intent(in) :: fpraddlt
                       ! Formal parameter of unique index of raddlt

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

! Output variable

      character(len=6), intent(out) :: radon
                       ! Control flag of mstranx radiation scheme

! Internal shared variables

      integer radopt   ! Option for turnig on mstranx radation scheme

      real raddlt(2)   ! Time interval of radiation scheme

! Internal private variables

      integer(kind=i8) :: imod1
      integer(kind=i8) :: imod2
                       ! Temporary

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpradopt,radopt)
      call getrname(fpraddlt  ,raddlt(1))
      call getrname(fpraddlt+1,raddlt(2))

! -----

! Get the control flag of radiation scheme.

      if(radopt==1) then

        ! ** mstrnX ** 

        call getrname(fpraddlt,raddlt(1))

        if(mod(ctime,10_i8*int(1.e2*(raddlt(1)+.001e0),i8)).eq.0_i8) then

          write(radon(1:6),'(a6)') 'motion'

        else

          write(radon(1:6),'(a6)') 'rest  '

        end if

      else if(radopt==10) then

        ! ** RRTM ** 

        call getrname(fpraddlt  ,raddlt(1))
        call getrname(fpraddlt+1,raddlt(2))

        imod1 = mod(ctime,10_i8*int(1.e2*(raddlt(1)+.001e0),i8))
        imod2 = mod(ctime,10_i8*int(1.e2*(raddlt(2)+.001e0),i8))

        if( (raddlt(1)>=0.01e0) .and. (imod1==0_i8) ) then

          if( (raddlt(2)>=0.01e0) .and. (imod2==0_i8) ) then

            write(radon(1:6),'(a6)') 'lsw   '

          else

            write(radon(1:6),'(a6)') 'sw    '

          end if

        else if( (raddlt(2)>=0.01e0) .and. (imod2==0_i8) ) then

          write(radon(1:6),'(a6)') 'lw    '

        else

          write(radon(1:6),'(a6)') 'rest  '

        end if

      else

        write(radon(1:6),'(a6)') 'rest  '

      end if

! -----

      end subroutine s_radstep

!-----7--------------------------------------------------------------7--

      end module m_radstep

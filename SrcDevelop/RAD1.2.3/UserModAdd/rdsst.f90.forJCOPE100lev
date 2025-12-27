!***********************************************************************
      module m_rdsst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2002/07/15, 2002/09/09, 2002/11/11, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2004/05/31,
!                   2005/02/10, 2007/01/20, 2007/08/24, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/02/27, 2009/03/31,
!                   2011/11/10, 2013/02/13, 2013/03/27
!
!     Author      : Mayumi Yoshioka
!     Modification: 2015/03/16
!
!      ** modifiied for reading MGDSST before reading JCOPE SST
!      ** (undef skip), 2011/01/20

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the 2 dimensional sea surface
!     temperature data file.

!-----7--------------------------------------------------------------7--

! ####################################################################
! ##                                                                ##
! ## You have to modify reading section to adjust your data set.    ##
! ##                                                                ##
! ####################################################################

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_comphy         ! added for JCOPE2 data
      use m_destroy
      use m_getname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_putunit

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdsst, s_rdsst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdsst

        module procedure s_rdsst

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
      subroutine s_rdsst(fpdatdir,fpncdat,fpwlngth,cdate,ctime,stat,    &
     &                   nid,njd,sstdat)
!***********************************************************************

! Input variables

      character(len=12), intent(in) :: cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: fpdatdir
                       ! Formal parameter of unique index of datdir

      integer, intent(in) :: fpncdat
                       ! Formal parameter of unique index of ncdat

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

! Output variables

      integer, intent(out) :: stat
                       ! Runtime status

      real, intent(out) :: sstdat(1:nid,1:njd)
                       ! Sea suface temperature in data

! Internal shared variables

      character(len=108) datdir
                       ! User specified directory for external data

      character(len=108) sstfl
                       ! Opened file name

      integer ncdat    ! Number of character of datdir

      integer wlngth   ! Word length of direct access file

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction

      integer ncfl     ! Number of character of sstfl

      integer iosst    ! Unit number of sea surface temperature file

      integer recsst   ! Record number of sea surface temperature file

      integer siz      ! Record length of sea surface temperature file

! For JCOPE2 data

      integer kd       ! Array index in z direction  !! for 3D ocean data

      integer sumc     ! counter to make zonal averaged SST

      real   sumt(njd) ! zonal averaged SST to fill undef grids

!-----7--------------------------------------------------------------7--

!! Open and read out the data from the 2 dimensional sea surface
!! temperature data file.

! Initialize the character variable.

      call inichar(datdir)

! -----

! Get the required namelist variables.

      call getname(fpdatdir,datdir)
      call getname(fpncdat,ncdat)
      call getname(fpwlngth,wlngth)

! -----

! Initialize the character variable.

      call inichar(sstfl)

! -----

! Get the unit number.

      call getunit(iosst)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Open the 2 dimensional sea surface temperature data file.

! ##### You will have to modify the following lines. #####

      siz=nid*njd*wlngth

      ncfl=24

      write(sstfl(1:ncfl),'(a8,a12,a4)') 'JCOPEsst',cdate(1:12),'.bin'
!!      write(sstfl(1:ncfl),'(a8,a12,a4)') 'data.sst',cdate(1:12),'.bin'

!!v3.2.1
!!      write(sstfl(1:ncfl),'(a12)') 'JCOPEsst.bin'
!!org
!!      write(sstfl(1:ncfl),'(a12)') 'data.sst.bin'

      open(iosst,iostat=stat,err=100,                                   &
     &     file=datdir(1:ncdat)//sstfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

! ####################

  100 if(stat.ne.0) then

        call destroy('rdsst   ',5,'cont',1,'              ',14,iosst,   &
     &               stat)

        return

      end if

      call outstd03('rdsst   ',5,sstfl,ncfl,iosst,1,1,ctime)

! -----

! Read out the data from the 2 dimensional sea surface temperature data
! file.

! ##### You will have to modify the following lines. #####

!! First,
!! Read 3-D JCOPE ocean data, kd=1(bottom), ..., kd=100(depth= -1m, SST)

      recsst=1

      do kd=1,100

      recsst=1

      read(iosst,rec=recsst,iostat=stat,err=110)                        &
     &    ((sstdat(id,jd),id=1,nid),jd=1,njd)

      recsst=recsst+1

      end do

!! Second,
!! Check undef grid and make zonal average SST(kd=100), sumt(jd)

      do jd = 1, njd

         sumc=0
         sumt(jd) = 0.0e0

         do id = 1, nid
!                                                   not for undef grid
            if( sstdat(id,jd) .gt. -100.e0 ) then

               sumc = sumc + 1
               sumt(jd) = sumt(jd) + sstdat(id,jd)

            endif

         end do

         if( sumc .gt. 0 ) then

            sumt(jd) = sumt(jd)/sumc

         else

            sumt(jd) = 27.0e0

         endif

      end do

!! Third,
!! Check undef grid and fill the grid with zonal averaged sumt(jd)

      do jd= 1, njd

         do id = 1, nid
!                                                   for undef grid
            if( sstdat(id,jd) .le. -100.e0 ) then

               sstdat(id,jd) = sumt(jd)

            endif
!                                                   convert Kelvin Deg.
            sstdat(id,jd) = sstdat(id,jd) + t0


         end do
         
      end do


!! Fourth,
!! Check undef grid again.

      do jd = 1, njd

         do id = 1, nid
!                                                   check undef grid
            if( sstdat(id,jd) .le. 0.e0 ) then

               write(6,*) 'undef SST grid !!', id, jd, sstdat(id,jd)
               stop

            endif
            
         end do

      end do

!! Fifth,
!! Check write, output gradsfile

      open( 60, file='jcopesst.grd', form='unformatted', &
    &         status='unknown', access='direct', recl=siz )

      recsst=1
      write( 60, rec=recsst )( ( sstdat(id,jd), id=1, nid), jd=1, njd)


! ####################

  110 if(stat.ne.0) then

        call destroy('rdsst   ',5,'cont',3,'              ',14,iosst,   &
     &               stat)

        return

      end if

      call outstd03('rdsst   ',5,sstfl,108,iosst,3,1,ctime)

! -----

! Close the 2 dimensional sea surface temperature data file.

      close(iosst,iostat=stat,err=120,status='keep')

  120 if(stat.ne.0) then

        call destroy('rdsst   ',5,'cont',2,'              ',14,iosst,   &
     &               stat)

        return

      end if

      call outstd03('rdsst   ',5,sstfl,108,iosst,2,1,ctime)

! -----

! Return the unit number.

      call putunit(iosst)

! -----

!! -----

      end subroutine s_rdsst

!-----7--------------------------------------------------------------7--

      end module m_rdsst

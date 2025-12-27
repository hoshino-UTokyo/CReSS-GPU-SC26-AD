!***********************************************************************
      module m_raddrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/20
!     Modification: 2010/12/17, 2010/12/21, 2013/02/13, 2013/03/27,
!                   2013/10/08

!     Author      : Hasegawa Koichi
!     Modification: 2013/11/28, 2014/06/10, 2014/07/01, 2014/08/25,
!                   2015/07/22

!     Author      : Oda Naotaka
!     Modification: 2013/11/12

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior subroutines for mstranx/RRTM radiation
!     scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_adjstlnd
      use m_atmosin
      use m_changept
      use m_chkerr
      use m_chkstd
      use m_comindx
      use m_commstrn
      use m_commpi
      use m_comrad
      use m_cpondpe
      use m_cress21d
      use m_destroy
      use m_getcname
      use m_getexner
      use m_getiname
      use m_getta3d
      use m_inichar
      use m_outstd03
      use m_radheat
      use m_rrtmdrv
      use m_zenith

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: raddrv, s_raddrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface raddrv

        module procedure s_raddrv

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
      subroutine s_raddrv(fpdatdir,fpncdat,fpadvopt,fpsfcopt,fpradopt,  &
     &                    radon,pdate,dtb,                              &
     &                    ni,nj,nk,nqw,nqi,nund,land,albe,beta,         &
     &                    lat,lon,zph,pbr,ptbr,rbr,ppp,ptpp,qvp,qwtrp,  &
     &                    qicep,tundp,ptpf,                             &
     &                    coseta,p,t,tmp1,tmp2)
!***********************************************************************

! Input variables

      character(len=6), intent(in) :: radon
                       ! Control flag of radiation scheme

      character(len=12), intent(in) :: pdate
                       ! Forecast date at 1 step past
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: fpdatdir
                       ! Formal parameter of unique index of datdir

      integer, intent(in) :: fpncdat
                       ! Formal parameter of unique index of ncdat

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpradopt
                       ! Formal parameter of unique index of radopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of water hydrometeor array

      integer, intent(in) :: nqi
                       ! Number of ice hydrometeor array

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: albe(0:ni+1,0:nj+1)
                       ! Albedo

      real, intent(in) :: beta(0:ni+1,0:nj+1)
                       ! Evapotranspiration efficiency

      real, intent(in) :: lat(0:ni+1,0:nj+1)
                       ! Latitude

      real, intent(in) :: lon(0:ni+1,0:nj+1)
                       ! Longitude

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(in) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(in) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Ground temperature at past

! Input and output variable

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

! Internal shared variables

      character(len=108) datdir
                       ! User specified directory for external data

      character(len=2) cbnd
                       ! Data file extension

      character(len=64) err
                       ! Error messages

      integer ncdat    ! Number of character of datdir

      integer advopt   ! Option for advection scheme
      integer sfcopt   ! Option for surface physics
      integer radopt   ! Option for turning on mstranx radiation scheme

      integer nln      ! number of layers in data

      real aland       ! Adjusted real land use category

      real dtb_sub     ! Substitute for dtb

      real, intent(inout) :: coseta(0:ni+1,0:nj+1)
                       ! cos (Zenith angle), use work array

      real, intent(inout) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(inout) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary (Exner function)

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary

! Internal private variables

      character(len=108) parafl1
      character(len=108) parafl2
      character(len=108) parafl3
                       ! User specified file name for parameter file

      integer broot    ! Broadcasting root

      integer io       ! Temporary
      integer stat
      integer stat1

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(datdir)

! -----

! Get the required namelist variables.

      call getcname(fpdatdir,datdir)
      call getiname(fpncdat ,ncdat)
      call getiname(fpadvopt,advopt)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpradopt,radopt)

! -----

! Reset the large time steps interval.

      if(advopt.le.3) then
        dtb_sub = 2.e0*dtb
      else
        dtb_sub = dtb
      end if

! -----

!!!! Perform radiation scheme.

      if( radopt==1 ) then

!!! Perform MstrnX radiation scheme at marked time.

        if(radon.eq.'motion') then

          if(mype==root) then
            write(*,*)
            write(*,*)' ### run MSTRNX radiation scheme ###' 
            write(*,*)
          end if

!! open read files.

          write(cbnd,'(i2)') kbnd

          if(mype.eq.root) then

!!         open(iug,file=datdir(1:ncdat)//'DataMSTRN/PARAG.'//cbnd,     &
!!   &          status='old')

!!         open(iup,file=datdir(1:ncdat)//'DataMSTRN/PARAPC.'//cbnd,    &
!!   &          status='old')

!!         open(iuv,file=datdir(1:ncdat)//'DataMSTRN/VARDATA.RM'//cbnd, &
!!   &          status='old')

            io   = iuv
            stat = 0

            parafl1(1:ncdat) = datdir(1:ncdat)
            parafl2(1:ncdat) = datdir(1:ncdat)
            parafl3(1:ncdat) = datdir(1:ncdat)

            write(parafl1(ncdat+1:ncdat+18),'(a16,a2)')                 &
     &                                     'DataMSTRN/PARAG.',cbnd
            open(iug,iostat=stat1,err=100,file=parafl1(1:ncdat+18),     &
     &               status='old',access='sequential',form='formatted', &
     &               blank='null',position='rewind',action='read')

  100       continue
            if(stat1/=0) then
              io   = iug
              stat = stat1
              go to 150
            end if

            write(parafl2(ncdat+1:ncdat+19),'(a17,a2)')                 &
     &                                     'DataMSTRN/PARAPC.',cbnd
            open(iup,iostat=stat1,err=110,file=parafl2(1:ncdat+19),     &
     &               status='old',access='sequential',form='formatted', &
     &               blank='null',position='rewind',action='read')

  110       continue
            if(stat1/=0) then
              io   = iup
              stat = stat1
              go to 150
            end if

            write(parafl3(ncdat+1:ncdat+22),'(a20,a2)')                 &
     &                                     'DataMSTRN/VARDATA.RM',cbnd
            open(iuv,iostat=stat1,err=120,file=parafl3(1:ncdat+22),     &
     &               status='old',access='sequential',form='formatted', &
     &               blank='null',position='rewind',action='read')

  120       continue
            if(stat1/=0) then
              io   = iuv
              stat = stat1
              go to 150
            end if

          else

            stat=0

          end if

  150     call chkerr(stat)

          if(mype.eq.root) then
            if(io>=iug) call outstd03('raddrv  ',6,parafl1,ncdat+18,iug,1,0,0_i8)
            if(io>=iup) call outstd03('raddrv  ',6,parafl2,ncdat+19,iup,1,0,0_i8)
            if(io>=iuv) call outstd03('raddrv  ',6,parafl3,ncdat+22,iuv,1,0,0_i8)
          end if

          if(stat<0) then

             if(mype.eq.-stat-1) then
               call destroy('raddrv  ',6,'cont',1,'              ',     &
     &                      14,io,stat)
             end if

             call cpondpe
             call destroy('raddrv  ',6,'stop',1001,'              ',    &
     &                    14,101,stat)

          end if

          broot=stat-1

          call chkstd(broot)

!! -----

!! reading standard atmospheric condition file only once.

          nln=nk-3

          if( iradcnt==0 ) then

! open file.

            if(mype==root) then

              parafl1(1:ncdat) = datdir(1:ncdat)
              write(parafl1(ncdat+1:ncdat+23),'(a23)')                  &
     &                                       'DataMSTRN/DATA.gas_ptcl'

              open(iud,iostat=stat,err=200,file=parafl1(1:ncdat+23),    &
     &             status='old',access='sequential',form='formatted',   &
     &             blank='null',position='rewind',action='read')

            else

              stat=0

            end if

  200       call chkerr(stat)

            if(mype.eq.root) then
              call outstd03('raddrv  ',6,parafl1,ncdat+23,iud,1,0,0_i8)
            end if

            if(stat<0) then

               if(mype.eq.-stat-1) then
                 call destroy('raddrv  ',6,'cont',1,'              ',   &
     &                        14,iud,stat)
               end if

               call cpondpe
               call destroy('raddrv  ',6,'stop',1001,'              ',  &
     &                      14,101,stat)

            end if

            broot=stat-1

            call chkstd(broot)

! -----

! read file.

            call atmosin(iud,kln,zalt_data,patm_data,tklv_data,         &
     &                   cpcl_data,cgas_data,ccfc)

! -----

! close file.

            if(mype==root) then
              close(iud)
            end if

          end if

!! -----

!! calculte the zenith angle.

          call zenith(pdate,ni,nj,lat,lon,coseta)

!! -----

!! calculate the total pressure variable and Exner function.

          call getexner(ni,nj,nk,pbr,ppp,tmp1,p)

!! -----

!! calculate the air temperature.

          call getta3d(ni,nj,nk,ptbr,tmp1,ptpp,t)

!! -----

!! perform radiative transfer.

          do j=1,nj-1
          do i=1,ni-1

! set solar zenith angle.

            ams=coseta(i,j)

! -----

! set ground surface parameters.

            if(sfcopt.eq.0) then

              prg(1)=4.1e0
              prg(2)=0.e0

              gtmp=t(i,j,2)

            else

              call adjstlnd(i,j,aland,ni,nj,land)

              prg(1)=aland
              prg(2)=beta(i,j)

              gtmp=tundp(i,j,1)

            end if

!DBG        gtmp = 260.94d0    !! debug for mstrnX-CAOS

! -----

! set the one dimentional variables for radiation scheme.

            call cress21d(idcphopt,idradobj,idcldlim,idreliqum,         &
     &                    idreiceum,iradcnt,                            &
     &                    i,j,ni,nj,nk,nqw,nqi,zph,rbr,p,t,qvp,         &
     &                    qwtrp,qicep,                                  &
     &                    kln,zalt_data,cpcl_data,cgas_data,            &
     &                    zl,pl,tl,pb,tb,cpcl,gdcfrc,cgas,ccfc)

! -----

! perform radiative transfer.

            err=''

            call dtrn3(iug,iup,iuv,ams,nln,pb,pl,tb,tl,gtmp,cpcl,gdcfrc,&
     &                 cgas,ccfc,prg,fd,fu,err)

            if(err.ne.'') then

               write(6,*) 'error: ',err

               stop

            end if

! -----

! get heating rate.

            call radheat(i,j,ni,nj,nk,pb,cgas,                          &
     &                   fd,fu,htrsd,htrsu,htrld,htrlu,rsnet,rlnet)

! -----

          end do
          end do

!! -----

!! close files.

          if(mype.eq.root) then

            close(iug,status='keep')
            close(iup,status='keep')
            close(iuv,status='keep')

          end if

!! -----

!! counter

          iradcnt = iradcnt+1

!! -----

        end if

!!! -----

      else if( radopt==10 ) then

!!! Perform RRTM radiation scheme at marked time.

        if( trim(radon).eq.'sw' .or. trim(radon).eq.'lw'                &
     &                          .or. trim(radon).eq.'lsw' ) then

!ORG      if(mype==root) then
!ORG        write(*,*)
!ORG        write(*,*)' ### run RRTM radiation scheme ###' 
!ORG        write(*,*)
!ORG      end if

!! reading standard atmospheric condition file only once.

          if( iradcnt==0 ) then

! open file.

            if(mype==root) then

              parafl1(1:ncdat) = datdir(1:ncdat)
              write(parafl1(ncdat+1:ncdat+23),'(a23)')                  &
     &                                       'DataMSTRN/DATA.gas_ptcl'

              open(iud,iostat=stat,err=300,file=parafl1(1:ncdat+23),    &
     &             status='old',access='sequential',form='formatted',   &
     &             blank='null',position='rewind',action='read')

            else

              stat=0

            end if

  300       call chkerr(stat)

            if(mype.eq.root) then
              call outstd03('raddrv  ',6,parafl1,ncdat+23,iud,1,0,0_i8)
            end if

            if(stat<0) then

               if(mype.eq.-stat-1) then
                 call destroy('raddrv  ',6,'cont',1,'              ',   &
     &                        14,iud,stat)
               end if

               call cpondpe
               call destroy('raddrv  ',6,'stop',1001,'              ',  &
     &                      14,101,stat)

            end if

            broot=stat-1

            call chkstd(broot)

! -----

! read file.

            call atmosin(iud,kln,zalt_data,patm_data,tklv_data,         &
     &                   cpcl_data,cgas_data,ccfc)

! -----

! close file.

            if(mype==root) then
              close(iud)
            end if

          end if

!! -----

!! calculte the zenith angle.

          call zenith(pdate,ni,nj,lat,lon,coseta)

!! -----

!! calculate the total pressure variable and Exner function.

          call getexner(ni,nj,nk,pbr,ppp,tmp1,p)

!! -----

!! calculate the air temperature.

          call getta3d(ni,nj,nk,ptbr,tmp1,ptpp,t)

!! -----

!! run radiation scheme

          call rrtmdrv(radon,idcphopt,idradobj,idmcaopt,idovpopt,       &
     &                 idcldlim,idreliqum,idreiceum,iradcnt,            &
     &                 ni,nj,nk,nqw,nqi,nund,htrsd,htrsu,htrld,htrlu,   &
     &                 rsnet,rlnet,p,t,qvp,tundp,qwtrp,qicep,zph,       &
     &                 tmp1,albe,coseta,tmp2,                           &
     &                 kln,zalt_data,cgas_data,ccfc)

!! -----

!! counter

          iradcnt = iradcnt+1

        end if 

!!! -----

      end if

!!!! -----

! Heating/Cooling atmosphere by radiation scheme.

      if( radopt==1 .or. radopt==10 ) then

        call changept(dtb_sub,ni,nj,nk,htrsd,htrsu,htrld,htrlu,ptpf)

      end if

! -----

      end subroutine s_raddrv

!-----7--------------------------------------------------------------7--

      end module m_raddrv

!***********************************************************************
      module m_outpbl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/12/10
!     Modification: 2002/04/02, 2002/07/03, 2002/12/02, 2003/01/20,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/09/01, 2003/10/31, 2003/12/12, 2004/02/01,
!                   2004/03/05, 2004/04/01, 2004/05/07, 2004/05/31,
!                   2004/06/10, 2004/08/20, 2004/09/01, 2004/09/10,
!                   2004/09/25, 2005/01/31, 2006/01/10, 2006/02/13,
!                   2006/05/12, 2006/09/21, 2006/11/06, 2007/01/05,
!                   2007/01/20, 2007/05/14, 2007/05/21, 2007/06/27,
!                   2007/07/30, 2007/09/04, 2007/10/19, 2008/01/11,
!                   2008/03/12, 2008/04/17, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2008/10/10, 2009/02/27, 2009/08/20,
!                   2011/08/18, 2011/09/22, 2013/01/28, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the surface variables to the dumped file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bulksfc
      use m_comprofile
      use m_dump_kernel
      use m_comcapt
      use m_comdmp
      use m_comindx
      use m_comphy
      use m_comtable
      use m_getcname
      use m_getiname
      use m_inichar
      use m_outdmp2d
      use m_rotuvm2s
      use m_setcst2d
      use m_setproj

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outpbl, s_outpbl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outpbl

        module procedure s_outpbl

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outpbl(fpdmpvar,fpdmplev,fmois,ni,nj,nk,nund,        &
     &                    za,lon,p,u,v,qv,ufrc,vfrc,ptfrc,qvfrc,        &
     &                    land,kai,z0m,z0h,ptv,qvsfc,rch,cm,ch,         &
     &                    tund,tice,hs,le,rgd,rsd,rld,rlu,cdl,cdm,cdh,  &
     &                    z10,z15,cm10,ch15,u10,v10,p15,pt15,qv15,      &
     &                    tsfc,cdave,usflx,vsflx,ptsflx,qvsflx)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: fpdmplev
                       ! Formal parameter of unique index of dmplev

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: za(0:ni+1,0:nj+1)
                       ! z physical coordinates at lowest plane

      real, intent(in) :: lon(0:ni+1,0:nj+1)
                       ! Longitude

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(in) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(in) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in potential temperature equation

      real, intent(in) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(in) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(in) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(in) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(in) :: qvsfc(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio on surface

      real, intent(in) :: rch(0:ni+1,0:nj+1)
                       ! Bulk Richardson number

      real, intent(in) :: cm(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity

      real, intent(in) :: ch(0:ni+1,0:nj+1)
                       ! Bulk coefficient for scalar

      real, intent(in) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature

      real, intent(in) :: tice(0:ni+1,0:nj+1)
                       ! Mixed ice surface temperature

      real, intent(in) :: hs(0:ni+1,0:nj+1)
                       ! Sensible heat

      real, intent(in) :: le(0:ni+1,0:nj+1)
                       ! Latent heat

      real, intent(in) :: rgd(0:ni+1,0:nj+1)
                       ! Global solar radiation

      real, intent(in) :: rsd(0:ni+1,0:nj+1)
                       ! Net downward short wave radiation

      real, intent(in) :: rld(0:ni+1,0:nj+1)
                       ! Downward long wave radiation

      real, intent(in) :: rlu(0:ni+1,0:nj+1)
                       ! Upward long wave radiation

      real, intent(in) :: cdl(0:ni+1,0:nj+1)
                       ! Cloud cover in lower layer

      real, intent(in) :: cdm(0:ni+1,0:nj+1)
                       ! Cloud cover in middle layer

      real, intent(in) :: cdh(0:ni+1,0:nj+1)
                       ! Cloud cover in upper layer

! Internal shared variables

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      integer dmplev   ! Option for z coordinates of dumped variables

      real rddwkp      ! ln(da0 / dv0) / wkappa

      real cpj(1:7)    ! Map projection parameters

      real, intent(inout) :: z10(0:ni+1,0:nj+1)
                       ! Constant height of 10m

      real, intent(inout) :: z15(0:ni+1,0:nj+1)
                       ! Constant height of 1.5m

      real, intent(inout) :: cm10(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity
                       ! at an altitude of 10m

      real, intent(inout) :: ch15(0:ni+1,0:nj+1)
                       ! Bulk coefficient for scalar
                       ! at an altitude of 1.5m

      real, intent(inout) :: u10(0:ni+1,0:nj+1)
                       ! x components of velocity at an altitude of 10m

      real, intent(inout) :: v10(0:ni+1,0:nj+1)
                       ! y components of velocity at an altitude of 10m

      real, intent(inout) :: p15(0:ni+1,0:nj+1)
                       ! Pressure at an altitude of 1.5m

      real, intent(inout) :: pt15(0:ni+1,0:nj+1)
                       ! Potential temperature at an altitude of 1.5m

      real, intent(inout) :: qv15(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio at an altitude of 1.5m

      real, intent(inout) :: tsfc(0:ni+1,0:nj+1)
                       ! Soil and sea surface temperature

      real, intent(inout) :: cdave(0:ni+1,0:nj+1)
                       ! Averaged cloud cover

      real, intent(inout) :: usflx(0:ni+1,0:nj+1)
                       ! Surface stress in x-z coordinates

      real, intent(inout) :: vsflx(0:ni+1,0:nj+1)
                       ! Surface stress in y-z coordinates

      real, intent(inout) :: ptsflx(0:ni+1,0:nj+1)
                       ! Surface heat flux

      real, intent(inout) :: qvsflx(0:ni+1,0:nj+1)
                       ! Surface moisture flux

! Internal private variables

      integer i        ! Array index in x drection
      integer j        ! Array index in y drection

      real a           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_outpbl = 0
      integer, parameter :: DUMP_TARGET_outpbl = 4
      logical, save :: dump_done_outpbl = .false.


!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(dmpvar)

! -----

! Get the required namelist variable.

      call getcname(fpdmpvar,dmpvar)

! -----

!! Read in the surface predicted variables to the dumped file.

      if(fmon(1:3).eq.'act'.and.dmpvar(13:13).eq.'o') then

! Get the required namelist variable.

        call getiname(fpdmplev,dmplev)

! -----

! Set the common used variable.

        rddwkp=log(da0/dv0)/wkappa

! -----

! Set the constant height of 10.0 and 1.5 m.

        call setcst2d(0,ni+1,0,nj+1,10.e0,z10)
        call setcst2d(0,ni+1,0,nj+1,1.5e0,z15)

! -----

! Calculate the bulk coefficients.

        call bulksfc(ni,nj,z10,land,kai,z0m,z0h,rch,cm10,u10)
        call bulksfc(ni,nj,z15,land,kai,z0m,z0h,rch,p15,ch15)

! -----

! Calculate the x and y components of velocity at an altitude of 10 m,
! the pressure, potential temperature and water vapor mixing ratio at
! an altitude of 1.5 m, the surface temperature, total cloud cover and
! surface fluxes.

!@llm start meta_info ----------------------------------------------------
! Location: outpbl.f90 :: s_outpbl
! Summary : Compute 10m wind, 1.5m pressure/temperature/humidity, surface temp,
!           cloud cover, and surface fluxes for output diagnostics.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region (pure arithmetic)
!   - No global/module variable writes (only intent(inout) arrays)
!   - No sync constructs (barrier, critical, atomic)
!   - Conditional branches based on fmois (dry vs moist) with separate do loops
! Next:
!   - Convert to OpenACC or OpenACC with data directives for arrays
!   - Collapse nested i,j loops for better GPU occupancy
! Runtime:
!   - Calls: 4
!   - AvgLoops: 806.4K
!   - TotalTime: 0.001s (0.00%)
!   - AvgTime: 0.215ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('outpbl.f90', 's_outpbl', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_outpbl = dump_call_count_outpbl + 1
if (dump_call_count_outpbl == DUMP_TARGET_outpbl .and. .not. dump_done_outpbl) then
  call dump_init('outpbl')
  call dump_scalar_c('dmpvar', dmpvar)
  call dump_scalar_i('dmplev', dmplev)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('nund', nund)
  call dump_scalar_r('epsav', epsav)
  call dump_array_2d('za.bin', za, 0, ni+1, 0, nj+1)
  call dump_array_2d('lon.bin', lon, 0, ni+1, 0, nj+1)
  call dump_array_3d('p.bin', p, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ufrc.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vfrc.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptfrc.bin', ptfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvfrc.bin', qvfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d_int('land.bin', land, 0, ni+1, 0, nj+1)
  call dump_array_2d('kai.bin', kai, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0m.bin', z0m, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0h.bin', z0h, 0, ni+1, 0, nj+1)
  call dump_array_3d('ptv.bin', ptv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('qvsfc.bin', qvsfc, 0, ni+1, 0, nj+1)
  call dump_array_2d('rch.bin', rch, 0, ni+1, 0, nj+1)
  call dump_array_2d('cm.bin', cm, 0, ni+1, 0, nj+1)
  call dump_array_2d('ch.bin', ch, 0, ni+1, 0, nj+1)
  call dump_array_3d('tund.bin', tund, 0, ni+1, 0, nj+1, 1, nund)
  call dump_array_2d('tice.bin', tice, 0, ni+1, 0, nj+1)
  call dump_array_2d('hs.bin', hs, 0, ni+1, 0, nj+1)
  call dump_array_2d('le.bin', le, 0, ni+1, 0, nj+1)
  call dump_array_2d('rgd.bin', rgd, 0, ni+1, 0, nj+1)
  call dump_array_2d('rsd.bin', rsd, 0, ni+1, 0, nj+1)
  call dump_array_2d('rld.bin', rld, 0, ni+1, 0, nj+1)
  call dump_array_2d('rlu.bin', rlu, 0, ni+1, 0, nj+1)
  call dump_array_2d('cdl.bin', cdl, 0, ni+1, 0, nj+1)
  call dump_array_2d('cdm.bin', cdm, 0, ni+1, 0, nj+1)
  call dump_array_2d('cdh.bin', cdh, 0, ni+1, 0, nj+1)
  call dump_array_2d('z10_in.bin', z10, 0, ni+1, 0, nj+1)
  call dump_array_2d('z15_in.bin', z15, 0, ni+1, 0, nj+1)
  call dump_array_2d('cm10_in.bin', cm10, 0, ni+1, 0, nj+1)
  call dump_array_2d('ch15_in.bin', ch15, 0, ni+1, 0, nj+1)
  call dump_array_2d('u10_in.bin', u10, 0, ni+1, 0, nj+1)
  call dump_array_2d('v10_in.bin', v10, 0, ni+1, 0, nj+1)
  call dump_array_2d('p15_in.bin', p15, 0, ni+1, 0, nj+1)
  call dump_array_2d('pt15_in.bin', pt15, 0, ni+1, 0, nj+1)
  call dump_array_2d('qv15_in.bin', qv15, 0, ni+1, 0, nj+1)
  call dump_array_2d('tsfc_in.bin', tsfc, 0, ni+1, 0, nj+1)
  call dump_array_2d('cdave_in.bin', cdave, 0, ni+1, 0, nj+1)
  call dump_array_2d('usflx_in.bin', usflx, 0, ni+1, 0, nj+1)
  call dump_array_2d('vsflx_in.bin', vsflx, 0, ni+1, 0, nj+1)
  call dump_array_2d('ptsflx_in.bin', ptsflx, 0, ni+1, 0, nj+1)
  call dump_array_2d('qvsflx_in.bin', qvsflx, 0, ni+1, 0, nj+1)
  ! FIXME: cdrat is array - call dump_scalar_r('cdrat', cdrat)
  call dump_scalar_c('fmois', fmois)
  call dump_scalar_r('rddwkp', rddwkp)
end if

!$omp parallel default(shared)

        if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j,a)

          do j=1,nj-1
          do i=1,ni-1
            a=.5e0*cm(i,j)/cm10(i,j)

            u10(i,j)=a*(u(i,j,2)+u(i+1,j,2))
            v10(i,j)=a*(v(i,j,2)+v(i,j+1,2))

            a=ch(i,j)/ch15(i,j)

            pt15(i,j)=ptv(i,j,1)+a*(ptv(i,j,2)-ptv(i,j,1))
            qv15(i,j)=0.e0

            p15(i,j)=.5e0*(p(i,j,1)*(za(i,j)-1.5e0)                     &
     &        +p(i,j,2)*(za(i,j)+1.5e0))/za(i,j)

            if(land(i,j).eq.1) then
              tsfc(i,j)=kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1)
            else
              tsfc(i,j)=tund(i,j,1)
            end if

            cdave(i,j)=0.e0

            usflx(i,j)=-.5e0*(ufrc(i,j,1)+ufrc(i+1,j,1))
            vsflx(i,j)=-.5e0*(vfrc(i,j,1)+vfrc(i,j+1,1))

            ptsflx(i,j)=-ptfrc(i,j,1)

            qvsflx(i,j)=0.e0

          end do
          end do

!$omp end do

        else if(fmois(1:5).eq.'moist') then

!$omp do schedule(runtime) private(i,j,a)

          do j=1,nj-1
          do i=1,ni-1
            a=.5e0*cm(i,j)/cm10(i,j)

            u10(i,j)=a*(u(i,j,2)+u(i+1,j,2))
            v10(i,j)=a*(v(i,j,2)+v(i,j+1,2))

            a=ch(i,j)/ch15(i,j)

            pt15(i,j)=ptv(i,j,1)+a*(ptv(i,j,2)-ptv(i,j,1))

            if(land(i,j).lt.0) then

              qv15(i,j)=qvsfc(i,j)+a*(qv(i,j,2)-qvsfc(i,j))             &
     &          *((1.e0+rddwkp*ch15(i,j))/(1.e0+rddwkp*ch(i,j)))

            else

              qv15(i,j)=qvsfc(i,j)+a*(qv(i,j,2)-qvsfc(i,j))

            end if

            pt15(i,j)=pt15(i,j)*(1.e0+qv15(i,j))/(1.e0+epsav*qv15(i,j))

            p15(i,j)=.5e0*(p(i,j,1)*(za(i,j)-1.5e0)                     &
     &        +p(i,j,2)*(za(i,j)+1.5e0))/za(i,j)

            if(land(i,j).eq.1) then
              tsfc(i,j)=kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1)
            else
              tsfc(i,j)=tund(i,j,1)
            end if

            cdave(i,j)                                                  &
     &        =cdrat(1)*cdl(i,j)+cdrat(2)*cdm(i,j)+cdrat(3)*cdh(i,j)

            usflx(i,j)=-.5e0*(ufrc(i,j,1)+ufrc(i+1,j,1))
            vsflx(i,j)=-.5e0*(vfrc(i,j,1)+vfrc(i,j+1,1))

            ptsflx(i,j)=-ptfrc(i,j,1)
            qvsflx(i,j)=-qvfrc(i,j,1)

          end do
          end do

!$omp end do

        end if

!$omp end parallel

! Dump output data at target call
if (dump_call_count_outpbl == DUMP_TARGET_outpbl .and. .not. dump_done_outpbl) then
  call dump_array_2d('z10_ref.bin', z10, 0, ni+1, 0, nj+1)
  call dump_array_2d('z15_ref.bin', z15, 0, ni+1, 0, nj+1)
  call dump_array_2d('cm10_ref.bin', cm10, 0, ni+1, 0, nj+1)
  call dump_array_2d('ch15_ref.bin', ch15, 0, ni+1, 0, nj+1)
  call dump_array_2d('u10_ref.bin', u10, 0, ni+1, 0, nj+1)
  call dump_array_2d('v10_ref.bin', v10, 0, ni+1, 0, nj+1)
  call dump_array_2d('p15_ref.bin', p15, 0, ni+1, 0, nj+1)
  call dump_array_2d('pt15_ref.bin', pt15, 0, ni+1, 0, nj+1)
  call dump_array_2d('qv15_ref.bin', qv15, 0, ni+1, 0, nj+1)
  call dump_array_2d('tsfc_ref.bin', tsfc, 0, ni+1, 0, nj+1)
  call dump_array_2d('cdave_ref.bin', cdave, 0, ni+1, 0, nj+1)
  call dump_array_2d('usflx_ref.bin', usflx, 0, ni+1, 0, nj+1)
  call dump_array_2d('vsflx_ref.bin', vsflx, 0, ni+1, 0, nj+1)
  call dump_array_2d('ptsflx_ref.bin', ptsflx, 0, ni+1, 0, nj+1)
  call dump_array_2d('qvsflx_ref.bin', qvsflx, 0, ni+1, 0, nj+1)
  call dump_finalize()
  dump_done_outpbl = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

! Dump the surface monitor variables.

        if(dmplev.lt.10) then

          call outdmp2d('us    ',2,capt(48),ncpt(48),'xx',ni,nj,u10)
          call outdmp2d('vs    ',2,capt(49),ncpt(49),'xx',ni,nj,v10)

        else

          call setproj(idmpopt,idnspol,idtlat1,idtlat2,cpj)

          call s_rotuvm2s(idmpopt,idnspol,idtlon,1,ni-1,1,nj-1,1,1,cpj, &
     &                    0,ni+1,0,nj+1,1,1,lon,u10,v10)

          call outdmp2d('us    ',2,capt(53),ncpt(53),'xx',ni,nj,u10)
          call outdmp2d('vs    ',2,capt(54),ncpt(54),'xx',ni,nj,v10)

        end if

        call outdmp2d('ps    ',2,capt(50),ncpt(50),'xx',ni,nj,p15)
        call outdmp2d('pts   ',3,capt(51),ncpt(51),'xx',ni,nj,pt15)
        call outdmp2d('qvs   ',3,capt(52),ncpt(52),'xx',ni,nj,qv15)

        call outdmp2d('tgs   ',3,capt(55),ncpt(55),'xx',ni,nj,tsfc)

        call outdmp2d('hs    ',2,capt(56),ncpt(56),'xx',ni,nj,hs)
        call outdmp2d('le    ',2,capt(57),ncpt(57),'xx',ni,nj,le)

        call outdmp2d('rgd   ',3,capt(78),ncpt(78),'xx',ni,nj,rgd)

        call outdmp2d('rsd   ',3,capt(58),ncpt(58),'xx',ni,nj,rsd)
        call outdmp2d('rld   ',3,capt(59),ncpt(59),'xx',ni,nj,rld)
        call outdmp2d('rlu   ',3,capt(60),ncpt(60),'xx',ni,nj,rlu)

        call outdmp2d('cdl   ',3,capt(79),ncpt(79),'xx',ni,nj,cdl)
        call outdmp2d('cdm   ',3,capt(80),ncpt(80),'xx',ni,nj,cdm)
        call outdmp2d('cdh   ',3,capt(81),ncpt(81),'xx',ni,nj,cdh)

        call outdmp2d('cdave ',5,capt(61),ncpt(61),'xx',ni,nj,cdave)

        call outdmp2d('usflx ',5,capt(62),ncpt(62),'xx',ni,nj,usflx)
        call outdmp2d('vsflx ',5,capt(63),ncpt(63),'xx',ni,nj,vsflx)
        call outdmp2d('ptsflx',6,capt(64),ncpt(64),'xx',ni,nj,ptsflx)
        call outdmp2d('qvsflx',6,capt(65),ncpt(65),'xx',ni,nj,qvsflx)

! -----

      end if

!! -----

      end subroutine s_outpbl

!-----7--------------------------------------------------------------7--

      end module m_outpbl

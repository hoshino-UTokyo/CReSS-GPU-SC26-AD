!***********************************************************************
      module m_forcept
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/20,
!                   1999/06/07, 1999/07/05, 1999/07/21, 1999/08/03,
!                   1999/08/09, 1999/08/18, 1999/08/23, 1999/09/16,
!                   1999/09/30, 1999/10/07, 1999/10/12, 1999/10/22,
!                   1999/11/01, 1999/11/19, 1999/11/24, 2000/01/17,
!                   2000/02/07, 2000/04/18, 2000/12/19, 2001/03/13,
!                   2001/04/15, 2001/05/29, 2001/06/06, 2001/07/13,
!                   2001/08/07, 2001/11/20, 2002/04/02, 2002/06/18,
!                   2002/07/23, 2002/08/15, 2002/09/09, 2002/10/31,
!                   2002/12/11, 2003/01/04, 2003/01/20, 2003/03/13,
!                   2003/03/21, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/09/01, 2003/10/10, 2003/11/28, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/15, 2004/05/31,
!                   2004/06/10, 2004/08/01, 2004/08/20, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/06/21,
!                   2006/09/21, 2006/11/06, 2007/05/07, 2007/07/30,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing term in the potential temperature equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advbspt
      use m_comprofile
      use m_dump_kernel
      use m_advs
      use m_comindx
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lsps
      use m_nlsms
      use m_s2gpv
      use m_smoo2s
      use m_smoo4s
      use m_turbflx
      use m_turbs
      use m_vsps

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forcept, s_forcept

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forcept

        module procedure s_forcept

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_forcept(fpnggvar,fplspvar,fpvspvar,                  &
     &                     fplspopt,fpvspopt,fpgwmopt,fpsmtopt,         &
     &                     fpadvopt,fptubopt,ksp0,nggdmp,gtinc,         &
     &                     ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,mf,rmf,     &
     &                     rmf8u,rmf8v,ptbr,rbr,rst,rstxu,rstxv,rstxwc, &
     &                     w,wp,ptp,ptpp,rkh8u,rkh8v,rkv8w,rbcxy,rbct,  &
     &                     ptpgpv,ptptd,ptfrc,pt,h3,tmp1,tmp2,          &
     &                     tmp3,tmp4,tmp5)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpnggvar
                       ! Formal parameter of unique index of nggvar

      integer, intent(in) :: fplspvar
                       ! Formal parameter of unique index of lspvar

      integer, intent(in) :: fpvspvar
                       ! Formal parameter of unique index of vspvar

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rstxu(0:ni+1,0:nj+1,1:nk)
                       ! u x base state density x Jacobian

      real, intent(in) :: rstxv(0:ni+1,0:nj+1,1:nk)
                       ! v x base state density x Jacobian

      real, intent(in) :: rstxwc(0:ni+1,0:nj+1,1:nk)
                       ! wc x base state density x Jacobian

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: rkh8u(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity / jcb
                       ! at u points

      real, intent(in) :: rkh8v(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity / jcb
                       ! at v points

      real, intent(in) :: rkv8w(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy diffusivity / jcb
                       ! at w points

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(in) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

! Input and output variable

      real, intent(inout) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in potential temperature equation

! Internal shared variables

      character(len=108) nggvar
                       ! Control flag of
                       ! analysis nudged variables to GPV

      character(len=108) lspvar
                       ! Control flag of
                       ! lateral sponge damped variables

      character(len=108) vspvar
                       ! Control flag of
                       ! vertical sponge damped variables

      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping
      integer gwmopt   ! Option for gravity wave mode integration

      integer smtopt   ! Option for numerical smoothing
      integer advopt   ! Option for advection scheme
      integer tubopt   ! Option for turbulent mixing

      real, intent(inout) :: pt(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature

      real, intent(inout) :: h3(0:ni+1,0:nj+1,1:nk)
                       ! z components of turbulent fluxes

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_forcept = 0
      integer, parameter :: DUMP_TARGET_forcept = 360
      logical, save :: dump_done_forcept = .false.


!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(nggvar)
      call inichar(lspvar)
      call inichar(vspvar)

! -----

! Get the required namelist variables.

      call getcname(fpnggvar,nggvar)
      call getcname(fplspvar,lspvar)
      call getcname(fpvspvar,vspvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpgwmopt,gwmopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fpadvopt,advopt)
      call getiname(fptubopt,tubopt)

! -----

! Get the potential temperature.

      if(tubopt.ge.1) then

!@llm start meta_info ----------------------------------------------------
! Location: forcept.f90 :: s_forcept
! Summary : Compute potential temperature pt by adding base state ptbr
!           and perturbation ptpp for turbulent mixing calculation.
! GPU diff: Easy
! Findings:
!   - Simple element-wise addition of two arrays
!   - No function calls inside parallel region
!   - No reductions or synchronization
!   - Only writes to pt array
! Next:
!   - Direct GPU kernel port with straightforward 3D mapping
!   - Consider fusing with subsequent turbulent mixing kernels
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 1.361s (0.05%)
!   - AvgTime: 3.780ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('forcept.f90', 's_forcept', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_forcept = dump_call_count_forcept + 1
if (dump_call_count_forcept == DUMP_TARGET_forcept .and. .not. dump_done_forcept) then
  call dump_init('forcept')
  call dump_scalar_c('nggvar', nggvar)
  call dump_scalar_c('lspvar', lspvar)
  call dump_scalar_c('vspvar', vspvar)
  call dump_scalar_i('lspopt', lspopt)
  call dump_scalar_i('vspopt', vspopt)
  call dump_scalar_i('gwmopt', gwmopt)
  call dump_scalar_i('smtopt', smtopt)
  call dump_scalar_i('advopt', advopt)
  call dump_scalar_i('tubopt', tubopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('nggdmp', nggdmp)
  call dump_scalar_r('gtinc', gtinc)
  call dump_array_3d('pt.bin', pt, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('jcb8u.bin', jcb8u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('jcb8v.bin', jcb8v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('mf.bin', mf, 0, ni+1, 0, nj+1)
  call dump_array_3d('rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call dump_array_3d('rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rstxu.bin', rstxu, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rstxv.bin', rstxv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rstxwc.bin', rstxwc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wp.bin', wp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptpp.bin', ptpp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rkh8u.bin', rkh8u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rkh8v.bin', rkh8v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rkv8w.bin', rkv8w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('rbcxy.bin', rbcxy, 1, ni, 1, nj)
  call dump_array_4d('rbct.bin', rbct, 1, ni, 1, nj, 1, nk, 1, 2)
  call dump_array_3d('ptpgpv.bin', ptpgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptptd.bin', ptptd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptfrc_in.bin', ptfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('h3_in.bin', h3, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_in.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp2_in.bin', tmp2, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp3_in.bin', tmp3, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp4_in.bin', tmp4, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp5_in.bin', tmp5, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            pt(i,j,k)=ptbr(i,j,k)+ptpp(i,j,k)
          end do
          end do

!$omp end do

        end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_forcept == DUMP_TARGET_forcept .and. .not. dump_done_forcept) then
  call dump_array_3d('ptfrc_ref.bin', ptfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('h3_ref.bin', h3, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp1_ref.bin', tmp1, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp2_ref.bin', tmp2, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp3_ref.bin', tmp3, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp4_ref.bin', tmp4, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tmp5_ref.bin', tmp5, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_forcept = .true.
end if


call profile_stop(prof_id1, loop_len)

      end if

! -----

! Calculate the advection.

      call advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,             &
     &          iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,       &
     &          ptp,ptfrc,tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

      if(mod(smtopt,10).eq.1) then

        call smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,ptpp,ptfrc,tmp1)

      end if

! -----

! Calculate the 4th order smoothing.

      if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

        call smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              idsmhcoe,idsmvcoe,ni,nj,nk,rbr,ptpp,ptfrc,          &
     &              tmp1,tmp2,tmp3,tmp4,tmp5)

      end if

! -----

! Calculate the non linear smoothing.

      if(smtopt.ge.11) then

        call nlsms(idnlhcoe,idnlvcoe,1.e0,ni,nj,nk,rbr,ptpp,ptfrc,      &
     &             tmp1,tmp2,tmp3,tmp4)

      end if

! -----

! Calculate the turbulent mixing.

      if(tubopt.ge.1) then

        call turbflx(idtrnopt,idsfcopt,iddxiv,iddyiv,iddziv,ni,nj,nk,   &
     &               j31,j32,jcb,pt,ptfrc,rkh8u,rkh8v,rkv8w,            &
     &               tmp1,tmp2,h3,tmp3,tmp4,tmp5)

        call turbs(idtrnopt,idmpopt,idmfcopt,iddxiv,iddyiv,iddziv,      &
     &             ni,nj,nk,j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,     &
     &             tmp1,tmp2,h3,ptfrc,tmp3,tmp4,tmp5)

      end if

! -----

! Perform the analysis nudging to GPV.

      if(nggdmp.gt.0.e0.and.nggvar(5:5).eq.'o') then

        call s2gpv(idgpvvar,9,nggdmp,gtinc,                             &
     &             ni,nj,nk,rst,ptpp,ptpgpv,ptptd,ptfrc)

      end if

! -----

! Calculate the lateral sponge damping.

      if(lspopt.ge.1.and.lspvar(5:5).eq.'o') then

        call lsps(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,9,gtinc, &
     &            ni,nj,nk,rst,ptpp,rbcxy,ptpgpv,ptptd,ptfrc,tmp1)

      end if

! -----

! Calculate the vertical sponge damping.

      if(vspopt.ge.1.and.vspvar(5:5).eq.'o') then

        call vsps(idgpvvar,idvspopt,9,ksp0,gtinc,ni,nj,nk,rst,ptpp,rbct,&
     &            ptpgpv,ptptd,ptfrc)

      end if

! -----

! Calculate the base state advection.

      if(gwmopt.eq.0) then

        if(advopt.le.3) then

          call advbspt(idgwmopt,iddziv,ni,nj,nk,ptbr,rbr,w,ptfrc,tmp1)

        else

          call advbspt(idgwmopt,iddziv,ni,nj,nk,ptbr,rbr,wp,ptfrc,tmp1)

        end if

      end if

! -----

      end subroutine s_forcept

!-----7--------------------------------------------------------------7--

      end module m_forcept

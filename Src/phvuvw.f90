!***********************************************************************
      module m_phvuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/09/13
!     Modification: 2001/12/11, 2002/04/02, 2002/06/06, 2002/08/15,
!                   2002/10/31, 2003/01/04, 2003/02/13, 2003/03/13,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/07/28,
!                   2003/10/31, 2003/11/05, 2003/11/28, 2003/12/12,
!                   2004/04/01, 2004/04/15, 2004/08/01, 2004/08/20,
!                   2004/09/01, 2006/04/03, 2006/09/21, 2006/11/06,
!                   2006/12/04, 2007/01/05, 2007/01/31, 2007/05/21,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2009/11/13, 2010/12/01, 2011/01/19,
!                   2011/08/09, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the velocity phase speed for the open boundary
!     conditions.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_dump_kernel
      use m_commpi
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: phvuvw, s_phvuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface phvuvw

        module procedure s_phvuvw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max
      intrinsic min
      intrinsic mod
      intrinsic real
      intrinsic sign

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_phvuvw(fpexbvar,fpexbopt,                            &
     &                    fpwbc,fpebc,fpsbc,fpnbc,fpmpopt,fpmfcopt,     &
     &                    fpdxiv,fpdyiv,fpgwave,dtb,dts,ni,nj,nk,       &
     &                    rmf,rmf8u,rmf8v,u,up,uf,v,vp,vf,w,wp,wf,      &
     &                    ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,u8v,v8u,        &
     &                    cpavex,cpavey)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: fpgwave
                       ! Formal parameter of unique index of gwave

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

! Output variables

      real, intent(out) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(out) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(out) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(out) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, intent(out) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(out) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer exbopt   ! Option for external boundary forcing

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer nim3     ! ni - 3

      integer njm1     ! nj - 1
      integer njm2     ! nj - 2
      integer njm3     ! nj - 3

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dx

      real gwave       ! Fastest gravity wave speed

      real gdxdt       ! gwave x dxiv x dtb
      real gdydt       ! gwave x dyiv x dtb

      real gdxdtn      ! - gwave x dxiv x dtb
      real gdydtn      ! - gwave x dyiv x dtb

      real dxdt        ! dxiv x dtb
      real dydt        ! dyiv x dtb

      real dxdt5       ! 0.5 x dxiv x dtb
      real dydt5       ! 0.5 x dyiv x dtb

      real dtsdb       ! dts / dtb

      real nkm2v       ! 1.0 / real(nk - 2)
      real nkm3v       ! 1.0 / real(nk - 3)

      real, intent(inout) :: u8v(0:nj+1,1:nk)
                       ! x components of velocity at v points

      real, intent(inout) :: v8u(0:ni+1,1:nk)
                       ! y components of velocity at u points

      real, intent(inout) :: cpavex(0:nj+1)
                       ! Vertically averaged phase speed in x direction

      real, intent(inout) :: cpavey(0:ni+1)
                       ! Vertically averaged phase speed in y direction

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_phvuvw = 0
      integer, parameter :: DUMP_TARGET_phvuvw = 360
      logical, save :: dump_done_phvuvw = .false.


!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpexbopt,exbopt)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpgwave,gwave)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      nim3=ni-3

      njm1=nj-1
      njm2=nj-2
      njm3=nj-3

      if(ebw.eq.1.and.isub.eq.0) then
        istr=2
      else
        istr=1
      end if

      if(ebe.eq.1.and.isub.eq.nisub-1) then
        iend=ni-1
      else
        iend=ni
      end if

      if(ebs.eq.1.and.jsub.eq.0) then
        jstr=2
      else
        jstr=1
      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then
        jend=nj-1
      else
        jend=nj
      end if

      gdxdt=gwave*dtb*dxiv
      gdydt=gwave*dtb*dyiv

      gdxdtn=-gwave*dtb*dxiv
      gdydtn=-gwave*dtb*dyiv

      dxdt=dtb*dxiv
      dydt=dtb*dyiv

      dxdt5=.5e0*dtb*dxiv
      dydt5=.5e0*dtb*dyiv

      dtsdb=dts/dtb

      nkm2v=1.e0/real(nk-2)
      nkm3v=1.e0/real(nk-3)

! -----

!!! Calculate the velocity phase speed for the open boundary conditions.

!@llm start meta_info ----------------------------------------------------
! Location: phvuvw.f90 :: s_phvuvw
! Summary : Calculate velocity phase speed (u,v,w) for open boundary
!           conditions on all four boundaries (west, east, south, north).
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region (only intrinsic abs, sign, min, max, mod)
!   - Writes to output arrays ucpx, ucpy, vcpx, vcpy, wcpx, wcpy
!   - Uses shared work arrays cpavex, cpavey, u8v, v8u for vertical averaging
!   - Complex conditional branching based on boundary condition options
!   - Multiple sequential k-loops with workshared inner j/i loops
!   - No explicit barriers but implicit at !$omp end do
! Next:
!   - Collapse nested loops where possible for better GPU occupancy
!   - Consider using OpenACC kernels directive with appropriate private clauses
!   - Boundary-only computation may benefit from separate small kernels
! Runtime:
!   - Calls: 360
!   - AvgLoops: 112.2K
!   - TotalTime: 5.198s (0.17%)
!   - AvgTime: 14.440ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('phvuvw.f90', 's_phvuvw', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) * int((nj-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_phvuvw = dump_call_count_phvuvw + 1
if (dump_call_count_phvuvw == DUMP_TARGET_phvuvw .and. .not. dump_done_phvuvw) then
  call dump_init('phvuvw')
  call dump_scalar_c('exbvar', exbvar)
  call dump_scalar_i('exbopt', exbopt)
  call dump_scalar_i('wbc', wbc)
  call dump_scalar_i('ebc', ebc)
  call dump_scalar_i('sbc', sbc)
  call dump_scalar_i('nbc', nbc)
  call dump_scalar_i('mpopt', mpopt)
  call dump_scalar_i('mfcopt', mfcopt)
  call dump_scalar_r('dxiv', dxiv)
  call dump_scalar_r('dyiv', dyiv)
  call dump_scalar_r('gwave', gwave)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('dtb', dtb)
  call dump_scalar_r('dts', dts)
  call dump_array_3d('rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call dump_array_3d('rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('up.bin', up, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('uf.bin', uf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vp.bin', vp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vf.bin', vf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wp.bin', wp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wf.bin', wf, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('u8v_in.bin', u8v, 0, nj+1, 1, nk)
  call dump_array_2d('v8u_in.bin', v8u, 0, ni+1, 1, nk)
  ! FIXME: cpavex is an array, not scalar
  ! ! FIXME: cpavex is array - call dump_scalar_r('cpavex', cpavex)
  ! FIXME: cpavey is an array, not scalar
  ! ! FIXME: cpavey is array - call dump_scalar_r('cpavey', cpavey)
  call dump_scalar_r('dtsdb', dtsdb)
  call dump_scalar_r('dxdt', dxdt)
  call dump_scalar_r('dxdt5', dxdt5)
  call dump_scalar_r('dydt', dydt)
  call dump_scalar_r('dydt5', dydt5)
  call dump_scalar_i('ebe', ebe)
  call dump_scalar_i('ebn', ebn)
  call dump_scalar_i('ebs', ebs)
  call dump_scalar_i('ebw', ebw)
  call dump_scalar_r('gdxdt', gdxdt)
  call dump_scalar_r('gdxdtn', gdxdtn)
  call dump_scalar_r('gdydt', gdydt)
  call dump_scalar_r('gdydtn', gdydtn)
  call dump_scalar_i('iend', iend)
  call dump_scalar_i('istr', istr)
  call dump_scalar_i('isub', isub)
  call dump_scalar_i('jend', jend)
  call dump_scalar_i('jstr', jstr)
  call dump_scalar_i('jsub', jsub)
  call dump_scalar_i('nim1', nim1)
  call dump_scalar_i('nim2', nim2)
  call dump_scalar_i('nim3', nim3)
  call dump_scalar_i('nisub', nisub)
  call dump_scalar_i('njm1', njm1)
  call dump_scalar_i('njm2', njm2)
  call dump_scalar_i('njm3', njm3)
  call dump_scalar_i('njsub', njsub)
  call dump_scalar_r('nkm2v', nkm2v)
  call dump_scalar_r('nkm3v', nkm3v)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_235)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    ! For wbc=7, set ucpx to constant gdxdtn on west boundary
    if (wbc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do j = 1, nj-1
          ucpx(j, k, 1) = gdxdtn
        end do
      end do
      !$acc end kernels
    end if

    ! For ebc=7, set ucpx to constant gdxdt on east boundary
    if (ebc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do j = 1, nj-1
          ucpx(j, k, 2) = gdxdt
        end do
      end do
      !$acc end kernels
    end if

    ! For sbc=7, set ucpy to constant gdydtn on south boundary
    if (sbc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do i = 1, ni-1
          ucpy(i, k, 1) = gdydtn
        end do
      end do
      !$acc end kernels
    end if

    ! For nbc=7, set ucpy to constant gdydt on north boundary
    if (nbc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do i = 1, ni-1
          ucpy(i, k, 2) = gdydt
        end do
      end do
      !$acc end kernels
    end if

    ! v component phase speeds
    if (wbc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do j = 1, nj
          vcpx(j, k, 1) = gdxdtn
        end do
      end do
      !$acc end kernels
    end if

    if (ebc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do j = 1, nj
          vcpx(j, k, 2) = gdxdt
        end do
      end do
      !$acc end kernels
    end if

    if (sbc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do i = 1, ni-1
          vcpy(i, k, 1) = gdydtn
        end do
      end do
      !$acc end kernels
    end if

    if (nbc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do i = 1, ni-1
          vcpy(i, k, 2) = gdydt
        end do
      end do
      !$acc end kernels
    end if

    ! w component phase speeds
    if (wbc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do j = 1, nj-1
          wcpx(j, k, 1) = gdxdtn
        end do
      end do
      !$acc end kernels
    end if

    if (ebc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do j = 1, nj-1
          wcpx(j, k, 2) = gdxdt
        end do
      end do
      !$acc end kernels
    end if

    if (sbc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do i = 1, ni-1
          wcpy(i, k, 1) = gdydtn
        end do
      end do
      !$acc end kernels
    end if

    if (nbc == 7) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do i = 1, ni-1
          wcpy(i, k, 2) = gdydt
        end do
      end do
      !$acc end kernels
    end if

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared) private(k)

!! Calculate the u phase speed.

! Calculate the u phase speed on the west boundary.

      if((ebw.eq.1.and.isub.eq.0).and.((wbc.ge.4.and.exbopt.eq.0)       &
     &  .or.(wbc.ge.4.and.exbopt.ge.1.and.exbvar(1:1).eq.'x'))) then

        if(mod(wbc,10).eq.4.or.mod(wbc,10).eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,1)=uf(2,j,k)+up(2,j,k)-2.e0*u(3,j,k)

              if(abs(ucpx(j,k,1)).lt.eps) then

                ucpx(j,k,1)=sign(eps,ucpx(j,k,1))

              end if

              ucpx(j,k,1)=min((uf(2,j,k)-up(2,j,k))/ucpx(j,k,1),gdxdtn)

            end do

!$omp end do

          end do

          if(mod(wbc,10).eq.5) then

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              cpavex(j)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=cpavex(j)+ucpx(j,k,1)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,1)=cpavex(j)
              end do

!$omp end do

            end do

          end if

        else if(mod(wbc,10).eq.6) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,1)=min(u(2,j,k)*dxdt,gdxdtn)
            end do

!$omp end do

          end do

        else if(wbc.eq.7) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,1)=gdxdtn
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,1)=max(ucpx(j,k,1),-rmf8u(2,j,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,1)=max(ucpx(j,k,1),-1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(j)

          do j=1,nj-1
            ucpx(j,k,1)=ucpx(j,k,1)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the u phase speed on the east boundary.

      if((ebe.eq.1.and.isub.eq.nisub-1).and.((ebc.ge.4.and.exbopt.eq.0) &
     &  .or.(ebc.ge.4.and.exbopt.ge.1.and.exbvar(1:1).eq.'x'))) then

        if(mod(ebc,10).eq.4.or.mod(ebc,10).eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,2)=2.e0*u(nim2,j,k)-uf(nim1,j,k)-up(nim1,j,k)

              if(abs(ucpx(j,k,2)).lt.eps) then

                ucpx(j,k,2)=sign(eps,ucpx(j,k,2))

              end if

              ucpx(j,k,2)                                               &
     &          =max((uf(nim1,j,k)-up(nim1,j,k))/ucpx(j,k,2),gdxdt)

            end do

!$omp end do

          end do

          if(mod(ebc,10).eq.5) then

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              cpavex(j)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=cpavex(j)+ucpx(j,k,2)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,2)=cpavex(j)
              end do

!$omp end do

            end do

          end if

        else if(mod(ebc,10).eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,2)=max(u(nim1,j,k)*dxdt,gdxdt)
            end do

!$omp end do

          end do

        else if(ebc.eq.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,2)=gdxdt
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,2)=min(ucpx(j,k,2),rmf8u(nim1,j,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,2)=min(ucpx(j,k,2),1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(j)

          do j=1,nj-1
            ucpx(j,k,2)=ucpx(j,k,2)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the u phase speed on the south boundary.

      if((ebs.eq.1.and.jsub.eq.0).and.((sbc.ge.4.and.exbopt.eq.0)       &
     &  .or.(sbc.ge.4.and.exbopt.ge.1.and.exbvar(1:1).eq.'x'))) then

        if(sbc.eq.4.or.sbc.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,1)=uf(i,2,k)+up(i,2,k)-2.e0*u(i,3,k)

              if(abs(ucpy(i,k,1)).lt.eps) then

                ucpy(i,k,1)=sign(eps,ucpy(i,k,1))

              end if

              ucpy(i,k,1)=min((uf(i,2,k)-up(i,2,k))/ucpy(i,k,1),gdydtn)

            end do

!$omp end do

          end do

          if(sbc.eq.5) then

!$omp do schedule(runtime) private(i)

            do i=1,ni
              cpavey(i)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                cpavey(i)=cpavey(i)+ucpy(i,k,1)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,1)=cpavey(i)
              end do

!$omp end do

            end do

          end if

        else if(sbc.eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=istr,iend
              v8u(i,k)=v(i-1,2,k)+v(i,2,k)
            end do

!$omp end do

          end do

          if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime)

            do k=2,nk-2
              v8u(1,k)=v8u(2,k)
            end do

!$omp end do

          else if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime)

            do k=2,nk-2
              v8u(ni,k)=v8u(ni-1,k)
            end do

!$omp end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,1)=min(v8u(i,k)*dydt5,gdydtn)
            end do

!$omp end do

          end do

        else if(sbc.ge.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,1)=gdydtn
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,1)=max(ucpy(i,k,1),-rmf8u(i,2,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,1)=max(ucpy(i,k,1),-1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(i)

          do i=1,ni
            ucpy(i,k,1)=ucpy(i,k,1)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the u phase speed on the north boundary.

      if((ebn.eq.1.and.jsub.eq.njsub-1).and.((nbc.ge.4.and.exbopt.eq.0) &
     &  .or.(nbc.ge.4.and.exbopt.ge.1.and.exbvar(1:1).eq.'x'))) then

        if(nbc.eq.4.or.nbc.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,2)=2.e0*u(i,njm3,k)-uf(i,njm2,k)-up(i,njm2,k)

              if(abs(ucpy(i,k,2)).lt.eps) then

                ucpy(i,k,2)=sign(eps,ucpy(i,k,2))

              end if

              ucpy(i,k,2)                                               &
     &          =max((uf(i,njm2,k)-up(i,njm2,k))/ucpy(i,k,2),gdydt)

            end do

!$omp end do

          end do

          if(nbc.eq.5) then

!$omp do schedule(runtime) private(i)

            do i=1,ni
              cpavey(i)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                cpavey(i)=cpavey(i)+ucpy(i,k,2)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,2)=cpavey(i)
              end do

!$omp end do

            end do

          end if

        else if(nbc.eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=istr,iend
              v8u(i,k)=v(i-1,njm1,k)+v(i,njm1,k)
            end do

!$omp end do

          end do

          if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime)

            do k=2,nk-2
              v8u(1,k)=v8u(2,k)
            end do

!$omp end do

          else if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime)

            do k=2,nk-2
              v8u(ni,k)=v8u(ni-1,k)
            end do

!$omp end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,2)=max(v8u(i,k)*dydt5,gdydt)
            end do

!$omp end do

          end do

        else if(nbc.ge.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,2)=gdydt
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,2)=min(ucpy(i,k,2),rmf8u(i,njm2,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,2)=min(ucpy(i,k,2),1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(i)

          do i=1,ni
            ucpy(i,k,2)=ucpy(i,k,2)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

!! -----

!! Calculate the v phase speed.

! Calculate the v phase speed on the west boundary.

      if((ebw.eq.1.and.isub.eq.0).and.((wbc.ge.4.and.exbopt.eq.0)       &
     &  .or.(wbc.ge.4.and.exbopt.ge.1.and.exbvar(2:2).eq.'x'))) then

        if(wbc.eq.4.or.wbc.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,1)=vf(2,j,k)+vp(2,j,k)-2.e0*v(3,j,k)

              if(abs(vcpx(j,k,1)).lt.eps) then

                vcpx(j,k,1)=sign(eps,vcpx(j,k,1))

              end if

              vcpx(j,k,1)=min((vf(2,j,k)-vp(2,j,k))/vcpx(j,k,1),gdxdtn)

            end do

!$omp end do

          end do

          if(wbc.eq.5) then

!$omp do schedule(runtime) private(j)

            do j=1,nj
              cpavex(j)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                cpavex(j)=cpavex(j)+vcpx(j,k,1)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,1)=cpavex(j)
              end do

!$omp end do

            end do

          end if

        else if(wbc.eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=jstr,jend
              u8v(j,k)=u(2,j-1,k)+u(2,j,k)
            end do

!$omp end do

          end do

          if(ebs.eq.1.and.jsub.eq.0) then

!$omp do schedule(runtime)

            do k=2,nk-2
              u8v(1,k)=u8v(2,k)
            end do

!$omp end do

          else if(ebn.eq.1.and.jsub.eq.njsub-1) then

!$omp do schedule(runtime)

            do k=2,nk-2
              u8v(nj,k)=u8v(nj-1,k)
            end do

!$omp end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,1)=min(u8v(j,k)*dxdt5,gdxdtn)
            end do

!$omp end do

          end do

        else if(wbc.ge.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,1)=gdxdtn
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,1)=max(vcpx(j,k,1),-rmf8v(2,j,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,1)=max(vcpx(j,k,1),-1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(j)

          do j=1,nj
            vcpx(j,k,1)=vcpx(j,k,1)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the v phase speed on the east boundary.

      if((ebe.eq.1.and.isub.eq.nisub-1).and.((ebc.ge.4.and.exbopt.eq.0) &
     &  .or.(ebc.ge.4.and.exbopt.ge.1.and.exbvar(2:2).eq.'x'))) then

        if(ebc.eq.4.or.ebc.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,2)=2.e0*v(nim3,j,k)-vf(nim2,j,k)-vp(nim2,j,k)

              if(abs(vcpx(j,k,2)).lt.eps) then

                vcpx(j,k,2)=sign(eps,vcpx(j,k,2))

              end if

              vcpx(j,k,2)                                               &
     &          =max((vf(nim2,j,k)-vp(nim2,j,k))/vcpx(j,k,2),gdxdt)

            end do

!$omp end do

          end do

          if(ebc.eq.5) then

!$omp do schedule(runtime) private(j)

            do j=1,nj
              cpavex(j)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                cpavex(j)=cpavex(j)+vcpx(j,k,2)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,2)=cpavex(j)
              end do

!$omp end do

            end do

          end if

        else if(ebc.eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=jstr,jend
              u8v(j,k)=u(nim1,j-1,k)+u(nim1,j,k)
            end do

!$omp end do

          end do

          if(ebs.eq.1.and.jsub.eq.0) then

!$omp do schedule(runtime)

            do k=2,nk-2
              u8v(1,k)=u8v(2,k)
            end do

!$omp end do

          else if(ebn.eq.1.and.jsub.eq.njsub-1) then

!$omp do schedule(runtime)

            do k=2,nk-2
              u8v(nj,k)=u8v(nj-1,k)
            end do

!$omp end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,2)=max(u8v(j,k)*dxdt5,gdxdt)
            end do

!$omp end do

          end do

        else if(ebc.ge.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,2)=gdxdt
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,2)=min(vcpx(j,k,2),rmf8v(nim2,j,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,2)=min(vcpx(j,k,2),1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(j)

          do j=1,nj
            vcpx(j,k,2)=vcpx(j,k,2)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the v phase speed on the south boundary.

      if((ebs.eq.1.and.jsub.eq.0).and.((sbc.ge.4.and.exbopt.eq.0)       &
     &  .or.(sbc.ge.4.and.exbopt.ge.1.and.exbvar(2:2).eq.'x'))) then

        if(mod(sbc,10).eq.4.or.mod(sbc,10).eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,1)=vf(i,2,k)+vp(i,2,k)-2.e0*v(i,3,k)

              if(abs(vcpy(i,k,1)).lt.eps) then

                vcpy(i,k,1)=sign(eps,vcpy(i,k,1))

              end if

              vcpy(i,k,1)=min((vf(i,2,k)-vp(i,2,k))/vcpy(i,k,1),gdydtn)

            end do

!$omp end do

          end do

          if(mod(sbc,10).eq.5) then

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              cpavey(i)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=cpavey(i)+vcpy(i,k,1)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,1)=cpavey(i)
              end do

!$omp end do

            end do

          end if

        else if(mod(sbc,10).eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,1)=min(v(i,2,k)*dydt,gdydtn)
            end do

!$omp end do

          end do

        else if(sbc.eq.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,1)=gdydtn
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,1)=max(vcpy(i,k,1),-rmf8v(i,2,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,1)=max(vcpy(i,k,1),-1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(i)

          do i=1,ni-1
            vcpy(i,k,1)=vcpy(i,k,1)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the v phase speed on the north boundary.

      if((ebn.eq.1.and.jsub.eq.njsub-1).and.((nbc.ge.4.and.exbopt.eq.0) &
     &  .or.(nbc.ge.4.and.exbopt.ge.1.and.exbvar(2:2).eq.'x'))) then

        if(mod(nbc,10).eq.4.or.mod(nbc,10).eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,2)=2.e0*v(i,njm2,k)-vf(i,njm1,k)-vp(i,njm1,k)

              if(abs(vcpy(i,k,2)).lt.eps) then

                vcpy(i,k,2)=sign(eps,vcpy(i,k,2))

              end if

              vcpy(i,k,2)                                               &
     &          =max((vf(i,njm1,k)-vp(i,njm1,k))/vcpy(i,k,2),gdydt)

            end do

!$omp end do

          end do

          if(mod(nbc,10).eq.5) then

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              cpavey(i)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=cpavey(i)+vcpy(i,k,2)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,2)=cpavey(i)
              end do

!$omp end do

            end do

          end if

        else if(mod(nbc,10).eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,2)=max(v(i,njm1,k)*dydt,gdydt)
            end do

!$omp end do

          end do

        else if(nbc.eq.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,2)=gdydt
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,2)=min(vcpy(i,k,2),rmf8v(i,njm1,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,2)=min(vcpy(i,k,2),1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(i)

          do i=1,ni-1
            vcpy(i,k,2)=vcpy(i,k,2)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

!! -----

!! Calculate the w phase speed.

! Calculate the w phase speed on the west boundary.

      if((ebw.eq.1.and.isub.eq.0).and.((wbc.ge.4.and.exbopt.eq.0)       &
     &  .or.(wbc.ge.4.and.exbopt.ge.1.and.exbvar(3:3).eq.'x'))) then

        if(wbc.eq.4.or.wbc.eq.5) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,1)=wf(2,j,k)+wp(2,j,k)-2.e0*w(3,j,k)

              if(abs(wcpx(j,k,1)).lt.eps) then

                wcpx(j,k,1)=sign(eps,wcpx(j,k,1))

              end if

              wcpx(j,k,1)=min((wf(2,j,k)-wp(2,j,k))/wcpx(j,k,1),gdxdtn)

            end do

!$omp end do

          end do

          if(wbc.eq.5) then

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              cpavex(j)=0.e0
            end do

!$omp end do

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=cpavex(j)+wcpx(j,k,1)*nkm2v
              end do

!$omp end do

            end do

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,1)=cpavex(j)
              end do

!$omp end do

            end do

          end if

        else if(wbc.eq.6) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,1)=min((u(2,j,k-1)+u(2,j,k))*dxdt5,gdxdtn)
            end do

!$omp end do

          end do

        else if(wbc.ge.7) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,1)=gdxdtn
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,1)=max(wcpx(j,k,1),-rmf(2,j,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,1)=max(wcpx(j,k,1),-1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-1

!$omp do schedule(runtime) private(j)

          do j=1,nj-1
            wcpx(j,k,1)=wcpx(j,k,1)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the w phase speed on the east boundary.

      if((ebe.eq.1.and.isub.eq.nisub-1).and.((ebc.ge.4.and.exbopt.eq.0) &
     &  .or.(ebc.ge.4.and.exbopt.ge.1.and.exbvar(3:3).eq.'x'))) then

        if(ebc.eq.4.or.ebc.eq.5) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,2)=2.e0*w(nim3,j,k)-wf(nim2,j,k)-wp(nim2,j,k)

              if(abs(wcpx(j,k,2)).lt.eps) then

                wcpx(j,k,2)=sign(eps,wcpx(j,k,2))

              end if

              wcpx(j,k,2)                                               &
     &          =max((wf(nim2,j,k)-wp(nim2,j,k))/wcpx(j,k,2),gdxdt)

            end do

!$omp end do

          end do

          if(ebc.eq.5) then

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              cpavex(j)=0.e0
            end do

!$omp end do

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=cpavex(j)+wcpx(j,k,2)*nkm2v
              end do

!$omp end do

            end do

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,2)=cpavex(j)
              end do

!$omp end do

            end do

          end if

        else if(ebc.eq.6) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,2)=max((u(nim1,j,k-1)+u(nim1,j,k))*dxdt5,gdxdt)
            end do

!$omp end do

          end do

        else if(ebc.ge.7) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,2)=gdxdt
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,2)=min(wcpx(j,k,2),rmf(nim2,j,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,2)=min(wcpx(j,k,2),1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-1

!$omp do schedule(runtime) private(j)

          do j=1,nj-1
            wcpx(j,k,2)=wcpx(j,k,2)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the w phase speed on the south boundary.

      if((ebs.eq.1.and.jsub.eq.0).and.((sbc.ge.4.and.exbopt.eq.0)       &
     &  .or.(sbc.ge.4.and.exbopt.ge.1.and.exbvar(3:3).eq.'x'))) then

        if(sbc.eq.4.or.sbc.eq.5) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,1)=wf(i,2,k)+wp(i,2,k)-2.e0*w(i,3,k)

              if(abs(wcpy(i,k,1)).lt.eps) then

                wcpy(i,k,1)=sign(eps,wcpy(i,k,1))

              end if

              wcpy(i,k,1)=min((wf(i,2,k)-wp(i,2,k))/wcpy(i,k,1),gdydtn)

            end do

!$omp end do

          end do

          if(sbc.eq.5) then

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              cpavey(i)=0.e0
            end do

!$omp end do

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=cpavey(i)+wcpy(i,k,1)*nkm2v
              end do

!$omp end do

            end do

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,1)=cpavey(i)
              end do

!$omp end do

            end do

          end if

        else if(sbc.eq.6) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,1)=min((v(i,2,k-1)+v(i,2,k))*dydt5,gdydtn)
            end do

!$omp end do

          end do

        else if(sbc.ge.7) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,1)=gdydtn
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,1)=max(wcpy(i,k,1),-rmf(i,2,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,1)=max(wcpy(i,k,1),-1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-1

!$omp do schedule(runtime) private(i)

          do i=1,ni-1
            wcpy(i,k,1)=wcpy(i,k,1)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the w phase speed on the north boundary.

      if((ebn.eq.1.and.jsub.eq.njsub-1).and.((nbc.ge.4.and.exbopt.eq.0) &
     &  .or.(nbc.ge.4.and.exbopt.ge.1.and.exbvar(3:3).eq.'x'))) then

        if(nbc.eq.4.or.nbc.eq.5) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,2)=2.e0*w(i,njm3,k)-wf(i,njm2,k)-wp(i,njm2,k)

              if(abs(wcpy(i,k,2)).lt.eps) then

                wcpy(i,k,2)=sign(eps,wcpy(i,k,2))

              end if

              wcpy(i,k,2)                                               &
     &          =max((wf(i,njm2,k)-wp(i,njm2,k))/wcpy(i,k,2),gdydt)

            end do

!$omp end do

          end do

          if(nbc.eq.5) then

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              cpavey(i)=0.e0
            end do

!$omp end do

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=cpavey(i)+wcpy(i,k,2)*nkm2v
              end do

!$omp end do

            end do

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,2)=cpavey(i)
              end do

!$omp end do

            end do

          end if

        else if(nbc.eq.6) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,2)=max((v(i,njm1,k-1)+v(i,njm1,k))*dydt5,gdydt)
            end do

!$omp end do

          end do

        else if(nbc.ge.7) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,2)=gdydt
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,2)=min(wcpy(i,k,2),rmf(i,njm2,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,2)=min(wcpy(i,k,2),1.e0)
            end do

!$omp end do

          end do

        end if

        do k=2,nk-1

!$omp do schedule(runtime) private(i)

          do i=1,ni-1
            wcpy(i,k,2)=wcpy(i,k,2)*dtsdb
          end do

!$omp end do

        end do

      end if

! -----

!! -----

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_phvuvw == DUMP_TARGET_phvuvw .and. .not. dump_done_phvuvw) then
  call dump_array_3d('ucpx_ref.bin', ucpx, 1, nj, 1, nk, 1, 2)
  call dump_array_3d('ucpy_ref.bin', ucpy, 1, ni, 1, nk, 1, 2)
  call dump_array_3d('vcpx_ref.bin', vcpx, 1, nj, 1, nk, 1, 2)
  call dump_array_3d('vcpy_ref.bin', vcpy, 1, ni, 1, nk, 1, 2)
  call dump_array_3d('wcpx_ref.bin', wcpx, 1, nj, 1, nk, 1, 2)
  call dump_array_3d('wcpy_ref.bin', wcpy, 1, ni, 1, nk, 1, 2)
  call dump_array_2d('u8v_ref.bin', u8v, 0, nj+1, 1, nk)
  call dump_array_2d('v8u_ref.bin', v8u, 0, ni+1, 1, nk)
  call dump_finalize()
  dump_done_phvuvw = .true.
end if


call profile_stop(prof_id1, loop_len)

!!! -----

      end subroutine s_phvuvw

!-----7--------------------------------------------------------------7--

      end module m_phvuvw

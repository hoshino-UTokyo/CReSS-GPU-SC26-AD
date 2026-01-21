!***********************************************************************
      module m_exbcu
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/07/05, 1999/08/03, 1999/08/09, 1999/09/06,
!                   1999/09/30, 1999/11/01, 2000/01/17, 2000/02/02,
!                   2000/04/18, 2001/01/15, 2001/03/13, 2001/06/06,
!                   2001/06/29, 2001/07/13, 2001/08/07, 2001/12/11,
!                   2002/04/02, 2002/06/06, 2002/07/23, 2002/08/15,
!                   2002/10/31, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2003/11/05, 2003/11/28, 2003/12/12, 2004/05/07,
!                   2004/08/01, 2004/08/20, 2005/01/31, 2005/02/10,
!                   2006/09/21, 2006/12/04, 2007/01/05, 2007/05/07,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the lateral boundary value to the external boundary value
!     for the x components of velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comprofile
      use m_dump_kernel
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

      public :: exbcu, s_exbcu

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface exbcu

        module procedure s_exbcu

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_exbcu(fpexbvar,fpwbc,fpebc,fpexnews,fpexnorm,        &
     &                   isstp,dts,gtinc,ni,nj,nk,ucpx,ucpy,ugpv,utd,u)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpexnews
                       ! Formal parameter of unique index of exnews

      integer, intent(in) :: fpexnorm
                       ! Formal parameter of unique index of exnorm

      integer, intent(in) :: isstp
                       ! Index of small time steps integration

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(in) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

! Input and output variable

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer nim1     ! ni - 1
      integer njm1     ! nj - 1
      integer njm2     ! nj - 2

      real exnews      ! Boundary damping coefficient

      real exnorm      ! Boundary damping coefficient
                       ! for u and v in normal

      real tdmpdt      ! exnews x dts
      real ndmpdt      ! exnorm x dts

      real tpdt        ! gtinc + real(isstp - 1) x dts

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real ub1         ! Temporary variable
      real ub2         ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_exbcu = 0
      integer, parameter :: DUMP_TARGET_exbcu = 14400
      logical, save :: dump_done_exbcu = .false.


!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getrname(fpexnews,exnews)
      call getrname(fpexnorm,exnorm)

! -----

! Set the common used variables.

      nim1=ni-1
      njm1=nj-1
      njm2=nj-2

      tdmpdt=exnews*dts
      ndmpdt=exnorm*dts

      tpdt=gtinc+real(isstp-1)*dts

! -----

!! Force the lateral boundary value to the external boundary value.

!@llm start meta_info ----------------------------------------------------
! Location: exbcu.f90 :: subroutine s_exbcu
! Summary : Forces lateral boundary values of u velocity to external
!           (GPV) boundary values with radiative/relaxation approach.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module variables (ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub).
!   - No synchronization constructs.
!   - Uses intrinsic abs() - GPU compatible.
!   - Boundary-position-dependent conditionals (W, E, S, N edges).
!   - Only boundary cells are updated - sparse computation.
! Next:
!   - Consider separate kernels for each boundary region.
!   - Boundary-only work has low arithmetic intensity on GPU.
!   - Domain decomposition flags need proper handling on GPU.
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 112.2K
!   - TotalTime: 3.626s (0.12%)
!   - AvgTime: 0.252ms
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('exbcu.f90', 's_exbcu', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) * int((nj-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_exbcu = dump_call_count_exbcu + 1
if (dump_call_count_exbcu == DUMP_TARGET_exbcu .and. .not. dump_done_exbcu) then
  call dump_init('exbcu')
  call dump_scalar_c('exbvar', exbvar)
  call dump_scalar_i('wbc', wbc)
  call dump_scalar_i('ebc', ebc)
  call dump_scalar_r('exnews', exnews)
  call dump_scalar_r('exnorm', exnorm)
  call dump_scalar_i('isstp', isstp)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('dts', dts)
  call dump_scalar_r('gtinc', gtinc)
  call dump_array_3d('u_in.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ucpx.bin', ucpx, 1, nj, 1, nk, 1, 2)
  call dump_array_3d('ucpy.bin', ucpy, 1, ni, 1, nk, 1, 2)
  call dump_array_3d('ugpv.bin', ugpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('utd.bin', utd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_i('ebe', ebe)
  call dump_scalar_i('ebn', ebn)
  call dump_scalar_i('ebs', ebs)
  call dump_scalar_i('ebw', ebw)
  call dump_scalar_i('isub', isub)
  call dump_scalar_i('jsub', jsub)
  call dump_scalar_r('ndmpdt', ndmpdt)
  call dump_scalar_i('nim1', nim1)
  call dump_scalar_i('nisub', nisub)
  call dump_scalar_i('njm1', njm1)
  call dump_scalar_i('njm2', njm2)
  call dump_scalar_i('njsub', njsub)
  call dump_scalar_r('tdmpdt', tdmpdt)
  call dump_scalar_r('tpdt', tpdt)
end if

!$omp parallel default(shared)

! Force the west boundary value to the external boundary value.

      if(ebw.eq.1.and.isub.eq.0) then

        if(abs(wbc).ne.1) then

          if(exbvar(1:1).eq.'-') then

!$omp do schedule(runtime) private(j,k,ub1,ub2)

            do k=2,nk-2
            do j=1,nj-1
              ub1=ugpv(1,j,k)+utd(1,j,k)*tpdt
              ub2=ugpv(2,j,k)+utd(2,j,k)*tpdt

              u(1,j,k)=u(1,j,k)+utd(1,j,k)*dts                          &
     &          -ucpx(j,k,1)*((u(2,j,k)-u(1,j,k))-(ub2-ub1))            &
     &          -ndmpdt*(u(1,j,k)-ub1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              u(1,j,k)=u(1,j,k)+utd(1,j,k)*dts
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Force the east boundary value to the external boundary value.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(abs(ebc).ne.1) then

          if(exbvar(1:1).eq.'-') then

!$omp do schedule(runtime) private(j,k,ub1,ub2)

            do k=2,nk-2
            do j=1,nj-1
              ub1=ugpv(ni,j,k)+utd(ni,j,k)*tpdt
              ub2=ugpv(nim1,j,k)+utd(nim1,j,k)*tpdt

              u(ni,j,k)=u(ni,j,k)+utd(ni,j,k)*dts                       &
     &          +ucpx(j,k,2)*((u(nim1,j,k)-u(ni,j,k))-(ub2-ub1))        &
     &          -ndmpdt*(u(ni,j,k)-ub1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              u(ni,j,k)=u(ni,j,k)+utd(ni,j,k)*dts
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Force the south boundary value to the external boundary value.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(exbvar(1:1).eq.'-'.or.exbvar(1:1).eq.'+') then

!$omp do schedule(runtime) private(i,k,ub1,ub2)

          do k=2,nk-2
          do i=2,ni-1
            ub1=ugpv(i,1,k)+utd(i,1,k)*tpdt
            ub2=ugpv(i,2,k)+utd(i,2,k)*tpdt

            u(i,1,k)=u(i,1,k)+utd(i,1,k)*dts                            &
     &        -ucpy(i,k,1)*((u(i,2,k)-u(i,1,k))-(ub2-ub1))              &
     &        -tdmpdt*(u(i,1,k)-ub1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-2
          do i=1,ni
            u(i,1,k)=u(i,1,k)+utd(i,1,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

! Force the north boundary value to the external boundary value.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(exbvar(1:1).eq.'-'.or.exbvar(1:1).eq.'+') then

!$omp do schedule(runtime) private(i,k,ub1,ub2)

          do k=2,nk-2
          do i=2,ni-1
            ub1=ugpv(i,njm1,k)+utd(i,njm1,k)*tpdt
            ub2=ugpv(i,njm2,k)+utd(i,njm2,k)*tpdt

            u(i,njm1,k)=u(i,njm1,k)+utd(i,njm1,k)*dts                   &
     &        +ucpy(i,k,2)*((u(i,njm2,k)-u(i,njm1,k))-(ub2-ub1))        &
     &        -tdmpdt*(u(i,njm1,k)-ub1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-2
          do i=1,ni
            u(i,njm1,k)=u(i,njm1,k)+utd(i,njm1,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_exbcu == DUMP_TARGET_exbcu .and. .not. dump_done_exbcu) then
  call dump_array_3d('u_ref.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_exbcu = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_exbcu

!-----7--------------------------------------------------------------7--

      end module m_exbcu

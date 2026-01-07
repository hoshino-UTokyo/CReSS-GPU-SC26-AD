!***********************************************************************
      module m_getarea
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/10/31
!     Modification: 2003/04/30, 2003/05/19, 2003/12/12, 2004/09/10,
!                   2005/04/04, 2006/11/06, 2006/12/04, 2007/01/05,
!                   2007/01/31, 2007/05/14, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/09, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the area of each boundary plane.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comprofile
      use m_dump_kernel
      use m_getiname
      use m_getrname
      use m_reducelb
      use m_reducevb

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getarea, s_getarea

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getarea

        module procedure s_getarea

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
      subroutine s_getarea(fpwbc,fpebc,fpmpopt,fpmfcopt,fpdx,fpdy,fpdz, &
     &                     area,ni,nj,nk,jcb8u,jcb8v,rmf,rmf8u,rmf8v)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

! Output variable

      real, intent(out) :: area(0:4)
                       ! Area of each boundary plane

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction

      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction
      real dz          ! Grid distance in z direction

      real dxdy        ! dx x dy
      real dxdz        ! dx x dz
      real dydz        ! dy x dz

      real area0       ! Area of top and bottom boundary plane

      real areaw       ! Area of west boundary plane
      real areae       ! Area of east boundary plane
      real areas       ! Area of south boundary plane
      real arean       ! Area of north boundary plane

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_getarea = 0
      integer, parameter :: DUMP_TARGET_getarea = 1
      logical, save :: dump_done_getarea = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpdz,dz)

! -----

! Set the common used variables.

      dxdy=dx*dy
      dxdz=dx*dz
      dydz=dy*dz

! -----

! Get the maximum and minimum do loops index.

      if(ebw.eq.1.and.isub.eq.0) then
        istr=1
      else
        istr=2
      end if

      if(ebe.eq.1.and.isub.eq.nisub-1) then
        iend=ni-1
      else
        iend=ni-2
      end if

      if(ebs.eq.1.and.jsub.eq.0) then
        jstr=1
      else
        jstr=2
      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then
        jend=nj-1
      else
        jend=nj-2
      end if

! -----

! Initialize the sumed variables.

      area0=0.e0

      areaw=0.e0
      areae=0.e0
      areas=0.e0
      arean=0.e0

! -----

! Get the area of each boundary plane.

!@llm start meta_info ----------------------------------------------------
! Location: getarea.f90 :: s_getarea
! Summary : Compute total area of each boundary plane (top/bottom and
!           lateral) by summing grid cell areas with map scale factors.
! GPU diff: Medium
! Findings:
!   - Multiple reduction operations (+ for area0, areaw, areae, areas, arean)
!   - Complex conditional logic based on mfcopt, mpopt, boundary flags
!   - Different loops for different boundary planes
!   - Uses MPI module variables (ebw, ebe, ebs, ebn, isub, jsub, etc.)
! Next:
!   - Use GPU reduction kernels for area summation
!   - May need separate kernels for each boundary plane
!   - Consider whether GPU overhead is worthwhile for boundary-only computation
! Runtime:
!   - Calls: 1
!   - AvgLoops: 806.4K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.120ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('getarea.f90', 's_getarea', &
   & 'OMP section 1')
end if
loop_len = int((jend)-(jstr)+1,8) * int((iend)-(istr)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_getarea = dump_call_count_getarea + 1
if (dump_call_count_getarea == DUMP_TARGET_getarea .and. .not. dump_done_getarea) then
  call dump_init('getarea')
  call dump_scalar_i('fpwbc', fpwbc)
  call dump_scalar_i('fpebc', fpebc)
  call dump_scalar_i('fpmpopt', fpmpopt)
  call dump_scalar_i('fpmfcopt', fpmfcopt)
  call dump_scalar_i('fpdx', fpdx)
  call dump_scalar_i('fpdy', fpdy)
  call dump_scalar_i('fpdz', fpdz)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('jcb8u.bin', jcb8u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('jcb8v.bin', jcb8v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call dump_array_3d('rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
end if

!$omp parallel default(shared)

      if(mfcopt.eq.0) then

!$omp do schedule(runtime) private(i,j) reduction(+: area0)

        do j=jstr,jend
        do i=istr,iend
          area0=area0+dxdy
        end do
        end do

!$omp end do

      else

        if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j) reduction(+: area0)

          do j=jstr,jend
          do i=istr,iend
            area0=area0+dxdy*rmf(i,j,2)
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j) reduction(+: area0)

          do j=jstr,jend
          do i=istr,iend
            area0=area0+dxdy*rmf(i,j,3)
          end do
          end do

!$omp end do

        end if

      end if

      if(ebw.eq.1.and.isub.eq.0.and.abs(wbc).ne.1) then

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

!$omp do schedule(runtime) private(j,k) reduction(+: areaw)

          do k=2,nk-2
          do j=jstr,jend
            areaw=areaw+dydz*rmf8u(1,j,2)*jcb8u(1,j,k)
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(j,k) reduction(+: areaw)

          do k=2,nk-2
          do j=jstr,jend
            areaw=areaw+dydz*jcb8u(1,j,k)
          end do
          end do

!$omp end do

        end if

      end if

      if(ebe.eq.1.and.isub.eq.nisub-1.and.abs(ebc).ne.1) then

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

!$omp do schedule(runtime) private(j,k) reduction(+: areae)

          do k=2,nk-2
          do j=jstr,jend
            areae=areae+dydz*rmf8u(ni,j,2)*jcb8u(ni,j,k)
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(j,k) reduction(+: areae)

          do k=2,nk-2
          do j=jstr,jend
            areae=areae+dydz*jcb8u(ni,j,k)
          end do
          end do

!$omp end do

        end if

      end if

      if(ebs.eq.1.and.jsub.eq.0) then

        if(mfcopt.eq.1.and.mpopt.ne.5) then

!$omp do schedule(runtime) private(i,k) reduction(+: areas)

          do k=2,nk-2
          do i=istr,iend
            areas=areas+dxdz*rmf8v(i,1,2)*jcb8v(i,1,k)
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k) reduction(+: areas)

          do k=2,nk-2
          do i=istr,iend
            areas=areas+dxdz*jcb8v(i,1,k)
          end do
          end do

!$omp end do

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(mfcopt.eq.1.and.mpopt.ne.5) then

!$omp do schedule(runtime) private(i,k) reduction(+: arean)

          do k=2,nk-2
          do i=istr,iend
            arean=arean+dxdz*rmf8v(i,nj,2)*jcb8v(i,nj,k)
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k) reduction(+: arean)

          do k=2,nk-2
          do i=istr,iend
            arean=arean+dxdz*jcb8v(i,nj,k)
          end do
          end do

!$omp end do

        end if

      end if

!$omp end parallel

! Dump output data at target call
if (dump_call_count_getarea == DUMP_TARGET_getarea .and. .not. dump_done_getarea) then
  call dump_finalize()
  dump_done_getarea = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

! Finally get the total area of each boundary plane.

      call reducevb(area0)

      call reducelb(areaw,areae,areas,arean)

      area(0)=area0

      area(1)=areaw
      area(2)=areae
      area(3)=areas
      area(4)=arean

! -----

      end subroutine s_getarea

!-----7--------------------------------------------------------------7--

      end module m_getarea

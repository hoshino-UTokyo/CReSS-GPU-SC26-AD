!***********************************************************************
      module m_jacobian
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/01/25, 1999/04/06, 1999/07/05,
!                   1999/08/03, 1999/08/18, 1999/08/23, 1999/09/30,
!                   1999/10/12, 1999/11/19, 2000/01/17, 2000/04/18,
!                   2000/12/18, 2001/04/15, 2001/06/06, 2001/06/29,
!                   2002/04/02, 2002/06/06, 2002/10/31, 2003/04/30,
!                   2003/05/19, 2003/07/15, 2003/11/05, 2004/03/05,
!                   2004/08/20, 2005/02/10, 2005/04/04, 2006/04/03,
!                   2006/09/30, 2006/11/06, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the transformation Jacobian.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc8u
      use m_comprofile
      use m_dump_kernel
      use m_bc8v
      use m_bcyclex
      use m_bcycley
      use m_combuf
      use m_comindx
      use m_getarea
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
      use m_var8uvw

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: jacobian, s_jacobian

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface jacobian

        module procedure s_jacobian

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_jacobian(fpwbc,fpebc,fpexbopt,fpadvopt,              &
     &                    fpsmtopt,fptubopt,area,ni,nj,nk,x,y,z,zph,    &
     &                    rmf,rmf8u,rmf8v,j31,j32,jcb,jcb8u,jcb8v,jcb8w)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: x(0:ni+1)
                       ! x coordinates

      real, intent(in) :: y(0:nj+1)
                       ! y coordinates

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

! Output variables

      real, intent(out) :: area(0:4)
                       ! Area of each boundary plane

      real, intent(out) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(out) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(out) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(out) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(out) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(out) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer exbopt   ! Option for external boundary forcing

      integer advopt   ! Option for advection scheme
      integer smtopt   ! Option for numerical smoothing
      integer tubopt   ! Option for turbulent mixing

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_jacobian = 0
      integer, parameter :: DUMP_TARGET_jacobian = 1
      logical, save :: dump_done_jacobian = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpexbopt,exbopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fptubopt,tubopt)

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

! Calculate the transformation Jacobian.

!@llm start meta_info ----------------------------------------------------
! Location: jacobian.f90 :: s_jacobian
! Summary : Calculate transformation Jacobian components (j31, j32, jcb) from
!           physical coordinates (x, y, z, zph)
! GPU diff: Easy
! Findings:
!   - Three separate loop nests computing j31, j32, and jcb arrays
!   - Simple arithmetic operations (subtraction, division)
!   - Private variables: k, i, j
!   - Reads from x, y, z, zph arrays
!   - Writes to j31, j32, jcb arrays
!   - No function calls within the parallel region
!   - No sync constructs or data dependencies between grid points
! Next:
!   - Can be directly ported to GPU with OpenACC parallel loop
!   - Consider fusing loops for better GPU memory access patterns
! Runtime:
!   - Calls: 1
!   - AvgLoops: 103.6M
!   - TotalTime: 0.008s (0.00%)
!   - AvgTime: 8.288ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('jacobian.f90', 's_jacobian', &
   & 'OMP section 1')
end if
loop_len = int((nk)-(1)+1,8) * int((nj)-(0)+1,8) * int((ni)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_jacobian = dump_call_count_jacobian + 1
if (dump_call_count_jacobian == DUMP_TARGET_jacobian .and. .not. dump_done_jacobian) then
  call dump_init('jacobian')
  call dump_scalar_i('wbc', wbc)
  call dump_scalar_i('ebc', ebc)
  call dump_scalar_i('exbopt', exbopt)
  call dump_scalar_i('advopt', advopt)
  call dump_scalar_i('smtopt', smtopt)
  call dump_scalar_i('tubopt', tubopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('zph.bin', zph, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call dump_array_3d('rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call dump_array_3d('rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  ! FIXME: x is an array, not scalar
  ! ! FIXME: x is array - call dump_scalar_r('x', x)
  ! FIXME: y is an array, not scalar
  ! ! FIXME: y is array - call dump_scalar_r('y', y)
  ! FIXME: z is an array, not scalar
  ! ! FIXME: z is array - call dump_scalar_r('z', z)
end if

!$omp parallel default(shared) private(k)

      do k=1,nk

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=1,ni
          j31(i,j,k)=2.e0*(zph(i-1,j,k)-zph(i,j,k))/(x(i+1)-x(i-1))
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=0,ni
          j32(i,j,k)=2.e0*(zph(i,j-1,k)-zph(i,j,k))/(y(j+1)-y(j-1))
        end do
        end do

!$omp end do

      end do

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          jcb(i,j,k)=(zph(i,j,k+1)-zph(i,j,k))/(z(k+1)-z(k))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_jacobian == DUMP_TARGET_jacobian .and. .not. dump_done_jacobian) then
  call dump_array_3d('j31_ref.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('j32_ref.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('jcb_ref.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('jcb8u_ref.bin', jcb8u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('jcb8v_ref.bin', jcb8v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('jcb8w_ref.bin', jcb8w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_jacobian = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

!! Set the lateral boundary conditions for j31 and j32.

      if(exbopt.eq.0) then

! Set the west and east boundary conditions.

        call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,j31,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,j31,1,1,    &
     &                  rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,j31,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,j31,1,1,    &
     &                  rbuf)

        call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,j31)

        call bc8u(idwbc,idebc,ni,nj,nk,j31)

! -----

! Set the south and north boundary conditions.

        call s_putbufsy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,j32,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,j32,1,1,    &
     &                  rbuf)

        call s_putbufgy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,j32,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,j32,1,1,    &
     &                  rbuf)

        call bcycley(idsbc,idnbc,3,1,nj-2,nj_sub,ni,nj,nk,j32)

        call bc8v(idsbc,idnbc,ni,nj,nk,j32)

! -----

      end if

!! -----

! Set the periodic boundary conditions for j31.

      if(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1) then

        call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,j31,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,j31,1,1,    &
     &                  rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,j31,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,j31,1,1,    &
     &                  rbuf)

        call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,j31)

      end if

! -----

! The Jacobian is averaged to the u, v and w points.

      call var8uvw(idwbc,idebc,idexbopt,ni,nj,nk,jcb,jcb8u,jcb8v,jcb8w)

! -----

! Get the area of each boundary plane.

      if(exbopt.ge.11) then

        call getarea(idwbc,idebc,idmpopt,idmfcopt,iddx,iddy,iddz,area,  &
     &               ni,nj,nk,jcb8u,jcb8v,rmf,rmf8u,rmf8v)

      end if

! -----

!! Exchange the value in the case the 4th order calculation is
!! performed.

      if(advopt.ge.2.or.tubopt.ge.1                                     &
     &  .or.mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

! Exchange the value in x direction.

        call s_putbufsx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,jcb8u,1,1,    &
     &                  sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,jcb8u,1,1,    &
     &                  rbuf)

        call s_putbufgx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,jcb8u,1,1,    &
     &                  sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,jcb8u,1,1,    &
     &                  rbuf)

        call bcyclex(idwbc,idebc,4,0,ni-3,ni+1,ni,nj,nk,jcb8u)

! -----

! Exchange the value in y direction.

        call s_putbufsy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,jcb8v,1,1,    &
     &                  sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,jcb8v,1,1,    &
     &                  rbuf)

        call s_putbufgy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,jcb8v,1,1,    &
     &                  sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,jcb8v,1,1,    &
     &                  rbuf)

        call bcycley(idsbc,idnbc,4,0,nj-3,nj+1,ni,nj,nk,jcb8v)

! -----

      end if

!! -----

      end subroutine s_jacobian

!-----7--------------------------------------------------------------7--

      end module m_jacobian

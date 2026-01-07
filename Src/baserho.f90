!***********************************************************************
      module m_baserho
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/01/25, 1999/04/06, 1999/07/05, 1999/08/23,
!                   1999/09/30, 1999/10/12, 1999/11/01, 2000/01/17,
!                   2000/12/18, 2002/04/02, 2002/06/06, 2003/04/30,
!                   2003/05/19, 2003/07/15, 2003/11/05, 2004/03/05,
!                   2004/08/20, 2005/02/10, 2006/04/03, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     the base state density is multiplyed by the Jacobian.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcyclex
      use m_comprofile
      use m_dump_kernel
      use m_bcycley
      use m_combuf
      use m_comindx
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

      public :: baserho, s_baserho

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface baserho

        module procedure s_baserho

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
      subroutine s_baserho(fpadvopt,fpsmtopt,ni,nj,nk,jcb,rbr,          &
     &                     rst,rst8u,rst8v,rst8w)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

! Output variables

      real, intent(out) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(out) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(out) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(out) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

! Internal shared variables

      integer advopt   ! Option for advection scheme
      integer smtopt   ! Option for numerical smoothing

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_baserho = 0
      integer, parameter :: DUMP_TARGET_baserho = 1
      logical, save :: dump_done_baserho = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpadvopt,advopt)
      call getiname(fpsmtopt,smtopt)

! -----

! The base state density is multiplyed by the Jacobian.

!@llm start meta_info ----------------------------------------------------
! Location: baserho.f90 :: s_baserho
! Summary : Multiply base state density by Jacobian to compute rst array
!           for use in atmospheric dynamics calculations
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to intent(out) array rst
!   - Simple element-wise computation: rst = abs(jcb) * rbr
!   - Outer k loop with inner parallel i,j loops
! Next:
!   - Straightforward GPU port with OpenACC parallel loops
!   - Collapse all three loops (k,j,i) for maximum parallelism
!   - Intrinsic abs function is GPU-compatible
! Runtime:
!   - Calls: 1
!   - AvgLoops: 102.9M
!   - TotalTime: 0.004s (0.00%)
!   - AvgTime: 3.750ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('baserho.f90', 's_baserho', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) * int((nj)-(0)+1,8) * int((ni)-(0)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_baserho = dump_call_count_baserho + 1
if (dump_call_count_baserho == DUMP_TARGET_baserho .and. .not. dump_done_baserho) then
  call dump_init('baserho')
  call dump_scalar_i('fpadvopt', fpadvopt)
  call dump_scalar_i('fpsmtopt', fpsmtopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          rst(i,j,k)=abs(jcb(i,j,k))*rbr(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_baserho == DUMP_TARGET_baserho .and. .not. dump_done_baserho) then
  call dump_array_3d('rst_ref.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rst8u_ref.bin', rst8u, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rst8v_ref.bin', rst8v, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rst8w_ref.bin', rst8w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_baserho = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

! The rst is averaged to the u, v and w points.

      call var8uvw(idwbc,idebc,idexbopt,ni,nj,nk,rst,rst8u,rst8v,rst8w)

! -----

!! Exchange the value in the case the 4th order calculation is
!! performed.

      if(advopt.ge.2.or.mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

! Exchange the value in x direction.

        call s_putbufsx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,rst8u,1,1,    &
     &                  sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,rst8u,1,1,    &
     &                  rbuf)

        call s_putbufgx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,rst8u,1,1,    &
     &                  sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,rst8u,1,1,    &
     &                  rbuf)

        call bcyclex(idwbc,idebc,4,0,ni-3,ni+1,ni,nj,nk,rst8u)

! -----

! Exchange the value in y direction.

        call s_putbufsy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,rst8v,1,1,    &
     &                  sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,rst8v,1,1,    &
     &                  rbuf)

        call s_putbufgy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,rst8v,1,1,    &
     &                  sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,rst8v,1,1,    &
     &                  rbuf)

        call bcycley(idsbc,idnbc,4,0,nj-3,nj+1,ni,nj,nk,rst8v)

! -----

      end if

!! -----

      end subroutine s_baserho

!-----7--------------------------------------------------------------7--

      end module m_baserho

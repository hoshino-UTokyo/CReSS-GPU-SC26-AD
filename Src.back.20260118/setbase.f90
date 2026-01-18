!***********************************************************************
      module m_setbase
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/05/20, 1999/06/28, 1999/07/05, 1999/08/03,
!                   1999/08/18, 1999/08/23, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2000/03/08, 2000/07/05,
!                   2001/01/15, 2001/04/15, 2001/05/29, 2001/06/29,
!                   2001/10/17, 2002/04/02, 2002/06/18, 2002/08/15,
!                   2002/09/09, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/01/09, 2004/04/15, 2004/09/10, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the base state variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcbase
      use m_comprofile
      use m_dump_kernel
      use m_comindx
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setbase, s_setbase

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setbase

        module procedure s_setbase

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_setbase(ni,nj,nk,zph,ubr,vbr,pbr,ptbr,qvbr,rbr,      &
     &                     zph8s,pibr,ptvbr)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Input and output variables

      real, intent(inout) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(inout) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(inout) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(inout) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(inout) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

! Output variable

      real, intent(out) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

! Internal shared variables

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

      real, intent(inout) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(inout) :: pibr(0:ni+1,0:nj+1,1:nk)
                       ! Base state Exner function

      real, intent(inout) :: ptvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state virtual potential temperature

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_setbase = 0
      integer, parameter :: DUMP_TARGET_setbase = 1
      logical, save :: dump_done_setbase = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rddvcp=rd/cp
      p0iv=1.e0/p0

! -----

!! Set the base state variables.

!@llm start meta_info ----------------------------------------------------
! Location: setbase.f90 :: s_setbase
! Summary : Set base state variables including z coordinates at scalar points,
!           Exner function, virtual potential temperature, and density
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No external function calls (only intrinsic exp, log)
!   - Two separate loop nests with different k ranges
!   - Writes to zph8s, ptvbr, pibr, rbr arrays
!   - Uses module constants from m_comphy (rd, cp, p0, epsav)
!   - No synchronization constructs
! Next:
!   - Convert to OpenACC with parallel loop collapse(3) for each loop nest
!   - Ensure module constants accessible on device
!   - Note: bcbase call after parallel region needs separate handling
! Runtime:
!   - Calls: 1
!   - AvgLoops: 102.9M
!   - TotalTime: 0.013s (0.00%)
!   - AvgTime: 13.291ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('setbase.f90', 's_setbase', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) * int((nj)-(0)+1,8) * int((ni)-(0)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_setbase = dump_call_count_setbase + 1
if (dump_call_count_setbase == DUMP_TARGET_setbase .and. .not. dump_done_setbase) then
  call dump_init('setbase')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('rd', rd)
  call dump_scalar_r('epsav', epsav)
  call dump_array_3d('zph.bin', zph, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ubr_in.bin', ubr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vbr_in.bin', vbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pbr_in.bin', pbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptbr_in.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvbr_in.bin', qvbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('zph8s_in.bin', zph8s, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pibr_in.bin', pibr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptvbr_in.bin', ptvbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('p0iv', p0iv)
  call dump_scalar_r('rddvcp', rddvcp)
end if

!$omp parallel default(shared) private(k)

! Calculate the z physical coordinates at the scalar, u and v points.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          zph8s(i,j,k)=.5e0*(zph(i,j,k+1)+zph(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Get the base state Exner function and the density.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          ptvbr(i,j,k)=ptbr(i,j,k)                                      &
     &      *(1.e0+epsav*qvbr(i,j,k))/(1.e0+qvbr(i,j,k))

          pibr(i,j,k)=exp(rddvcp*log(p0iv*pbr(i,j,k)))

          rbr(i,j,k)=pbr(i,j,k)/(rd*ptvbr(i,j,k)*pibr(i,j,k))

        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_setbase == DUMP_TARGET_setbase .and. .not. dump_done_setbase) then
  call dump_array_3d('rbr_ref.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ubr_ref.bin', ubr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vbr_ref.bin', vbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pbr_ref.bin', pbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptbr_ref.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvbr_ref.bin', qvbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('zph8s_ref.bin', zph8s, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pibr_ref.bin', pibr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptvbr_ref.bin', ptvbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_setbase = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

! Set the boundary conditions.

      call bcbase(idsmtopt,ni,nj,nk,                                    &
     &            zph8s,ubr,vbr,pbr,ptbr,qvbr,rbr,pibr,ptvbr)

! -----

      end subroutine s_setbase

!-----7--------------------------------------------------------------7--

      end module m_setbase

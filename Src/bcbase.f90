!***********************************************************************
      module m_bcbase
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/04/06, 1999/07/05, 1999/08/03, 1999/08/18,
!                   1999/08/23, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/03/08, 2000/07/05, 2001/04/15,
!                   2001/05/29, 2001/12/11, 2002/04/02, 2002/06/18,
!                   2002/08/15, 2003/04/30, 2003/05/19, 2003/11/05,
!                   2003/12/12, 2004/01/09, 2004/03/05, 2004/04/15,
!                   2004/08/20, 2005/02/10, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the bottom and the top boundary conditions for the base state
!     variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcyclex
      use m_comprofile
      use m_dump_kernel
      use m_bcycley
      use m_combuf
      use m_comindx
      use m_comphy
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

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bcbase, s_bcbase

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcbase

        module procedure s_bcbase

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_bcbase(fpsmtopt,ni,nj,nk,zph8s,                      &
     &                    ubr,vbr,pbr,ptbr,qvbr,rbr,pibr,ptvbr)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

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

      real, intent(inout) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(inout) :: pibr(0:ni+1,0:nj+1,1:nk)
                       ! Base state Exner function

      real, intent(inout) :: ptvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state virtual potential temperature

! Internal shared variables

      integer smtopt   ! Option for numerical smoothing

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      real gdvcp2      ! 2.0 x g / cp

      real cpdvrd      ! cp / rd

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_bcbase = 0
      integer, parameter :: DUMP_TARGET_bcbase = 1
      logical, save :: dump_done_bcbase = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpsmtopt,smtopt)

! -----

! Set the common used variables

      nkm1=nk-1
      nkm2=nk-2

      gdvcp2=2.e0*g/cp

      cpdvrd=cp/rd

! -----

!! Set the bottom and the top boundary conditions for the base state
!! variables.

!@llm start meta_info ----------------------------------------------------
! Location: bcbase.f90 :: s_bcbase
! Summary : Sets bottom and top boundary conditions for base state variables
!           (ubr, vbr, ptbr, qvbr, ptvbr, pibr, pbr, rbr) using extrapolation.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Multiple 2D loops over i,j for different variables
!   - Uses physical constants from m_comphy (g, cp, rd, p0)
!   - Exner function BC requires exp/log calculations
!   - Pressure and density BCs depend on previously computed pibr and ptvbr
! Next:
!   - Convert to OpenACC with collapsed i,j loops
!   - Ensure data dependencies between loops are respected (ptvbr before pibr, pibr before pbr/rbr)
!   - Consider fusing independent loops for better kernel efficiency
! Runtime:
!   - Calls: 1
!   - AvgLoops: 809.1K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.447ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('bcbase.f90', 's_bcbase', &
   & 'OMP section 1')
end if
loop_len = int((nj)-(0)+1,8) * int((ni)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_bcbase = dump_call_count_bcbase + 1
if (dump_call_count_bcbase == DUMP_TARGET_bcbase .and. .not. dump_done_bcbase) then
  call dump_init('bcbase')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nkm1', nkm1)
  call dump_scalar_i('nkm2', nkm2)
  call dump_scalar_r('gdvcp2', gdvcp2)
  call dump_scalar_r('cpdvrd', cpdvrd)
  call dump_array_3d('zph8s.bin', zph8s, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ubr_in.bin', ubr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vbr_in.bin', vbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pbr_in.bin', pbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptbr_in.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvbr_in.bin', qvbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rbr_in.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pibr_in.bin', pibr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptvbr_in.bin', ptvbr, 0, ni+1, 0, nj+1, 1, nk)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_031)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    !$acc kernels
    !$acc loop independent collapse(2)
    do j = 0, nj
      do i = 1, ni
        ubr(i,j,1) = ubr(i,j,2)
        ubr(i,j,nkm1) = ubr(i,j,nkm2)
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(2)
    do j = 1, nj
      do i = 0, ni
        vbr(i,j,1) = vbr(i,j,2)
        vbr(i,j,nkm1) = vbr(i,j,nkm2)
      end do
    end do
    !$acc end kernels

    ! Set the bottom and the top boundary conditions for the base state
    ! potential temperature and water vapor mixing ratio.
    !$acc kernels
    !$acc loop independent collapse(2)
    do j = 0, nj
      do i = 0, ni
        ptbr(i,j,1) = ptbr(i,j,2)
        ptbr(i,j,nkm1) = ptbr(i,j,nkm2)
        qvbr(i,j,1) = qvbr(i,j,2)
        qvbr(i,j,nkm1) = qvbr(i,j,nkm2)
      end do
    end do
    !$acc end kernels

    ! Set the bottom and the top boundary conditions for the base state
    ! virtual potential temperature.
    !$acc kernels
    !$acc loop independent collapse(2)
    do j = 0, nj
      do i = 0, ni
        ptvbr(i,j,1) = ptvbr(i,j,2)
        ptvbr(i,j,nkm1) = ptvbr(i,j,nkm2)
      end do
    end do
    !$acc end kernels

    ! Set the bottom and the top boundary conditions for the base state
    ! Exner function.
    !$acc kernels
    !$acc loop independent collapse(2)
    do j = 0, nj
      do i = 0, ni
        pibr(i,j,1) = pibr(i,j,2) &
          + gdvcp2 * (zph8s(i,j,2) - zph8s(i,j,1)) / (ptvbr(i,j,1) + ptvbr(i,j,2))
        pibr(i,j,nkm1) = pibr(i,j,nkm2) &
          - gdvcp2 * (zph8s(i,j,nkm1) - zph8s(i,j,nkm2)) / (ptvbr(i,j,nkm2) + ptvbr(i,j,nkm1))
      end do
    end do
    !$acc end kernels

    ! Set the bottom and the top boundary conditions for the base state
    ! pressure and the base state density.
    !$acc kernels
    !$acc loop independent collapse(2)
    do j = 0, nj
      do i = 0, ni
        pbr(i,j,1) = p0 * exp(cpdvrd * log(pibr(i,j,1)))
        pbr(i,j,nkm1) = p0 * exp(cpdvrd * log(pibr(i,j,nkm1)))
        rbr(i,j,1) = pbr(i,j,1) / (rd * ptvbr(i,j,1) * pibr(i,j,1))
        rbr(i,j,nkm1) = pbr(i,j,nkm1) / (rd * ptvbr(i,j,nkm1) * pibr(i,j,nkm1))
      end do
    end do
    !$acc end kernels

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared)

! Set the bottom and the top boundary conditions for the base state
! velocity.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=1,ni
        ubr(i,j,1)=ubr(i,j,2)
        ubr(i,j,nkm1)=ubr(i,j,nkm2)
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj
      do i=0,ni
        vbr(i,j,1)=vbr(i,j,2)
        vbr(i,j,nkm1)=vbr(i,j,nkm2)
      end do
      end do

!$omp end do

! -----

! Set the bottom and the top boundary conditions for the base state
! potential temperature and water vapor mixing ratio.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        ptbr(i,j,1)=ptbr(i,j,2)
        ptbr(i,j,nkm1)=ptbr(i,j,nkm2)

        qvbr(i,j,1)=qvbr(i,j,2)
        qvbr(i,j,nkm1)=qvbr(i,j,nkm2)

      end do
      end do

!$omp end do

! -----

! Set the bottom and the top boundary conditions for the base state
! virtual potential temperature.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        ptvbr(i,j,1)=ptvbr(i,j,2)
        ptvbr(i,j,nkm1)=ptvbr(i,j,nkm2)
      end do
      end do

!$omp end do

! -----

! Set the bottom and the top boundary conditions for the base state
! exnar function.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        pibr(i,j,1)=pibr(i,j,2)                                         &
     &   +gdvcp2*(zph8s(i,j,2)-zph8s(i,j,1))/(ptvbr(i,j,1)+ptvbr(i,j,2))

        pibr(i,j,nkm1)=pibr(i,j,nkm2)-gdvcp2*(zph8s(i,j,nkm1)           &
     &   -zph8s(i,j,nkm2))/(ptvbr(i,j,nkm2)+ptvbr(i,j,nkm1))

      end do
      end do

!$omp end do

! -----

! Set the bottom and the top boundary conditions for the base state
! pressure and the base state density.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        pbr(i,j,1)=p0*exp(cpdvrd*log(pibr(i,j,1)))
        pbr(i,j,nkm1)=p0*exp(cpdvrd*log(pibr(i,j,nkm1)))

        rbr(i,j,1)=pbr(i,j,1)/(rd*ptvbr(i,j,1)*pibr(i,j,1))
        rbr(i,j,nkm1)=pbr(i,j,nkm1)/(rd*ptvbr(i,j,nkm1)*pibr(i,j,nkm1))

      end do
      end do

!$omp end do

! -----

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_bcbase == DUMP_TARGET_bcbase .and. .not. dump_done_bcbase) then
  call dump_array_3d('ubr_ref.bin', ubr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vbr_ref.bin', vbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pbr_ref.bin', pbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptbr_ref.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvbr_ref.bin', qvbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rbr_ref.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pibr_ref.bin', pibr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptvbr_ref.bin', ptvbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_bcbase = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

!! Exchange the value in the case the 4th order calculation is
!! performed.

      if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

! Exchange the value in x direction.

        call s_putbufsx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,ubr,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,ubr,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,ubr,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,ubr,1,1,rbuf)

        call bcyclex(idwbc,idebc,4,0,ni-3,ni+1,ni,nj,nk,ubr)

! -----

! Exchange the value in y direction.

        call s_putbufsy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,vbr,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,vbr,1,1,rbuf)

        call s_putbufgy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,vbr,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,vbr,1,1,rbuf)

        call bcycley(idsbc,idnbc,4,0,nj-3,nj+1,ni,nj,nk,vbr)

! -----

      end if

!! -----

      end subroutine s_bcbase

!-----7--------------------------------------------------------------7--

      end module m_bcbase

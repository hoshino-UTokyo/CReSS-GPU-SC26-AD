!***********************************************************************
      module m_sfcflx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/11/14, 2001/12/03, 2001/12/10, 2002/04/02,
!                   2002/12/02, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/10/31, 2003/12/12, 2004/02/01, 2004/03/05,
!                   2004/04/01, 2004/05/07, 2004/07/01, 2004/08/20,
!                   2004/09/01, 2004/09/10, 2006/05/12, 2007/09/04,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/06/16, 2011/06/01, 2011/09/22, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the bulk and the exchange coefficients of surface flux.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bulksfc
      use m_comprofile
      use m_dump_kernel
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: sfcflx, s_sfcflx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sfcflx

        module procedure s_sfcflx

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
      subroutine s_sfcflx(ni,nj,nk,za,rbr,land,kai,z0m,z0h,va,rch,      &
     &                    cm,ch,ce,ct,cq)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: za(0:ni+1,0:nj+1)
                       ! z physical coordinates at lowest plane

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(in) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(in) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(in) :: va(0:ni+1,0:nj+1)
                       ! Magnitude of velocity at lowest plane

      real, intent(in) :: rch(0:ni+1,0:nj+1)
                       ! Bulk Richardson number

! Output variables

      real, intent(out) :: cm(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity

      real, intent(out) :: ch(0:ni+1,0:nj+1)
                       ! Bulk coefficient for scalar

      real, intent(out) :: ce(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface momentum flux

      real, intent(out) :: ct(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface heat flux

      real, intent(out) :: cq(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface moisture flux

! Internal shared variable

      real rddwkp      ! ln(da0 / dv0) / wkappa

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real a           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_sfcflx = 0
      integer, parameter :: DUMP_TARGET_sfcflx = 361
      logical, save :: dump_done_sfcflx = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variable.

      rddwkp=log(da0/dv0)/wkappa

! -----

! Calculate the bulk coefficients of surface flux.

      call bulksfc(ni,nj,za,land,kai,z0m,z0h,rch,cm,ch)

! -----

! Finally get the exchange coefficients of surface flux.

!@llm start meta_info ----------------------------------------------------
! Location: sfcflx.f90 :: subroutine s_sfcflx
! Summary : Calculates exchange coefficients for surface momentum, heat,
!           and moisture fluxes from bulk coefficients.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Reads module constants (da0, dv0, wkappa) from comphy.
!   - No synchronization constructs.
!   - Simple 2D loop over (i,j) - surface arrays only.
!   - Land/sea conditional for moisture exchange coefficient.
!   - All grid points are independent (embarrassingly parallel).
!   - Uses intrinsic log() outside parallel region.
! Next:
!   - Direct OpenACC kernels should work well.
!   - 2D arrays only - good memory access pattern.
!   - Consider fusing with bulksfc call if beneficial.
! Runtime:
!   - Calls: 361
!   - AvgLoops: 806.4K
!   - TotalTime: 0.015s (0.00%)
!   - AvgTime: 0.043ms
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('sfcflx.f90', 's_sfcflx', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_sfcflx = dump_call_count_sfcflx + 1
if (dump_call_count_sfcflx == DUMP_TARGET_sfcflx .and. .not. dump_done_sfcflx) then
  call dump_init('sfcflx')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_2d('za.bin', za, 0, ni+1, 0, nj+1)
  call dump_array_3d('rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d_int('land.bin', land, 0, ni+1, 0, nj+1)
  call dump_array_2d('kai.bin', kai, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0m.bin', z0m, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0h.bin', z0h, 0, ni+1, 0, nj+1)
  call dump_array_2d('va.bin', va, 0, ni+1, 0, nj+1)
  call dump_array_2d('rch.bin', rch, 0, ni+1, 0, nj+1)
end if

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,a)

      do j=1,nj-1
      do i=1,ni-1

        a=rbr(i,j,2)*cm(i,j)*va(i,j)

        ce(i,j)=a*cm(i,j)
        ct(i,j)=a*ch(i,j)

        if(land(i,j).lt.0) then

          cq(i,j)=ct(i,j)/(1.e0+rddwkp*ch(i,j))

        else

          cq(i,j)=ct(i,j)

        end if

      end do
      end do

!$omp end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_sfcflx == DUMP_TARGET_sfcflx .and. .not. dump_done_sfcflx) then
  call dump_array_2d('cm_ref.bin', cm, 0, ni+1, 0, nj+1)
  call dump_array_2d('ch_ref.bin', ch, 0, ni+1, 0, nj+1)
  call dump_array_2d('ce_ref.bin', ce, 0, ni+1, 0, nj+1)
  call dump_array_2d('ct_ref.bin', ct, 0, ni+1, 0, nj+1)
  call dump_array_2d('cq_ref.bin', cq, 0, ni+1, 0, nj+1)
  call dump_finalize()
  dump_done_sfcflx = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_sfcflx

!-----7--------------------------------------------------------------7--

      end module m_sfcflx

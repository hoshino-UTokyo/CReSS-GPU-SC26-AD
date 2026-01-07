!***********************************************************************
      module m_roughitr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/15
!     Modification: 2001/10/18, 2001/11/14, 2001/12/20, 2002/04/02,
!                   2002/12/02, 2003/02/05, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/10/31, 2003/12/12, 2004/02/01,
!                   2003/03/05, 2004/04/01, 2004/04/15, 2004/05/07,
!                   2004/08/20, 2004/09/01, 2004/11/10, 2006/01/10,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/06/16, 2009/11/13, 2011/06/01,
!                   2011/11/10, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the roughness length on the sea surface by iteration at
!     forecast start time.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bulksfc
      use m_comprofile
      use m_dump_kernel
      use m_chkitr
      use m_comphy
      use m_getrich
      use m_setcst2d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: roughitr, s_roughitr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface roughitr

        module procedure s_roughitr

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_roughitr(ni,nj,nk,za,land,kai,ptv,va,z0m,z0h,        &
     &                      rch,cm,ch,dz0m)
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

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(in) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(in) :: va(0:ni+1,0:nj+1)
                       ! Magnitude of velocity at lowest plane

! Input and output variables

      real, intent(inout) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(inout) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

! Internal shared variables

      real itcon       ! Control flag of continuation of iteration

      real z0meps      ! Minimum value of convergence of iteration

      real, intent(inout) :: rch(0:ni+1,0:nj+1)
                       ! Bulk Richardson number

      real, intent(inout) :: cm(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity

      real, intent(inout) :: ch(0:ni+1,0:nj+1)
                       ! Bulk coefficient for scalar

      real, intent(inout) :: dz0m(0:ni+1,0:nj+1)
                       ! Variations of z0m

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real ust         ! Friction velocity

      real z0itr       ! Roughness length at current loop


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_roughitr = 0
      integer, parameter :: DUMP_TARGET_roughitr = 16
      logical, save :: dump_done_roughitr = .false.


!-----7--------------------------------------------------------------7--

!!! Calculate the roughness length on the sea surface by iteration at
!!! forecast start time.

! Set the common used variable.

      z0meps=1.e-5

! -----

! Fill in the array with constants.

      call setcst2d(0,ni+1,0,nj+1,1.e0,dz0m)

! -----

!! Iterate the roughness parameter on the sea surface.

      iterate: do

! Calculate the bulk Richardson number on the surface.

        call getrich(ni,nj,nk,za,land,kai,z0m,z0h,ptv,va,rch)

! -----

! Calculate the bulk coefficients of surface flux.

        call bulksfc(ni,nj,za,land,kai,z0m,z0h,rch,cm,ch)

! -----

! If the variations of the roughness parameter is greater than z0meps,
! perform the iteration.

!@llm start meta_info ----------------------------------------------------
! Location: roughitr.f90 :: s_roughitr
! Summary : Calculate sea surface roughness length by iteration using bulk coefficient
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - Located inside iterate: do loop (outer iteration loop is serial)
!   - Uses intrinsic functions: abs, max
!   - Conditional logic based on land use (land < 3) and ust threshold
!   - Updates z0m, z0h, dz0m arrays
!   - External subroutine calls getrich, bulksfc, chkitr before/after parallel region
!   - No synchronization constructs inside parallel region
! Next:
!   - Convert to OpenACC with teams distribute parallel for
!   - Keep iteration control on host, only offload inner 2D loop
!   - May need to manage z0m, z0h data between iterations on GPU
! Runtime:
!   - Calls: 16
!   - AvgLoops: 806.4K
!   - TotalTime: 0.001s (0.00%)
!   - AvgTime: 0.033ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('roughitr.f90', 's_roughitr', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_roughitr = dump_call_count_roughitr + 1
if (dump_call_count_roughitr == DUMP_TARGET_roughitr .and. .not. dump_done_roughitr) then
  call dump_init('roughitr')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_2d('za.bin', za, 0, ni+1, 0, nj+1)
  call dump_array_2d_int('land.bin', land, 0, ni+1, 0, nj+1)
  call dump_array_2d('kai.bin', kai, 0, ni+1, 0, nj+1)
  call dump_array_3d('ptv.bin', ptv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('va.bin', va, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0m_in.bin', z0m, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0h_in.bin', z0h, 0, ni+1, 0, nj+1)
  call dump_array_2d('rch_in.bin', rch, 0, ni+1, 0, nj+1)
  call dump_array_2d('cm_in.bin', cm, 0, ni+1, 0, nj+1)
  call dump_array_2d('ch_in.bin', ch, 0, ni+1, 0, nj+1)
  call dump_array_2d('dz0m_in.bin', dz0m, 0, ni+1, 0, nj+1)
end if

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,ust,z0itr)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).lt.3) then

            ust=cm(i,j)*va(i,j)

            if(ust.lt.1.08e0) then

              z0itr=max(-34.7e-6+8.28e-4*ust,z0min)

            else

              z0itr=max(-.277e-2+3.39e-3*ust,z0min)

            end if

            dz0m(i,j)=abs(z0m(i,j)/z0itr-1.e0)

            z0m(i,j)=z0itr
            z0h(i,j)=z0m(i,j)

          else

            dz0m(i,j)=0.e0

          end if

        end do
        end do

!$omp end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_roughitr == DUMP_TARGET_roughitr .and. .not. dump_done_roughitr) then
  call dump_array_2d('z0m_ref.bin', z0m, 0, ni+1, 0, nj+1)
  call dump_array_2d('z0h_ref.bin', z0h, 0, ni+1, 0, nj+1)
  call dump_array_2d('rch_ref.bin', rch, 0, ni+1, 0, nj+1)
  call dump_array_2d('cm_ref.bin', cm, 0, ni+1, 0, nj+1)
  call dump_array_2d('ch_ref.bin', ch, 0, ni+1, 0, nj+1)
  call dump_array_2d('dz0m_ref.bin', dz0m, 0, ni+1, 0, nj+1)
  call dump_finalize()
  dump_done_roughitr = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

! If the variations of the roughness parameter is greater than z0meps,
! perform the iteration again.

        call s_chkitr('common',1,ni-1,1,nj-1,1,1,itcon,ni,nj,1,dz0m)

        if(itcon.lt.z0meps) then

          exit iterate

        end if

! -----

      end do iterate

!! -----

!!! -----

      end subroutine s_roughitr

!-----7--------------------------------------------------------------7--

      end module m_roughitr

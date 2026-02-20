!***********************************************************************
      module m_getvdens
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/04/01
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the inverse of base state density.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing
      use m_comprofile
      use m_dump_kernel

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getvdens, s_getvdens

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getvdens

        module procedure s_getvdens

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getvdens(ni,nj,nk,rbr,rbv)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

! Output variable

      real, intent(out) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_getvdens = 0
      integer, parameter :: DUMP_TARGET_getvdens = 360
      logical, save :: dump_done_getvdens = .false.


!-----7--------------------------------------------------------------7--

! Calculate the inverse of base state density.

!@llm start meta_info ----------------------------------------------------
! Location: getvdens.f90 :: s_getvdens
! Summary : Compute inverse of base state density (1/rbr) for all 3D
!           grid points for use in momentum equations.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Simple element-wise division: rbv = 1/rbr
!   - No global writes, only output array rbv is modified
!   - No synchronization constructs other than implicit barrier
! Next:
!   - Direct translation to OpenACC with collapsed loops
!   - Consider loop collapse for k,j,i dimensions
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 1.102s (0.04%)
!   - AvgTime: 3.061ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('getvdens.f90', 's_getvdens', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_getvdens = dump_call_count_getvdens + 1
if (dump_call_count_getvdens == DUMP_TARGET_getvdens .and. .not. dump_done_getvdens) then
  call dump_init('getvdens')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_142)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    !$acc kernels
    do k = 1, nk-1
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          rbv(i,j,k) = 1.e0 / rbr(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          rbv(i,j,k)=1.e0/rbr(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

#endif

! Dump output data at target call
if (dump_call_count_getvdens == DUMP_TARGET_getvdens .and. .not. dump_done_getvdens) then
  call dump_array_3d('rbv_ref.bin', rbv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_getvdens = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_getvdens

!-----7--------------------------------------------------------------7--

      end module m_getvdens

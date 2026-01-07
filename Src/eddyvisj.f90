!***********************************************************************
      module m_eddyvisj
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/02/13
!     Modification: 2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     the eddy viscosity is devided by Jacobian.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: eddyvisj, s_eddyvisj

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface eddyvisj

        module procedure s_eddyvisj

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
      subroutine s_eddyvisj(fpmfcopt,ni,nj,nk,jcb,mf,rkh,rkv)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

! Input and output variables

      real, intent(inout) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy diffusivity / Jacobian

      real, intent(inout) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy diffusivity / Jacobian

! Internal shared variable

      integer mfcopt   ! Option for map scale factor

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real jcbiv       ! Inverse of Jacobian


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_eddyvisj = 0
      integer, parameter :: DUMP_TARGET_eddyvisj = 360
      logical, save :: dump_done_eddyvisj = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpmfcopt,mfcopt)

! -----

! The eddy viscosity is devided by Jacobian.

!@llm start meta_info ----------------------------------------------------
! Location: eddyvisj.f90 :: s_eddyvisj
! Summary : Divide eddy viscosity by Jacobian, with optional map scale
!           factor multiplication for horizontal component
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - Simple element-wise division operations
!   - Conditional branch based on mfcopt option (map scale factor)
!   - Writes to rkh, rkv arrays (no race conditions)
!   - No synchronization constructs besides implicit barrier
! Next:
!   - Straightforward conversion to OpenACC or OpenACC kernels
!   - Collapse k,j,i loops for maximum parallelism
!   - Consider using a single kernel with conditional inside for both paths
! Runtime:
!   - Calls: 360
!   - AvgLoops: 102.4M
!   - TotalTime: 2.137s (0.07%)
!   - AvgTime: 5.937ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('eddyvisj.f90', 's_eddyvisj', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_eddyvisj = dump_call_count_eddyvisj + 1
if (dump_call_count_eddyvisj == DUMP_TARGET_eddyvisj .and. .not. dump_done_eddyvisj) then
  call dump_init('eddyvisj')
  call dump_scalar_i('fpmfcopt', fpmfcopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_2d('mf.bin', mf, 0, ni+1, 0, nj+1)
  call dump_array_3d('rkh_in.bin', rkh, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rkv_in.bin', rkv, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

      if(mfcopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv)

          do j=1,nj-1
          do i=1,ni-1
            jcbiv=1.e0/jcb(i,j,k)

            rkh(i,j,k)=jcbiv*rkh(i,j,k)
            rkv(i,j,k)=jcbiv*rkv(i,j,k)

          end do
          end do

!$omp end do

        end do

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv)

          do j=1,nj-1
          do i=1,ni-1
            jcbiv=1.e0/jcb(i,j,k)

            rkh(i,j,k)=jcbiv*mf(i,j)*rkh(i,j,k)
            rkv(i,j,k)=jcbiv*rkv(i,j,k)

          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! Dump output data at target call
if (dump_call_count_eddyvisj == DUMP_TARGET_eddyvisj .and. .not. dump_done_eddyvisj) then
  call dump_array_3d('rkh_ref.bin', rkh, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rkv_ref.bin', rkv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_eddyvisj = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_eddyvisj

!-----7--------------------------------------------------------------7--

      end module m_eddyvisj

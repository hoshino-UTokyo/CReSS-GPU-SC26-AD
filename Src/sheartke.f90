!***********************************************************************
      module m_sheartke
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/10/12
!     Modification: 1999/11/01, 2000/01/17, 2000/07/05, 2000/12/19,
!                   2001/05/29, 2001/11/20, 2002/04/02, 2003/01/04,
!                   2003/04/30, 2003/05/19, 2003/11/28, 2003/12/12,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/11/10, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the shear production in the turbulent kinetic energy
!     equation.

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

      public :: sheartke, s_sheartke

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sheartke

        module procedure s_sheartke

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
      subroutine s_sheartke(ni,nj,nk,jcb,ssq,rkv,tkefrc)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: ssq(0:ni+1,0:nj+1,1:nk)
                       ! Magnitude of deformation squared

      real, intent(in) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x virtical eddy viscosity

! Input and output variable

      real, intent(inout) :: tkefrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in turbulent kinetic energy equation

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_sheartke = 0
      integer, parameter :: DUMP_TARGET_sheartke = 360
      logical, save :: dump_done_sheartke = .false.


!-----7--------------------------------------------------------------7--

! Calculate the shear production in the turbulent kinetic energy
! equation.

!@llm start meta_info ----------------------------------------------------
! Location: sheartke.f90 :: s_sheartke
! Summary : Calculates shear production term in TKE equation by adding
!           Jacobian times eddy viscosity times deformation squared to forcing.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Simple 3D loop with single multiply-add operation per point
!   - All loops independent with private i,j,k indices
!   - No synchronization constructs
!   - Minimal computation per grid point
! Next:
!   - Straightforward GPU port with 3D kernel
!   - Use OpenACC/OpenACC with collapse(3)
!   - Good candidate for kernel fusion with other TKE terms
! Runtime:
!   - Calls: 360
!   - AvgLoops: 100.4M
!   - TotalTime: 1.960s (0.07%)
!   - AvgTime: 5.446ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('sheartke.f90', 's_sheartke', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) &
     & * int((nj-2)-(2)+1,8) &
     & * int((ni-2)-(2)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_sheartke = dump_call_count_sheartke + 1
if (dump_call_count_sheartke == DUMP_TARGET_sheartke .and. .not. dump_done_sheartke) then
  call dump_init('sheartke')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ssq.bin', ssq, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('rkv.bin', rkv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tkefrc_in.bin', tkefrc, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tkefrc(i,j,k)=tkefrc(i,j,k)+jcb(i,j,k)*rkv(i,j,k)*ssq(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_sheartke == DUMP_TARGET_sheartke .and. .not. dump_done_sheartke) then
  call dump_array_3d('tkefrc_ref.bin', tkefrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_sheartke = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_sheartke

!-----7--------------------------------------------------------------7--

      end module m_sheartke

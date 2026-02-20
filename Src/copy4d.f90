!***********************************************************************
      module m_copy4d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/08/23, 1999/10/12,
!                   2000/01/17, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2006/01/10, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     copy 4 dimensional invar to outvar.

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

      public :: copy4d, s_copy4d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface copy4d

        module procedure s_copy4d

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
      subroutine s_copy4d(imin,imax,jmin,jmax,kmin,kmax,nmin,nmax,      &
     &                    invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      integer, intent(in) :: nmin
                       ! Minimum array index in 4th direction

      integer, intent(in) :: nmax
                       ! Maximum array index in 4th direction

      real, intent(in) :: invar(imin:imax,jmin:jmax,kmin:kmax,nmin:nmax)
                       ! Coping variable

! Output variable

      real, intent(out) ::                                              &
     &                   outvar(imin:imax,jmin:jmax,kmin:kmax,nmin:nmax)
                       ! Copied variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_copy4d = 0
      integer, parameter :: DUMP_TARGET_copy4d = 3
      logical, save :: dump_done_copy4d = .false.


!-----7--------------------------------------------------------------7--

! Copy the invar to the outvar.

!@llm start meta_info ----------------------------------------------------
! Location: copy4d.f90 :: s_copy4d
! Summary : Simple 4D array copy from invar to outvar.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls
!   - Trivial memory copy operation
!   - Independent element-wise operations
!   - Outer n,k loops are sequential in current OpenMP structure
! Next:
!   - Collapse all four loops (n,k,j,i) for better GPU occupancy
!   - Consider using device-to-device memcpy for efficiency
!   - May be better to keep data resident on GPU and avoid copy calls
! Runtime:
!   - Calls: 3
!   - AvgLoops: 277.1M
!   - TotalTime: 0.026s (0.00%)
!   - AvgTime: 8.662ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('copy4d.f90', 's_copy4d', &
   & 'OMP section 1')
end if
loop_len = int((nmax)-(nmin)+1,8) &
     & * int((kmax)-(kmin)+1,8) &
     & * int((jmax)-(jmin)+1,8) &
     & * int((imax)-(imin)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_copy4d = dump_call_count_copy4d + 1
if (dump_call_count_copy4d == DUMP_TARGET_copy4d .and. .not. dump_done_copy4d) then
  call dump_init('copy4d')
  call dump_scalar_i('imin', imin)
  call dump_scalar_i('imax', imax)
  call dump_scalar_i('jmin', jmin)
  call dump_scalar_i('jmax', jmax)
  call dump_scalar_i('kmin', kmin)
  call dump_scalar_i('kmax', kmax)
  call dump_scalar_i('nmin', nmin)
  call dump_scalar_i('nmax', nmax)
  call dump_array_4d('invar.bin', invar, imin, imax, jmin, jmax, kmin, kmax, nmin, nmax)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_068)
! GPU version (OpenACC)
!$acc kernels
!$acc loop independent
      do n=nmin,nmax
!$acc loop independent
        do k=kmin,kmax
!$acc loop independent
          do j=jmin,jmax
!$acc loop independent
          do i=imin,imax
            outvar(i,j,k,n)=invar(i,j,k,n)
          end do
          end do
        end do
      end do
!$acc end kernels
#else
! CPU version (OpenMP) - Original code preserved
!$omp parallel default(shared) private(k,n)

      do n=nmin,nmax

        do k=kmin,kmax

!$omp do schedule(runtime) private(i,j)

          do j=jmin,jmax
          do i=imin,imax
            outvar(i,j,k,n)=invar(i,j,k,n)
          end do
          end do

!$omp end do

        end do

      end do

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_copy4d == DUMP_TARGET_copy4d .and. .not. dump_done_copy4d) then
  call dump_array_4d('outvar_ref.bin', outvar, imin, imax, jmin, jmax, kmin, kmax, nmin, nmax)
  call dump_finalize()
  dump_done_copy4d = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_copy4d

!-----7--------------------------------------------------------------7--

      end module m_copy4d

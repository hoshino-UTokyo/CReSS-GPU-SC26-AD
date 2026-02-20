!***********************************************************************
      module m_trilat
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/05/10, 1999/08/23,
!                   1999/09/30, 1999/10/12, 1999/11/01, 2000/01/17,
!                   2000/07/05, 2000/12/18, 2001/03/13, 2002/04/02,
!                   2003/01/04, 2003/04/30, 2003/05/19, 2003/10/31,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the Coriolis parameters x 0.25.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_dump_kernel
      use m_comphy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: trilat, s_trilat

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface trilat

        module procedure s_trilat

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic sin
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_trilat(fpcoropt,ni,nj,lat,fc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcoropt
                       ! Formal parameter of unique index of coropt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: lat(0:ni+1,0:nj+1)
                       ! Latitude

! Output variables

      real, intent(out) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

! Internal shared variables

      integer coropt   ! Option for Coriolis force

      real omega5      ! 0.5 x omega

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real sinlat      ! sin(latitude)


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_trilat = 0
      integer, parameter :: DUMP_TARGET_trilat = 1
      logical, save :: dump_done_trilat = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpcoropt,coropt)

! -----

! Calculate 0.5 omega.

      omega5=.5e0*omega

! -----

! Calculate the Coriolis parameters x 0.25.

!@llm start meta_info ----------------------------------------------------
! Location: trilat.f90 :: s_trilat
! Summary : Calculate Coriolis parameters (fc) from latitude using
!           sin/sqrt for vertical and horizontal components
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Calls intrinsic functions only (sin, sqrt)
!   - Writes to fc array (2 components)
!   - Conditional branches based on coropt (1 or 2)
!   - No synchronization constructs
! Next:
!   - Straightforward GPU port
!   - Trigonometric functions available on GPU
!   - Can be computed once and cached if lat doesn't change
! Runtime:
!   - Calls: 1
!   - AvgLoops: 806.4K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.042ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('trilat.f90', 's_trilat', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_trilat = dump_call_count_trilat + 1
if (dump_call_count_trilat == DUMP_TARGET_trilat .and. .not. dump_done_trilat) then
  call dump_init('trilat')
  ! FIXME: coropt is array - call dump_scalar_i('coropt', coropt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_array_2d('lat.bin', lat, 0, ni+1, 0, nj+1)
  call dump_scalar_r('omega5', omega5)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_333)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    if (coropt == 1) then
      !$acc kernels
      !$acc loop independent collapse(2) private(sinlat)
      do j = 1, nj-1
        do i = 1, ni-1
          sinlat = sin(lat(i,j) * d2r)
          fc(i,j,1) = omega5 * sinlat
          fc(i,j,2) = 0.e0
        end do
      end do
      !$acc end kernels
    else if (coropt == 2) then
      !$acc kernels
      !$acc loop independent collapse(2) private(sinlat)
      do j = 1, nj-1
        do i = 1, ni-1
          sinlat = sin(lat(i,j) * d2r)
          fc(i,j,1) = omega5 * sinlat
          fc(i,j,2) = omega5 * sqrt(1.e0 - sinlat * sinlat)
        end do
      end do
      !$acc end kernels
    end if

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared)

      if(coropt.eq.1) then

!$omp do schedule(runtime) private(i,j,sinlat)

        do j=1,nj-1
        do i=1,ni-1
          sinlat=sin(lat(i,j)*d2r)

          fc(i,j,1)=omega5*sinlat
          fc(i,j,2)=0.e0

        end do
        end do

!$omp end do

      else if(coropt.eq.2) then

!$omp do schedule(runtime) private(i,j,sinlat)

        do j=1,nj-1
        do i=1,ni-1
          sinlat=sin(lat(i,j)*d2r)

          fc(i,j,1)=omega5*sinlat
          fc(i,j,2)=omega5*sqrt(1.e0-sinlat*sinlat)

        end do
        end do

!$omp end do

      end if

!$omp end parallel

#endif

! Dump output data at target call
if (dump_call_count_trilat == DUMP_TARGET_trilat .and. .not. dump_done_trilat) then
  call dump_array_3d('fc_ref.bin', fc, 0, ni+1, 0, nj+1, 1, 2)
  call dump_finalize()
  dump_done_trilat = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_trilat

!-----7--------------------------------------------------------------7--

      end module m_trilat

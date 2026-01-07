!***********************************************************************
      module m_advbspi
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/03, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/07/05, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2003/12/26, 2004/08/20, 2005/08/05,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2010/02/01, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the base state pressure advection for the horizontally
!     explicit and vertically implicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_getrname
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: advbspi, s_advbspi

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface advbspi

        module procedure s_advbspi

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
      subroutine s_advbspi(fpweicoe,fproc,ni,nj,nk,rst,w,pfrc,          &
     &                     phdiv,pvdiv,fp)
!***********************************************************************

! Input variables

      character(len=4), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpweicoe
                       ! Formal parameter of unique index of weicoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure x Jacobian

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: pfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation
                       ! in large time steps

      real, intent(in) :: phdiv(0:ni+1,0:nj+1,1:nk)
                       ! Horizontal divergence value
                       ! in pressure equation

      real, intent(in) :: pvdiv(0:ni+1,0:nj+1,1:nk)
                       ! Vertical divergence value
                       ! in pressure equation

! Input and output variable

      real, intent(inout) :: fp(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

! Internal shared variables

      real weicoe      ! Weighting coefficient for implicit method

      real weic1m      ! 1.0 - weicoe

      real g05         ! 0.5 x g

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_advbspi = 0
      integer, parameter :: DUMP_TARGET_advbspi = 28800
      logical, save :: dump_done_advbspi = .false.

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpweicoe,weicoe)

! -----

! Set the common used variables.

      weic1m=1.e0-weicoe

      g05=.5e0*g

! -----

! Calculate the base state pressure advection.

!@llm start meta_info ----------------------------------------------------
! Location: advbspi.f90 :: subroutine s_advbspi
! Summary : Calculates base state pressure advection for horizontally
!           explicit and vertically implicit method (forward/backward).
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - Calls getrname() before parallel region (not inside).
!   - Uses module constant g from comphy.
!   - Conditional on fproc determines backward vs forward calculation.
!   - Pure arithmetic, all GPU compatible.
!   - All grid points independent (embarrassingly parallel).
! Next:
!   - Direct OpenACC kernels with collapse(3) for (k,j,i).
!   - May split into two kernels for back/fore branches.
! Runtime:
!   - Calls: 28800
!   - AvgLoops: 100.4M
!   - TotalTime: 172.825s (5.80%)
!   - AvgTime: 6.001ms
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('advbspi.f90', 's_advbspi', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8) &
     & * int((nj-2)-(2)+1,8) &
     & * int((ni-2)-(2)+1,8)

! Dump input data at target call
dump_call_count_advbspi = dump_call_count_advbspi + 1
if (dump_call_count_advbspi == DUMP_TARGET_advbspi .and. .not. dump_done_advbspi) then
  call dump_init('advbspi')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('weicoe', weicoe)
  call dump_scalar_r('g', g)
  call dump_scalar_c('fproc', fproc)
  call dump_array_3d('rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pfrc.bin', pfrc, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('phdiv.bin', phdiv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pvdiv.bin', pvdiv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('fp_in.bin', fp, 0, ni+1, 0, nj+1, 1, nk)
end if

call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

      if(fproc(1:4).eq.'back') then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
           fp(i,j,k)=pfrc(i,j,k)+phdiv(i,j,k)                           &
     &       +weic1m*(g05*(w(i,j,k)+w(i,j,k+1))*rst(i,j,k)+pvdiv(i,j,k))
          end do
          end do

!$omp end do

        end do

      else if(fproc(1:4).eq.'fore') then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
           fp(i,j,k)=fp(i,j,k)                                          &
     &       +weicoe*(g05*(w(i,j,k)+w(i,j,k+1))*rst(i,j,k)+pvdiv(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! Dump output data at target call
if (dump_call_count_advbspi == DUMP_TARGET_advbspi .and. .not. dump_done_advbspi) then
  call dump_array_3d('fp_ref.bin', fp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_advbspi = .true.
end if

! -----

      end subroutine s_advbspi

!-----7--------------------------------------------------------------7--

      end module m_advbspi

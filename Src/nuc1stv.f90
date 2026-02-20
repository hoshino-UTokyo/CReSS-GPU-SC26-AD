!***********************************************************************
      module m_nuc1stv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2001/06/29, 2001/10/18,
!                   2001/11/20, 2001/12/11, 2002/01/07, 2002/01/15,
!                   2002/04/02, 2002/12/06, 2003/01/04, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/11/05, 2003/12/12,
!                   2004/03/22, 2004/08/01, 2004/09/01, 2004/09/10,
!                   2004/10/12, 2004/12/17, 2005/04/04, 2006/02/13,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the nucleation rate of the deposition or sorption.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_dump_kernel

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: nuc1stv, s_nuc1stv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface nuc1stv

        module procedure s_nuc1stv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_nuc1stv(thresq,ni,nj,nk,rbv,qv,qi,t,qvsi,nuvi)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: qvsi(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio for ice

! Output variable

      real, intent(out) :: nuvi(0:ni+1,0:nj+1,1:nk)
                       ! Nucleation rate of deposition or sorption

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real qvssi       ! Super saturation mixing ratio for ice

      real a           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_nuc1stv = 0
      integer, parameter :: DUMP_TARGET_nuc1stv = 45720
      logical, save :: dump_done_nuc1stv = .false.


!-----7--------------------------------------------------------------7--

!! Calculate the nucleation rate of the deposition or sorption.

!@llm start meta_info ----------------------------------------------------
! Location: nuc1stv.f90 :: s_nuc1stv
! Summary : Calculate nucleation rate of deposition/sorption for ice
!           based on supersaturation and temperature conditions.
! GPU diff: Easy
! Findings:
!   - Conditional branch for nk.eq.1 vs nk.gt.1 cases
!   - Simple conditionals on qv, qvsi, t values
!   - Uses exp, max, min intrinsics - GPU compatible
!   - Output array nuvi written independently per grid point
!   - No inter-thread dependencies; fully parallel
! Next:
!   - Direct port to GPU kernel with minimal changes
!   - Branch divergence from conditionals is manageable
!   - Consider using predication for conditional assignments
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 1.910s (0.06%)
!   - AvgTime: 0.042ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('nuc1stv.f90', 's_nuc1stv', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_nuc1stv = dump_call_count_nuc1stv + 1
if (dump_call_count_nuc1stv == DUMP_TARGET_nuc1stv .and. .not. dump_done_nuc1stv) then
  call dump_init('nuc1stv')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('thresq', thresq)
  call dump_array_3d('t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('t0', t0)
  call dump_array_3d('rbv.bin', rbv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qi.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvsi.bin', qvsi, 0, ni+1, 0, nj+1, 1, nk)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_213)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    if (nk == 1) then
      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent private(qvssi, a)
        do i = 1, ni-1
          if (qv(i,j,1) > thresq) then
            if (qv(i,j,1) > qvsi(i,j,1) .and. t(i,j,1) > tlow .and. t(i,j,1) < t0) then
              qvssi = qv(i,j,1) - qvsi(i,j,1)
              a = 15.25 * qv(i,j,1) / qvsi(i,j,1) - 10.08
              if (a < 40.0) then
                nuvi(i,j,1) = max(min(mi0 * exp(a) * rbv(i,j,1) - qi(i,j,1), qvssi), 0.0)
              else
                nuvi(i,j,1) = qvssi
              end if
            else
              nuvi(i,j,1) = 0.0
            end if
          else
            nuvi(i,j,1) = 0.0
          end if
        end do
      end do
      !$acc end kernels

    else
      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(qvssi, a)
          do i = 1, ni-1
            if (qv(i,j,k) > thresq) then
              if (qv(i,j,k) > qvsi(i,j,k) .and. t(i,j,k) > tlow .and. t(i,j,k) < t0) then
                qvssi = qv(i,j,k) - qvsi(i,j,k)
                a = 15.25 * qv(i,j,k) / qvsi(i,j,k) - 10.08
                if (a < 40.0) then
                  nuvi(i,j,k) = max(min(mi0 * exp(a) * rbv(i,j,k) - qi(i,j,k), qvssi), 0.0)
                else
                  nuvi(i,j,k) = qvssi
                end if
              else
                nuvi(i,j,k) = 0.0
              end if
            else
              nuvi(i,j,k) = 0.0
            end if
          end do
        end do
      end do
      !$acc end kernels
    end if

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared) private(k)

! In the case nk = 1.

      if(nk.eq.1) then

!$omp do schedule(runtime) private(i,j,qvssi,a)

        do j=1,nj-1
        do i=1,ni-1

          if(qv(i,j,1).gt.thresq) then

            if(qv(i,j,1).gt.qvsi(i,j,1)                                 &
     &        .and.t(i,j,1).gt.tlow.and.t(i,j,1).lt.t0) then

              qvssi=qv(i,j,1)-qvsi(i,j,1)

              a=15.25e0*qv(i,j,1)/qvsi(i,j,1)-10.08e0

              if(a.lt.40.e0) then

                nuvi(i,j,1)                                             &
     &            =max(min(mi0*exp(a)*rbv(i,j,1)-qi(i,j,1),qvssi),0.e0)

              else

                nuvi(i,j,1)=qvssi

              end if

            else

              nuvi(i,j,1)=0.e0

            end if

          else

            nuvi(i,j,1)=0.e0

          end if

        end do
        end do

!$omp end do

! -----

! In the case nk > 1.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,qvssi,a)

          do j=1,nj-1
          do i=1,ni-1

            if(qv(i,j,k).gt.thresq) then

              if(qv(i,j,k).gt.qvsi(i,j,k)                               &
     &          .and.t(i,j,k).gt.tlow.and.t(i,j,k).lt.t0) then

                qvssi=qv(i,j,k)-qvsi(i,j,k)

                a=15.25e0*qv(i,j,k)/qvsi(i,j,k)-10.08e0

                if(a.lt.40.e0) then

                  nuvi(i,j,k)=max(0.e0,                                 &
     &              min(mi0*exp(a)*rbv(i,j,k)-qi(i,j,k),qvssi))

                else

                  nuvi(i,j,k)=qvssi

                end if

              else

                nuvi(i,j,k)=0.e0

              end if

            else

              nuvi(i,j,k)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_nuc1stv == DUMP_TARGET_nuc1stv .and. .not. dump_done_nuc1stv) then
  call dump_array_3d('nuvi_ref.bin', nuvi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_nuc1stv = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_nuc1stv

!-----7--------------------------------------------------------------7--

      end module m_nuc1stv

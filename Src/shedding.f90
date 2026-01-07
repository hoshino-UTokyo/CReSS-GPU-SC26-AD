!***********************************************************************
      module m_shedding
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2001/06/29, 2001/10/18, 2001/11/20,
!                   2002/01/15, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the shedding rate from the snow and graupel to the rain
!     water.

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

      public :: shedding, s_shedding

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface shedding

        module procedure s_shedding

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
      subroutine s_shedding(thresq,ni,nj,nk,qs,qg,t,clcs,clcg,clrs,clrg,&
     &                      clig,clsg,pgwet,shsr,shgr)
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

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and snow

      real, intent(in) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and graupel

      real, intent(in) :: clrs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between rain water and snow

      real, intent(in) :: clrg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between rain water and graupel

      real, intent(in) :: clig(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud ice and graupel

      real, intent(in) :: clsg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between snow and graupel

      real, intent(in) :: pgwet(0:ni+1,0:nj+1,1:nk)
                       ! Graupel production rate for moist process

! Output variables

      real, intent(out) :: shsr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from snow

      real, intent(out) :: shgr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from graupel

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_shedding = 0
      integer, parameter :: DUMP_TARGET_shedding = 45720
      logical, save :: dump_done_shedding = .false.


!-----7--------------------------------------------------------------7--

!!! Calculate the shedding rate.

!@llm start meta_info ----------------------------------------------------
! Location: shedding.f90 :: s_shedding
! Summary : Calculates shedding rates of liquid water from snow and graupel to rain,
!           based on temperature and collection/production rates.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses module constant t0 from m_comphy for temperature threshold
!   - Conditional branches based on mixing ratio thresholds and temperature
!   - Handles nk=1 case separately (2D) vs nk>1 case (3D)
!   - All loops independent with private i,j,k indices
!   - No synchronization constructs
! Next:
!   - Straightforward GPU port with conditional logic preserved
!   - Use OpenACC/OpenACC with collapse for nested loops
!   - Consider single kernel handling both nk cases with runtime check
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 1.185s (0.04%)
!   - AvgTime: 0.026ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('shedding.f90', 's_shedding', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_shedding = dump_call_count_shedding + 1
if (dump_call_count_shedding == DUMP_TARGET_shedding .and. .not. dump_done_shedding) then
  call dump_init('shedding')
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('thresq', thresq)
  call dump_array_3d('t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('t0', t0)
  call dump_array_3d('qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qg.bin', qg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcs.bin', clcs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clcg.bin', clcg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrs.bin', clrs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrg.bin', clrg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clig.bin', clig, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsg.bin', clsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pgwet.bin', pgwet, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

!! In the case nk = 1.

      if(nk.eq.1) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

! Calculate the shedding rate from the snow to the rain water.

          if(qs(i,j,1).gt.thresq) then

            if(t(i,j,1).ge.t0) then

              shsr(i,j,1)=clcs(i,j,1)+clrs(i,j,1)

            else

              shsr(i,j,1)=0.e0

            end if

          else

            shsr(i,j,1)=0.e0

          end if

! -----

! Calculate the shedding rate from the graupel to the rain water.

          if(qg(i,j,1).gt.thresq) then

            if(t(i,j,1).ge.t0) then

              shgr(i,j,1)=clcg(i,j,1)+clrg(i,j,1)

            else

              if(pgwet(i,j,1).gt.0.e0) then

                shgr(i,j,1)=(clcg(i,j,1)+clrg(i,j,1)                    &
     &            +clig(i,j,1)+clsg(i,j,1))-pgwet(i,j,1)

              else

                shgr(i,j,1)=0.e0

              end if

            end if

          else

            shgr(i,j,1)=0.e0

          end if

! -----

        end do
        end do

!$omp end do

!! -----

!! In the case nk > 1.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the shedding rate from the snow to the rain water.

            if(qs(i,j,k).gt.thresq) then

              if(t(i,j,k).ge.t0) then

                shsr(i,j,k)=clcs(i,j,k)+clrs(i,j,k)

              else

                shsr(i,j,k)=0.e0

              end if

            else

              shsr(i,j,k)=0.e0

            end if

! -----

! Calculate the shedding rate from the graupel to the rain water.

            if(qg(i,j,k).gt.thresq) then

              if(t(i,j,k).ge.t0) then

                shgr(i,j,k)=clcg(i,j,k)+clrg(i,j,k)

              else

                if(pgwet(i,j,k).gt.0.e0) then

                  shgr(i,j,k)=(clcg(i,j,k)+clrg(i,j,k)                  &
     &              +clig(i,j,k)+clsg(i,j,k))-pgwet(i,j,k)

                else

                  shgr(i,j,k)=0.e0

                end if

              end if

            else

              shgr(i,j,k)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

        end do

      end if

!! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_shedding == DUMP_TARGET_shedding .and. .not. dump_done_shedding) then
  call dump_array_3d('shsr_ref.bin', shsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('shgr_ref.bin', shgr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_shedding = .true.
end if


call profile_stop(prof_id1, loop_len)

!!! -----

      end subroutine s_shedding

!-----7--------------------------------------------------------------7--

      end module m_shedding

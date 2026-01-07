!***********************************************************************
      module m_distrpg
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2001/05/29, 2001/06/29, 2001/11/20,
!                   2002/04/02, 2002/09/09, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2004/06/10, 2004/09/01, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/04/04, 2005/09/30,
!                   2007/10/19, 2007/11/26, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/03/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the distribution ratio at which the collisions between
!     rain water and snow and reset the collection rate.

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

      public :: distrpg, s_distrpg

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface distrpg

        module procedure s_distrpg

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_distrpg(cphopt,thresq,ni,nj,nk,qs,t,diaqr,diaqs,     &
     &                     clrs,clsr,clrsn,clsrn,clrsg)
!***********************************************************************

! Input variables

      integer, intent(in) :: cphopt
                       ! Option for cloud micro physics

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

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: diaqr(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of rain water

      real, intent(in) :: diaqs(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of snow

! Input and output variables

      real, intent(inout) :: clrs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate from rain water to snow

      real, intent(inout) :: clsr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate from snow to rain water

      real, intent(inout) :: clrsn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, intent(inout) :: clsrn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

! Output variable

      real, intent(out) :: clrsg(0:ni+1,0:nj+1,1:nk)
                       ! Production rate of graupel
                       ! from collection rate form rain to snow

! Internal shared variables

      real rhos2       ! rhos x rhos
      real rhow2       ! rhow x rhow

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real alpha       ! Distribution ratio between snow and graupel

      real alpha1      ! 1.0 - alpha

      real a           ! Temporary variable
      real b           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_distrpg = 0
      integer, parameter :: DUMP_TARGET_distrpg = 45720
      logical, save :: dump_done_distrpg = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rhos2=rhos*rhos
      rhow2=rhow*rhow

! -----

!!! Calculate the distribution ratio at which the collisions between
!!! rain water and snow and reset the collection rate.

!@llm start meta_info ----------------------------------------------------
! Location: distrpg.f90 :: s_distrpg
! Summary : Distribute collision rates between rain and snow to graupel,
!           based on diameter ratios and temperature thresholds.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region (only intrinsic: abs)
!   - No global/module variable writes
!   - No synchronization constructs
!   - Multiple branches (nk, cphopt) but all loops are data-parallel
!   - Conditional updates per grid point (temperature check, threshold)
! Next:
!   - Direct OpenACC with collapse(2) on j-i loops
!   - Conditionals inside loop are fine for GPU (divergent but manageable)
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 1.393s (0.05%)
!   - AvgTime: 0.030ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('distrpg.f90', 's_distrpg', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_distrpg = dump_call_count_distrpg + 1
if (dump_call_count_distrpg == DUMP_TARGET_distrpg .and. .not. dump_done_distrpg) then
  call dump_init('distrpg')
  call dump_scalar_i('cphopt', cphopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('thresq', thresq)
  call dump_array_3d('t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('t0', t0)
  call dump_array_3d('qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('diaqr.bin', diaqr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('diaqs.bin', diaqs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrs_in.bin', clrs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsr_in.bin', clsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrsn_in.bin', clrsn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsrn_in.bin', clsrn, 0, ni+1, 0, nj+1, 1, nk)
end if

!$omp parallel default(shared) private(k)

!! In the case nk = 1.

      if(nk.eq.1) then

! Perform calculating in the case the option abs(cphopt) is equal to 2.

        if(abs(cphopt).eq.2) then

!$omp do schedule(runtime) private(i,j,alpha,alpha1,a,b)

          do j=1,nj-1
          do i=1,ni-1

            if(qs(i,j,1).gt.thresq) then

              if(t(i,j,1).lt.t0) then

                a=diaqs(i,j,1)*diaqs(i,j,1)
                b=diaqr(i,j,1)*diaqr(i,j,1)

                a=rhos2*a*a*a
                b=rhow2*b*b*b

                alpha=a/(a+b)
                alpha1=1.e0-alpha

                clrsg(i,j,1)=alpha1*clrs(i,j,1)

                clrs(i,j,1)=alpha*clrs(i,j,1)

                clsr(i,j,1)=alpha1*clsr(i,j,1)

              else

                clrsg(i,j,1)=clrs(i,j,1)

                clrs(i,j,1)=0.e0

              end if

            else

              clrsg(i,j,1)=0.e0

            end if

          end do
          end do

!$omp end do

! -----

! Perform calculating in the case the option abs(cphopt) is greater
! than 2.

        else if(abs(cphopt).ge.3) then

!$omp do schedule(runtime) private(i,j,alpha,alpha1,a,b)

          do j=1,nj-1
          do i=1,ni-1

            if(qs(i,j,1).gt.thresq) then

              if(t(i,j,1).lt.t0) then

                a=diaqs(i,j,1)*diaqs(i,j,1)
                b=diaqr(i,j,1)*diaqr(i,j,1)

                a=rhos2*a*a*a
                b=rhow2*b*b*b

                alpha=a/(a+b)
                alpha1=1.e0-alpha

                clrsg(i,j,1)=alpha1*clrs(i,j,1)

                clrs(i,j,1)=alpha*clrs(i,j,1)

                clsr(i,j,1)=alpha1*clsr(i,j,1)

                clrsn(i,j,1)=alpha1*clrsn(i,j,1)
                clsrn(i,j,1)=alpha1*clsrn(i,j,1)

              else

                clrsg(i,j,1)=clrs(i,j,1)

                clrs(i,j,1)=0.e0

              end if

            else

              clrsg(i,j,1)=0.e0

            end if

          end do
          end do

!$omp end do

        end if

! -----

!! -----

!! In the case nk > 1.

      else

! Perform calculating in the case the option abs(cphopt) is equal to 2.

        if(abs(cphopt).eq.2) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,alpha,alpha1,a,b)

            do j=1,nj-1
            do i=1,ni-1

              if(qs(i,j,k).gt.thresq) then

                if(t(i,j,k).lt.t0) then

                  a=diaqs(i,j,k)*diaqs(i,j,k)
                  b=diaqr(i,j,k)*diaqr(i,j,k)

                  a=rhos2*a*a*a
                  b=rhow2*b*b*b

                  alpha=a/(a+b)
                  alpha1=1.e0-alpha

                  clrsg(i,j,k)=alpha1*clrs(i,j,k)

                  clrs(i,j,k)=alpha*clrs(i,j,k)

                  clsr(i,j,k)=alpha1*clsr(i,j,k)

                else

                  clrsg(i,j,k)=clrs(i,j,k)

                  clrs(i,j,k)=0.e0

                end if

              else

                clrsg(i,j,k)=0.e0

              end if

            end do
            end do

!$omp end do

          end do

! -----

! Perform calculating in the case the option abs(cphopt) is greater
! than 2.

        else if(abs(cphopt).ge.3) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,alpha,alpha1,a,b)

            do j=1,nj-1
            do i=1,ni-1

              if(qs(i,j,k).gt.thresq) then

                if(t(i,j,k).lt.t0) then

                  a=diaqs(i,j,k)*diaqs(i,j,k)
                  b=diaqr(i,j,k)*diaqr(i,j,k)

                  a=rhos2*a*a*a
                  b=rhow2*b*b*b

                  alpha=a/(a+b)
                  alpha1=1.e0-alpha

                  clrsg(i,j,k)=alpha1*clrs(i,j,k)

                  clrs(i,j,k)=alpha*clrs(i,j,k)

                  clsr(i,j,k)=alpha1*clsr(i,j,k)

                  clrsn(i,j,k)=alpha1*clrsn(i,j,k)
                  clsrn(i,j,k)=alpha1*clsrn(i,j,k)

                else

                  clrsg(i,j,k)=clrs(i,j,k)

                  clrs(i,j,k)=0.e0

                end if

              else

                clrsg(i,j,k)=0.e0

              end if

            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!$omp end parallel

! Dump output data at target call
if (dump_call_count_distrpg == DUMP_TARGET_distrpg .and. .not. dump_done_distrpg) then
  call dump_array_3d('clrsg_ref.bin', clrsg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrs_ref.bin', clrs, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsr_ref.bin', clsr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clrsn_ref.bin', clrsn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('clsrn_ref.bin', clsrn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_distrpg = .true.
end if


call profile_stop(prof_id1, loop_len)

!!! -----

      end subroutine s_distrpg

!-----7--------------------------------------------------------------7--

      end module m_distrpg

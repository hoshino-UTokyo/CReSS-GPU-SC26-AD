!***********************************************************************
      module m_freezing
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2000/11/17, 2001/06/29,
!                   2001/10/18, 2001/11/20, 2001/12/11, 2002/01/15,
!                   2002/04/02, 2002/12/02, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/10/31, 2003/11/05, 2003/12/12,
!                   2004/04/01, 2004/05/31, 2004/06/10, 2004/07/10,
!                   2004/09/01, 2004/09/10, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2005/01/31, 2005/04/04, 2006/04/03,
!                   2007/10/19, 2007/11/26, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/11/13, 2011/03/18, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the freezing rate from the rain water to the graupel.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_dump_kernel
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: freezing, s_freezing

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface freezing

        module procedure s_freezing

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic exp
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_freezing(cphopt,dtb,thresq,ni,nj,nk,qr,ncr,tcel,     &
     &                      diaqr,frrg,frrgn)
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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

      real, intent(in) :: ncr(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water

      real, intent(in) :: tcel(0:ni+1,0:nj+1,1:nk)
                       ! Ambient air temperature

      real, intent(in) :: diaqr(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of rain water

! Output variables

      real, intent(out) :: frrg(0:ni+1,0:nj+1,1:nk)
                       ! Freezing rate from rain water to graupel

      real, intent(out) :: frrgn(0:ni+1,0:nj+1,1:nk)
                       ! Freezing rate for concentrations
                       ! from rain water to graupel

! Internal shared variables

      real tclow       ! tlow - t0

      real cfrrg       ! Coefficient of freezing rate

      real cfrrgn      ! Coefficient of freezing rate
                       ! for concentrations

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real diaqr3      ! diaqr^3

      real a           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_freezing = 0
      integer, parameter :: DUMP_TARGET_freezing = 45720
      logical, save :: dump_done_freezing = .false.


!-----7--------------------------------------------------------------7--

! Set the common used variables.

      tclow=tlow-t0

      cfrrg=20.e0*100.e0*cc*cc*rhow*dtb
      cfrrgn=100.e0*cc*dtb

! -----

!!! Calculate the freezing rate from the rain water to the graupel.

!@llm start meta_info ----------------------------------------------------
! Location: freezing.f90 :: s_freezing
! Summary : Calculate freezing rate from rain water to graupel based on
!           temperature and rain water content with different cphopt modes.
! GPU diff: Medium
! Findings:
!   - Complex conditional logic (nk==1 vs nk>1, cphopt==2 vs >=3)
!   - Uses exp intrinsic function
!   - Multiple threshold checks (thresq, t0cel, tclow)
!   - Writes to frrg and frrgn arrays
!   - No reductions or synchronization
! Next:
!   - Port with GPU kernels handling conditionals via masks or separate kernels
!   - Consider separating cphopt==2 and cphopt>=3 cases for clarity
! Runtime:
!   - Calls: 45720
!   - AvgLoops: 806.4K
!   - TotalTime: 1.296s (0.04%)
!   - AvgTime: 0.028ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('freezing.f90', 's_freezing', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_freezing = dump_call_count_freezing + 1
if (dump_call_count_freezing == DUMP_TARGET_freezing .and. .not. dump_done_freezing) then
  call dump_init('freezing')
  call dump_scalar_i('cphopt', cphopt)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('dtb', dtb)
  call dump_scalar_r('thresq', thresq)
  call dump_array_3d('qr.bin', qr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ncr.bin', ncr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('tcel.bin', tcel, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('diaqr.bin', diaqr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('cfrrg', cfrrg)
  call dump_scalar_r('cfrrgn', cfrrgn)
  call dump_scalar_r('tclow', tclow)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_114)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------
    cfrrgn = 100.e0 * cc * dtb

    if (nk == 1) then
      if (abs(cphopt) == 2) then
        !$acc kernels
        !$acc loop independent collapse(2) private(diaqr3)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qr(i,j,1) > thresq) then
              if (tcel(i,j,1) < t0cel) then
                if (tcel(i,j,1) <= tclow) then
                  frrg(i,j,1) = qr(i,j,1)
                else
                  diaqr3 = diaqr(i,j,1) * diaqr(i,j,1) * diaqr(i,j,1)
                  frrg(i,j,1) = min(cfrrg * diaqr3 * diaqr3 * &
                    (exp(-0.66e0 * tcel(i,j,1)) - 1.e0) * ncr(i,j,1), qr(i,j,1))
                end if
              else
                frrg(i,j,1) = 0.e0
              end if
            else
              frrg(i,j,1) = 0.e0
            end if
          end do
        end do
        !$acc end kernels
      else if (abs(cphopt) >= 3) then
        !$acc kernels
        !$acc loop independent collapse(2) private(diaqr3, a)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qr(i,j,1) > thresq) then
              if (tcel(i,j,1) < t0cel) then
                if (tcel(i,j,1) <= tclow) then
                  frrg(i,j,1) = qr(i,j,1)
                  frrgn(i,j,1) = ncr(i,j,1)
                else
                  a = (exp(-0.66e0 * tcel(i,j,1)) - 1.e0) * ncr(i,j,1)
                  diaqr3 = diaqr(i,j,1) * diaqr(i,j,1) * diaqr(i,j,1)
                  frrg(i,j,1) = min(cfrrg * diaqr3 * diaqr3 * a, qr(i,j,1))
                  frrgn(i,j,1) = min(cfrrgn * diaqr3 * a, ncr(i,j,1))
                end if
              else
                frrg(i,j,1) = 0.e0
                frrgn(i,j,1) = 0.e0
              end if
            else
              frrg(i,j,1) = 0.e0
              frrgn(i,j,1) = 0.e0
            end if
          end do
        end do
        !$acc end kernels
      end if
    else
      if (abs(cphopt) == 2) then
        !$acc kernels
        !$acc loop independent collapse(3) private(diaqr3)
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni-1
              if (qr(i,j,k) > thresq) then
                if (tcel(i,j,k) < t0cel) then
                  if (tcel(i,j,k) <= tclow) then
                    frrg(i,j,k) = qr(i,j,k)
                  else
                    diaqr3 = diaqr(i,j,k) * diaqr(i,j,k) * diaqr(i,j,k)
                    frrg(i,j,k) = min(qr(i,j,k), cfrrg * diaqr3 * diaqr3 * &
                      (exp(-0.66e0 * tcel(i,j,k)) - 1.e0) * ncr(i,j,k))
                  end if
                else
                  frrg(i,j,k) = 0.e0
                end if
              else
                frrg(i,j,k) = 0.e0
              end if
            end do
          end do
        end do
        !$acc end kernels
      else if (abs(cphopt) >= 3) then
        !$acc kernels
        !$acc loop independent collapse(3) private(diaqr3, a)
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni-1
              if (qr(i,j,k) > thresq) then
                if (tcel(i,j,k) < t0cel) then
                  if (tcel(i,j,k) <= tclow) then
                    frrg(i,j,k) = qr(i,j,k)
                    frrgn(i,j,k) = ncr(i,j,k)
                  else
                    a = (exp(-0.66e0 * tcel(i,j,k)) - 1.e0) * ncr(i,j,k)
                    diaqr3 = diaqr(i,j,k) * diaqr(i,j,k) * diaqr(i,j,k)
                    frrg(i,j,k) = min(cfrrg * diaqr3 * diaqr3 * a, qr(i,j,k))
                    frrgn(i,j,k) = min(cfrrgn * diaqr3 * a, ncr(i,j,k))
                  end if
                else
                  frrg(i,j,k) = 0.e0
                  frrgn(i,j,k) = 0.e0
                end if
              else
                frrg(i,j,k) = 0.e0
                frrgn(i,j,k) = 0.e0
              end if
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared) private(k)

!! In the case nk = 1.

      if(nk.eq.1) then

! Perform calculating in the case the option abs(cphopt) is equal to 2.

        if(abs(cphopt).eq.2) then

!$omp do schedule(runtime) private(i,j,diaqr3)

          do j=1,nj-1
          do i=1,ni-1

            if(qr(i,j,1).gt.thresq) then

              if(tcel(i,j,1).lt.t0cel) then

                if(tcel(i,j,1).le.tclow) then

                  frrg(i,j,1)=qr(i,j,1)

                else

                  diaqr3=diaqr(i,j,1)*diaqr(i,j,1)*diaqr(i,j,1)

                  frrg(i,j,1)=min(cfrrg*diaqr3*diaqr3                   &
     &             *(exp(-.66e0*tcel(i,j,1))-1.e0)*ncr(i,j,1),qr(i,j,1))

                end if

              else

                frrg(i,j,1)=0.e0

              end if

            else

              frrg(i,j,1)=0.e0

            end if

          end do
          end do

!$omp end do

! -----

! Perform calculating in the case the option abs(cphopt) is greater
! than 2.

        else if(abs(cphopt).ge.3) then

!$omp do schedule(runtime) private(i,j,diaqr3,a)

          do j=1,nj-1
          do i=1,ni-1

            if(qr(i,j,1).gt.thresq) then

              if(tcel(i,j,1).lt.t0cel) then

                if(tcel(i,j,1).le.tclow) then

                  frrg(i,j,1)=qr(i,j,1)
                  frrgn(i,j,1)=ncr(i,j,1)

                else

                  a=(exp(-.66e0*tcel(i,j,1))-1.e0)*ncr(i,j,1)

                  diaqr3=diaqr(i,j,1)*diaqr(i,j,1)*diaqr(i,j,1)

                  frrg(i,j,1)=min(cfrrg*diaqr3*diaqr3*a,qr(i,j,1))
                  frrgn(i,j,1)=min(cfrrgn*diaqr3*a,ncr(i,j,1))

                end if

              else

                frrg(i,j,1)=0.e0
                frrgn(i,j,1)=0.e0

              end if

            else

              frrg(i,j,1)=0.e0
              frrgn(i,j,1)=0.e0

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

!$omp do schedule(runtime) private(i,j,diaqr3)

            do j=1,nj-1
            do i=1,ni-1

              if(qr(i,j,k).gt.thresq) then

                if(tcel(i,j,k).lt.t0cel) then

                  if(tcel(i,j,k).le.tclow) then

                    frrg(i,j,k)=qr(i,j,k)

                  else

                    diaqr3=diaqr(i,j,k)*diaqr(i,j,k)*diaqr(i,j,k)

                    frrg(i,j,k)=min(qr(i,j,k),cfrrg*diaqr3*diaqr3       &
     &                *(exp(-.66e0*tcel(i,j,k))-1.e0)*ncr(i,j,k))

                  end if

                else

                  frrg(i,j,k)=0.e0

                end if

              else

                frrg(i,j,k)=0.e0

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

!$omp do schedule(runtime) private(i,j,diaqr3,a)

            do j=1,nj-1
            do i=1,ni-1

              if(qr(i,j,k).gt.thresq) then

                if(tcel(i,j,k).lt.t0cel) then

                  if(tcel(i,j,k).le.tclow) then

                    frrg(i,j,k)=qr(i,j,k)
                    frrgn(i,j,k)=ncr(i,j,k)

                  else

                    a=(exp(-.66e0*tcel(i,j,k))-1.e0)*ncr(i,j,k)

                    diaqr3=diaqr(i,j,k)*diaqr(i,j,k)*diaqr(i,j,k)

                    frrg(i,j,k)=min(cfrrg*diaqr3*diaqr3*a,qr(i,j,k))
                    frrgn(i,j,k)=min(cfrrgn*diaqr3*a,ncr(i,j,k))

                  end if

                else

                  frrg(i,j,k)=0.e0
                  frrgn(i,j,k)=0.e0

                end if

              else

                frrg(i,j,k)=0.e0
                frrgn(i,j,k)=0.e0

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

#endif

! Dump output data at target call
if (dump_call_count_freezing == DUMP_TARGET_freezing .and. .not. dump_done_freezing) then
  call dump_array_3d('frrg_ref.bin', frrg, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('frrgn_ref.bin', frrgn, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_freezing = .true.
end if


call profile_stop(prof_id1, loop_len)

!!! -----

      end subroutine s_freezing

!-----7--------------------------------------------------------------7--

      end module m_freezing

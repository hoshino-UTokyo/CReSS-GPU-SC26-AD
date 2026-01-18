!***********************************************************************
      module m_siadjst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/05/31
!     Modification: 2004/06/10, 2004/09/01, 2004/09/10, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/01/07, 2005/01/31,
!                   2005/04/04, 2005/10/05, 2006/02/13, 2006/09/30,
!                   2007/10/19, 2007/11/26, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the saturation adjustment for ice.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_dump_kernel
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: siadjst, s_siadjst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface siadjst

        module procedure s_siadjst

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_siadjst(fpthresq,ni,nj,nk,ptbr,pi,p,ptp,qv,qi,nci)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpthresq
                       ! Formal parameter of unique index of thresq

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

! Input and output variables

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(inout) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(inout) :: nci(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice

! Internal shared variables

      real thresq      ! Minimum threshold value of mixing ratio

      real cwmci       ! cw - ci

      real mi0iv       ! Inverse of mi0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real t           ! Air temperature

      real tcel        ! Ambient air temperature

      real esi         ! Saturation vapor pressure for ice
      real qvsi        ! Saturation mixing ratio for ice

      real lscpi       ! Latent heat of sublimation / (cp x pi)

      real dqi         ! Variation of cloud ice mixing ratio

      real a           ! Temporary variable
      real b           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_siadjst = 0
      integer, parameter :: DUMP_TARGET_siadjst = 720
      logical, save :: dump_done_siadjst = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpthresq,thresq)

! -----

! Set the common used variables.

      cwmci=cw-ci

      mi0iv=1.e0/mi0

! -----

! Perform the saturation adjustment.

!@llm start meta_info ----------------------------------------------------
! Location: siadjst.f90 :: s_siadjst
! Summary : Performs saturation adjustment for ice, converting between water vapor
!           and cloud ice based on saturation conditions at low temperatures.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses intrinsic functions: exp, log
!   - Uses module constants from m_comphy (t0, tlow, es0, epsva, lv0, lf0, etc.)
!   - Complex thermodynamic calculations with multiple conditional branches
!   - Two-iteration adjustment loop structure within each grid point
!   - Updates ptp, qv, qi, nci arrays (multiple output variables)
!   - All loops independent with private i,j,k and local scalar variables
!   - No synchronization constructs
! Next:
!   - Port exp/log intrinsics directly (GPU-compatible)
!   - May need to handle thread divergence from nested conditionals
!   - Consider data regions for the 4 updated 3D arrays
! Runtime:
!   - Calls: 720
!   - AvgLoops: 102.4M
!   - TotalTime: 4.520s (0.15%)
!   - AvgTime: 6.277ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('siadjst.f90', 's_siadjst', &
   & 'OMP section 1')
end if
loop_len = int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_siadjst = dump_call_count_siadjst + 1
if (dump_call_count_siadjst == DUMP_TARGET_siadjst .and. .not. dump_done_siadjst) then
  call dump_init('siadjst')
  call dump_scalar_r('thresq', thresq)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_array_3d('p.bin', p, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('cp', cp)
  call dump_scalar_r('t0', t0)
  call dump_scalar_r('epsva', epsva)
  call dump_scalar_r('es0', es0)
  call dump_scalar_r('lv0', lv0)
  call dump_scalar_r('lf0', lf0)
  call dump_array_3d('ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pi.bin', pi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptp_in.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qv_in.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qi_in.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nci_in.bin', nci, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('cwmci', cwmci)
  call dump_scalar_r('mi0iv', mi0iv)
end if

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,t,tcel,esi,qvsi,lscpi,dqi,a,b)

        do j=1,nj-1
        do i=1,ni-1
          t=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

          if(t.le.tlow) then

            tcel=t-t0

            a=1.e0/(t-7.66e0)
            b=a*tcel

            esi=es0*exp(21.875e0*b)

            qvsi=epsva*esi/(p(i,j,k)-esi)

            if(qi(i,j,k).gt.thresq.or.qv(i,j,k).gt.qvsi) then

              lscpi=(lv0*exp((.167e0+3.67e-4*t)*log(t0/t))              &
     &          +(lf0+cwmci*tcel))/(cp*pi(i,j,k))

              dqi=(qvsi-qv(i,j,k))                                      &
     &          /(1.e0+21.875e0*a*(1.e0-b)*qvsi*lscpi*pi(i,j,k))

              if(qi(i,j,k).gt.dqi) then

                if(qi(i,j,k).gt.thresq) then

                  nci(i,j,k)=nci(i,j,k)-dqi*nci(i,j,k)/qi(i,j,k)

                else

                  nci(i,j,k)=nci(i,j,k)-dqi*mi0iv

                end if

                ptp(i,j,k)=ptp(i,j,k)-dqi*lscpi

                qv(i,j,k)=qv(i,j,k)+dqi
                qi(i,j,k)=qi(i,j,k)-dqi

              else

                nci(i,j,k)=0.e0

                ptp(i,j,k)=ptp(i,j,k)-qi(i,j,k)*lscpi

                qv(i,j,k)=qv(i,j,k)+qi(i,j,k)
                qi(i,j,k)=0.e0

              end if

            end if

            t=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)

            if(t.le.tlow) then

              tcel=t-t0

              a=1.e0/(t-7.66e0)
              b=a*tcel

              esi=es0*exp(21.875e0*b)

              qvsi=epsva*esi/(p(i,j,k)-esi)

              if(qi(i,j,k).gt.thresq.or.qv(i,j,k).gt.qvsi) then

                lscpi=(lv0*exp((.167e0+3.67e-4*t)*log(t0/t))            &
     &            +(lf0+cwmci*tcel))/(cp*pi(i,j,k))

                dqi=(qvsi-qv(i,j,k))                                    &
     &            /(1.e0+21.875e0*a*(1.e0-b)*qvsi*lscpi*pi(i,j,k))

                if(qi(i,j,k).gt.dqi) then

                  if(qi(i,j,k).gt.thresq) then

                    nci(i,j,k)=nci(i,j,k)-dqi*nci(i,j,k)/qi(i,j,k)

                  else

                    nci(i,j,k)=nci(i,j,k)-dqi*mi0iv

                  end if

                  ptp(i,j,k)=ptp(i,j,k)-dqi*lscpi

                  qv(i,j,k)=qv(i,j,k)+dqi
                  qi(i,j,k)=qi(i,j,k)-dqi

                else

                  nci(i,j,k)=0.e0

                  ptp(i,j,k)=ptp(i,j,k)-qi(i,j,k)*lscpi

                  qv(i,j,k)=qv(i,j,k)+qi(i,j,k)
                  qi(i,j,k)=0.e0

                end if

              end if

            end if

          end if

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_siadjst == DUMP_TARGET_siadjst .and. .not. dump_done_siadjst) then
  call dump_array_3d('ptp_ref.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qv_ref.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qi_ref.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('nci_ref.bin', nci, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_siadjst = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_siadjst

!-----7--------------------------------------------------------------7--

      end module m_siadjst

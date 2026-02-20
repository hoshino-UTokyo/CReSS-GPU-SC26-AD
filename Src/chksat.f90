!***********************************************************************
      module m_chksat
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2008/07/01
!     Modification: 2008/08/25, 2009/02/27, 2009/11/13, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check and avoid the super saturation mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile
      use m_dump_kernel
      use m_getindx

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chksat, s_chksat

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chksat

        module procedure s_chksat

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chksat(fproc,xo,imin,imax,jmin,jmax,kmin,kmax,       &
     &                    pbr,ptbr,pp,ptp,qv)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fproc
                       ! Control flag of processing type

      character(len=3), intent(in) :: xo
                       ! Control flag of variable arrangement

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

      real, intent(in) :: pbr(imin:imax,jmin:jmax,kmin:kmax)
                       ! Base state pressure

      real, intent(in) :: ptbr(imin:imax,jmin:jmax,kmin:kmax)
                       ! Base state potential temperature

      real, intent(in) :: pp(imin:imax,jmin:jmax,kmin:kmax)
                       ! Pressure perturbation

      real, intent(in) :: ptp(imin:imax,jmin:jmax,kmin:kmax)
                       ! Potential temperature perturbation

! Input and output variable

      real, intent(inout) :: qv(imin:imax,jmin:jmax,kmin:kmax)
                       ! Water vapor mixing ratio

! Internal shared variables

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction
      integer kstr     ! Minimum do loops index in z direction
      integer kend     ! Maximum do loops index in z direction

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real t           ! Temperature

      real es          ! Saturation vapor pressure

      real qvs         ! Saturation mixing ratio

      real pres        ! Full pressure

      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_chksat = 0
      integer, parameter :: DUMP_TARGET_chksat = 3
      logical, save :: dump_done_chksat = .false.


!-----7--------------------------------------------------------------7--

! Get the maximum and minimim indices of do loops.

      call getindx(xo,imin,imax,jmin,jmax,istr,iend,jstr,jend)

      if(xo(3:3).eq.'o') then
        kstr=1
        kend=kmax
      else if(xo(3:3).eq.'x') then
        kstr=2
        kend=kmax-2
      end if

! -----

! Set the common used variables.

      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Check and avoid the super saturation mixing ratio.

!@llm start meta_info ----------------------------------------------------
! Location: chksat.f90 :: s_chksat
! Summary : Limit water vapor mixing ratio (qv) to saturation value computed
!           from pressure and temperature fields
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function/subroutine calls inside parallel region
!   - No reductions or synchronization constructs
!   - Uses intrinsic functions (exp, log, min)
!   - Modifies output array qv in-place
!   - Conditional branches based on fproc flag and temperature threshold (tlow)
!   - Saturation vapor pressure computed using Clausius-Clapeyron approximation
! Next:
!   - Direct OpenACC with collapse(2) for inner loops
!   - exp/log functions have GPU intrinsic support
! Runtime:
!   - Calls: 3
!   - AvgLoops: 101.2M
!   - TotalTime: 0.040s (0.00%)
!   - AvgTime: 13.411ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('chksat.f90', 's_chksat', &
   & 'OMP section 1')
end if
loop_len = int((kend)-(kstr)+1,8) &
     & * int((jend)-(jstr)+1,8) &
     & * int((iend)-(istr)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_chksat = dump_call_count_chksat + 1
if (dump_call_count_chksat == DUMP_TARGET_chksat .and. .not. dump_done_chksat) then
  call dump_init('chksat')
  call dump_scalar_i('imin', imin)
  call dump_scalar_i('imax', imax)
  call dump_scalar_i('jmin', jmin)
  call dump_scalar_i('jmax', jmax)
  call dump_scalar_i('kmin', kmin)
  call dump_scalar_i('kmax', kmax)
  call dump_scalar_r('t0', t0)
  call dump_scalar_r('epsva', epsva)
  call dump_scalar_r('es0', es0)
  call dump_array_3d('pbr.bin', pbr, imin, imax, jmin, jmax, kmin, kmax)
  call dump_array_3d('ptbr.bin', ptbr, imin, imax, jmin, jmax, kmin, kmax)
  call dump_array_3d('pp.bin', pp, imin, imax, jmin, jmax, kmin, kmax)
  call dump_array_3d('ptp.bin', ptp, imin, imax, jmin, jmax, kmin, kmax)
  call dump_array_3d('qv_in.bin', qv, imin, imax, jmin, jmax, kmin, kmax)
  call dump_scalar_c('fproc', fproc)
  call dump_scalar_i('iend', iend)
  call dump_scalar_i('istr', istr)
  call dump_scalar_i('jend', jend)
  call dump_scalar_i('jstr', jstr)
  call dump_scalar_i('kend', kend)
  call dump_scalar_i('kstr', kstr)
  call dump_scalar_r('p0iv', p0iv)
  call dump_scalar_r('rddvcp', rddvcp)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_055)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------

      if(fproc(1:3).eq.'bar') then

    !$acc kernels
    !$acc loop independent
        do k=kstr,kend
          !$acc loop independent
          do j=jstr,jend
            !$acc loop independent
            do i=istr,iend

              t=ptbr(i,j,k)*exp(rddvcp*log(p0iv*pbr(i,j,k)))

              if(t.gt.tlow) then
                es=es0*exp(17.269e0*(t-t0)/(t-35.86e0))
                qvs=epsva*es/(pbr(i,j,k)-es)
              else
                es=es0*exp(21.875e0*(t-t0)/(t-7.66e0))
                qvs=epsva*es/(pbr(i,j,k)-es)
              end if

              qv(i,j,k)=min(qv(i,j,k),qvs)

            end do
          end do
        end do
    !$acc end kernels

      else if(fproc(1:5).eq.'total') then

    !$acc kernels
    !$acc loop independent
        do k=kstr,kend
          !$acc loop independent
          do j=jstr,jend
            !$acc loop independent
            do i=istr,iend

              pres=pbr(i,j,k)+pp(i,j,k)
              t=(ptbr(i,j,k)+ptp(i,j,k))                                &
     &          *exp(rddvcp*log(p0iv*pres))

              if(t.gt.tlow) then
                es=es0*exp(17.269e0*(t-t0)/(t-35.86e0))
                qvs=epsva*es/(pres-es)
              else
                es=es0*exp(21.875e0*(t-t0)/(t-7.66e0))
                qvs=epsva*es/(pres-es)
              end if

              qv(i,j,k)=min(qv(i,j,k),qvs)

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

      if(fproc(1:3).eq.'bar') then

        do k=kstr,kend

!$omp do schedule(runtime) private(i,j,t,es,qvs)

          do j=jstr,jend
          do i=istr,iend

            t=ptbr(i,j,k)*exp(rddvcp*log(p0iv*pbr(i,j,k)))

            if(t.gt.tlow) then

              es=es0*exp(17.269e0*(t-t0)/(t-35.86e0))

              qvs=epsva*es/(pbr(i,j,k)-es)

            else

              es=es0*exp(21.875e0*(t-t0)/(t-7.66e0))

              qvs=epsva*es/(pbr(i,j,k)-es)

            end if

            qv(i,j,k)=min(qv(i,j,k),qvs)

          end do
          end do

!$omp end do

        end do

      else if(fproc(1:5).eq.'total') then

        do k=kstr,kend

!$omp do schedule(runtime) private(i,j,t,es,qvs)

          do j=jstr,jend
          do i=istr,iend

            t=(ptbr(i,j,k)+ptp(i,j,k))                                  &
     &        *exp(rddvcp*log(p0iv*(pbr(i,j,k)+pp(i,j,k))))

            if(t.gt.tlow) then

              es=es0*exp(17.269e0*(t-t0)/(t-35.86e0))

              qvs=epsva*es/((pbr(i,j,k)+pp(i,j,k))-es)

            else

              es=es0*exp(21.875e0*(t-t0)/(t-7.66e0))

              qvs=epsva*es/((pbr(i,j,k)+pp(i,j,k))-es)

            end if

            qv(i,j,k)=min(qv(i,j,k),qvs)

          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_chksat == DUMP_TARGET_chksat .and. .not. dump_done_chksat) then
  call dump_array_3d('qv_ref.bin', qv, imin, imax, jmin, jmax, kmin, kmax)
  call dump_finalize()
  dump_done_chksat = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_chksat

!-----7--------------------------------------------------------------7--

      end module m_chksat

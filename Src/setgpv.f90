!***********************************************************************
      module m_setgpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/12/12
!     Modification: 2006/01/10, 2006/09/21, 2007/09/25, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2008/12/11, 2009/02/27,
!                   2009/03/23, 2011/08/18, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the interpolated GPV variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_comprofile
      use m_dump_kernel
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setgpv, s_setgpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setgpv

        module procedure s_setgpv

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
      subroutine s_setgpv(fpgpvvar,fpcphopt,fphaiopt,fpgpvitv,ird,      &
     &                    ni,nj,nk,nqw,nqi,ugpv,utd,vgpv,vtd,wgpv,wtd,  &
     &                    ppgpv,pptd,ptpgpv,ptptd,qvgpv,qvtd,           &
     &                    qwgpv,qwtd,qigpv,qitd)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpgpvitv
                       ! Formal parameter of unique index of gpvitv

      integer, intent(in) :: ird
                       ! Index of count to read out in rdgpvnxt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

! Input and output variables

      real, intent(inout) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(inout) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(inout) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

      real, intent(inout) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: pptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

      real, intent(inout) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

      real, intent(inout) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(inout) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

      real, intent(inout) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(inout) :: qwtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of water hydrometeor of GPV data

      real, intent(inout) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

      real, intent(inout) :: qitd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of ice hydrometeor of GPV data

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

      real gpvitv      ! Time interval of GPV data

      real gpviv       ! Inverse of gpvitv

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_setgpv = 0
      integer, parameter :: DUMP_TARGET_setgpv = 2
      logical, save :: dump_done_setgpv = .false.


!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(gpvvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getrname(fpgpvitv,gpvitv)

! -----

! Set the common used variable.

      gpviv=1.e0/gpvitv

! -----

!! Set the interpolated GPV variables.

!@llm start meta_info ----------------------------------------------------
! Location: setgpv.f90 :: s_setgpv
! Summary : Sets interpolated GPV (Grid Point Value) variables including velocity, pressure,
!           temperature, and hydrometeor time tendencies or values depending on read index.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Multiple conditional branches based on gpvvar flags and cphopt/haiopt options
!   - All loops are independent with private i,j,k indices
!   - Multiple separate do-omp do blocks for different variable categories
!   - No synchronization constructs beyond implicit barriers at omp end do
! Next:
!   - Consider collapsing nested conditionals into unified kernels
!   - Use OpenACC/OpenACC with data regions for array transfers
!   - May benefit from kernel fusion for related variable updates
! Runtime:
!   - Calls: 2
!   - AvgLoops: 103.4M
!   - TotalTime: 0.035s (0.00%)
!   - AvgTime: 17.427ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('setgpv.f90', 's_setgpv', &
   & 'OMP section 1')
end if
loop_len = int((nk)-(1)+1,8) * int((nj)-(1)+1,8) * int((ni)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_setgpv = dump_call_count_setgpv + 1
if (dump_call_count_setgpv == DUMP_TARGET_setgpv .and. .not. dump_done_setgpv) then
  call dump_init('setgpv')
  call dump_scalar_c('gpvvar', gpvvar)
  call dump_scalar_i('cphopt', cphopt)
  call dump_scalar_i('haiopt', haiopt)
  call dump_scalar_r('gpvitv', gpvitv)
  call dump_scalar_i('ird', ird)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_i('nqw', nqw)
  call dump_scalar_i('nqi', nqi)
  call dump_array_3d('ugpv_in.bin', ugpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('utd_in.bin', utd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vgpv_in.bin', vgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vtd_in.bin', vtd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wgpv_in.bin', wgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wtd_in.bin', wtd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ppgpv_in.bin', ppgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pptd_in.bin', pptd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptpgpv_in.bin', ptpgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptptd_in.bin', ptptd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvgpv_in.bin', qvgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvtd_in.bin', qvtd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_4d('qwgpv_in.bin', qwgpv, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
  call dump_array_4d('qwtd_in.bin', qwtd, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
  call dump_array_4d('qigpv_in.bin', qigpv, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)
  call dump_array_4d('qitd_in.bin', qitd, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)
  call dump_scalar_r('gpviv', gpviv)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_288)
! GPU version (OpenACC)

! Set the time tendency of variables at current marked time.

      if(ird.eq.1) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk
          !$acc loop independent
          do j=1,nj
            !$acc loop independent
            do i=1,ni
              utd(i,j,k)=(utd(i,j,k)-ugpv(i,j,k))*gpviv
              vtd(i,j,k)=(vtd(i,j,k)-vgpv(i,j,k))*gpviv
              pptd(i,j,k)=(pptd(i,j,k)-ppgpv(i,j,k))*gpviv
              ptptd(i,j,k)=(ptptd(i,j,k)-ptpgpv(i,j,k))*gpviv
            end do
          end do
        end do
        !$acc end kernels

        if(gpvvar(1:1).eq.'o') then
          !$acc kernels
          !$acc loop independent
          do k=1,nk
            !$acc loop independent
            do j=1,nj
              !$acc loop independent
              do i=1,ni
                wtd(i,j,k)=(wtd(i,j,k)-wgpv(i,j,k))*gpviv
              end do
            end do
          end do
          !$acc end kernels
        end if

        if(gpvvar(2:2).eq.'o') then
          !$acc kernels
          !$acc loop independent
          do k=1,nk
            !$acc loop independent
            do j=1,nj
              !$acc loop independent
              do i=1,ni
                qvtd(i,j,k)=(qvtd(i,j,k)-qvgpv(i,j,k))*gpviv
              end do
            end do
          end do
          !$acc end kernels
        end if

        if(abs(cphopt).lt.10) then
          if(abs(cphopt).ge.1) then
            if(gpvvar(3:3).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qwtd(i,j,k,1)=(qwtd(i,j,k,1)-qwgpv(i,j,k,1))*gpviv
                  end do
                end do
              end do
              !$acc end kernels
            end if
            if(gpvvar(4:4).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qwtd(i,j,k,2)=(qwtd(i,j,k,2)-qwgpv(i,j,k,2))*gpviv
                  end do
                end do
              end do
              !$acc end kernels
            end if
          end if
          if(abs(cphopt).ge.2) then
            if(gpvvar(5:5).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qitd(i,j,k,1)=(qitd(i,j,k,1)-qigpv(i,j,k,1))*gpviv
                  end do
                end do
              end do
              !$acc end kernels
            end if
            if(gpvvar(6:6).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qitd(i,j,k,2)=(qitd(i,j,k,2)-qigpv(i,j,k,2))*gpviv
                  end do
                end do
              end do
              !$acc end kernels
            end if
            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qitd(i,j,k,3)=(qitd(i,j,k,3)-qigpv(i,j,k,3))*gpviv
                  end do
                end do
              end do
              !$acc end kernels
            end if
            if(gpvvar(8:8).eq.'o') then
              if(haiopt.eq.1) then
                !$acc kernels
                !$acc loop independent
                do k=1,nk
                  !$acc loop independent
                  do j=1,nj
                    !$acc loop independent
                    do i=1,ni
                      qitd(i,j,k,4)=(qitd(i,j,k,4)-qigpv(i,j,k,4))*gpviv
                    end do
                  end do
                end do
                !$acc end kernels
              end if
            end if
          end if
        end if

      end if

! Set the variables at current marked time.

      if(ird.eq.2) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk
          !$acc loop independent
          do j=1,nj
            !$acc loop independent
            do i=1,ni
              ugpv(i,j,k)=utd(i,j,k)
              vgpv(i,j,k)=vtd(i,j,k)
              ppgpv(i,j,k)=pptd(i,j,k)
              ptpgpv(i,j,k)=ptptd(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

        if(gpvvar(1:1).eq.'o') then
          !$acc kernels
          !$acc loop independent
          do k=1,nk
            !$acc loop independent
            do j=1,nj
              !$acc loop independent
              do i=1,ni
                wgpv(i,j,k)=wtd(i,j,k)
              end do
            end do
          end do
          !$acc end kernels
        end if

        if(gpvvar(2:2).eq.'o') then
          !$acc kernels
          !$acc loop independent
          do k=1,nk
            !$acc loop independent
            do j=1,nj
              !$acc loop independent
              do i=1,ni
                qvgpv(i,j,k)=qvtd(i,j,k)
              end do
            end do
          end do
          !$acc end kernels
        end if

        if(abs(cphopt).lt.10) then
          if(abs(cphopt).ge.1) then
            if(gpvvar(3:3).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qwgpv(i,j,k,1)=qwtd(i,j,k,1)
                  end do
                end do
              end do
              !$acc end kernels
            end if
            if(gpvvar(4:4).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qwgpv(i,j,k,2)=qwtd(i,j,k,2)
                  end do
                end do
              end do
              !$acc end kernels
            end if
          end if
          if(abs(cphopt).ge.2) then
            if(gpvvar(5:5).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qigpv(i,j,k,1)=qitd(i,j,k,1)
                  end do
                end do
              end do
              !$acc end kernels
            end if
            if(gpvvar(6:6).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qigpv(i,j,k,2)=qitd(i,j,k,2)
                  end do
                end do
              end do
              !$acc end kernels
            end if
            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then
              !$acc kernels
              !$acc loop independent
              do k=1,nk
                !$acc loop independent
                do j=1,nj
                  !$acc loop independent
                  do i=1,ni
                    qigpv(i,j,k,3)=qitd(i,j,k,3)
                  end do
                end do
              end do
              !$acc end kernels
            end if
            if(gpvvar(8:8).eq.'o') then
              if(haiopt.eq.1) then
                !$acc kernels
                !$acc loop independent
                do k=1,nk
                  !$acc loop independent
                  do j=1,nj
                    !$acc loop independent
                    do i=1,ni
                      qigpv(i,j,k,4)=qitd(i,j,k,4)
                    end do
                  end do
                end do
                !$acc end kernels
              end if
            end if
          end if
        end if

      end if

#else
! CPU version (OpenMP) - Original code preserved

!$omp parallel default(shared) private(k)

! Set the time tendency of variables at current marked time.

      if(ird.eq.1) then

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=1,nj
          do i=1,ni
            utd(i,j,k)=(utd(i,j,k)-ugpv(i,j,k))*gpviv
            vtd(i,j,k)=(vtd(i,j,k)-vgpv(i,j,k))*gpviv

            pptd(i,j,k)=(pptd(i,j,k)-ppgpv(i,j,k))*gpviv
            ptptd(i,j,k)=(ptptd(i,j,k)-ptpgpv(i,j,k))*gpviv

          end do
          end do

!$omp end do

        end do

        if(gpvvar(1:1).eq.'o') then

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni
              wtd(i,j,k)=(wtd(i,j,k)-wgpv(i,j,k))*gpviv
            end do
            end do

!$omp end do

          end do

        end if

        if(gpvvar(2:2).eq.'o') then

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni
              qvtd(i,j,k)=(qvtd(i,j,k)-qvgpv(i,j,k))*gpviv
            end do
            end do

!$omp end do

          end do

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qwtd(i,j,k,1)=(qwtd(i,j,k,1)-qwgpv(i,j,k,1))*gpviv
                end do
                end do

!$omp end do

              end do

            end if

            if(gpvvar(4:4).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qwtd(i,j,k,2)=(qwtd(i,j,k,2)-qwgpv(i,j,k,2))*gpviv
                end do
                end do

!$omp end do

              end do

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qitd(i,j,k,1)=(qitd(i,j,k,1)-qigpv(i,j,k,1))*gpviv
                end do
                end do

!$omp end do

              end do

            end if

            if(gpvvar(6:6).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qitd(i,j,k,2)=(qitd(i,j,k,2)-qigpv(i,j,k,2))*gpviv
                end do
                end do

!$omp end do

              end do

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qitd(i,j,k,3)=(qitd(i,j,k,3)-qigpv(i,j,k,3))*gpviv
                end do
                end do

!$omp end do

              end do

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                do k=1,nk

!$omp do schedule(runtime) private(i,j)

                  do j=1,nj
                  do i=1,ni
                    qitd(i,j,k,4)=(qitd(i,j,k,4)-qigpv(i,j,k,4))*gpviv
                  end do
                  end do

!$omp end do

                end do

              end if

            end if

          end if

        end if

      end if

! -----

! Set the variables at current marked time.

      if(ird.eq.2) then

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=1,nj
          do i=1,ni
            ugpv(i,j,k)=utd(i,j,k)
            vgpv(i,j,k)=vtd(i,j,k)

            ppgpv(i,j,k)=pptd(i,j,k)
            ptpgpv(i,j,k)=ptptd(i,j,k)

          end do
          end do

!$omp end do

        end do

        if(gpvvar(1:1).eq.'o') then

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni
              wgpv(i,j,k)=wtd(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

        if(gpvvar(2:2).eq.'o') then

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni
              qvgpv(i,j,k)=qvtd(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qwgpv(i,j,k,1)=qwtd(i,j,k,1)
                end do
                end do

!$omp end do

              end do

            end if

            if(gpvvar(4:4).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qwgpv(i,j,k,2)=qwtd(i,j,k,2)
                end do
                end do

!$omp end do

              end do

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qigpv(i,j,k,1)=qitd(i,j,k,1)
                end do
                end do

!$omp end do

              end do

            end if

            if(gpvvar(6:6).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qigpv(i,j,k,2)=qitd(i,j,k,2)
                end do
                end do

!$omp end do

              end do

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni
                  qigpv(i,j,k,3)=qitd(i,j,k,3)
                end do
                end do

!$omp end do

              end do

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                do k=1,nk

!$omp do schedule(runtime) private(i,j)

                  do j=1,nj
                  do i=1,ni
                    qigpv(i,j,k,4)=qitd(i,j,k,4)
                  end do
                  end do

!$omp end do

                end do

              end if

            end if

          end if

        end if

      end if

! -----

!$omp end parallel
#endif

! Dump output data at target call
if (dump_call_count_setgpv == DUMP_TARGET_setgpv .and. .not. dump_done_setgpv) then
  call dump_array_3d('ugpv_ref.bin', ugpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('utd_ref.bin', utd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vgpv_ref.bin', vgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('vtd_ref.bin', vtd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wgpv_ref.bin', wgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('wtd_ref.bin', wtd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ppgpv_ref.bin', ppgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('pptd_ref.bin', pptd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptpgpv_ref.bin', ptpgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('ptptd_ref.bin', ptptd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvgpv_ref.bin', qvgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('qvtd_ref.bin', qvtd, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_4d('qwgpv_ref.bin', qwgpv, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
  call dump_array_4d('qwtd_ref.bin', qwtd, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
  call dump_array_4d('qigpv_ref.bin', qigpv, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)
  call dump_array_4d('qitd_ref.bin', qitd, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)
  call dump_finalize()
  dump_done_setgpv = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_setgpv

!-----7--------------------------------------------------------------7--

      end module m_setgpv

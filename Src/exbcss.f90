!***********************************************************************
      module m_exbcss
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/07/05, 1999/08/03, 1999/08/09, 1999/09/30,
!                   1999/11/01, 2000/01/17, 2000/02/02, 2000/04/18,
!                   2001/01/15, 2001/03/13, 2001/04/15, 2001/05/29,
!                   2001/06/06, 2001/06/29, 2001/07/13, 2001/08/07,
!                   2001/11/20, 2001/12/11, 2002/04/02, 2002/06/06,
!                   2002/07/23, 2002/08/15, 2002/10/31, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/11/05,
!                   2003/11/28, 2003/12/12, 2004/04/15, 2004/08/20,
!                   2005/01/31, 2005/02/10, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/31, 2007/05/07, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2008/12/11, 2009/02/27,
!                   2009/03/23, 2011/09/22, 2013/01/28, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the lateral boundary value to the external boundary value
!     for scalar variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comprofile
      use m_dump_kernel
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: exbcss, s_exbcss

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface exbcss

        module procedure s_exbcss

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_exbcss(fpexbvar,fpwbc,fpebc,fpexnews,ape,            &
     &                    isstp,dts,gtinc,ni,nj,nk,scpx,scpy,sgpv,std,s)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpexnews
                       ! Formal parameter of unique index of exnews

      integer, intent(in) :: ape
                       ! Pointer of exbvar

      integer, intent(in) :: isstp
                       ! Index of small time steps integration

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: scpx(1:nj,1:nk,1:2)
                       ! Phase speed of optional scalar variable
                       ! on west and east boundary

      real, intent(in) :: scpy(1:ni,1:nk,1:2)
                       ! Phase speed of optional scalar variable
                       ! on south and north boundary

      real, intent(in) :: sgpv(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable of GPV data
                       ! at marked time

      real, intent(in) :: std(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! optional scalar variable of GPV data

! Input and output variable

      real, intent(inout) :: s(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer njm1     ! nj - 1
      integer njm2     ! nj - 2

      real exnews      ! Boundary damping coefficient

      real dmpdt       ! exnews x dts

      real tpdt        ! gtinc + real(isstp - 1) x dts

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real sb1         ! Temporary variable
      real sb2         ! Temporary variable
      real sb2i        ! Temporary variable
      real sb2j        ! Temporary variable

      real radwe       ! Temporary variable
      real radsn       ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_exbcss = 0
      integer, parameter :: DUMP_TARGET_exbcss = 14400
      logical, save :: dump_done_exbcss = .false.


!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getrname(fpexnews,exnews)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      njm1=nj-1
      njm2=nj-2

      dmpdt=exnews*dts

      tpdt=gtinc+real(isstp-1)*dts

! -----

!! Force the lateral boundary value to the external boundary value.

!@llm start meta_info ----------------------------------------------------
! Location: exbcss.f90 :: s_exbcss
! Summary : Force lateral boundary values to external GPV boundary values
!           for scalar variables using radiation boundary conditions
! GPU diff: Hard
! Findings:
!   - No omp_get_thread_num usage
!   - Uses MPI domain decomposition variables (ebw, ebe, ebs, ebn, isub, jsub)
!   - Many conditional branches based on boundary location and options
!   - Processes corners, west, east, south, north boundaries separately
!   - Updates s array at domain boundaries only
!   - Small time step integration (dts) for acoustic mode
!   - exbvar character flags control which boundaries are active
! Next:
!   - Boundary-only operations may not benefit much from GPU
!   - Consider keeping boundary conditions on CPU if main computation on GPU
!   - If porting, need separate small kernels for each boundary section
!   - MPI communication patterns need careful handling with GPU buffers
! Runtime:
!   - Calls: 14400
!   - AvgLoops: 125
!   - TotalTime: 3.244s (0.11%)
!   - AvgTime: 0.225ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('exbcss.f90', 's_exbcss', &
   & 'OMP section 1')
end if
loop_len = int((nk-2)-(2)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_exbcss = dump_call_count_exbcss + 1
if (dump_call_count_exbcss == DUMP_TARGET_exbcss .and. .not. dump_done_exbcss) then
  call dump_init('exbcss')
  call dump_scalar_c('exbvar', exbvar)
  call dump_scalar_i('wbc', wbc)
  call dump_scalar_i('ebc', ebc)
  call dump_scalar_r('exnews', exnews)
  call dump_scalar_i('ape', ape)
  call dump_scalar_i('isstp', isstp)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_i('nk', nk)
  call dump_scalar_r('dts', dts)
  call dump_scalar_r('gtinc', gtinc)
  call dump_array_3d('s_in.bin', s, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('scpx.bin', scpx, 1, nj, 1, nk, 1, 2)
  call dump_array_3d('scpy.bin', scpy, 1, ni, 1, nk, 1, 2)
  call dump_array_3d('sgpv.bin', sgpv, 0, ni+1, 0, nj+1, 1, nk)
  call dump_array_3d('std.bin', std, 0, ni+1, 0, nj+1, 1, nk)
  call dump_scalar_r('dmpdt', dmpdt)
  call dump_scalar_i('ebe', ebe)
  call dump_scalar_i('ebn', ebn)
  call dump_scalar_i('ebs', ebs)
  call dump_scalar_i('ebw', ebw)
  call dump_scalar_i('isub', isub)
  call dump_scalar_i('jsub', jsub)
  call dump_scalar_i('nim1', nim1)
  call dump_scalar_i('nim2', nim2)
  call dump_scalar_i('nisub', nisub)
  call dump_scalar_i('njm1', njm1)
  call dump_scalar_i('njm2', njm2)
  call dump_scalar_i('njsub', njsub)
  call dump_scalar_r('tpdt', tpdt)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_104)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------

    ! Force boundary at four corners
    if (abs(wbc) /= 1 .and. abs(ebc) /= 1) then
      if (exbvar(ape:ape) == '-') then
        if (ebs == 1 .and. jsub == 0) then
          if (ebw == 1 .and. isub == 0) then
            !$acc kernels
            !$acc loop independent private(sb1,sb2i,sb2j,radwe,radsn)
            do k = 2, nk-2
              sb1 = sgpv(1,1,k) + std(1,1,k)*tpdt
              sb2i = sgpv(2,1,k) + std(2,1,k)*tpdt
              sb2j = sgpv(1,2,k) + std(1,2,k)*tpdt
              radwe = scpx(1,k,1)*((s(2,1,k)-s(1,1,k))-(sb2i-sb1))
              radsn = scpy(1,k,1)*((s(1,2,k)-s(1,1,k))-(sb2j-sb1))
              s(1,1,k) = s(1,1,k) + std(1,1,k)*dts &
                - (radwe+radsn) - dmpdt*(s(1,1,k)-sb1)
            end do
            !$acc end kernels
          end if
          if (ebe == 1 .and. isub == nisub-1) then
            !$acc kernels
            !$acc loop independent private(sb1,sb2i,sb2j,radwe,radsn)
            do k = 2, nk-2
              sb1 = sgpv(nim1,1,k) + std(nim1,1,k)*tpdt
              sb2i = sgpv(nim2,1,k) + std(nim2,1,k)*tpdt
              sb2j = sgpv(nim1,2,k) + std(nim1,2,k)*tpdt
              radwe = scpx(1,k,2)*((s(nim2,1,k)-s(nim1,1,k))-(sb2i-sb1))
              radsn = scpy(nim1,k,1)*((s(nim1,2,k)-s(nim1,1,k))-(sb2j-sb1))
              s(nim1,1,k) = s(nim1,1,k) + std(nim1,1,k)*dts &
                + (radwe-radsn) - dmpdt*(s(nim1,1,k)-sb1)
            end do
            !$acc end kernels
          end if
        end if
        if (ebn == 1 .and. jsub == njsub-1) then
          if (ebw == 1 .and. isub == 0) then
            !$acc kernels
            !$acc loop independent private(sb1,sb2i,sb2j,radwe,radsn)
            do k = 2, nk-2
              sb1 = sgpv(1,njm1,k) + std(1,njm1,k)*tpdt
              sb2i = sgpv(2,njm1,k) + std(2,njm1,k)*tpdt
              sb2j = sgpv(1,njm2,k) + std(1,njm2,k)*tpdt
              radwe = scpx(njm1,k,1)*((s(2,njm1,k)-s(1,njm1,k))-(sb2i-sb1))
              radsn = scpy(1,k,2)*((s(1,njm2,k)-s(1,njm1,k))-(sb2j-sb1))
              s(1,njm1,k) = s(1,njm1,k) + std(1,njm1,k)*dts &
                - (radwe-radsn) - dmpdt*(s(1,njm1,k)-sb1)
            end do
            !$acc end kernels
          end if
          if (ebe == 1 .and. isub == nisub-1) then
            !$acc kernels
            !$acc loop independent private(sb1,sb2i,sb2j,radwe,radsn)
            do k = 2, nk-2
              sb1 = sgpv(nim1,njm1,k) + std(nim1,njm1,k)*tpdt
              sb2i = sgpv(nim2,njm1,k) + std(nim2,njm1,k)*tpdt
              sb2j = sgpv(nim1,njm2,k) + std(nim1,njm2,k)*tpdt
              radwe = scpx(njm1,k,2)*((s(nim2,njm1,k)-s(nim1,njm1,k))-(sb2i-sb1))
              radsn = scpy(nim1,k,2)*((s(nim1,njm2,k)-s(nim1,njm1,k))-(sb2j-sb1))
              s(nim1,njm1,k) = s(nim1,njm1,k) + std(nim1,njm1,k)*dts &
                + (radwe+radsn) - dmpdt*(s(nim1,njm1,k)-sb1)
            end do
            !$acc end kernels
          end if
        end if
      end if
    end if

    ! Force west boundary
    if (ebw == 1 .and. isub == 0) then
      if (abs(wbc) /= 1) then
        if (exbvar(ape:ape) == '-') then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent private(sb1,sb2)
            do j = 2, nj-2
              sb1 = sgpv(1,j,k) + std(1,j,k)*tpdt
              sb2 = sgpv(2,j,k) + std(2,j,k)*tpdt
              s(1,j,k) = s(1,j,k) + std(1,j,k)*dts &
                - scpx(j,k,1)*((s(2,j,k)-s(1,j,k))-(sb2-sb1)) &
                - dmpdt*(s(1,j,k)-sb1)
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 2, nj-2
              s(1,j,k) = s(1,j,k) + std(1,j,k)*dts
            end do
          end do
          !$acc end kernels
        end if
      end if
    end if

    ! Force east boundary
    if (ebe == 1 .and. isub == nisub-1) then
      if (abs(ebc) /= 1) then
        if (exbvar(ape:ape) == '-') then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent private(sb1,sb2)
            do j = 2, nj-2
              sb1 = sgpv(nim1,j,k) + std(nim1,j,k)*tpdt
              sb2 = sgpv(nim2,j,k) + std(nim2,j,k)*tpdt
              s(nim1,j,k) = s(nim1,j,k) + std(nim1,j,k)*dts &
                + scpx(j,k,2)*((s(nim2,j,k)-s(nim1,j,k))-(sb2-sb1)) &
                - dmpdt*(s(nim1,j,k)-sb1)
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 2, nj-2
              s(nim1,j,k) = s(nim1,j,k) + std(nim1,j,k)*dts
            end do
          end do
          !$acc end kernels
        end if
      end if
    end if

    ! Force south boundary
    if (ebs == 1 .and. jsub == 0) then
      if (exbvar(ape:ape) == '-') then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent private(sb1,sb2)
          do i = 2, ni-2
            sb1 = sgpv(i,1,k) + std(i,1,k)*tpdt
            sb2 = sgpv(i,2,k) + std(i,2,k)*tpdt
            s(i,1,k) = s(i,1,k) + std(i,1,k)*dts &
              - scpy(i,k,1)*((s(i,2,k)-s(i,1,k))-(sb2-sb1)) &
              - dmpdt*(s(i,1,k)-sb1)
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            s(i,1,k) = s(i,1,k) + std(i,1,k)*dts
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! Force north boundary
    if (ebn == 1 .and. jsub == njsub-1) then
      if (exbvar(ape:ape) == '-') then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent private(sb1,sb2)
          do i = 2, ni-2
            sb1 = sgpv(i,njm1,k) + std(i,njm1,k)*tpdt
            sb2 = sgpv(i,njm2,k) + std(i,njm2,k)*tpdt
            s(i,njm1,k) = s(i,njm1,k) + std(i,njm1,k)*dts &
              + scpy(i,k,2)*((s(i,njm2,k)-s(i,njm1,k))-(sb2-sb1)) &
              - dmpdt*(s(i,njm1,k)-sb1)
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            s(i,njm1,k) = s(i,njm1,k) + std(i,njm1,k)*dts
          end do
        end do
        !$acc end kernels
      end if
    end if

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared)

! Force the boundary value to the external boundary value at the four
! corners.

      if(abs(wbc).ne.1.and.abs(ebc).ne.1) then

        if(exbvar(ape:ape).eq.'-') then

          if(ebs.eq.1.and.jsub.eq.0) then

            if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,sb1,sb2i,sb2j,radwe,radsn)

              do k=2,nk-2
                sb1=sgpv(1,1,k)+std(1,1,k)*tpdt
                sb2i=sgpv(2,1,k)+std(2,1,k)*tpdt
                sb2j=sgpv(1,2,k)+std(1,2,k)*tpdt

                radwe=scpx(1,k,1)*((s(2,1,k)-s(1,1,k))-(sb2i-sb1))
                radsn=scpy(1,k,1)*((s(1,2,k)-s(1,1,k))-(sb2j-sb1))

                s(1,1,k)=s(1,1,k)+std(1,1,k)*dts                        &
     &            -(radwe+radsn)-dmpdt*(s(1,1,k)-sb1)

              end do

!$omp end do

            end if

            if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,sb1,sb2i,sb2j,radwe,radsn)

              do k=2,nk-2
                sb1=sgpv(nim1,1,k)+std(nim1,1,k)*tpdt
                sb2i=sgpv(nim2,1,k)+std(nim2,1,k)*tpdt
                sb2j=sgpv(nim1,2,k)+std(nim1,2,k)*tpdt

                radwe=scpx(1,k,2)                                       &
     &            *((s(nim2,1,k)-s(nim1,1,k))-(sb2i-sb1))

                radsn=scpy(nim1,k,1)                                    &
     &            *((s(nim1,2,k)-s(nim1,1,k))-(sb2j-sb1))

                s(nim1,1,k)=s(nim1,1,k)+std(nim1,1,k)*dts               &
     &            +(radwe-radsn)-dmpdt*(s(nim1,1,k)-sb1)

              end do

!$omp end do

            end if

          end if

          if(ebn.eq.1.and.jsub.eq.njsub-1) then

            if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,sb1,sb2i,sb2j,radwe,radsn)

              do k=2,nk-2
                sb1=sgpv(1,njm1,k)+std(1,njm1,k)*tpdt
                sb2i=sgpv(2,njm1,k)+std(2,njm1,k)*tpdt
                sb2j=sgpv(1,njm2,k)+std(1,njm2,k)*tpdt

                radwe=scpx(njm1,k,1)                                    &
     &            *((s(2,njm1,k)-s(1,njm1,k))-(sb2i-sb1))

                radsn=scpy(1,k,2)                                       &
     &            *((s(1,njm2,k)-s(1,njm1,k))-(sb2j-sb1))

                s(1,njm1,k)=s(1,njm1,k)+std(1,njm1,k)*dts               &
     &            -(radwe-radsn)-dmpdt*(s(1,njm1,k)-sb1)

              end do

!$omp end do

            end if

            if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,sb1,sb2i,sb2j,radwe,radsn)

              do k=2,nk-2
                sb1=sgpv(nim1,njm1,k)+std(nim1,njm1,k)*tpdt
                sb2i=sgpv(nim2,njm1,k)+std(nim2,njm1,k)*tpdt
                sb2j=sgpv(nim1,njm2,k)+std(nim1,njm2,k)*tpdt

                radwe=scpx(njm1,k,2)                                    &
     &            *((s(nim2,njm1,k)-s(nim1,njm1,k))-(sb2i-sb1))

                radsn=scpy(nim1,k,2)                                    &
     &            *((s(nim1,njm2,k)-s(nim1,njm1,k))-(sb2j-sb1))

                s(nim1,njm1,k)=s(nim1,njm1,k)+std(nim1,njm1,k)*dts      &
     &            +(radwe+radsn)-dmpdt*(s(nim1,njm1,k)-sb1)

              end do

!$omp end do

            end if

          end if

        end if

      end if

! -----

! Force the west boundary value to the external boundary value.

      if(ebw.eq.1.and.isub.eq.0) then

        if(abs(wbc).ne.1) then

          if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(j,k,sb1,sb2)

            do k=2,nk-2
            do j=2,nj-2
              sb1=sgpv(1,j,k)+std(1,j,k)*tpdt
              sb2=sgpv(2,j,k)+std(2,j,k)*tpdt

              s(1,j,k)=s(1,j,k)+std(1,j,k)*dts                          &
     &          -scpx(j,k,1)*((s(2,j,k)-s(1,j,k))-(sb2-sb1))            &
     &          -dmpdt*(s(1,j,k)-sb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              s(1,j,k)=s(1,j,k)+std(1,j,k)*dts
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Force the east boundary value to the external boundary value.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(abs(ebc).ne.1) then

          if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(j,k,sb1,sb2)

            do k=2,nk-2
            do j=2,nj-2
              sb1=sgpv(nim1,j,k)+std(nim1,j,k)*tpdt
              sb2=sgpv(nim2,j,k)+std(nim2,j,k)*tpdt

              s(nim1,j,k)=s(nim1,j,k)+std(nim1,j,k)*dts                 &
     &          +scpx(j,k,2)*((s(nim2,j,k)-s(nim1,j,k))-(sb2-sb1))      &
     &          -dmpdt*(s(nim1,j,k)-sb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              s(nim1,j,k)=s(nim1,j,k)+std(nim1,j,k)*dts
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Force the south boundary value to the external boundary value.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(i,k,sb1,sb2)

          do k=2,nk-2
          do i=2,ni-2
            sb1=sgpv(i,1,k)+std(i,1,k)*tpdt
            sb2=sgpv(i,2,k)+std(i,2,k)*tpdt

            s(i,1,k)=s(i,1,k)+std(i,1,k)*dts                            &
     &        -scpy(i,k,1)*((s(i,2,k)-s(i,1,k))-(sb2-sb1))              &
     &        -dmpdt*(s(i,1,k)-sb1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-2
          do i=1,ni-1
            s(i,1,k)=s(i,1,k)+std(i,1,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

! Force the north boundary value to the external boundary value.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(i,k,sb1,sb2)

          do k=2,nk-2
          do i=2,ni-2
            sb1=sgpv(i,njm1,k)+std(i,njm1,k)*tpdt
            sb2=sgpv(i,njm2,k)+std(i,njm2,k)*tpdt

            s(i,njm1,k)=s(i,njm1,k)+std(i,njm1,k)*dts                   &
     &        +scpy(i,k,2)*((s(i,njm2,k)-s(i,njm1,k))-(sb2-sb1))        &
     &        -dmpdt*(s(i,njm1,k)-sb1)

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k)

          do k=2,nk-2
          do i=1,ni-1
            s(i,njm1,k)=s(i,njm1,k)+std(i,njm1,k)*dts
          end do
          end do

!$omp end do

        end if

      end if

! -----

!$omp end parallel

#endif

! Dump output data at target call
if (dump_call_count_exbcss == DUMP_TARGET_exbcss .and. .not. dump_done_exbcss) then
  call dump_array_3d('s_ref.bin', s, 0, ni+1, 0, nj+1, 1, nk)
  call dump_finalize()
  dump_done_exbcss = .true.
end if


call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_exbcss

!-----7--------------------------------------------------------------7--

      end module m_exbcss

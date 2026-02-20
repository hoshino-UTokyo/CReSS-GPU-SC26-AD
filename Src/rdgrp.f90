!***********************************************************************
      module m_rdgrp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/20, 2007/04/11, 2007/08/24,
!                   2007/10/19, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2008/10/10, 2009/01/05, 2009/02/27, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the group domain arrangement file
!     and check them.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_castgrp
      use m_comprofile
      use m_dump_kernel
      use m_chkerr
      use m_chkstd
      use m_comgrp
      use m_comkind
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_putunit

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdgrp, s_rdgrp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdgrp

        module procedure s_rdgrp

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
      subroutine s_rdgrp(fpexprim,fpcrsdir,fpncexp,fpnccrs)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpncexp
                       ! Formal parameter of unique index of ncexp

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) grpfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer igc      ! Array index in x direction
      integer jgc      ! Array index in y direction

      integer iered    ! Index of group domain at east boundary
                       ! in reductional entire domain
      integer jnred    ! Index of group domain at north boundary
                       ! in reductional entire domain

      integer nsrtmp   ! Temporary used maximum serial number
                       ! of group domain in entire domain

      integer iogrp    ! Unit number of group domain arrangement file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

! Internal private variables

      integer igc_sub  ! Substitute for igc
      integer jgc_sub  ! Substitute for jgc


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_rdgrp = 0
      integer, parameter :: DUMP_TARGET_rdgrp = 1
      logical, save :: dump_done_rdgrp = .false.


!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(exprim)
      call inichar(crsdir)

! -----

! Get the required namelist variables.

      call getcname(fpexprim,exprim)
      call getcname(fpcrsdir,crsdir)
      call getiname(fpncexp,ncexp)
      call getiname(fpnccrs,nccrs)

! -----

!! Open and read out the data from the group domain arrangement file.

      if(ngrp.ge.3.and.ngrp.ne.nsrl) then

! Initialize the character variable.

        if(mype.eq.root) then

          call inichar(grpfl)

        end if

! -----

! Get the unit number.

        if(mype.eq.root) then

          call getunit(iogrp)

        end if

! -----

! Open the group domain arrangement file.

        if(mype.eq.root) then

          grpfl(1:ncexp)=exprim(1:ncexp)

          write(grpfl(ncexp+1:ncexp+11),'(a11)') 'arrange.txt'

          open(iogrp,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//grpfl(1:ncexp+11),                 &
     &         status='old',access='sequential',form='formatted',       &
     &         blank='null',position='rewind',action='read')

        else

          stat=0

        end if

  100   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('rdgrp   ',5,'cont',1,'              ',14,     &
     &                   iogrp,stat)

          end if

          call cpondpe

          call destroy('rdgrp   ',5,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('rdgrp   ',5,grpfl,ncexp+11,iogrp,1,0,0_i8)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Read out the data from the group domain arrangement file.

        if(mype.eq.root) then

          do_head: do

            read(iogrp,'(a1)',iostat=stat,end=110,err=110)              &
     &        chrgrp(1,1)(1:1)

            if(chrgrp(1,1)(1:1).ne.'#') then

              backspace(iogrp,iostat=stat,err=110)

              exit do_head

            end if

          end do do_head

          do jgc=njgrp,1,-1

            read(iogrp,'(4999a1)',iostat=stat,end=110,err=110)          &
     &          (chrgrp(igc,jgc)(1:1),igc=1,nigrp)

          end do

        else

          stat=0

        end if

  110   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('rdgrp   ',5,'cont',3,'              ',14,     &
     &                   iogrp,stat)

          end if

          call cpondpe

          call destroy('rdgrp   ',5,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('rdgrp   ',5,grpfl,108,iogrp,3,0,0_i8)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Close the group domain arrangement file.

        if(mype.eq.root) then

          close(iogrp,iostat=stat,err=120,status='keep')

        else

          stat=0

        end if

  120   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('rdgrp   ',5,'cont',2,'              ',14,     &
     &                   iogrp,stat)

          end if

          call cpondpe

          call destroy('rdgrp   ',5,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('rdgrp   ',5,grpfl,108,iogrp,2,0,0_i8)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Return the unit number.

        if(mype.eq.root) then

          call putunit(iogrp)

        end if

! -----

! Broadcast the group domain arrangemnet table.

        call castgrp

! -----

!! -----

! Initialize the group domain arrangement.

      else

        do jgc=1,njgrp
        do igc=1,nigrp

          write(chrgrp(igc,jgc)(1:1),'(a1)') 'o'

        end do
        end do

      end if

! -----

! Check the group domain arrangement.

      nsrtmp=-1

      do_jgc: do jgc=1,njgrp
      do_igc: do igc=1,nigrp

        if(chrgrp(igc,jgc)(1:1).eq.'x') then

          grpxy(igc,jgc)=-1

        else if(chrgrp(igc,jgc)(1:1).eq.'o') then

          nsrtmp=nsrtmp+1

          if(nsrtmp.ge.nsrl) then

            exit do_jgc

          end if

          grpxy(igc,jgc)=nsrtmp

          xgrp(nsrtmp)=igc
          ygrp(nsrtmp)=jgc

        else

          nsrtmp=-1

          exit do_jgc

        end if

      end do do_igc
      end do do_jgc

      if(nsrtmp.eq.nsrl-1) then
        stat=0
      else
        stat=1
      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdgrp   ',5,'cont',11,'              ',14,101,  &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdgrp   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! Set the parameters of reductional entire domain.

! Initialize the parameters of reductional entire domain.

      iwred=nigrp+1
      jsred=njgrp+1

      iered=0
      jnred=0

! -----

! Get the parameters of reductional entire domain.

!@llm start meta_info ----------------------------------------------------
! Location: rdgrp.f90 :: s_rdgrp
! Summary : Find min/max indices of active group domains to determine reductional domain bounds
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Uses min/max reductions on iwred, jsred, iered, jnred
!   - No synchronization constructs
!   - Simple 2D loop with conditional bounds checking
! Next:
!   - Convert to OpenACC with teams distribute parallel for and reduction clause
!   - Alternatively use OpenACC with parallel loop reduction
! Runtime:
!   - Calls: 1
!   - AvgLoops: 1
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.014ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('rdgrp.f90', 's_rdgrp', &
   & 'OMP section 1')
end if
loop_len = int((njgrp)-(1)+1,8) * int((nigrp)-(1)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_rdgrp = dump_call_count_rdgrp + 1
if (dump_call_count_rdgrp == DUMP_TARGET_rdgrp .and. .not. dump_done_rdgrp) then
  call dump_init('rdgrp')
  call dump_scalar_i('nigrp', nigrp)
  call dump_scalar_i('njgrp', njgrp)
  call dump_array_2d_int('grpxy.bin', grpxy, 1, nigrp, 1, njgrp)
  call dump_scalar_i('iwred_in', iwred)
  call dump_scalar_i('jsred_in', jsred)
  call dump_scalar_i('iered_in', iered)
  call dump_scalar_i('jnred_in', jnred)
end if

!$omp parallel default(shared)

!$omp do schedule(runtime) private(igc_sub,jgc_sub)                     &
!$omp&   reduction(min: iwred,jsred) reduction(max: iered,jnred)

      do jgc_sub=1,njgrp
      do igc_sub=1,nigrp

        if(grpxy(igc_sub,jgc_sub).ge.0) then

          if(igc_sub.lt.iwred) then
            iwred=igc_sub
          end if

          if(igc_sub.gt.iered) then
            iered=igc_sub
          end if

          if(jgc_sub.lt.jsred) then
            jsred=jgc_sub
          end if

          if(jgc_sub.gt.jnred) then
            jnred=jgc_sub
          end if

        end if

      end do
      end do

!$omp end do

!$omp end parallel

! Dump output data at target call
if (dump_call_count_rdgrp == DUMP_TARGET_rdgrp .and. .not. dump_done_rdgrp) then
  call dump_scalar_i('iwred_ref', iwred)
  call dump_scalar_i('jsred_ref', jsred)
  call dump_scalar_i('iered_ref', iered)
  call dump_scalar_i('jnred_ref', jnred)
  call dump_finalize()
  dump_done_rdgrp = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

! Finally reset the parameters of reductional entire domain.

      iwred=iwred-1
      jsred=jsred-1

      iered=iered-1
      jnred=jnred-1

      nired=iered-iwred+1
      njred=jnred-jsred+1

      nred=nired*njred

! -----

!! -----

      end subroutine s_rdgrp

!-----7--------------------------------------------------------------7--

      end module m_rdgrp

!***********************************************************************
      module m_allociot
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 2000/01/17,
!                   2003/04/30, 2003/05/19, 2004/01/09, 2004/08/31,
!                   2004/09/25, 2005/01/31, 2005/02/10, 2006/09/21,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/01/31,
!                   2007/05/14, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/01/05, 2009/02/27, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the table of unit numbers.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comprofile
      use m_comfile
      use m_comionum
      use m_commpi
      use m_cpondpe
      use m_destroy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allociot, s_allociot

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allociot

        module procedure s_allociot

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
      subroutine s_allociot
!***********************************************************************

! Internal shared variables

      character(len=26) chkfl
                       ! Opened file name to check

      integer iio      ! Index of unit numbers table

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

! Internal private variable

      integer iio_sub  ! Substitute for iio


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

!! Get the maximum number of unit.

! Set the maximum number of unit.

      if(mype.eq.root) then

        nio=-1

        do_iio_sub: do iio=11,nsub+11

          write(chkfl(1:26),'(a17,i5.5,a4)')                            &
     &                      'simultaneous.unit',iio,'.txt'

          open(iio,iostat=stat,err=100,                                 &
     &         file=curdir(1:2)//chkfl(1:26),                           &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='readwrite')

  100     if(stat.ne.0) then

            nio=iio-11

            exit do_iio_sub

          end if

        end do do_iio_sub

        if(nio.eq.0) then

          stat=1

        else

          if(nio.lt.0) then
            nio=nsub+1
          end if

          do iio=11,nio+10

            close(iio,iostat=stat,err=110,status='delete')

          end do

          if(nio.eq.1) then
            stat=1
          end if

        end if

      else

        nio=nsub+1

        stat=0

      end if

! -----

! If error occured, call the procedure destroy.

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allociot',8,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allociot',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

!! Allocate the table of unit numbers.

! Perform allocate.

      stat=0

      allocate(iolst(1:nio),stat=cstat)

      stat=stat+abs(cstat)

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allociot',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allociot',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

! Initialize the table of unit numbers.

!@llm start meta_info ----------------------------------------------------
! Location: allociot.f90 :: s_allociot
! Summary : Initialize the I/O unit number table (iolst) with sequential
!           unit numbers starting from 11
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to module-level array iolst from m_comionum
!   - Simple 1D loop with no data dependencies
!   - Very small loop size (nio typically small)
! Next:
!   - Trivial GPU port but likely not worth offloading
!   - Small array size means CPU execution is faster
!   - Keep on CPU; initialization is one-time cost
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('allociot.f90', 's_allociot', &
   & 'OMP section 1')
end if
loop_len = int((nio)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared)

!$omp do schedule(runtime) private(iio_sub)

      do iio_sub=1,nio
        iolst(iio_sub)=iio_sub+10
      end do

!$omp end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_allociot

!-----7--------------------------------------------------------------7--

      end module m_allociot

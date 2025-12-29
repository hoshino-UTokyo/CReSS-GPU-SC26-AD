!***********************************************************************
      module m_setasl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the interpolated aerosol variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname
      use m_comprofile

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setasl, s_setasl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setasl

        module procedure s_setasl

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
      subroutine s_setasl(fpaslitv,ird,ni,nj,nk,nqa,qagpv,qatd)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpaslitv
                       ! Formal parameter of unique index of aslitv

      integer, intent(in) :: ird
                       ! Index of count to read out in rdaslnxt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

! Input and output variables

      real, intent(inout) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

      real, intent(inout) :: qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Time tendency of mixing ratio of aerosol data

! Internal shared variables

      real aslitv      ! Time interval of aerosol data

      real asliv       ! Inverse of aslitv

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpaslitv,aslitv)

! -----

! Set the common used variable.

      asliv=1.e0/aslitv

! -----

!! Set the interpolated aerosol variables.

!@llm start meta_info ----------------------------------------------------
! Location: setasl.f90 :: s_setasl
! Summary : Set interpolated aerosol variables by computing time tendency
!           or copying values based on read index (ird)
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Simple 4D loops with direct array operations
!   - Conditional execution based on ird (if ird==1 or ird==2)
!   - Writes to qatd and qagpv arrays
!   - No synchronization constructs
! Next:
!   - Convert to OpenACC with parallel loop collapse(4)
!   - Single parallel region can be offloaded with appropriate data clauses
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('setasl.f90', 's_setasl', &
   & 'OMP section 1')
end if
loop_len = int((nqa(0))-(1)+1,8) &
     & * int((nk)-(1)+1,8) &
     & * int((nj)-(1)+1,8) &
     & * int((ni)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k,n)

! Set the time tendency of variables at current marked time.

      if(ird.eq.1) then

        do n=1,nqa(0)

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni
              qatd(i,j,k,n)=(qatd(i,j,k,n)-qagpv(i,j,k,n))*asliv
            end do
            end do

!$omp end do

          end do

        end do

      end if

! -----

! Set the variables at current marked time.

      if(ird.eq.2) then

        do n=1,nqa(0)

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni
              qagpv(i,j,k,n)=qatd(i,j,k,n)
            end do
            end do

!$omp end do

          end do

        end do

      end if

! -----

!$omp end parallel

call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_setasl

!-----7--------------------------------------------------------------7--

      end module m_setasl

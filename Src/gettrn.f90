!***********************************************************************
      module m_gettrn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/06/07, 1999/06/21,
!                   1999/08/18, 1999/08/23, 1999/09/01, 1999/09/06,
!                   1999/09/30, 1999/10/12, 2000/01/17, 2000/04/18,
!                   2000/12/18, 2001/01/15, 2001/05/29, 2001/11/20,
!                   2002/04/02, 2002/06/18, 2002/07/15, 2002/10/31,
!                   2003/04/30, 2003/05/19, 2003/10/31, 2005/08/05,
!                   2005/12/13, 2006/09/21, 2007/01/05, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the terrain height.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comprofile
      use m_dump_kernel
      use m_getiname
      use m_getrname
      use m_rdtrn

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: gettrn, s_gettrn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface gettrn

        module procedure s_gettrn

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_gettrn(fptrnopt,fpzsfc,fpmnthgh,fpmntwx,fpmntwy,     &
     &                    fpmntcx,fpmntcy,dvname,ncdvn,fmsg,            &
     &                    ni,nj,xs,ys,ht)
!***********************************************************************

! Input variables

      character(len=12), intent(in) :: dvname
                       ! Optional data variable name

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpzsfc
                       ! Formal parameter of unique index of zsfc

      integer, intent(in) :: fpmnthgh
                       ! Formal parameter of unique index of mnthgh

      integer, intent(in) :: fpmntwx
                       ! Formal parameter of unique index of mntwx

      integer, intent(in) :: fpmntwy
                       ! Formal parameter of unique index of mntwy

      integer, intent(in) :: fpmntcx
                       ! Formal parameter of unique index of mntcx

      integer, intent(in) :: fpmntcy
                       ! Formal parameter of unique index of mntcy

      integer, intent(in) :: ncdvn
                       ! Number of character of dvname

      integer, intent(in) :: fmsg
                       ! Control flag of message type in outstd03

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: xs(0:ni+1)
                       ! x coordinates at scalar points

      real, intent(in) :: ys(0:nj+1)
                       ! y coordinates at scalar points

! Output variable

      real, intent(out) :: ht(0:ni+1,0:nj+1)
                       ! Terrain height

! Internal shared variables

      integer trnopt   ! Option for terrain height setting

      real zsfc        ! Sea surface terrain height

      real mnthgh(1:2) ! Flat or bell shaped mountain height
                       ! and base level height

      real mntwx       ! Bell shaped mountain width in x direction
      real mntwy       ! Bell shaped mountain width in y direction

      real mntcx       ! Center in x coordinates of
                       ! bell shaped mountain

      real mntcy       ! Center in y coordinates of
                       ! bell shaped mountain

      real wxiv        ! 1.0 / mntwx
      real wyiv        ! 1.0 / mntwy

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real a           ! Temporary variable
      real b           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer, save :: prof_id2 = -1
      integer(8) :: loop_len

      ! Dump variables
      integer, save :: dump_call_count_gettrn = 0
      integer, parameter :: DUMP_TARGET_gettrn = 1
      logical, save :: dump_done_gettrn = .false.


!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fptrnopt,trnopt)

! -----

!! Fill in the array ht with the surface height in the case the flat
!! terrain is applied.

      if(trnopt.eq.0) then

! Get the required namelist variables.

        call getrname(fpzsfc,zsfc)
        call getrname(fpmnthgh,mnthgh(1))
        call getrname(fpmnthgh+1,mnthgh(2))

! -----

! Set the flat terrain.

!@llm start meta_info ----------------------------------------------------
! Location: gettrn.f90 :: s_gettrn (trnopt=0 branch)
! Summary : Initialize terrain height array with constant flat value
!           using max of mountain height+base or sea surface height.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Uses intrinsic max function (GPU compatible)
!   - Simple constant assignment to all grid points
!   - No global writes, only output array ht is modified
! Next:
!   - Direct translation to OpenACC with teams distribute
!   - Consider using GPU memset-like operation for constant fill
! Runtime:
!   - Calls: 1
!   - AvgLoops: 810.0K
!   - TotalTime: 0.000s (0.00%)
!   - AvgTime: 0.015ms
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('gettrn.f90', 's_gettrn', &
   & 'OMP section 1')
end if
loop_len = int((nj)-(0)+1,8) * int((ni)-(0)+1,8)
call profile_start(prof_id1)


! Dump input data at target call
dump_call_count_gettrn = dump_call_count_gettrn + 1
if (dump_call_count_gettrn == DUMP_TARGET_gettrn .and. .not. dump_done_gettrn) then
  call dump_init('gettrn')
  call dump_scalar_i('trnopt', trnopt)
  call dump_scalar_r('zsfc', zsfc)
  call dump_scalar_r('mnthgh1', mnthgh(1))
  call dump_scalar_r('mnthgh2', mnthgh(2))
  call dump_scalar_r('mntwx', mntwx)
  call dump_scalar_r('mntwy', mntwy)
  call dump_scalar_r('mntcx', mntcx)
  call dump_scalar_r('mntcy', mntcy)
  call dump_scalar_i('ncdvn', ncdvn)
  call dump_scalar_i('fmsg', fmsg)
  call dump_scalar_i('ni', ni)
  call dump_scalar_i('nj', nj)
  call dump_scalar_r('wxiv', wxiv)
  call dump_scalar_r('wyiv', wyiv)
  ! FIXME: xs is an array, not scalar
  ! ! FIXME: xs is array - call dump_scalar_r('xs', xs)
  ! FIXME: ys is an array, not scalar
  ! ! FIXME: ys is array - call dump_scalar_r('ys', ys)
end if

#if defined(USE_GPU) && !defined(DISABLE_GPU_140)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------

    !$acc kernels
    !$acc loop independent
    do j = 0, nj
      !$acc loop independent
      do i = 0, ni
        ht(i,j) = max(mnthgh(1)+mnthgh(2),zsfc)
      end do
    end do
    !$acc end kernels

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          ht(i,j)=max(mnthgh(1)+mnthgh(2),zsfc)
        end do
        end do

!$omp end do

!$omp end parallel

#endif

! Dump output data at target call
if (dump_call_count_gettrn == DUMP_TARGET_gettrn .and. .not. dump_done_gettrn) then
  call dump_array_2d('ht_ref.bin', ht, 0, ni+1, 0, nj+1)
  call dump_finalize()
  dump_done_gettrn = .true.
end if


call profile_stop(prof_id1, loop_len)

! -----

!! -----

!! Set the bell shaped mountain.

      else if(trnopt.eq.1) then

! Get the required namelist variables.

        call getrname(fpzsfc,zsfc)
        call getrname(fpmnthgh,mnthgh(1))
        call getrname(fpmnthgh+1,mnthgh(2))
        call getrname(fpmntwx,mntwx)
        call getrname(fpmntwy,mntwy)
        call getrname(fpmntcx,mntcx)
        call getrname(fpmntcy,mntcy)

! -----

! Set the common used variables.

        wxiv=1.e0/mntwx
        wyiv=1.e0/mntwy

! -----

! Set the bell shaped mountain.

!@llm start meta_info ----------------------------------------------------
! Location: gettrn.f90 :: s_gettrn (trnopt=1 branch)
! Summary : Generate bell-shaped mountain terrain using Gaussian-like
!           formula with configurable height, width, and center position.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Uses intrinsic max function (GPU compatible)
!   - Reads from 1D arrays xs(i), ys(j) - need to ensure GPU accessible
!   - Simple arithmetic with division and max
!   - No global writes, only output array ht is modified
! Next:
!   - Direct translation to OpenACC with teams distribute
!   - Ensure xs and ys arrays are mapped to device
!@llm end meta_info ------------------------------------------------------
#if defined(USE_GPU) && !defined(DISABLE_GPU_141)
!----------------------------------------------------------------------
! GPU version (OpenACC)
!----------------------------------------------------------------------

    ! NOTE: gettrn_sec2 benchmark uses max reduction for htmax,
    ! but the actual section 2 in source is the bell-shaped mountain computation.
    ! We use acc kernels for the bell-shaped terrain calculation.
    !$acc kernels
    !$acc loop independent
    do j = 0, nj
      !$acc loop independent private(a,b)
      do i = 0, ni
        a = wxiv*(xs(i)-mntcx)
        b = wyiv*(ys(j)-mntcy)
        ht(i,j) = max(mnthgh(1)/(1.e0+(a*a+b*b))+mnthgh(2),zsfc)
      end do
    end do
    !$acc end kernels

#else
!----------------------------------------------------------------------
! CPU version (OpenMP) - Original code preserved
!----------------------------------------------------------------------
!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,a,b)

        do j=0,nj
        do i=0,ni
          a=wxiv*(xs(i)-mntcx)
          b=wyiv*(ys(j)-mntcy)

          ht(i,j)=max(mnthgh(1)/(1.e0+(a*a+b*b))+mnthgh(2),zsfc)

        end do
        end do

!$omp end do

!$omp end parallel

#endif

! -----

!! -----

! Read out the data from the terrain file.

      else if(trnopt.eq.2) then

        call rdtrn(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,          &
     &             dvname,ncdvn,fmsg,ni,nj,ht)

      end if

! -----

      end subroutine s_gettrn

!-----7--------------------------------------------------------------7--

      end module m_gettrn

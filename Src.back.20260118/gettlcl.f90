!***********************************************************************
      module m_gettlcl
!***********************************************************************

!     Author      : Monoe Daisuke, Sakakibara Atsushi
!     Date        : 2010/09/16, 2010/09/22
!     Modification: 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the Lifting Condensation Level.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comprofile

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: gettlcl, s_gettlcl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface gettlcl

        module procedure s_gettlcl

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
      subroutine s_gettlcl(datype,k,ni,nj,nk,pbr,ptbr,pp,ptp,qv,tlcl)
!***********************************************************************

! Input variables

      character(len=1), intent(in) :: datype
                       ! Input data type

      integer, intent(in) :: k
                       ! Input array index in z direction

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

! Output variable

      real, intent(out) :: tlcl(0:ni+1,0:nj+1)
                       ! Air Temperature at Lifting Condensation Level

! Internal shared variables

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real t           ! Air temperature
      real p           ! Pressure
      real ea          ! Pertial vapor pressure at specified plane

      real a           ! Temporary variable
      real b           ! Temporary variable


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Calculate the air temperature at Lifting Condensation Level.

!@llm start meta_info ----------------------------------------------------
! Location: gettlcl.f90 :: s_gettlcl
! Summary : Calculate temperature at Lifting Condensation Level using
!           Bolton's formula based on mixing ratio or relative humidity.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls within parallel region
!   - Uses intrinsic exp and log functions (GPU compatible)
!   - Conditional branch (if/else if) selects calculation method
!   - No global writes, only output array tlcl is modified
!   - No synchronization constructs other than implicit barriers
! Next:
!   - Direct translation to OpenACC with teams distribute
!   - Both branches are simple arithmetic, suitable for GPU
!   - Consider unifying branches or using separate kernels per datype
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('gettlcl.f90', 's_gettlcl', &
   & 'OMP section 1')
end if
loop_len = int((nj-1)-(1)+1,8) * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared)

      if(datype(1:1).eq.'m') then

!$omp do schedule(runtime) private(i,j,p,t,ea,a)

        do j=1,nj-1
        do i=1,ni-1

          p=pbr(i,j,k)+pp(i,j,k)
          t=(ptbr(i,j,k)+ptp(i,j,k))*exp(rddvcp*log(p0iv*p))

          ea=p*qv(i,j,k)/(epsva+qv(i,j,k))

          a=3.5e0*log(t)-log(1.e-2*ea)-4.805e0

          tlcl(i,j)=2.84e3/a+5.5e1

        end do
        end do

!$omp end do

      else if(datype(1:1).eq.'r') then

!$omp do schedule(runtime) private(i,j,p,t,ea,a,b)

        do j=1,nj-1
        do i=1,ni-1

          p=pbr(i,j,k)+pp(i,j,k)
          t=(ptbr(i,j,k)+ptp(i,j,k))*exp(rddvcp*log(p0iv*p))

          a=t-5.5e1
          b=1.e0/a

          a=log(1.e-2*qv(i,j,k))/2.84e3

          tlcl(i,j)=1.e0/(b-a)+5.5e1

        end do
        end do

!$omp end do

      end if

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_gettlcl

!-----7--------------------------------------------------------------7--

      end module m_gettlcl

!***********************************************************************
      module m_bc4news
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/31, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions at the four corners.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comprofile
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bc4news, s_bc4news

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bc4news

        module procedure s_bc4news

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
      subroutine s_bc4news(fpwbc,fpebc,fpsbc,fpnbc,isw,ise,jss,jsn,     &
     &                     ni,nj,kmax,var)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: isw
                       ! Specified index of west corner in x direction

      integer, intent(in) :: ise
                       ! Specified index of east corner in x direction

      integer, intent(in) :: jss
                       ! Specified index of south corner in y direction

      integer, intent(in) :: jsn
                       ! Specified index of north corner in y direction

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

! Input and output variable

      real, intent(inout) :: var(0:ni+1,0:nj+1,1:kmax)
                       ! Optional variable

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer iswp1    ! isw + 1
      integer isem1    ! ise - 1

      integer jssp1    ! jss + 1
      integer jsnm1    ! jsn - 1

! Internal private variable

      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

! Set the common used variables.

      iswp1=isw+1
      isem1=ise-1

      jssp1=jss+1
      jsnm1=jsn-1

! -----

! Set the boundary conditions at the four corners.

!@llm start meta_info ----------------------------------------------------
! Location: bc4news.f90 :: s_bc4news
! Summary : Sets boundary conditions at the four corners (SW, SE, NW, NE) by
!           averaging adjacent boundary values for optional 3D variable.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Simple 1D loops over k dimension with direct array assignments
!   - Uses module variables from m_commpi (ebsw, ebse, ebnw, ebne, isub, jsub, nisub, njsub)
!   - Multiple conditionally executed small loops based on domain decomposition position
! Next:
!   - Convert to OpenACC with collapsed loops
!   - Consider merging the four conditional loops into a single kernel with conditional logic
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('bc4news.f90', 's_bc4news', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared)

      if(abs(wbc).ne.1.or.abs(ebc).ne.1                                 &
     &  .or.abs(sbc).ne.1.or.abs(nbc).ne.1) then

        if(ebsw.eq.1.and.isub.eq.0.and.jsub.eq.0) then

!$omp do schedule(runtime) private(k)

          do k=1,kmax
            var(isw,jss,k)=.5e0*(var(iswp1,jss,k)+var(isw,jssp1,k))
          end do

!$omp end do

        end if

        if(ebse.eq.1.and.isub.eq.nisub-1.and.jsub.eq.0) then

!$omp do schedule(runtime) private(k)

          do k=1,kmax
            var(ise,jss,k)=.5e0*(var(isem1,jss,k)+var(ise,jssp1,k))
          end do

!$omp end do

        end if

        if(ebnw.eq.1.and.isub.eq.0.and.jsub.eq.njsub-1) then

!$omp do schedule(runtime) private(k)

          do k=1,kmax
            var(isw,jsn,k)=.5e0*(var(iswp1,jsn,k)+var(isw,jsnm1,k))
          end do

!$omp end do

        end if

        if(ebne.eq.1.and.isub.eq.nisub-1.and.jsub.eq.njsub-1) then

!$omp do schedule(runtime) private(k)

          do k=1,kmax
            var(ise,jsn,k)=.5e0*(var(isem1,jsn,k)+var(ise,jsnm1,k))
          end do

!$omp end do

        end if

      end if

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_bc4news

!-----7--------------------------------------------------------------7--

      end module m_bc4news

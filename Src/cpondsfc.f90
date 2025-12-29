!***********************************************************************
      module m_cpondsfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2003/04/30, 2003/05/19, 2004/05/07, 2004/08/01,
!                   2004/09/10, 2005/02/10, 2006/01/10, 2006/02/13,
!                   2006/09/21, 2007/01/20, 2007/05/14, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/01/05,
!                   2009/02/27, 2009/11/13, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     correspond the surface data to the land use categories.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comprofile
      use m_commath
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_getindx
      use m_setcst2d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: cpondsfc, s_cpondsfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface cpondsfc

        module procedure s_cpondsfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_cpondsfc(fpnumctg_lnd,xo,lnduse_lnd,sfcvar_lnd,      &
     &                      imin,imax,jmin,jmax,land,sfcvar)
!***********************************************************************

! Input variables

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: fpnumctg_lnd
                       ! Formal parameter of unique index of numctg_lnd

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      integer, intent(in) :: land(imin:imax,jmin:jmax)
                       ! Land use on surface

      integer, intent(in) :: lnduse_lnd(1:100)
                       ! User specified land use namelist table

      real, intent(in) :: sfcvar_lnd(1:100)
                       ! User specified surface namelist table

! Output variable

      real, intent(out) :: sfcvar(imin:imax,jmin:jmax)
                       ! Optional surface variable

! Internal shared variables

      integer numctg_lnd
                       ! Number of land use data categories

      integer stat     ! Runtime status

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      real rstat       ! Real runtime status

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer ic       ! Index of user specified land use namelist table


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpnumctg_lnd,numctg_lnd)

! -----

! Get the maximum and minimim indices of do loops.

      call getindx(xo,imin,imax,jmin,jmax,istr,iend,jstr,jend)

! -----

! Fill in the array sfcvar with - 1.0 x 10^35.

      call setcst2d(imin,imax,jmin,jmax,lim35n,sfcvar)

! -----

! Initialize the processed variable, rstat.

      rstat=lim36

! -----

!! Correspond the surface data to the land use categories and check
!! errors.

!@llm start meta_info ----------------------------------------------------
! Location: cpondsfc.f90 :: s_cpondsfc
! Summary : Maps surface variables to land use categories from lookup table,
!           then performs error checking with reduction.
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - No external function calls inside parallel region
!   - First loop: lookup table matching with conditional assignment
!   - Second loop: reduction(min: rstat) for error checking
!   - Uses lnduse_lnd and sfcvar_lnd lookup tables (size 100)
!   - Potential race condition if multiple categories match same grid point
! Next:
!   - Lookup tables should be placed in constant memory on GPU
!   - Reduction can use GPU reduction primitives
!   - Consider restructuring first loop to avoid potential race condition
!   - May need atomic operations or different algorithm for category matching
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('cpondsfc.f90', 's_cpondsfc', &
   & 'OMP section 1')
end if
loop_len = int((numctg_lnd)-(1)+1,8) &
     & * int((jend)-(jstr)+1,8) &
     & * int((iend)-(istr)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared)

! Correspond the surface data to the land use categories.

!$omp do schedule(runtime) private(i,j,ic)

      do ic=1,numctg_lnd

        do j=jstr,jend
        do i=istr,iend

          if(land(i,j).eq.lnduse_lnd(ic)) then

            sfcvar(i,j)=sfcvar_lnd(ic)

          end if

        end do
        end do

      end do

!$omp end do

! -----

! Check errors.

!$omp do schedule(runtime) private(i,j) reduction(min: rstat)

      do j=jstr,jend
      do i=istr,iend
        rstat=min(sfcvar(i,j),rstat)
      end do
      end do

!$omp end do

! -----

!$omp end parallel

call profile_stop(prof_id1, loop_len)

!! -----

! If error occured, call the procedure destroy.

      if(rstat.lt.lim34n) then
        stat=1
      else
        stat=0
      end if

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('cpondsfc',8,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('cpondsfc',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

      end subroutine s_cpondsfc

!-----7--------------------------------------------------------------7--

      end module m_cpondsfc

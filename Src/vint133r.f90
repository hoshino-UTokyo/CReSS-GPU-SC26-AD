!***********************************************************************
      module m_vint133r
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2003/04/30, 2003/05/19, 2003/12/12, 2004/09/10,
!                   2006/09/21, 2007/10/19, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the variable to the model grid vertically.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_getindx

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vint133r, s_vint133r

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vint133r

        module procedure s_vint133r

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
      subroutine s_vint133r(xo,ni,nj,nk,zph,outvar,nlev,z1d,invar)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nlev
                       ! Horizontally averaged vertical dimension

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: z1d(1:nlev)
                       ! Horizontally averaged z physical coordinates

      real, intent(in) :: invar(0:ni+1,0:nj+1,1:nlev)
                       ! Optional variable
                       ! at horizontally averraged plane

! Output variable

      real, intent(out) :: outvar(0:ni+1,0:nj+1,1:nk)
                       ! Optional interpolated variable

! Internal shared variables

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction
      integer kend     ! Maximum do loops index in z direction

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer kl       ! Averaged array index in z direction

      real dk          ! Distance in z direction
                       ! between model and averaged points


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the maximum and minimim indices of do loops.

      call getindx(xo,1,ni,1,nj,istr,iend,jstr,jend)

      if(xo(3:3).eq.'o') then
        kend=nk-1
      else if(xo(3:3).eq.'x') then
        kend=nk-2
      end if

! -----

!! Interpolate the variable to the model grid vertically.

!@llm start meta_info ----------------------------------------------------
! Location: vint133r.f90 :: s_vint133r
! Summary : Interpolate variable to model grid with undefined value handling (radar data)
! GPU diff: Medium
! Findings:
!   - No omp_get_thread_num usage
!   - getindx() called before parallel region (safe)
!   - Uses lim35n, lim34n from m_commath for undefined value markers
!   - Reads from zph, invar (3D), z1d (1D); writes to outvar (3D)
!   - First section: fill undefined values outside flat plane range
!   - Second section: interpolate with validity check on input values
!   - Conditional branches for level selection and validity checking
! Next:
!   - Collapse k,j,i loops for GPU parallelism
!   - Use OpenACC teams distribute parallel do collapse(3)
!   - Ensure lim35n, lim34n constants are accessible on device
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('vint133r.f90', 's_vint133r', &
   & 'OMP section 1')
end if
loop_len = int((kend)-(2)+1,8) &
     & * int((jend)-(jstr)+1,8) &
     & * int((iend)-(istr)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k,kl)

! Fill in the undifined value outside of the flat plane.

      do k=2,kend

!$omp do schedule(runtime) private(i,j)

        do j=jstr,jend
        do i=istr,iend

          if(z1d(1).gt.zph(i,j,k).or.z1d(nlev).le.zph(i,j,k)) then

            outvar(i,j,k)=lim35n

          end if

        end do
        end do

!$omp end do

      end do

! -----

! Interpolate the variable.

      do kl=1,nlev-1

        do k=2,kend

!$omp do schedule(runtime) private(i,j,dk)

          do j=jstr,jend
          do i=istr,iend

            if(z1d(kl).le.zph(i,j,k).and.z1d(kl+1).gt.zph(i,j,k)) then

              if(invar(i,j,kl).gt.lim34n                                &
     &          .and.invar(i,j,kl+1).gt.lim34n) then

                dk=(zph(i,j,k)-z1d(kl))/(z1d(kl+1)-z1d(kl))

                outvar(i,j,k)=(1.e0-dk)*invar(i,j,kl)+dk*invar(i,j,kl+1)

              else

                outvar(i,j,k)=lim35n

              end if

            end if

          end do
          end do

!$omp end do

        end do

      end do

! -----

!$omp end parallel

call profile_stop(prof_id1, loop_len)

!! -----

      end subroutine s_vint133r

!-----7--------------------------------------------------------------7--

      end module m_vint133r

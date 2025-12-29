!***********************************************************************
      module m_adjstbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/08/08
!     Modification: 2006/09/30, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     adjust the mean water mass to be between their boundaries.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing
      use m_comprofile

      implicit none

! Default access control

      private

! Exceptional access control

      public :: adjstbw, s_adjstbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstbw

        module procedure s_adjstbw

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
      subroutine s_adjstbw(ni,nj,nk,nqw,nnw,bmw,mwbin,nwbin)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      real, intent(in) :: bmw(1:nqw+1,1:3)
                       ! Mass at water bin boundary [g]

! Input and output variables

      real, intent(inout) :: mwbin(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Total water mass [g/cm^3]

      real, intent(inout) :: nwbin(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations [1/cm^3]

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in water bin categories

      real mwbr        ! Mean water mass


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Adjust the mean water mass to be between their boundaries.

!@llm start meta_info ----------------------------------------------------
! Location: adjstbw.f90 :: subroutine s_adjstbw
! Summary : Adjusts mean water mass for bin microphysics to stay within
!           bin boundaries, resetting concentration if mass is invalid.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_* usage.
!   - No function calls inside parallel region.
!   - Pure arithmetic and conditional logic only.
!   - All grid points independent (embarrassingly parallel).
!   - 4D array access with bin categories in outer loop.
!   - Uses intrinsic conditionals, no sync constructs.
! Next:
!   - Direct OpenACC kernels with collapse for (n,k,j,i).
!   - bmw array is small and can be copied to device.
!@llm end meta_info ------------------------------------------------------


! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('adjstbw.f90', 's_adjstbw', &
   & 'OMP section 1')
end if
loop_len = int((nqw)-(1)+1,8) &
     & * int((nk-1)-(1)+1,8) &
     & * int((nj-1)-(1)+1,8) &
     & * int((ni-1)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k,n)

      do n=1,nqw

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,mwbr)

          do j=1,nj-1
          do i=1,ni-1

            if(mwbin(i,j,k,n).gt.0.e0.and.nwbin(i,j,k,n).gt.0.e0) then

              mwbr=mwbin(i,j,k,n)/nwbin(i,j,k,n)

              if(mwbr.lt.bmw(n,3)) then

                nwbin(i,j,k,n)=mwbin(i,j,k,n)/bmw(n,3)

              end if

              if(mwbr.gt.bmw(n+1,2)) then

                nwbin(i,j,k,n)=mwbin(i,j,k,n)/bmw(n+1,2)

              end if

            else

              mwbin(i,j,k,n)=0.e0
              nwbin(i,j,k,n)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_adjstbw

!-----7--------------------------------------------------------------7--

      end module m_adjstbw

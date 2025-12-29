!***********************************************************************
      module m_getbufsx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 1999/03/25, 1999/04/06, 1999/07/05, 1999/08/18,
!                   1999/08/23, 1999/10/07, 1999/10/12, 1999/11/01,
!                   1999/11/19, 2000/01/17, 2000/04/18, 2000/07/05,
!                   2001/05/29, 2002/04/02, 2002/06/06, 2003/04/30,
!                   2003/05/19, 2004/08/20, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the receiving buffer in x direction for sub domain.

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

      public :: getbufsx, s_getbufsx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getbufsx

        module procedure s_getbufsx

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
      subroutine s_getbufsx(fpwbc,fpebc,fproc,iwrcv,iercv,ni,nj,kmax,   &
     &                      var,ib,nb,rbufx)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: iwrcv
                       ! Array index of west received value

      integer, intent(in) :: iercv
                       ! Array index of east received value

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      integer, intent(in) :: ib
                       ! Exchanging variables number

      integer, intent(in) :: nb
                       ! Number of exchanging variables

      real, intent(in) :: rbufx(0:nj+1,1:kmax,1:2*nb)
                       ! Receiving buffer in x direction

! Output variable

      real, intent(out) :: var(0:ni+1,0:nj+1,1:kmax)
                       ! Optional exchanged variable

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundaty conditions

! Internal private variables

      integer j        ! Array index in y direction
      integer k        ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)

! -----

!! Get the receiving buffer in x direction.

      if(nisub.ge.2) then

        if(fproc(1:3).eq.'all'.or.(wbc.eq.-1.and.ebc.eq.-1)) then

!@llm start meta_info ----------------------------------------------------
! Location: getbufsx.f90 :: s_getbufsx
! Summary : Fills west/east halo regions from MPI receive buffer in x direction
!           for subdomain boundary exchange
! GPU diff: Easy
! Findings:
!   - No omp_get_thread usage
!   - No function calls inside parallel region
!   - Simple array copy from rbufx to var at boundary indices
!   - Multiple conditional branches based on boundary conditions (wbc, ebc)
!   - No synchronization constructs beyond implicit barriers
! Next:
!   - Convert to OpenACC or OpenACC data region
!   - Ensure rbufx and var are mapped appropriately
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('getbufsx.f90', 's_getbufsx', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(1)+1,8) * int((nj+1)-(0)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

! Fill in the west halo regions with the received value.

          if(wbc.eq.-1) then

            if(fproc(1:3).eq.'all') then

              do k=1,kmax

!$omp do schedule(runtime) private(j)

                do j=0,nj+1
                  var(iwrcv,j,k)=rbufx(j,k,ib)
                end do

!$omp end do

              end do

            else if(fproc(1:3).eq.'bnd') then

              if(isub.eq.0) then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    var(iwrcv,j,k)=rbufx(j,k,ib)
                  end do

!$omp end do

                end do

              end if

            end if

          else if(wbc.ne.-1) then

            if(isub.ne.0) then

              do k=1,kmax

!$omp do schedule(runtime) private(j)

                do j=0,nj+1
                  var(iwrcv,j,k)=rbufx(j,k,ib)
                end do

!$omp end do

              end do

            end if

          end if

! -----

! Fill in the east halo regions with the received value.

          if(ebc.eq.-1) then

            if(fproc(1:3).eq.'all') then

              do k=1,kmax

!$omp do schedule(runtime) private(j)

                do j=0,nj+1
                  var(iercv,j,k)=rbufx(j,k,nb+ib)
                end do

!$omp end do

              end do

            else if(fproc(1:3).eq.'bnd') then

              if(isub.eq.nisub-1) then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    var(iercv,j,k)=rbufx(j,k,nb+ib)
                  end do

!$omp end do

                end do

              end if

            end if

          else if(ebc.ne.-1) then

            if(isub.ne.nisub-1) then

              do k=1,kmax

!$omp do schedule(runtime) private(j)

                do j=0,nj+1
                  var(iercv,j,k)=rbufx(j,k,nb+ib)
                end do

!$omp end do

              end do

            end if

          end if

! -----

!$omp end parallel

call profile_stop(prof_id1, loop_len)

        end if

      end if

!! -----

      end subroutine s_getbufsx

!-----7--------------------------------------------------------------7--

      end module m_getbufsx

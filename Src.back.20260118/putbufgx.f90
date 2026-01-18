!***********************************************************************
      module m_putbufgx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/20, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     put the sending buffer in x direction for group domain.

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

      public :: putbufgx, s_putbufgx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface putbufgx

        module procedure s_putbufgx

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
      subroutine s_putbufgx(fpwbc,fpebc,fproc,iwsnd,iesnd,ni,nj,kmax,   &
     &                      var,ib,nb,sbufx)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: iwsnd
                       ! Array index of west sended value

      integer, intent(in) :: iesnd
                       ! Array index of east sended value

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

      real, intent(in) :: var(0:ni+1,0:nj+1,1:kmax)
                       ! Optional exchanged variable

! Output variable

      real, intent(out) :: sbufx(0:nj+1,1:kmax,1:2*nb)
                       ! Sending buffer in x direction

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

!! Put the sending buffer in x direction.

      if(nigrp.ge.2.and.(wbc.ne.-1.and.ebc.ne.-1)) then

        if(fproc(1:3).eq.'all'.or.(wbc.eq.1.and.ebc.eq.1)) then

!@llm start meta_info ----------------------------------------------------
! Location: putbufgx.f90 :: s_putbufgx
! Summary : Fill sending buffer in x direction for group domain boundary
!           exchange (west and east halo regions).
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - No function calls inside parallel region
!   - Writes to output array sbufx from input var
!   - Simple 2D copy operations (j,k loops)
!   - Conditional execution based on subdomain position (isub, igrp)
!   - No explicit barriers but implicit at !$omp end do
! Next:
!   - Use OpenACC parallel loop for buffer packing
!   - Consider async data transfers for overlap with computation
!   - May keep on host if buffer sizes are small relative to transfer cost
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('putbufgx.f90', 's_putbufgx', &
   & 'OMP section 1')
end if
loop_len = int((kmax)-(1)+1,8) * int((nj+1)-(0)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(k)

! Fill in the sending buffer with the value in the west halo regions.

          if(isub.eq.0) then

            if(wbc.eq.1) then

              if(fproc(1:3).eq.'all') then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    sbufx(j,k,ib)=var(iwsnd,j,k)
                  end do

!$omp end do

                end do

              else if(fproc(1:3).eq.'bnd') then

                if(igrp.eq.0) then

                  do k=1,kmax

!$omp do schedule(runtime) private(j)

                    do j=0,nj+1
                      sbufx(j,k,ib)=var(iwsnd,j,k)
                    end do

!$omp end do

                  end do

                end if

              end if

            else if(wbc.ne.1) then

              if(ebw.eq.0) then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    sbufx(j,k,ib)=var(iwsnd,j,k)
                  end do

!$omp end do

                end do

              end if

            end if

          end if

! -----

! Fill in the sending buffer with the value in the east halo regions.

          if(isub.eq.nisub-1) then

            if(ebc.eq.1) then

              if(fproc(1:3).eq.'all') then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    sbufx(j,k,nb+ib)=var(iesnd,j,k)
                  end do

!$omp end do

                end do

              else if(fproc(1:3).eq.'bnd') then

                if(igrp.eq.nigrp-1) then

                  do k=1,kmax

!$omp do schedule(runtime) private(j)

                    do j=0,nj+1
                      sbufx(j,k,nb+ib)=var(iesnd,j,k)
                    end do

!$omp end do

                  end do

                end if

              end if

            else if(ebc.ne.1) then

              if(ebe.eq.0) then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    sbufx(j,k,nb+ib)=var(iesnd,j,k)
                  end do

!$omp end do

                end do

              end if

            end if

          end if

! -----

!$omp end parallel

call profile_stop(prof_id1, loop_len)

        end if

      end if

!! -----

      end subroutine s_putbufgx

!-----7--------------------------------------------------------------7--

      end module m_putbufgx

!***********************************************************************
      module m_atmosin
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/20
!     Modification: 2010/12/17

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read standard atomospheric profile for mstranx radiation scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commpi
      use m_commstrn
      use m_defmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: atmosin, s_atmosin

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface atmosin

        module procedure s_atmosin

      end interface

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_atmosin(iud,nln,cpcl_data,gdcfrc_data,cgas_data,ccfc)
!***********************************************************************

! Input variable

      integer, intent(in) :: iud
                       ! Unit number

! Output variables

      integer, intent(out) :: nln
                       ! Maximum number of layers

      real(kind=r8), intent(out) :: cpcl_data(1:kln,1:kpcl,1:kpclc)
                       ! Parameter packet for particulates in sublayers

      real(kind=r8), intent(out) :: gdcfrc_data(1:kln)
                       ! Cloud cover late in sublayers

      real(kind=r8), intent(out) :: cgas_data(1:kln,1:kmol)
                       ! Gas concentrations in sublayers

      real(kind=r8), intent(out) :: ccfc(1:kcfc)
                       ! CFCs concentrations

! Internal shared variables

      integer ilev     ! Do loop index for layers
      integer ipcl     ! Do loop index for parameter packets

      integer ierr     ! Temporary variable

      real(kind=r8) tmp
                       ! Temporary variable

!-----7--------------------------------------------------------------7--

! Read standard atomospheric profile for mstranx radiation scheme.

      if(mype.eq.root) then

        read(iud,*) nln,tmp

        do ilev=1,nln

          read(iud,'(1x)')
          read(iud,'(1x)')

          do ipcl=1,kpcl

            read(iud,*) cpcl_data(ilev,ipcl,1:kpclc)

          end do

          read(iud,'(1x)')

          read(iud,*) gdcfrc_data(ilev)

          read(iud,'(1x)')

          read(iud,*) cgas_data(ilev,1:kmol)

        end do

        read(iud,'(1x)')
        read(iud,'(1x)')
        read(iud,'(1x)')

      end if

      ccfc(1:kcfc)=1.d-10

! -----

! Broadcast reading data.

      call mpi_bcast(nln,1,mpi_integer,root,mpi_comm_cress,ierr)

      call mpi_bcast(cpcl_data,                                         &
     &               kln*kpcl*kpclc,mpi_real8,root,mpi_comm_cress,ierr)

      call mpi_bcast(gdcfrc_data,                                       &
     &               kln,mpi_real8,root,mpi_comm_cress,ierr)

      call mpi_bcast(cgas_data,                                         &
     &               kln*kmol,mpi_real8,root,mpi_comm_cress,ierr)

! -----

      end subroutine s_atmosin

!-----7--------------------------------------------------------------7--

      end module m_atmosin

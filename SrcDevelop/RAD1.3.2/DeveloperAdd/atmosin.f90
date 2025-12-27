!***********************************************************************
      module m_atmosin
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2010/10/20
!     Modification: 2010/12/17

!     Author      : Hasegawa Koichi
!     Modification: 2014/07/01

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read standard atomospheric profile for radiation scheme.

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
      subroutine s_atmosin(iud,nln,zalt_data,patm_data,tklv_data,       &
     &                     cpcl_data,cgas_data,ccfc)
!***********************************************************************

! Input variable

      integer, intent(in) :: iud
                       ! Unit number

      integer, intent(in) :: nln
                       ! Layers number of external data

! Output variables

      real(kind=r8), intent(out) :: zalt_data(1:nln)
                       ! external data (altitude[km])

      real(kind=r8), intent(out) :: patm_data(1:nln)
                       ! external data (pressure[atm])

      real(kind=r8), intent(out) :: tklv_data(1:nln)
                       ! external data (temperature[K])

      real(kind=r8), intent(out) :: cpcl_data(1:nln,1:kpcl,1:kpclc)
                       ! external data (particulates concentrations)

      real(kind=r8), intent(out) :: cgas_data(1:nln,1:kmol)
                       ! external data (gas concentrations)

      real(kind=r8), intent(out) :: ccfc(1:kcfc)
                       ! CFCs concentrations

! Internal private variables

      integer ierr     ! Temporary variable
      integer ilev     ! Do loop index for layers

!-----7--------------------------------------------------------------7--

! Read standard atomospheric profile for radiation scheme.

      if(mype.eq.root) then

        read(iud,*)

        do ilev=nln,1,-1

          read(iud,*) zalt_data(ilev),patm_data(ilev),tklv_data(ilev),  &
     &                cgas_data(ilev,1:kmol),cpcl_data(ilev,1:kpcl,1),  &
     &                cpcl_data(ilev,1:kpcl,2)

        end do

      end if

      ccfc(1:kcfc)=1.d-10

! -----

! Broadcast reading data.

      call mpi_bcast(zalt_data,nln,mpi_real8,root,mpi_comm_cress,ierr)
      call mpi_bcast(patm_data,nln,mpi_real8,root,mpi_comm_cress,ierr)
      call mpi_bcast(tklv_data,nln,mpi_real8,root,mpi_comm_cress,ierr)

      call mpi_bcast(cpcl_data,                                         &
     &               nln*kpcl*kpclc,mpi_real8,root,mpi_comm_cress,ierr)

      call mpi_bcast(cgas_data,                                         &
     &               nln*kmol,mpi_real8,root,mpi_comm_cress,ierr)

! -----

      end subroutine s_atmosin

!-----7--------------------------------------------------------------7--

      end module m_atmosin

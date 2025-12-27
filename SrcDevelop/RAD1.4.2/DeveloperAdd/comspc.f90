!***********************************************************************
      module m_comspc
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2016/12/23, 2017/02/01, 2017/06/09

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the array for spectral nudging for ISPACK library.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      public

! Exceptional access control

!     none

!-----7--------------------------------------------------------------7--

! Module variables

      logical, save :: spn_mul_flag
                       ! If nk-1-nsl < numpe, flag to classify whether 
                       ! multiple processing or not in each horizontal plane.
                       ! .true. = multiple, .false. = single.

      integer, save, dimension(5) :: LITI
                       ! x component of numbers of prime factors (2,3,5,7)

      integer, save, dimension(5) :: LITJ
                       ! y component of numbers of prime factors (2,3,5,7)

      integer, allocatable, dimension(:), save :: spnk2pe
                       ! Array to classify SPN processor number from k

      logical, allocatable, dimension(:), save :: spnpe2k
                       ! Array to classify k from SPN processor number

      real, allocatable, dimension(:,:,:), save :: spnufrc
                       ! Forcing term in u equation due to spectral nudging

      real, allocatable, dimension(:,:,:), save :: spnvfrc
                       ! Forcing term in v equation due to spectral nudging

      real, allocatable, dimension(:,:,:), save :: spnsfrc
                       ! Forcing term in scalar equation due to spectral nudging

      double precision, allocatable, dimension(:), save :: LTI
                       ! working array for x direction of FFT

      double precision, allocatable, dimension(:), save :: LTJ
                       ! working array for y direction of FFT

      double precision, allocatable, dimension(:,:), save :: WGL
                       ! working array of FFT

! Module procedure

!     none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

!     none

!-----7--------------------------------------------------------------7--

      end module m_comspc

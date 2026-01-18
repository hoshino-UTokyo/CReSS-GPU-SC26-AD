!***********************************************************************
      module m_dbz2kg
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/11/20, 2003/04/30, 2003/05/19, 2004/04/15,
!                   2004/05/07, 2004/09/10, 2007/10/19, 2008/04/17,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     convert the hydrometeor mesuremement from [dBZe] to [kg/m^3]

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comprofile
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: dbz2kg, s_dbz2kg

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface dbz2kg

        module procedure s_dbz2kg

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
      subroutine s_dbz2kg(fprdrcoe_rdr,fprdrexp_rdr,nid,njd,nkd,qpdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fprdrcoe_rdr
                       ! Formal parameter of unique index of rdrcoe_rdr

      integer, intent(in) :: fprdrexp_rdr
                       ! Formal parameter of unique index of rdrexp_rdr

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

! Input and output variable

      real, intent(inout) :: qpdat(1:nid,1:njd,1:nkd)
                       ! Precipitation mixing ratio in data

! Internal shared variables

      real rdrcoe_rdr  ! Coeffient of converter of [dBZe] to [kg/m^3]
      real rdrexp_rdr  ! Exponent of converter of [dBZe] to [kg/m^3]

! Internal private variables

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction
      integer kd       ! Array index in z direction


      ! Profiling variables
      integer, save :: prof_id1 = -1
      integer(8) :: loop_len

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fprdrcoe_rdr,rdrcoe_rdr)
      call getrname(fprdrexp_rdr,rdrexp_rdr)

! -----

! Convert the hydrometeor mesuremement from [dBZe] to [kg/m^3]

!@llm start meta_info ----------------------------------------------------
! Location: dbz2kg.f90 :: s_dbz2kg
! Summary : Convert radar reflectivity from dBZe to precipitation mixing
!           ratio in kg/m^3 using exponential/logarithmic transformation.
! GPU diff: Easy
! Findings:
!   - No omp_get_thread_num usage
!   - Uses intrinsic exp and log functions (GPU-compatible)
!   - Simple 3D loop with independent point-wise operations
!   - Writes to qpdat array (in-place modification)
!   - Conditional check on lim34n threshold
! Next:
!   - Direct conversion to OpenACC or OpenACC with collapsed loops
!   - Data managed automatically via Unified Memory
!   - No synchronization needed between iterations
!@llm end meta_info ------------------------------------------------------

! Register profiling section (first call only)
if (prof_id1 < 0) then
  prof_id1 = profile_register('dbz2kg.f90', 's_dbz2kg', &
   & 'OMP section 1')
end if
loop_len = int((nkd)-(1)+1,8) * int((njd)-(1)+1,8) * int((nid)-(1)+1,8)
call profile_start(prof_id1)

!$omp parallel default(shared) private(kd)

      do kd=1,nkd

!$omp do schedule(runtime) private(id,jd)

        do jd=1,njd
        do id=1,nid

          if(qpdat(id,jd,kd).gt.lim34n) then

            qpdat(id,jd,kd)=rdrcoe_rdr                                  &
     &        *exp(rdrexp_rdr*log(qpdat(id,jd,kd)))

          end if

        end do
        end do

!$omp end do

      end do

!$omp end parallel

call profile_stop(prof_id1, loop_len)

! -----

      end subroutine s_dbz2kg

!-----7--------------------------------------------------------------7--

      end module m_dbz2kg

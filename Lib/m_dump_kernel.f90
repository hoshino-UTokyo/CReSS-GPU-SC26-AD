!***********************************************************************
! Backward Compatibility Module for m_dump_kernel
!***********************************************************************
!
! This module provides backward compatibility with existing code that
! uses the m_dump_kernel module. All functionality is provided by the
! unified m_kernel_io module.
!
! Usage (same as before):
!   use m_dump_kernel
!   call dump_init('kernel_name')
!   call dump_scalar_i('ni', ni)
!   call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
!   call dump_finalize()
!
!***********************************************************************
module m_dump_kernel
  use m_kernel_io
  implicit none
  private

  ! Re-export dump routines with original names
  public :: dump_init, dump_finalize
  public :: dump_scalar_i, dump_scalar_i8, dump_scalar_r, dump_scalar_d, dump_scalar_c
  public :: dump_scalar_s  ! Alias for dump_scalar_c
  public :: dump_array_1d, dump_array_2d, dump_array_3d, dump_array_4d
  public :: dump_array_1d_int, dump_array_2d_int, dump_array_3d_int
  public :: dump_array_1d_i, dump_array_2d_i, dump_array_3d_i  ! Aliases

  ! Aliases
  interface dump_scalar_s
    module procedure dump_scalar_c_wrapper
  end interface

  interface dump_array_1d_i
    module procedure dump_array_1d_int_wrapper
  end interface

  interface dump_array_2d_i
    module procedure dump_array_2d_int_wrapper
  end interface

  interface dump_array_3d_i
    module procedure dump_array_3d_int_wrapper
  end interface

contains

  !=====================================================================
  ! Initialization and finalization
  !=====================================================================
  subroutine dump_init(kernel_name)
    character(len=*), intent(in) :: kernel_name
    call kio_dump_init(kernel_name)
  end subroutine dump_init

  subroutine dump_finalize()
    call kio_dump_finalize()
  end subroutine dump_finalize

  !=====================================================================
  ! Scalar dump routines
  !=====================================================================
  subroutine dump_scalar_i(name, val)
    character(len=*), intent(in) :: name
    integer, intent(in) :: val
    call kio_dump_scalar_i(name, val)
  end subroutine dump_scalar_i

  subroutine dump_scalar_i8(name, val)
    character(len=*), intent(in) :: name
    integer(8), intent(in) :: val
    call kio_dump_scalar_i8(name, val)
  end subroutine dump_scalar_i8

  subroutine dump_scalar_r(name, val)
    character(len=*), intent(in) :: name
    real, intent(in) :: val
    call kio_dump_scalar_r(name, val)
  end subroutine dump_scalar_r

  subroutine dump_scalar_d(name, val)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: val
    call kio_dump_scalar_d(name, val)
  end subroutine dump_scalar_d

  subroutine dump_scalar_c(name, val)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    call kio_dump_scalar_c(name, val)
  end subroutine dump_scalar_c

  subroutine dump_scalar_c_wrapper(name, val)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    call kio_dump_scalar_c(name, val)
  end subroutine dump_scalar_c_wrapper

  !=====================================================================
  ! Real array dump routines
  !=====================================================================
  subroutine dump_array_1d(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    real, intent(in) :: arr(i1:i2)
    call kio_dump_array_1d(filename, arr, i1, i2)
  end subroutine dump_array_1d

  subroutine dump_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(in) :: arr(i1:i2, j1:j2)
    call kio_dump_array_2d(filename, arr, i1, i2, j1, j2)
  end subroutine dump_array_2d

  subroutine dump_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2)
    call kio_dump_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
  end subroutine dump_array_3d

  subroutine dump_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)
    call kio_dump_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
  end subroutine dump_array_4d

  !=====================================================================
  ! Integer array dump routines
  !=====================================================================
  subroutine dump_array_1d_int(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    integer, intent(in) :: arr(i1:i2)
    call kio_dump_array_1d_int(filename, arr, i1, i2)
  end subroutine dump_array_1d_int

  subroutine dump_array_2d_int(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    integer, intent(in) :: arr(i1:i2, j1:j2)
    call kio_dump_array_2d_int(filename, arr, i1, i2, j1, j2)
  end subroutine dump_array_2d_int

  subroutine dump_array_3d_int(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    integer, intent(in) :: arr(i1:i2, j1:j2, k1:k2)
    call kio_dump_array_3d_int(filename, arr, i1, i2, j1, j2, k1, k2)
  end subroutine dump_array_3d_int

  subroutine dump_array_1d_int_wrapper(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    integer, intent(in) :: arr(i1:i2)
    call kio_dump_array_1d_int(filename, arr, i1, i2)
  end subroutine dump_array_1d_int_wrapper

  subroutine dump_array_2d_int_wrapper(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    integer, intent(in) :: arr(i1:i2, j1:j2)
    call kio_dump_array_2d_int(filename, arr, i1, i2, j1, j2)
  end subroutine dump_array_2d_int_wrapper

  subroutine dump_array_3d_int_wrapper(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    integer, intent(in) :: arr(i1:i2, j1:j2, k1:k2)
    call kio_dump_array_3d_int(filename, arr, i1, i2, j1, j2, k1, k2)
  end subroutine dump_array_3d_int_wrapper

end module m_dump_kernel

!***********************************************************************
! Module for dumping kernel input/output data from running simulation
!***********************************************************************
!
! Usage:
!   1. Add "use m_dump_kernel" to the target source file
!   2. Insert dump calls before/after the OpenMP parallel section
!   3. Run simulation once to generate data files
!   4. Use generated files for standalone benchmark
!
! Example:
!   use m_dump_kernel
!   ...
!   integer, save :: dump_call_count = 0
!   integer, parameter :: DUMP_TARGET_CALL = 14400  ! From profiling
!   logical, save :: dump_done = .false.
!   ...
!   dump_call_count = dump_call_count + 1
!   if (dump_call_count == DUMP_TARGET_CALL .and. .not. dump_done) then
!     call dump_init('kernel_name')
!     call dump_scalar_i('ni', ni)
!     call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
!     ! ... more dumps
!   end if
!   !$omp parallel ...
!   !$omp end parallel
!   if (dump_call_count == DUMP_TARGET_CALL .and. .not. dump_done) then
!     call dump_array_3d('output.bin', output, ...)
!     call dump_finalize()
!     dump_done = .true.
!   end if
!
!***********************************************************************
module m_dump_kernel
  implicit none
  private

  ! Public subroutines
  public :: dump_init, dump_finalize
  public :: dump_scalar_i, dump_scalar_i8, dump_scalar_r, dump_scalar_d, dump_scalar_c
  public :: dump_scalar_s  ! Alias for dump_scalar_c (string)
  public :: dump_array_1d, dump_array_2d, dump_array_3d, dump_array_4d
  public :: dump_array_1d_int, dump_array_2d_int, dump_array_3d_int
  public :: dump_array_1d_i, dump_array_2d_i, dump_array_3d_i  ! Aliases

  ! Module variables
  character(len=256), save :: dump_dir = './kernel_dump'
  integer, save :: dump_unit = 100
  integer, save :: param_unit = 101
  logical, save :: param_file_open = .false.

  ! Interfaces for aliases
  interface dump_scalar_s
    module procedure dump_scalar_c
  end interface

  interface dump_array_1d_i
    module procedure dump_array_1d_int
  end interface

  interface dump_array_2d_i
    module procedure dump_array_2d_int
  end interface

  interface dump_array_3d_i
    module procedure dump_array_3d_int
  end interface

contains

  !=====================================================================
  ! Initialize dump directory and open params file
  !=====================================================================
  subroutine dump_init(kernel_name)
    character(len=*), intent(in) :: kernel_name
    integer :: ios

    dump_dir = './kernel_dump/'//trim(kernel_name)
    call execute_command_line('mkdir -p '//trim(dump_dir), wait=.true.)

    ! Open params file for appending scalar values
    open(unit=param_unit, file=trim(dump_dir)//'/params.txt', &
         status='replace', iostat=ios)
    if (ios == 0) then
      param_file_open = .true.
    end if

    write(*,'(A)') '[DUMP] ========================================'
    write(*,'(A)') '[DUMP] Kernel: '//trim(kernel_name)
    write(*,'(A)') '[DUMP] Output: '//trim(dump_dir)
    write(*,'(A)') '[DUMP] ========================================'

  end subroutine dump_init

  !=====================================================================
  ! Finalize (close files and print summary)
  !=====================================================================
  subroutine dump_finalize()
    if (param_file_open) then
      close(param_unit)
      param_file_open = .false.
    end if
    write(*,'(A)') '[DUMP] ========================================'
    write(*,'(A)') '[DUMP] Complete. Data saved to: '//trim(dump_dir)
    write(*,'(A)') '[DUMP] ========================================'
  end subroutine dump_finalize

  !=====================================================================
  ! Dump scalar integer
  !=====================================================================
  subroutine dump_scalar_i(name, val)
    character(len=*), intent(in) :: name
    integer, intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,I12)') trim(name), ' = ', val
    end if
    write(*,'(A,A,I12)') '[DUMP] ', trim(name)//' = ', val

  end subroutine dump_scalar_i

  !=====================================================================
  ! Dump scalar integer (8-byte / 64-bit)
  !=====================================================================
  subroutine dump_scalar_i8(name, val)
    character(len=*), intent(in) :: name
    integer(8), intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,I20)') trim(name), ' = ', val
    end if
    write(*,'(A,A,I20)') '[DUMP] ', trim(name)//' = ', val

  end subroutine dump_scalar_i8

  !=====================================================================
  ! Dump scalar real (single precision)
  !=====================================================================
  subroutine dump_scalar_r(name, val)
    character(len=*), intent(in) :: name
    real, intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,ES20.12)') trim(name), ' = ', val
    end if
    write(*,'(A,A,ES15.7)') '[DUMP] ', trim(name)//' = ', val

  end subroutine dump_scalar_r

  !=====================================================================
  ! Dump scalar real (double precision)
  !=====================================================================
  subroutine dump_scalar_d(name, val)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,ES24.16)') trim(name), ' = ', val
    end if
    write(*,'(A,A,ES20.12)') '[DUMP] ', trim(name)//' = ', val

  end subroutine dump_scalar_d

  !=====================================================================
  ! Dump scalar character string
  !=====================================================================
  subroutine dump_scalar_c(name, val)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,A)') trim(name), ' = ', trim(val)
    end if
    write(*,'(A,A,A)') '[DUMP] ', trim(name)//' = ', trim(val)

  end subroutine dump_scalar_c

  !=====================================================================
  ! Dump 1D real array to binary file
  !=====================================================================
  subroutine dump_array_1d(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    real, intent(in) :: arr(i1:i2)

    character(len=512) :: filepath
    integer :: nelements

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = i2 - i1 + 1
    write(*,'(A,I10,A)') '[DUMP] 1D: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_1d

  !=====================================================================
  ! Dump 2D real array to binary file
  !=====================================================================
  subroutine dump_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(in) :: arr(i1:i2, j1:j2)

    character(len=512) :: filepath
    integer(8) :: nelements

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = int(i2-i1+1,8) * int(j2-j1+1,8)
    write(*,'(A,I12,A)') '[DUMP] 2D: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_2d

  !=====================================================================
  ! Dump 3D real array to binary file
  !=====================================================================
  subroutine dump_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2)

    character(len=512) :: filepath
    integer(8) :: nelements

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = int(i2-i1+1,8) * int(j2-j1+1,8) * int(k2-k1+1,8)
    write(*,'(A,I12,A)') '[DUMP] 3D: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_3d

  !=====================================================================
  ! Dump 4D real array to binary file
  !=====================================================================
  subroutine dump_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)

    character(len=512) :: filepath
    integer(8) :: nelements

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = int(i2-i1+1,8) * int(j2-j1+1,8) * int(k2-k1+1,8) * int(l2-l1+1,8)
    write(*,'(A,I14,A)') '[DUMP] 4D: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_4d

  !=====================================================================
  ! Dump 1D integer array to binary file
  !=====================================================================
  subroutine dump_array_1d_int(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    integer, intent(in) :: arr(i1:i2)

    character(len=512) :: filepath
    integer :: nelements

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = i2 - i1 + 1
    write(*,'(A,I10,A)') '[DUMP] 1D int: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_1d_int

  !=====================================================================
  ! Dump 2D integer array to binary file
  !=====================================================================
  subroutine dump_array_2d_int(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    integer, intent(in) :: arr(i1:i2, j1:j2)

    character(len=512) :: filepath
    integer(8) :: nelements

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = int(i2-i1+1,8) * int(j2-j1+1,8)
    write(*,'(A,I12,A)') '[DUMP] 2D int: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_2d_int

  !=====================================================================
  ! Dump 3D integer array to binary file
  !=====================================================================
  subroutine dump_array_3d_int(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    integer, intent(in) :: arr(i1:i2, j1:j2, k1:k2)

    character(len=512) :: filepath
    integer(8) :: nelements

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = int(i2-i1+1,8) * int(j2-j1+1,8) * int(k2-k1+1,8)
    write(*,'(A,I12,A)') '[DUMP] 3D int: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_3d_int

end module m_dump_kernel

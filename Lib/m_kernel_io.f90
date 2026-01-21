!***********************************************************************
! Unified Kernel I/O Library Module
!***********************************************************************
!
! This module provides unified routines for:
!   - Dumping kernel data to binary files (for data generation)
!   - Loading kernel data from binary files (for benchmarks)
!
! IMPORTANT: Fortran runtime has an internal I/O buffer limit of ~1GB.
! This module checks array sizes before I/O operations and raises errors
! if the limit would be exceeded. For 4D arrays, data is written/read
! in 2D slices to avoid this limitation.
!
! Usage Example (Dumping):
!   use m_kernel_io
!   call kio_dump_init('kernel_name')
!   call kio_dump_scalar_i('ni', ni)
!   call kio_dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
!   call kio_dump_finalize()
!
! Usage Example (Loading):
!   use m_kernel_io
!   call kio_load_params('data/params.txt')
!   call kio_get_param_i('ni', ni)
!   call kio_load_array_3d('data/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
!
!***********************************************************************
module m_kernel_io
  implicit none
  private

  !=====================================================================
  ! Module constants and configurable limits
  !=====================================================================
  ! Default buffer limit: 1GB (conservative estimate)
  ! This can be overridden by calling kio_set_buffer_limit() or
  ! setting environment variable KIO_BUFFER_LIMIT_MB
  integer(8), parameter :: KIO_DEFAULT_BUFFER_BYTES = 1073741824_8  ! 1GB

  ! Current buffer limit (can be modified at runtime)
  integer(8), save, public :: KIO_BUFFER_LIMIT_BYTES = 1073741824_8
  integer(8), save, public :: KIO_MAX_ELEMENTS_4BYTE = 268435456_8   ! 1GB / 4
  integer(8), save, public :: KIO_MAX_ELEMENTS_8BYTE = 134217728_8   ! 1GB / 8

  ! Flag to track if limits have been initialized from environment
  logical, save :: limits_initialized = .false.

  !=====================================================================
  ! Public interfaces - Dump routines
  !=====================================================================
  public :: kio_dump_init, kio_dump_finalize
  public :: kio_dump_scalar_i, kio_dump_scalar_i8
  public :: kio_dump_scalar_r, kio_dump_scalar_d, kio_dump_scalar_c
  public :: kio_dump_array_1d, kio_dump_array_2d, kio_dump_array_3d, kio_dump_array_4d
  public :: kio_dump_array_1d_int, kio_dump_array_2d_int, kio_dump_array_3d_int

  !=====================================================================
  ! Public interfaces - Load routines
  !=====================================================================
  public :: kio_load_params, kio_params_loaded
  public :: kio_get_param_i, kio_get_param_r, kio_get_param_d, kio_get_param_c
  public :: kio_load_array_1d, kio_load_array_2d, kio_load_array_3d, kio_load_array_4d
  public :: kio_load_array_1d_int, kio_load_array_2d_int, kio_load_array_3d_int

  !=====================================================================
  ! Public interfaces - Utility routines
  !=====================================================================
  public :: kio_check_size, kio_check_size_bytes
  public :: kio_error
  public :: kio_set_buffer_limit, kio_get_buffer_limit
  public :: kio_init_from_env

  !=====================================================================
  ! Module variables - Dump state
  !=====================================================================
  character(len=512), save :: dump_dir = './kernel_dump'
  integer, save :: dump_unit = 100
  integer, save :: param_unit = 101
  logical, save :: param_file_open = .false.

  !=====================================================================
  ! Module variables - Load state (parameter storage)
  !=====================================================================
  integer, parameter :: MAX_PARAMS = 256
  integer, parameter :: MAX_KEY_LEN = 64
  integer, parameter :: MAX_VAL_LEN = 256

  type :: param_entry
    character(len=MAX_KEY_LEN) :: key = ''
    character(len=MAX_VAL_LEN) :: value = ''
  end type param_entry

  type(param_entry), save :: params(MAX_PARAMS)
  integer, save :: num_params = 0
  logical, save :: params_loaded = .false.

contains

  !=====================================================================
  !=====================================================================
  !  ERROR HANDLING
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Raise fatal error with message
  !---------------------------------------------------------------------
  subroutine kio_error(msg)
    character(len=*), intent(in) :: msg

    write(*,'(A)') ''
    write(*,'(A)') '***** KERNEL I/O ERROR *****'
    write(*,'(A)') trim(msg)
    write(*,'(A)') '****************************'
    write(*,'(A)') ''
    error stop 1

  end subroutine kio_error

  !=====================================================================
  !=====================================================================
  !  SIZE CHECKING ROUTINES
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Check if number of elements exceeds 1GB buffer limit for 4-byte data
  ! Returns .true. if size is OK, .false. if exceeds limit
  !---------------------------------------------------------------------
  function kio_check_size(nelements) result(ok)
    integer(8), intent(in) :: nelements
    logical :: ok

    ok = (nelements <= KIO_MAX_ELEMENTS_4BYTE)

  end function kio_check_size

  !---------------------------------------------------------------------
  ! Check if byte size exceeds 1GB buffer limit
  ! Returns .true. if size is OK, .false. if exceeds limit
  !---------------------------------------------------------------------
  function kio_check_size_bytes(nbytes) result(ok)
    integer(8), intent(in) :: nbytes
    logical :: ok

    ok = (nbytes <= KIO_BUFFER_LIMIT_BYTES)

  end function kio_check_size_bytes

  !---------------------------------------------------------------------
  ! Internal: Check size and raise error if exceeded
  !---------------------------------------------------------------------
  subroutine check_size_or_error(nelements, element_bytes, array_name, dims_info)
    integer(8), intent(in) :: nelements
    integer, intent(in) :: element_bytes
    character(len=*), intent(in) :: array_name
    character(len=*), intent(in) :: dims_info

    integer(8) :: total_bytes, max_elements
    character(len=1024) :: errmsg

    total_bytes = nelements * int(element_bytes, 8)
    max_elements = KIO_BUFFER_LIMIT_BYTES / int(element_bytes, 8)

    if (total_bytes > KIO_BUFFER_LIMIT_BYTES) then
      write(errmsg, '(A,A,A,I0,A,I0,A,A,A,I0,A)') &
        'Array size exceeds 1GB I/O buffer limit!', char(10), &
        '  Array: ', trim(array_name), char(10), &
        '  Dimensions: ', trim(dims_info), char(10), &
        '  Elements: ', nelements, ' (max: ', max_elements, ')', char(10), &
        '  Bytes: ', total_bytes, ' (max: ', KIO_BUFFER_LIMIT_BYTES, ')'
      call kio_error(trim(errmsg))
    end if

  end subroutine check_size_or_error

  !=====================================================================
  !=====================================================================
  !  BUFFER LIMIT CONFIGURATION
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Set buffer limit in megabytes
  ! Call this before any I/O operations to customize the limit
  !---------------------------------------------------------------------
  subroutine kio_set_buffer_limit(limit_mb)
    integer, intent(in) :: limit_mb

    if (limit_mb <= 0) then
      call kio_error('Buffer limit must be positive')
    end if

    KIO_BUFFER_LIMIT_BYTES = int(limit_mb, 8) * 1048576_8  ! MB to bytes
    KIO_MAX_ELEMENTS_4BYTE = KIO_BUFFER_LIMIT_BYTES / 4_8
    KIO_MAX_ELEMENTS_8BYTE = KIO_BUFFER_LIMIT_BYTES / 8_8
    limits_initialized = .true.

    write(*,'(A,I0,A)') '[KIO] Buffer limit set to ', limit_mb, ' MB'

  end subroutine kio_set_buffer_limit

  !---------------------------------------------------------------------
  ! Get current buffer limit in bytes
  !---------------------------------------------------------------------
  function kio_get_buffer_limit() result(limit_bytes)
    integer(8) :: limit_bytes
    call ensure_limits_initialized()
    limit_bytes = KIO_BUFFER_LIMIT_BYTES
  end function kio_get_buffer_limit

  !---------------------------------------------------------------------
  ! Initialize buffer limit from environment variable KIO_BUFFER_LIMIT_MB
  ! Returns .true. if environment variable was found and applied
  !---------------------------------------------------------------------
  function kio_init_from_env() result(found)
    logical :: found
    character(len=32) :: env_value
    integer :: limit_mb, ios

    found = .false.
    call get_environment_variable('KIO_BUFFER_LIMIT_MB', env_value, status=ios)

    if (ios == 0 .and. len_trim(env_value) > 0) then
      read(env_value, *, iostat=ios) limit_mb
      if (ios == 0 .and. limit_mb > 0) then
        call kio_set_buffer_limit(limit_mb)
        found = .true.
      end if
    end if

    limits_initialized = .true.

  end function kio_init_from_env

  !---------------------------------------------------------------------
  ! Internal: Ensure limits are initialized (check env on first use)
  !---------------------------------------------------------------------
  subroutine ensure_limits_initialized()
    logical :: dummy

    if (.not. limits_initialized) then
      dummy = kio_init_from_env()
      ! If env not set, keep defaults (already initialized)
      limits_initialized = .true.
    end if

  end subroutine ensure_limits_initialized

  !=====================================================================
  !=====================================================================
  !  DUMP ROUTINES - INITIALIZATION
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Initialize dump directory and open params file
  !---------------------------------------------------------------------
  subroutine kio_dump_init(kernel_name)
    character(len=*), intent(in) :: kernel_name
    integer :: ios

    ! Ensure buffer limits are initialized (check environment variable)
    call ensure_limits_initialized()

    dump_dir = './kernel_dump/'//trim(kernel_name)
    call execute_command_line('mkdir -p '//trim(dump_dir), wait=.true.)

    open(unit=param_unit, file=trim(dump_dir)//'/params.txt', &
         status='replace', iostat=ios)
    if (ios == 0) then
      param_file_open = .true.
    end if

    write(*,'(A)') '[KIO] ========================================'
    write(*,'(A)') '[KIO] Kernel: '//trim(kernel_name)
    write(*,'(A)') '[KIO] Output: '//trim(dump_dir)
    write(*,'(A)') '[KIO] ========================================'

  end subroutine kio_dump_init

  !---------------------------------------------------------------------
  ! Finalize dump (close files and print summary)
  !---------------------------------------------------------------------
  subroutine kio_dump_finalize()
    if (param_file_open) then
      close(param_unit)
      param_file_open = .false.
    end if
    write(*,'(A)') '[KIO] ========================================'
    write(*,'(A)') '[KIO] Dump complete. Data saved to: '//trim(dump_dir)
    write(*,'(A)') '[KIO] ========================================'
  end subroutine kio_dump_finalize

  !=====================================================================
  !=====================================================================
  !  DUMP ROUTINES - SCALARS
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Dump scalar integer (4-byte)
  !---------------------------------------------------------------------
  subroutine kio_dump_scalar_i(name, val)
    character(len=*), intent(in) :: name
    integer, intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,I12)') trim(name), ' = ', val
      flush(param_unit)
    end if
    write(*,'(A,A,I12)') '[KIO] ', trim(name)//' = ', val

  end subroutine kio_dump_scalar_i

  !---------------------------------------------------------------------
  ! Dump scalar integer (8-byte)
  !---------------------------------------------------------------------
  subroutine kio_dump_scalar_i8(name, val)
    character(len=*), intent(in) :: name
    integer(8), intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,I20)') trim(name), ' = ', val
      flush(param_unit)
    end if
    write(*,'(A,A,I20)') '[KIO] ', trim(name)//' = ', val

  end subroutine kio_dump_scalar_i8

  !---------------------------------------------------------------------
  ! Dump scalar real (4-byte single precision)
  !---------------------------------------------------------------------
  subroutine kio_dump_scalar_r(name, val)
    character(len=*), intent(in) :: name
    real, intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,ES20.12)') trim(name), ' = ', val
      flush(param_unit)
    end if
    write(*,'(A,A,ES15.7)') '[KIO] ', trim(name)//' = ', val

  end subroutine kio_dump_scalar_r

  !---------------------------------------------------------------------
  ! Dump scalar real (8-byte double precision)
  !---------------------------------------------------------------------
  subroutine kio_dump_scalar_d(name, val)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,ES24.16)') trim(name), ' = ', val
      flush(param_unit)
    end if
    write(*,'(A,A,ES20.12)') '[KIO] ', trim(name)//' = ', val

  end subroutine kio_dump_scalar_d

  !---------------------------------------------------------------------
  ! Dump scalar character string
  !---------------------------------------------------------------------
  subroutine kio_dump_scalar_c(name, val)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val

    if (param_file_open) then
      write(param_unit, '(A,A,A)') trim(name), ' = ', trim(val)
      flush(param_unit)
    end if
    write(*,'(A,A,A)') '[KIO] ', trim(name)//' = ', trim(val)

  end subroutine kio_dump_scalar_c

  !=====================================================================
  !=====================================================================
  !  DUMP ROUTINES - REAL ARRAYS
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Dump 1D real array to binary file
  ! Auto-splits into chunks if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_dump_array_1d(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    real, intent(in) :: arr(i1:i2)

    character(len=512) :: filepath
    integer(8) :: nelements, chunk_size, chunk_start, chunk_end, remaining
    integer :: ichunk

    nelements = int(i2-i1+1, 8)

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      ! Direct write if within limit
      write(dump_unit) arr
      write(*,'(A,I12,A)') '[KIO] 1D: '//trim(filename)//' (', nelements, ')'
    else
      ! Chunked write
      chunk_size = KIO_MAX_ELEMENTS_4BYTE
      ichunk = 0
      chunk_start = int(i1, 8)
      remaining = nelements

      do while (remaining > 0)
        chunk_end = min(chunk_start + chunk_size - 1, int(i2, 8))
        write(dump_unit) arr(int(chunk_start):int(chunk_end))
        ichunk = ichunk + 1
        remaining = remaining - (chunk_end - chunk_start + 1)
        chunk_start = chunk_end + 1
      end do

      write(*,'(A,I12,A,I0,A)') '[KIO] 1D: '//trim(filename)//' (', nelements, ') [', ichunk, ' chunks]'
    end if

    close(dump_unit)

  end subroutine kio_dump_array_1d

  !---------------------------------------------------------------------
  ! Dump 2D real array to binary file
  ! Auto-splits into 1D slices if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_dump_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(in) :: arr(i1:i2, j1:j2)

    character(len=512) :: filepath
    integer(8) :: nelements, slice_elements
    integer :: j

    nelements = int(i2-i1+1, 8) * int(j2-j1+1, 8)
    slice_elements = int(i2-i1+1, 8)

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      ! Direct write if within limit
      write(dump_unit) arr
      write(*,'(A,I12,A)') '[KIO] 2D: '//trim(filename)//' (', nelements, ')'
    else
      ! Write in 1D slices (row by row)
      do j = j1, j2
        write(dump_unit) arr(i1:i2, j)
      end do
      write(*,'(A,I12,A)') '[KIO] 2D: '//trim(filename)//' (', nelements, ') [1D slices]'
    end if

    close(dump_unit)

  end subroutine kio_dump_array_2d

  !---------------------------------------------------------------------
  ! Dump 3D real array to binary file
  ! Auto-splits into 2D slices if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_dump_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2)

    character(len=512) :: filepath
    integer(8) :: nelements, slice_elements
    integer :: k

    nelements = int(i2-i1+1, 8) * int(j2-j1+1, 8) * int(k2-k1+1, 8)
    slice_elements = int(i2-i1+1, 8) * int(j2-j1+1, 8)

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      ! Direct write if within limit
      write(dump_unit) arr
      write(*,'(A,I12,A)') '[KIO] 3D: '//trim(filename)//' (', nelements, ')'
    else
      ! Write in 2D slices (k-loop)
      ! Also check if single 2D slice exceeds limit
      if (slice_elements > KIO_MAX_ELEMENTS_4BYTE) then
        call kio_error('2D slice exceeds 1GB limit. Cannot auto-split further.')
      end if
      do k = k1, k2
        write(dump_unit) arr(i1:i2, j1:j2, k)
      end do
      write(*,'(A,I12,A)') '[KIO] 3D: '//trim(filename)//' (', nelements, ') [2D slices]'
    end if

    close(dump_unit)

  end subroutine kio_dump_array_3d

  !---------------------------------------------------------------------
  ! Dump 4D real array to binary file
  ! Always writes in 2D slices to avoid 1GB buffer overflow
  !---------------------------------------------------------------------
  subroutine kio_dump_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)

    character(len=512) :: filepath
    integer(8) :: nelements, slice_elements
    integer :: k, l

    nelements = int(i2-i1+1, 8) * int(j2-j1+1, 8) * int(k2-k1+1, 8) * int(l2-l1+1, 8)
    slice_elements = int(i2-i1+1, 8) * int(j2-j1+1, 8)

    ! Check 2D slice size
    if (slice_elements > KIO_MAX_ELEMENTS_4BYTE) then
      call kio_error('2D slice exceeds 1GB limit. Cannot auto-split further.')
    end if

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')

    ! Always write in 2D slices
    do l = l1, l2
      do k = k1, k2
        write(dump_unit) arr(i1:i2, j1:j2, k, l)
      end do
    end do

    close(dump_unit)

    write(*,'(A,I14,A)') '[KIO] 4D: '//trim(filename)//' (', nelements, ') [2D slices]'

  end subroutine kio_dump_array_4d

  !=====================================================================
  !=====================================================================
  !  DUMP ROUTINES - INTEGER ARRAYS
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Dump 1D integer array to binary file
  ! Auto-splits into chunks if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_dump_array_1d_int(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    integer, intent(in) :: arr(i1:i2)

    character(len=512) :: filepath
    integer(8) :: nelements, chunk_size, chunk_start, chunk_end, remaining
    integer :: ichunk

    nelements = int(i2-i1+1, 8)

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      write(dump_unit) arr
      write(*,'(A,I12,A)') '[KIO] 1D int: '//trim(filename)//' (', nelements, ')'
    else
      chunk_size = KIO_MAX_ELEMENTS_4BYTE
      ichunk = 0
      chunk_start = int(i1, 8)
      remaining = nelements

      do while (remaining > 0)
        chunk_end = min(chunk_start + chunk_size - 1, int(i2, 8))
        write(dump_unit) arr(int(chunk_start):int(chunk_end))
        ichunk = ichunk + 1
        remaining = remaining - (chunk_end - chunk_start + 1)
        chunk_start = chunk_end + 1
      end do

      write(*,'(A,I12,A,I0,A)') '[KIO] 1D int: '//trim(filename)//' (', nelements, ') [', ichunk, ' chunks]'
    end if

    close(dump_unit)

  end subroutine kio_dump_array_1d_int

  !---------------------------------------------------------------------
  ! Dump 2D integer array to binary file
  ! Auto-splits into 1D slices if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_dump_array_2d_int(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    integer, intent(in) :: arr(i1:i2, j1:j2)

    character(len=512) :: filepath
    integer(8) :: nelements
    integer :: j

    nelements = int(i2-i1+1, 8) * int(j2-j1+1, 8)

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      write(dump_unit) arr
      write(*,'(A,I12,A)') '[KIO] 2D int: '//trim(filename)//' (', nelements, ')'
    else
      do j = j1, j2
        write(dump_unit) arr(i1:i2, j)
      end do
      write(*,'(A,I12,A)') '[KIO] 2D int: '//trim(filename)//' (', nelements, ') [1D slices]'
    end if

    close(dump_unit)

  end subroutine kio_dump_array_2d_int

  !---------------------------------------------------------------------
  ! Dump 3D integer array to binary file
  ! Auto-splits into 2D slices if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_dump_array_3d_int(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    integer, intent(in) :: arr(i1:i2, j1:j2, k1:k2)

    character(len=512) :: filepath
    integer(8) :: nelements, slice_elements
    integer :: k

    nelements = int(i2-i1+1, 8) * int(j2-j1+1, 8) * int(k2-k1+1, 8)
    slice_elements = int(i2-i1+1, 8) * int(j2-j1+1, 8)

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      write(dump_unit) arr
      write(*,'(A,I12,A)') '[KIO] 3D int: '//trim(filename)//' (', nelements, ')'
    else
      if (slice_elements > KIO_MAX_ELEMENTS_4BYTE) then
        call kio_error('2D slice exceeds 1GB limit. Cannot auto-split further.')
      end if
      do k = k1, k2
        write(dump_unit) arr(i1:i2, j1:j2, k)
      end do
      write(*,'(A,I12,A)') '[KIO] 3D int: '//trim(filename)//' (', nelements, ') [2D slices]'
    end if

    close(dump_unit)

  end subroutine kio_dump_array_3d_int

  !=====================================================================
  !=====================================================================
  !  LOAD ROUTINES - PARAMETERS
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Load parameters from text file (key = value format)
  !---------------------------------------------------------------------
  subroutine kio_load_params(filename)
    character(len=*), intent(in) :: filename

    integer :: ios, eq_pos
    character(len=512) :: line, key, val

    ! Ensure buffer limits are initialized (check environment variable)
    call ensure_limits_initialized()

    ! Reset state
    num_params = 0
    params_loaded = .false.

    open(unit=20, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      call kio_error('Cannot open parameter file: '//trim(filename))
    end if

    do while (.true.)
      read(20, '(A)', iostat=ios) line
      if (ios /= 0) exit

      ! Skip empty lines and comments
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))

        if (num_params < MAX_PARAMS) then
          num_params = num_params + 1
          params(num_params)%key = trim(key)
          params(num_params)%value = trim(val)
        else
          call kio_error('Too many parameters (max 256)')
        end if
      end if
    end do

    close(20)
    params_loaded = .true.

  end subroutine kio_load_params

  !---------------------------------------------------------------------
  ! Check if parameters are loaded
  !---------------------------------------------------------------------
  function kio_params_loaded() result(loaded)
    logical :: loaded
    loaded = params_loaded
  end function kio_params_loaded

  !---------------------------------------------------------------------
  ! Get integer parameter by name
  !---------------------------------------------------------------------
  subroutine kio_get_param_i(name, val, default_val)
    character(len=*), intent(in) :: name
    integer, intent(out) :: val
    integer, intent(in), optional :: default_val

    integer :: i, ios

    do i = 1, num_params
      if (trim(params(i)%key) == trim(name)) then
        read(params(i)%value, *, iostat=ios) val
        if (ios == 0) return
      end if
    end do

    if (present(default_val)) then
      val = default_val
    else
      call kio_error('Parameter not found: '//trim(name))
    end if

  end subroutine kio_get_param_i

  !---------------------------------------------------------------------
  ! Get real (4-byte) parameter by name
  !---------------------------------------------------------------------
  subroutine kio_get_param_r(name, val, default_val)
    character(len=*), intent(in) :: name
    real, intent(out) :: val
    real, intent(in), optional :: default_val

    integer :: i, ios

    do i = 1, num_params
      if (trim(params(i)%key) == trim(name)) then
        read(params(i)%value, *, iostat=ios) val
        if (ios == 0) return
      end if
    end do

    if (present(default_val)) then
      val = default_val
    else
      call kio_error('Parameter not found: '//trim(name))
    end if

  end subroutine kio_get_param_r

  !---------------------------------------------------------------------
  ! Get real (8-byte) parameter by name
  !---------------------------------------------------------------------
  subroutine kio_get_param_d(name, val, default_val)
    character(len=*), intent(in) :: name
    real(8), intent(out) :: val
    real(8), intent(in), optional :: default_val

    integer :: i, ios

    do i = 1, num_params
      if (trim(params(i)%key) == trim(name)) then
        read(params(i)%value, *, iostat=ios) val
        if (ios == 0) return
      end if
    end do

    if (present(default_val)) then
      val = default_val
    else
      call kio_error('Parameter not found: '//trim(name))
    end if

  end subroutine kio_get_param_d

  !---------------------------------------------------------------------
  ! Get character parameter by name
  !---------------------------------------------------------------------
  subroutine kio_get_param_c(name, val, default_val)
    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: val
    character(len=*), intent(in), optional :: default_val

    integer :: i

    do i = 1, num_params
      if (trim(params(i)%key) == trim(name)) then
        val = trim(params(i)%value)
        return
      end if
    end do

    if (present(default_val)) then
      val = default_val
    else
      call kio_error('Parameter not found: '//trim(name))
    end if

  end subroutine kio_get_param_c

  !=====================================================================
  !=====================================================================
  !  LOAD ROUTINES - REAL ARRAYS
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Load 1D real array from binary file
  ! Auto-splits into chunks if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_load_array_1d(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    real, intent(out) :: arr(i1:i2)

    integer(8) :: nelements, chunk_size, chunk_start, chunk_end, remaining
    integer :: ios

    nelements = int(i2-i1+1, 8)

    open(unit=20, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      call kio_error('Cannot open file: '//trim(filename))
    end if

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      read(20, iostat=ios) arr
      if (ios /= 0) then
        call kio_error('Error reading file: '//trim(filename))
      end if
    else
      chunk_size = KIO_MAX_ELEMENTS_4BYTE
      chunk_start = int(i1, 8)
      remaining = nelements

      do while (remaining > 0)
        chunk_end = min(chunk_start + chunk_size - 1, int(i2, 8))
        read(20, iostat=ios) arr(int(chunk_start):int(chunk_end))
        if (ios /= 0) then
          call kio_error('Error reading file (chunk): '//trim(filename))
        end if
        remaining = remaining - (chunk_end - chunk_start + 1)
        chunk_start = chunk_end + 1
      end do
    end if

    close(20)

  end subroutine kio_load_array_1d

  !---------------------------------------------------------------------
  ! Load 2D real array from binary file
  ! Auto-splits into 1D slices if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_load_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: arr(i1:i2, j1:j2)

    integer(8) :: nelements
    integer :: ios, j

    nelements = int(i2-i1+1, 8) * int(j2-j1+1, 8)

    open(unit=20, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      call kio_error('Cannot open file: '//trim(filename))
    end if

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      read(20, iostat=ios) arr
      if (ios /= 0) then
        call kio_error('Error reading file: '//trim(filename))
      end if
    else
      do j = j1, j2
        read(20, iostat=ios) arr(i1:i2, j)
        if (ios /= 0) then
          call kio_error('Error reading file (slice): '//trim(filename))
        end if
      end do
    end if

    close(20)

  end subroutine kio_load_array_2d

  !---------------------------------------------------------------------
  ! Load 3D real array from binary file
  ! Auto-splits into 2D slices if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_load_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)

    integer(8) :: nelements, slice_elements
    integer :: ios, k

    nelements = int(i2-i1+1, 8) * int(j2-j1+1, 8) * int(k2-k1+1, 8)
    slice_elements = int(i2-i1+1, 8) * int(j2-j1+1, 8)

    open(unit=20, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      call kio_error('Cannot open file: '//trim(filename))
    end if

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      read(20, iostat=ios) arr
      if (ios /= 0) then
        call kio_error('Error reading file: '//trim(filename))
      end if
    else
      if (slice_elements > KIO_MAX_ELEMENTS_4BYTE) then
        call kio_error('2D slice exceeds 1GB limit. Cannot auto-split further.')
      end if
      do k = k1, k2
        read(20, iostat=ios) arr(i1:i2, j1:j2, k)
        if (ios /= 0) then
          call kio_error('Error reading file (2D slice): '//trim(filename))
        end if
      end do
    end if

    close(20)

  end subroutine kio_load_array_3d

  !---------------------------------------------------------------------
  ! Load 4D real array from binary file
  ! Always reads in 2D slices to avoid 1GB buffer overflow
  !---------------------------------------------------------------------
  subroutine kio_load_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)

    integer(8) :: slice_elements
    integer :: ios, k, l

    slice_elements = int(i2-i1+1, 8) * int(j2-j1+1, 8)

    if (slice_elements > KIO_MAX_ELEMENTS_4BYTE) then
      call kio_error('2D slice exceeds 1GB limit. Cannot auto-split further.')
    end if

    open(unit=20, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      call kio_error('Cannot open file: '//trim(filename))
    end if

    ! Always read in 2D slices
    do l = l1, l2
      do k = k1, k2
        read(20, iostat=ios) arr(i1:i2, j1:j2, k, l)
        if (ios /= 0) then
          call kio_error('Error reading file (4D slice): '//trim(filename))
        end if
      end do
    end do

    close(20)

  end subroutine kio_load_array_4d

  !=====================================================================
  !=====================================================================
  !  LOAD ROUTINES - INTEGER ARRAYS
  !=====================================================================
  !=====================================================================

  !---------------------------------------------------------------------
  ! Load 1D integer array from binary file
  ! Auto-splits into chunks if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_load_array_1d_int(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    integer, intent(out) :: arr(i1:i2)

    integer(8) :: nelements, chunk_size, chunk_start, chunk_end, remaining
    integer :: ios

    nelements = int(i2-i1+1, 8)

    open(unit=20, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      call kio_error('Cannot open file: '//trim(filename))
    end if

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      read(20, iostat=ios) arr
      if (ios /= 0) then
        call kio_error('Error reading file: '//trim(filename))
      end if
    else
      chunk_size = KIO_MAX_ELEMENTS_4BYTE
      chunk_start = int(i1, 8)
      remaining = nelements

      do while (remaining > 0)
        chunk_end = min(chunk_start + chunk_size - 1, int(i2, 8))
        read(20, iostat=ios) arr(int(chunk_start):int(chunk_end))
        if (ios /= 0) then
          call kio_error('Error reading file (chunk): '//trim(filename))
        end if
        remaining = remaining - (chunk_end - chunk_start + 1)
        chunk_start = chunk_end + 1
      end do
    end if

    close(20)

  end subroutine kio_load_array_1d_int

  !---------------------------------------------------------------------
  ! Load 2D integer array from binary file
  ! Auto-splits into 1D slices if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_load_array_2d_int(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    integer, intent(out) :: arr(i1:i2, j1:j2)

    integer(8) :: nelements
    integer :: ios, j

    nelements = int(i2-i1+1, 8) * int(j2-j1+1, 8)

    open(unit=20, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      call kio_error('Cannot open file: '//trim(filename))
    end if

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      read(20, iostat=ios) arr
      if (ios /= 0) then
        call kio_error('Error reading file: '//trim(filename))
      end if
    else
      do j = j1, j2
        read(20, iostat=ios) arr(i1:i2, j)
        if (ios /= 0) then
          call kio_error('Error reading file (slice): '//trim(filename))
        end if
      end do
    end if

    close(20)

  end subroutine kio_load_array_2d_int

  !---------------------------------------------------------------------
  ! Load 3D integer array from binary file
  ! Auto-splits into 2D slices if array exceeds 1GB buffer limit
  !---------------------------------------------------------------------
  subroutine kio_load_array_3d_int(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    integer, intent(out) :: arr(i1:i2, j1:j2, k1:k2)

    integer(8) :: nelements, slice_elements
    integer :: ios, k

    nelements = int(i2-i1+1, 8) * int(j2-j1+1, 8) * int(k2-k1+1, 8)
    slice_elements = int(i2-i1+1, 8) * int(j2-j1+1, 8)

    open(unit=20, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      call kio_error('Cannot open file: '//trim(filename))
    end if

    if (nelements <= KIO_MAX_ELEMENTS_4BYTE) then
      read(20, iostat=ios) arr
      if (ios /= 0) then
        call kio_error('Error reading file: '//trim(filename))
      end if
    else
      if (slice_elements > KIO_MAX_ELEMENTS_4BYTE) then
        call kio_error('2D slice exceeds 1GB limit. Cannot auto-split further.')
      end if
      do k = k1, k2
        read(20, iostat=ios) arr(i1:i2, j1:j2, k)
        if (ios /= 0) then
          call kio_error('Error reading file (2D slice): '//trim(filename))
        end if
      end do
    end if

    close(20)

  end subroutine kio_load_array_3d_int

end module m_kernel_io

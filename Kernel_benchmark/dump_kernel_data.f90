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
! Environment Variables:
!   DUMP_TARGETS - Comma-separated list of kernel names to dump
!                  e.g., "steps,diverpih,advp"
!                  If not set or empty, all kernels are dumped
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
!     if (dump_is_active()) then  ! Check if this kernel should be dumped
!       call dump_scalar_i('ni', ni)
!       call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
!     end if
!   end if
!   !$omp parallel ...
!   !$omp end parallel
!   if (dump_call_count == DUMP_TARGET_CALL .and. .not. dump_done) then
!     if (dump_is_active()) then
!       call dump_array_3d('output.bin', output, ...)
!       call dump_finalize()
!     end if
!     dump_done = .true.
!   end if
!
!***********************************************************************
module m_dump_kernel
  implicit none
  private

  ! Public subroutines
  public :: dump_init, dump_finalize, dump_is_active
  public :: dump_scalar_i, dump_scalar_i8, dump_scalar_r, dump_scalar_d, dump_scalar_c
  public :: dump_scalar_s  ! Alias for dump_scalar_c (string)
  public :: dump_array_1d, dump_array_1d_r8, dump_array_2d, dump_array_3d, dump_array_4d
  public :: dump_array_1d_int, dump_array_2d_int, dump_array_3d_int
  public :: dump_array_1d_i, dump_array_2d_i, dump_array_3d_i  ! Aliases

  ! Module variables
  character(len=256), save :: dump_dir = './kernel_dump'
  integer, save :: dump_unit = 100
  integer, save :: param_unit = 101
  logical, save :: param_file_open = .false.

  ! Selective dump control
  character(len=4096), save :: dump_targets = ''     ! Comma-separated target list
  logical, save :: dump_targets_loaded = .false.
  logical, save :: dump_active = .false.             ! Is current kernel in target list?
  integer, save :: dump_count = 0                    ! Total dumps completed

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
  ! Load dump targets from environment variable (called once)
  !=====================================================================
  subroutine load_dump_targets()
    integer :: env_len

    if (dump_targets_loaded) return

    call get_environment_variable('DUMP_TARGETS', dump_targets, length=env_len)

    if (env_len > 0) then
      write(*,'(A)') '[DUMP] ========================================'
      write(*,'(A)') '[DUMP] Selective dump mode enabled'
      write(*,'(A)') '[DUMP] DUMP_TARGETS: '//trim(dump_targets)
      write(*,'(A)') '[DUMP] ========================================'
    end if

    dump_targets_loaded = .true.
  end subroutine load_dump_targets

  !=====================================================================
  ! Check if kernel name is in the target list
  !=====================================================================
  function is_kernel_in_targets(kernel_name) result(found)
    character(len=*), intent(in) :: kernel_name
    logical :: found
    integer :: pos

    found = .false.

    ! If no targets specified, dump all
    if (len_trim(dump_targets) == 0) then
      found = .true.
      return
    end if

    ! Search for kernel name in comma-separated list
    ! Add commas to handle partial matches: ",name," or "name," or ",name"
    pos = index(','//trim(dump_targets)//',', ','//trim(kernel_name)//',')
    found = (pos > 0)

  end function is_kernel_in_targets

  !=====================================================================
  ! Check if dump data already exists for a kernel
  !=====================================================================
  function dump_exists(kernel_name) result(exists)
    character(len=*), intent(in) :: kernel_name
    logical :: exists
    character(len=512) :: check_path

    check_path = './kernel_dump/'//trim(kernel_name)//'/params.txt'
    inquire(file=trim(check_path), exist=exists)

  end function dump_exists

  !=====================================================================
  ! Initialize dump directory and open params file
  !=====================================================================
  subroutine dump_init(kernel_name)
    character(len=*), intent(in) :: kernel_name
    integer :: ios

    ! Load targets on first call
    call load_dump_targets()

    ! Check if this kernel should be dumped
    dump_active = is_kernel_in_targets(kernel_name)

    if (.not. dump_active) then
      write(*,'(A)') '[DUMP] SKIP: '//trim(kernel_name)//' (not in DUMP_TARGETS)'
      return
    end if

    ! Check if dump data already exists
    if (dump_exists(kernel_name)) then
      write(*,'(A)') '[DUMP] SKIP: '//trim(kernel_name)//' (data already exists)'
      dump_active = .false.
      return
    end if

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
  ! Check if current dump is active (kernel is in target list)
  !=====================================================================
  function dump_is_active() result(active)
    logical :: active
    active = dump_active
  end function dump_is_active

  !=====================================================================
  ! Finalize (close files and print summary)
  !=====================================================================
  subroutine dump_finalize()
    if (.not. dump_active) return

    if (param_file_open) then
      close(param_unit)
      param_file_open = .false.
    end if

    dump_count = dump_count + 1

    write(*,'(A)') '[DUMP] ========================================'
    write(*,'(A)') '[DUMP] Complete. Data saved to: '//trim(dump_dir)
    write(*,'(A,I6)') '[DUMP] Total dumps so far: ', dump_count
    write(*,'(A)') '[DUMP] ========================================'

    dump_active = .false.
  end subroutine dump_finalize

  !=====================================================================
  ! Dump scalar integer
  !=====================================================================
  subroutine dump_scalar_i(name, val)
    character(len=*), intent(in) :: name
    integer, intent(in) :: val

    if (.not. dump_active) return

    if (param_file_open) then
      write(param_unit, '(A,A,I12)') trim(name), ' = ', val
      flush(param_unit)
    end if
    write(*,'(A,A,I12)') '[DUMP] ', trim(name)//' = ', val

  end subroutine dump_scalar_i

  !=====================================================================
  ! Dump scalar integer (8-byte / 64-bit)
  !=====================================================================
  subroutine dump_scalar_i8(name, val)
    character(len=*), intent(in) :: name
    integer(8), intent(in) :: val

    if (.not. dump_active) return

    if (param_file_open) then
      write(param_unit, '(A,A,I20)') trim(name), ' = ', val
      flush(param_unit)
    end if
    write(*,'(A,A,I20)') '[DUMP] ', trim(name)//' = ', val

  end subroutine dump_scalar_i8

  !=====================================================================
  ! Dump scalar real (single precision)
  !=====================================================================
  subroutine dump_scalar_r(name, val)
    character(len=*), intent(in) :: name
    real, intent(in) :: val

    if (.not. dump_active) return

    if (param_file_open) then
      write(param_unit, '(A,A,ES20.12)') trim(name), ' = ', val
      flush(param_unit)
    end if
    write(*,'(A,A,ES15.7)') '[DUMP] ', trim(name)//' = ', val

  end subroutine dump_scalar_r

  !=====================================================================
  ! Dump scalar real (double precision)
  !=====================================================================
  subroutine dump_scalar_d(name, val)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: val

    if (.not. dump_active) return

    if (param_file_open) then
      write(param_unit, '(A,A,ES24.16)') trim(name), ' = ', val
      flush(param_unit)
    end if
    write(*,'(A,A,ES20.12)') '[DUMP] ', trim(name)//' = ', val

  end subroutine dump_scalar_d

  !=====================================================================
  ! Dump scalar character string
  !=====================================================================
  subroutine dump_scalar_c(name, val)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val

    if (.not. dump_active) return

    if (param_file_open) then
      write(param_unit, '(A,A,A)') trim(name), ' = ', trim(val)
      flush(param_unit)
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

    if (.not. dump_active) return

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = i2 - i1 + 1
    write(*,'(A,I10,A)') '[DUMP] 1D: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_1d

  !=====================================================================
  ! Dump 1D real*8 array to binary file
  !=====================================================================
  subroutine dump_array_1d_r8(filename, arr, i1, i2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2
    real(8), intent(in) :: arr(i1:i2)

    character(len=512) :: filepath
    integer :: nelements

    if (.not. dump_active) return

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = i2 - i1 + 1
    write(*,'(A,I10,A)') '[DUMP] 1D_r8: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_1d_r8

  !=====================================================================
  ! Dump 2D real array to binary file
  !=====================================================================
  subroutine dump_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(in) :: arr(i1:i2, j1:j2)

    character(len=512) :: filepath
    integer(8) :: nelements

    if (.not. dump_active) return

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
  ! Note: Write row-by-row to minimize temporary array size
  !=====================================================================
  subroutine dump_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2)

    character(len=512) :: filepath
    integer(8) :: nelements
    integer :: j, k

    if (.not. dump_active) return

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')

    ! Write row-by-row to minimize temporary array size
    do k = k1, k2
      do j = j1, j2
        write(dump_unit) arr(i1:i2, j, k)
      end do
      flush(dump_unit)
    end do

    close(dump_unit)

    nelements = int(i2-i1+1,8) * int(j2-j1+1,8) * int(k2-k1+1,8)
    write(*,'(A,I12,A)') '[DUMP] 3D: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_3d

  !=====================================================================
  ! Dump 4D real array to binary file
  ! Note: Write row-by-row to minimize temporary array size
  !=====================================================================
  subroutine dump_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)

    character(len=512) :: filepath
    integer(8) :: nelements
    integer :: j, k, l

    if (.not. dump_active) return

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')

    ! Write row-by-row to minimize temporary array size
    do l = l1, l2
      do k = k1, k2
        do j = j1, j2
          write(dump_unit) arr(i1:i2, j, k, l)
        end do
        flush(dump_unit)
      end do
    end do

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

    if (.not. dump_active) return

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

    if (.not. dump_active) return

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

    if (.not. dump_active) return

    filepath = trim(dump_dir)//'/'//trim(filename)
    open(unit=dump_unit, file=filepath, status='replace', &
         access='stream', form='unformatted')
    write(dump_unit) arr
    close(dump_unit)

    nelements = int(i2-i1+1,8) * int(j2-j1+1,8) * int(k2-k1+1,8)
    write(*,'(A,I12,A)') '[DUMP] 3D int: '//trim(filename)//' (', nelements, ')'

  end subroutine dump_array_3d_int

end module m_dump_kernel

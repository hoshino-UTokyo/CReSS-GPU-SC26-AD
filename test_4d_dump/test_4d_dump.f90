!***********************************************************************
! Minimal reproduction test for 4D array dump segfault
!***********************************************************************
! Problem: dump_array_4d crashes with segfault in steps.f90
! Array size: (0:900, 0:900, 1:128, 1:2) = 207,821,056 elements
!
! Hypothesis: Array slice arr(i1:i2, j1:j2, k, l) creates a temp array
!             that exceeds memory limits
!***********************************************************************
program test_4d_dump
  implicit none

  ! Parameters matching actual solver dimensions
  integer, parameter :: ni = 899
  integer, parameter :: nj = 899
  integer, parameter :: nk = 128
  integer, parameter :: nqw = 2

  ! Bounds
  integer, parameter :: i1 = 0, i2 = ni+1   ! 901 elements
  integer, parameter :: j1 = 0, j2 = nj+1   ! 901 elements
  integer, parameter :: k1 = 1, k2 = nk     ! 128 elements
  integer, parameter :: l1 = 1, l2 = nqw    ! 2 elements

  ! 4D array - allocatable to avoid stack overflow
  real, allocatable :: arr(:,:,:,:)

  ! I/O variables
  integer :: dump_unit = 100
  integer :: k, l, ios
  integer(8) :: nelements, slice_size
  character(len=256) :: filepath

  ! Calculate sizes
  slice_size = int(i2-i1+1, 8) * int(j2-j1+1, 8)  ! 2D slice size
  nelements = slice_size * int(k2-k1+1, 8) * int(l2-l1+1, 8)  ! Total elements

  print '(A)', '=============================================='
  print '(A)', ' 4D Array Dump Segfault Reproduction Test'
  print '(A)', '=============================================='
  print '(A,I5,A,I5,A,I5,A,I5,A,I5,A,I5,A,I5,A,I5,A)', ' Array bounds: (', i1, ':', i2, ', ', &
        j1, ':', j2, ', ', k1, ':', k2, ', ', l1, ':', l2, ')'
  print '(A,I15)', ' Total elements: ', nelements
  print '(A,F10.2,A)', ' Memory required: ', real(nelements * 4) / 1024.0 / 1024.0, ' MB'
  print '(A,I15)', ' 2D slice size: ', slice_size
  print '(A,F10.2,A)', ' Slice memory: ', real(slice_size * 4) / 1024.0 / 1024.0, ' MB'
  print '(A)', '=============================================='

  ! Allocate array
  print '(A)', 'Allocating array...'
  allocate(arr(i1:i2, j1:j2, k1:k2, l1:l2), stat=ios)
  if (ios /= 0) then
    print '(A)', 'ERROR: Failed to allocate array!'
    stop 1
  end if
  print '(A)', 'Allocation successful.'

  ! Initialize with test data
  print '(A)', 'Initializing array with test data...'
  arr = 1.0
  print '(A)', 'Initialization done.'

  ! Create output directory
  call execute_command_line('mkdir -p ./test_output', wait=.true.)
  filepath = './test_output/test_4d.bin'

  !-----------------------------------------------------------------
  ! Test 1: Current method (array slice notation)
  !-----------------------------------------------------------------
  print '(A)', ''
  print '(A)', '=== Test 1: Array slice notation (current method) ==='
  print '(A)', 'This is the method that crashes in solver...'

  open(unit=dump_unit, file=trim(filepath), status='replace', &
       access='stream', form='unformatted', iostat=ios)
  if (ios /= 0) then
    print '(A)', 'ERROR: Failed to open file!'
    stop 1
  end if

  do l = l1, l2
    do k = k1, k2
      ! This line may crash due to temporary array creation
      write(dump_unit, iostat=ios) arr(i1:i2, j1:j2, k, l)
      if (ios /= 0) then
        print '(A,I0,A,I0)', 'ERROR: Write failed at k=', k, ', l=', l
        stop 1
      end if
      flush(dump_unit)

      ! Progress indicator
      if (mod(k, 32) == 0) then
        print '(A,I0,A,I0)', '  Written l=', l, ', k=', k
      end if
    end do
  end do

  close(dump_unit)
  print '(A)', 'Test 1 PASSED: Array slice notation works.'

  ! Clean up test file
  call execute_command_line('rm -f '//trim(filepath), wait=.true.)

  !-----------------------------------------------------------------
  ! Test 2: Element-by-element method (alternative)
  !-----------------------------------------------------------------
  print '(A)', ''
  print '(A)', '=== Test 2: Element-by-element write (alternative) ==='

  open(unit=dump_unit, file=trim(filepath), status='replace', &
       access='stream', form='unformatted', iostat=ios)
  if (ios /= 0) then
    print '(A)', 'ERROR: Failed to open file!'
    stop 1
  end if

  block
    integer :: i, j
    do l = l1, l2
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            write(dump_unit) arr(i, j, k, l)
          end do
        end do
        flush(dump_unit)

        if (mod(k, 32) == 0) then
          print '(A,I0,A,I0)', '  Written l=', l, ', k=', k
        end if
      end do
    end do
  end block

  close(dump_unit)
  print '(A)', 'Test 2 PASSED: Element-by-element write works.'

  ! Clean up test file
  call execute_command_line('rm -f '//trim(filepath), wait=.true.)

  !-----------------------------------------------------------------
  ! Test 3: Row-by-row write (efficient alternative)
  !-----------------------------------------------------------------
  print '(A)', ''
  print '(A)', '=== Test 3: Row-by-row write (more efficient alternative) ==='

  open(unit=dump_unit, file=trim(filepath), status='replace', &
       access='stream', form='unformatted', iostat=ios)
  if (ios /= 0) then
    print '(A)', 'ERROR: Failed to open file!'
    stop 1
  end if

  block
    integer :: j
    do l = l1, l2
      do k = k1, k2
        do j = j1, j2
          ! Write one row at a time - minimal temp array (901 elements)
          write(dump_unit) arr(i1:i2, j, k, l)
        end do
        flush(dump_unit)

        if (mod(k, 32) == 0) then
          print '(A,I0,A,I0)', '  Written l=', l, ', k=', k
        end if
      end do
    end do
  end block

  close(dump_unit)
  print '(A)', 'Test 3 PASSED: Row-by-row write works.'

  ! Clean up
  call execute_command_line('rm -f '//trim(filepath), wait=.true.)
  call execute_command_line('rmdir ./test_output 2>/dev/null', wait=.true.)
  deallocate(arr)

  print '(A)', ''
  print '(A)', '=============================================='
  print '(A)', ' All tests completed!'
  print '(A)', ' If Test 1 fails but Test 2/3 pass,'
  print '(A)', ' use row-by-row write as the fix.'
  print '(A)', '=============================================='

end program test_4d_dump

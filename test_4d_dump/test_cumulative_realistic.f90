!***********************************************************************
! Realistic cumulative memory test simulating actual solver dump
!***********************************************************************
! Simulates the solver's dump behavior:
! - Many 3D arrays dumped before 4D arrays
! - Tests both 2D slice method (original) and 1D row method (fixed)
! - Command line: method (1=2D, 2=1D), num_3d, num_4d
!***********************************************************************
program test_cumulative_realistic
  implicit none

  ! Actual solver dimensions
  integer, parameter :: ni = 899
  integer, parameter :: nj = 899
  integer, parameter :: nk = 128
  integer, parameter :: nqw = 2

  ! Bounds
  integer, parameter :: i1 = 0, i2 = ni+1   ! 901 elements
  integer, parameter :: j1 = 0, j2 = nj+1   ! 901 elements
  integer, parameter :: k1 = 1, k2 = nk     ! 128 elements
  integer, parameter :: l1 = 1, l2 = nqw    ! 2 elements

  ! Arrays
  real, allocatable :: arr3d(:,:,:)
  real, allocatable :: arr4d(:,:,:,:)

  ! Parameters
  integer :: method          ! 1=2D slice, 2=1D row
  integer :: num_3d_dumps    ! Number of 3D dumps
  integer :: num_4d_dumps    ! Number of 4D dumps
  character(len=32) :: arg

  integer :: n, ios
  integer(8) :: total_bytes
  real :: start_time, end_time

  ! Get command line arguments
  method = 1
  num_3d_dumps = 100
  num_4d_dumps = 6

  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg, *) method
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg)
    read(arg, *) num_3d_dumps
  end if
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg)
    read(arg, *) num_4d_dumps
  end if

  print '(A)', '====================================================='
  print '(A)', ' Realistic Cumulative Dump Test'
  print '(A)', '====================================================='
  if (method == 1) then
    print '(A)', ' Method: 2D slice (ORIGINAL - may fail)'
  else
    print '(A)', ' Method: 1D row (FIXED)'
  end if
  print '(A,I10)', ' Number of 3D dumps: ', num_3d_dumps
  print '(A,I10)', ' Number of 4D dumps: ', num_4d_dumps
  print '(A,F12.2,A)', ' 3D array size: ', real(i2-i1+1)*real(j2-j1+1)*real(k2-k1+1)*4/1024/1024, ' MB'
  print '(A,F12.2,A)', ' 4D array size: ', real(i2-i1+1)*real(j2-j1+1)*real(k2-k1+1)*real(l2-l1+1)*4/1024/1024, ' MB'
  print '(A)', '====================================================='
  print '(A)', ''

  call execute_command_line('mkdir -p ./test_output', wait=.true.)

  ! Allocate arrays
  print '(A)', 'Allocating 3D array...'
  allocate(arr3d(i1:i2, j1:j2, k1:k2), stat=ios)
  if (ios /= 0) stop 'ERROR: Failed to allocate 3D array!'
  arr3d = 1.0

  print '(A)', 'Allocating 4D array...'
  allocate(arr4d(i1:i2, j1:j2, k1:k2, l1:l2), stat=ios)
  if (ios /= 0) stop 'ERROR: Failed to allocate 4D array!'
  arr4d = 2.0

  print '(A)', 'Allocation done.'
  print '(A)', ''

  call cpu_time(start_time)
  total_bytes = 0

  !-----------------------------------------------------------------
  ! Phase 1: Multiple 3D array dumps
  !-----------------------------------------------------------------
  print '(A)', '=== Phase 1: 3D array dumps ==='
  do n = 1, num_3d_dumps
    if (method == 1) then
      call dump_3d_2dslice(arr3d, i1, i2, j1, j2, k1, k2, n)
    else
      call dump_3d_1drow(arr3d, i1, i2, j1, j2, k1, k2, n)
    end if
    total_bytes = total_bytes + int(i2-i1+1,8)*int(j2-j1+1,8)*int(k2-k1+1,8)*4

    if (mod(n, 20) == 0) then
      print '(A,I6,A,F10.2,A)', '  Completed ', n, ' 3D dumps (', real(total_bytes)/1024/1024/1024, ' GB total)'
    end if
  end do
  print '(A,I6,A)', 'Phase 1 done: ', num_3d_dumps, ' 3D arrays dumped'
  print '(A)', ''

  !-----------------------------------------------------------------
  ! Phase 2: 4D array dumps (this is where solver crashes)
  !-----------------------------------------------------------------
  print '(A)', '=== Phase 2: 4D array dumps (crash point) ==='
  do n = 1, num_4d_dumps
    print '(A,I3,A)', '  Starting 4D dump ', n, '...'
    if (method == 1) then
      call dump_4d_2dslice(arr4d, i1, i2, j1, j2, k1, k2, l1, l2, n)
    else
      call dump_4d_1drow(arr4d, i1, i2, j1, j2, k1, k2, l1, l2, n)
    end if
    total_bytes = total_bytes + int(i2-i1+1,8)*int(j2-j1+1,8)*int(k2-k1+1,8)*int(l2-l1+1,8)*4
    print '(A,I3,A,F10.2,A)', '  Completed 4D dump ', n, ' (', real(total_bytes)/1024/1024/1024, ' GB total)'
  end do
  print '(A)', 'Phase 2 done!'

  call cpu_time(end_time)

  ! Cleanup
  call execute_command_line('rm -rf ./test_output', wait=.true.)
  deallocate(arr3d, arr4d)

  print '(A)', ''
  print '(A)', '====================================================='
  print '(A)', ' TEST PASSED!'
  print '(A,F10.2,A)', ' Total data written: ', real(total_bytes)/1024/1024/1024, ' GB'
  print '(A,F10.2,A)', ' Time elapsed: ', end_time - start_time, ' seconds'
  print '(A)', '====================================================='

contains

  !-----------------------------------------------------------------
  ! 3D dump using 2D slices (original method)
  !-----------------------------------------------------------------
  subroutine dump_3d_2dslice(arr, i1, i2, j1, j2, k1, k2, iter)
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, iter
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2)
    character(len=256) :: filepath
    integer :: k, ios

    write(filepath, '(A,I4.4,A)') './test_output/arr3d_', iter, '.bin'
    open(unit=100, file=trim(filepath), status='replace', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) stop 'ERROR: Failed to open 3D file!'

    do k = k1, k2
      write(100) arr(i1:i2, j1:j2, k)   ! 2D slice ~3.1MB
      flush(100)
    end do

    close(100)
  end subroutine dump_3d_2dslice

  !-----------------------------------------------------------------
  ! 3D dump using 1D rows (fixed method)
  !-----------------------------------------------------------------
  subroutine dump_3d_1drow(arr, i1, i2, j1, j2, k1, k2, iter)
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, iter
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2)
    character(len=256) :: filepath
    integer :: j, k, ios

    write(filepath, '(A,I4.4,A)') './test_output/arr3d_', iter, '.bin'
    open(unit=100, file=trim(filepath), status='replace', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) stop 'ERROR: Failed to open 3D file!'

    do k = k1, k2
      do j = j1, j2
        write(100) arr(i1:i2, j, k)     ! 1D row ~3.5KB
      end do
      flush(100)
    end do

    close(100)
  end subroutine dump_3d_1drow

  !-----------------------------------------------------------------
  ! 4D dump using 2D slices (original method - causes crash)
  !-----------------------------------------------------------------
  subroutine dump_4d_2dslice(arr, i1, i2, j1, j2, k1, k2, l1, l2, iter)
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2, iter
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)
    character(len=256) :: filepath
    integer :: k, l, ios

    write(filepath, '(A,I4.4,A)') './test_output/arr4d_', iter, '.bin'
    open(unit=100, file=trim(filepath), status='replace', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) stop 'ERROR: Failed to open 4D file!'

    do l = l1, l2
      do k = k1, k2
        write(100) arr(i1:i2, j1:j2, k, l)  ! 2D slice ~3.1MB
        flush(100)
      end do
    end do

    close(100)
  end subroutine dump_4d_2dslice

  !-----------------------------------------------------------------
  ! 4D dump using 1D rows (fixed method)
  !-----------------------------------------------------------------
  subroutine dump_4d_1drow(arr, i1, i2, j1, j2, k1, k2, l1, l2, iter)
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2, iter
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)
    character(len=256) :: filepath
    integer :: j, k, l, ios

    write(filepath, '(A,I4.4,A)') './test_output/arr4d_', iter, '.bin'
    open(unit=100, file=trim(filepath), status='replace', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) stop 'ERROR: Failed to open 4D file!'

    do l = l1, l2
      do k = k1, k2
        do j = j1, j2
          write(100) arr(i1:i2, j, k, l)    ! 1D row ~3.5KB
        end do
        flush(100)
      end do
    end do

    close(100)
  end subroutine dump_4d_1drow

end program test_cumulative_realistic

!***********************************************************************
! Test to investigate temporary array creation in array slicing
!***********************************************************************
! This test checks if passing array slices to write() creates temp arrays
! and if this causes memory accumulation
!***********************************************************************
program test_temp_array
  implicit none

  integer, parameter :: ni = 899
  integer, parameter :: nj = 899
  integer, parameter :: nk = 128
  integer, parameter :: nqw = 2

  integer, parameter :: i1 = 0, i2 = ni+1
  integer, parameter :: j1 = 0, j2 = nj+1
  integer, parameter :: k1 = 1, k2 = nk
  integer, parameter :: l1 = 1, l2 = nqw

  real, allocatable :: arr(:,:,:,:)
  integer :: dump_unit = 100
  integer :: j, k, l, ios, iter
  character(len=256) :: filepath

  ! Number of 4D dumps to simulate
  integer :: num_dumps
  character(len=32) :: arg

  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg, *) num_dumps
  else
    num_dumps = 10
  end if

  print '(A)', '================================================='
  print '(A)', ' Temporary Array Memory Test'
  print '(A)', '================================================='
  print '(A,I10)', ' Number of 4D dumps: ', num_dumps
  print '(A)', '================================================='

  call execute_command_line('mkdir -p ./test_output', wait=.true.)

  ! Allocate 4D array
  print '(A)', 'Allocating 4D array...'
  allocate(arr(i1:i2, j1:j2, k1:k2, l1:l2), stat=ios)
  if (ios /= 0) then
    print '(A)', 'ERROR: Failed to allocate!'
    stop 1
  end if
  arr = 1.0
  print '(A)', 'Allocation done.'
  print '(A)', ''

  ! Call dump subroutine multiple times
  do iter = 1, num_dumps
    print '(A,I5,A)', '=== Dump iteration ', iter, ' ==='
    call dump_4d_rowwise(arr, i1, i2, j1, j2, k1, k2, l1, l2, iter)
  end do

  ! Cleanup
  call execute_command_line('rm -rf ./test_output', wait=.true.)
  deallocate(arr)

  print '(A)', ''
  print '(A)', '================================================='
  print '(A)', ' TEST PASSED!'
  print '(A)', '================================================='

contains

  subroutine dump_4d_rowwise(arr, i1, i2, j1, j2, k1, k2, l1, l2, iter)
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2, iter
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)

    character(len=256) :: filepath
    integer :: j, k, l, ios
    integer :: dump_unit = 100

    write(filepath, '(A,I4.4,A)') './test_output/dump4d_', iter, '.bin'
    open(unit=dump_unit, file=trim(filepath), status='replace', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      print '(A)', 'ERROR: Failed to open file!'
      stop 1
    end if

    ! Row-by-row write to minimize temporary array size
    do l = l1, l2
      do k = k1, k2
        do j = j1, j2
          write(dump_unit) arr(i1:i2, j, k, l)
        end do
        flush(dump_unit)
      end do
    end do

    close(dump_unit)
    print '(A,I5,A)', '  Completed dump ', iter, ' (row-by-row)'

  end subroutine dump_4d_rowwise

end program test_temp_array

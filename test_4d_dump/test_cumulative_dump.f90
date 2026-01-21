!***********************************************************************
! Test to simulate cumulative memory usage during solver dump
!***********************************************************************
! In actual solver run:
! - 753 3D arrays (each ~395MB) are dumped before 4D array
! - 4D array dump crashes with segfault after ~650 dumps
! - Error: "Cgroup memsw limit exceeded"
!
! This test simulates the cumulative memory effect
!***********************************************************************
program test_cumulative_dump
  implicit none

  integer, parameter :: ni = 899
  integer, parameter :: nj = 899
  integer, parameter :: nk = 128
  integer, parameter :: nqw = 2

  ! Bounds for 3D
  integer, parameter :: i1 = 0, i2 = ni+1
  integer, parameter :: j1 = 0, j2 = nj+1
  integer, parameter :: k1 = 1, k2 = nk
  integer, parameter :: l1 = 1, l2 = nqw

  ! Arrays
  real, allocatable :: arr3d(:,:,:)
  real, allocatable :: arr4d(:,:,:,:)

  ! I/O
  integer :: dump_unit = 100
  integer :: k, l, n, ios
  integer(8) :: nelements_3d, nelements_4d
  character(len=256) :: filepath

  ! Number of 3D dumps to simulate before 4D dump
  integer :: num_3d_dumps
  character(len=32) :: arg

  nelements_3d = int(i2-i1+1, 8) * int(j2-j1+1, 8) * int(k2-k1+1, 8)
  nelements_4d = nelements_3d * int(l2-l1+1, 8)

  ! Get number of 3D dumps from command line
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg, *) num_3d_dumps
  else
    num_3d_dumps = 100  ! Default
  end if

  print '(A)', '================================================='
  print '(A)', ' Cumulative Memory Test for 4D Dump Crash'
  print '(A)', '================================================='
  print '(A,I10)', ' Number of 3D dumps before 4D: ', num_3d_dumps
  print '(A,I18)', ' 3D array elements: ', nelements_3d
  print '(A,F12.2,A)', ' 3D array size: ', real(nelements_3d * 4) / 1024.0 / 1024.0, ' MB'
  print '(A,I18)', ' 4D array elements: ', nelements_4d
  print '(A,F12.2,A)', ' 4D array size: ', real(nelements_4d * 4) / 1024.0 / 1024.0, ' MB'
  print '(A)', '================================================='

  call execute_command_line('mkdir -p ./test_output', wait=.true.)

  ! Allocate arrays
  print '(A)', 'Allocating 3D array...'
  allocate(arr3d(i1:i2, j1:j2, k1:k2), stat=ios)
  if (ios /= 0) then
    print '(A)', 'ERROR: Failed to allocate 3D array!'
    stop 1
  end if
  arr3d = 1.0

  print '(A)', 'Allocating 4D array...'
  allocate(arr4d(i1:i2, j1:j2, k1:k2, l1:l2), stat=ios)
  if (ios /= 0) then
    print '(A)', 'ERROR: Failed to allocate 4D array!'
    stop 1
  end if
  arr4d = 2.0
  print '(A)', 'Allocation done.'
  print '(A)', ''

  !-----------------------------------------------------------------
  ! Simulate multiple 3D dumps (like in solver)
  !-----------------------------------------------------------------
  print '(A)', '=== Phase 1: Multiple 3D array dumps ==='

  do n = 1, num_3d_dumps
    write(filepath, '(A,I4.4,A)') './test_output/test3d_', n, '.bin'

    open(unit=dump_unit, file=trim(filepath), status='replace', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      print '(A,I10)', 'ERROR: Failed to open file at n=', n
      stop 1
    end if

    ! Write in 2D slices (like current dump_array_3d)
    do k = k1, k2
      write(dump_unit, iostat=ios) arr3d(i1:i2, j1:j2, k)
      if (ios /= 0) then
        print '(A,I10,A,I10)', 'ERROR: Write failed at n=', n, ', k=', k
        stop 1
      end if
      flush(dump_unit)
    end do

    close(dump_unit)

    if (mod(n, 50) == 0) then
      print '(A,I10,A)', '  Completed ', n, ' 3D dumps'
    end if
  end do

  print '(A,I10,A)', 'Phase 1 done: ', num_3d_dumps, ' 3D arrays dumped'
  print '(A)', ''

  !-----------------------------------------------------------------
  ! Now try 4D dump (this is where solver crashes)
  !-----------------------------------------------------------------
  print '(A)', '=== Phase 2: 4D array dump (crash point) ==='

  filepath = './test_output/test4d.bin'
  open(unit=dump_unit, file=trim(filepath), status='replace', &
       access='stream', form='unformatted', iostat=ios)
  if (ios /= 0) then
    print '(A)', 'ERROR: Failed to open 4D output file!'
    stop 1
  end if

  print '(A)', 'Writing 4D array in 2D slices...'

  do l = l1, l2
    do k = k1, k2
      ! This is the line that crashes in solver
      write(dump_unit, iostat=ios) arr4d(i1:i2, j1:j2, k, l)
      if (ios /= 0) then
        print '(A,I10,A,I10)', 'ERROR: 4D write failed at k=', k, ', l=', l
        stop 1
      end if
      flush(dump_unit)

      if (mod(k, 32) == 0) then
        print '(A,I5,A,I5)', '  Written l=', l, ', k=', k
      end if
    end do
  end do

  close(dump_unit)
  print '(A)', '4D array dump completed successfully!'

  ! Cleanup
  print '(A)', ''
  print '(A)', 'Cleaning up...'
  call execute_command_line('rm -rf ./test_output', wait=.true.)
  deallocate(arr3d, arr4d)

  print '(A)', '================================================='
  print '(A)', ' TEST PASSED!'
  print '(A,I10,A)', ' Successfully dumped ', num_3d_dumps, ' 3D arrays'
  print '(A)', ' and 1 4D array without crash.'
  print '(A)', '================================================='

end program test_cumulative_dump

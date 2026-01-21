! Memory comparison test: 2D slice vs 1D row writes
program test_memory
  implicit none
  integer, parameter :: ni = 899, nj = 899, nk = 128
  integer, parameter :: i1 = 0, i2 = ni+1, j1 = 0, j2 = nj+1, k1 = 1, k2 = nk
  real, allocatable :: arr(:,:,:)
  integer :: j, k, ios, method
  character(len=32) :: arg

  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    read(arg, *) method
  else
    method = 1  ! 1=2D slice, 2=1D row
  end if

  allocate(arr(i1:i2, j1:j2, k1:k2))
  arr = 1.0
  call execute_command_line('mkdir -p ./test_output', wait=.true.)
  open(unit=100, file='./test_output/test.bin', status='replace', &
       access='stream', form='unformatted')

  print '(A,I1)', 'Method: ', method

  if (method == 1) then
    ! 2D slice write (original method)
    print '(A)', 'Writing with 2D slices...'
    do k = k1, k2
      write(100) arr(i1:i2, j1:j2, k)
      flush(100)
      if (mod(k,32)==0) print '(A,I5)', '  k=', k
    end do
  else
    ! 1D row write (new method)
    print '(A)', 'Writing with 1D rows...'
    do k = k1, k2
      do j = j1, j2
        write(100) arr(i1:i2, j, k)
      end do
      flush(100)
      if (mod(k,32)==0) print '(A,I5)', '  k=', k
    end do
  end if

  close(100)
  call execute_command_line('rm -rf ./test_output', wait=.true.)
  deallocate(arr)
  print '(A)', 'Done!'
end program

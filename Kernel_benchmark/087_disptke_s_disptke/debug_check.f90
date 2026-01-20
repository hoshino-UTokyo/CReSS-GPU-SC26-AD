program debug_check
  implicit none
  integer :: ni, nj, nk, ios
  real, allocatable :: tkefrc_ref(:,:,:), tkefrc(:,:,:)
  real, allocatable :: priv(:,:,:), tke(:,:,:), rst(:,:,:)
  integer :: i, j, k
  real :: diff, max_diff, ref_val, comp_val

  ni = 899; nj = 899; nk = 128

  allocate(tkefrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(priv(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))

  open(10, file='data/tkefrc_ref.bin', access='stream', form='unformatted')
  read(10) tkefrc_ref
  close(10)

  open(10, file='data/priv.bin', access='stream', form='unformatted')
  read(10) priv
  close(10)

  open(10, file='data/tke.bin', access='stream', form='unformatted')
  read(10) tke
  close(10)

  open(10, file='data/rst.bin', access='stream', form='unformatted')
  read(10) rst
  close(10)

  write(*,*) 'Sample tkefrc_ref values at (100,100,k):'
  do k = 10, 15
    write(*,'(A,I3,A,ES15.7)') '  k=', k, ': ', tkefrc_ref(100,100,k)
  end do

  write(*,*) 'Sample priv values at (100,100,k):'
  do k = 10, 15
    write(*,'(A,I3,A,ES15.7)') '  k=', k, ': ', priv(100,100,k)
  end do

  write(*,*) 'Sample tke values at (100,100,k):'
  do k = 10, 15
    write(*,'(A,I3,A,ES15.7)') '  k=', k, ': ', tke(100,100,k)
  end do

  write(*,*) 'Sample rst values at (100,100,k):'
  do k = 10, 15
    write(*,'(A,I3,A,ES15.7)') '  k=', k, ': ', rst(100,100,k)
  end do

  write(*,*) 'Checking for zeros in priv:'
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        if (priv(i,j,k) == 0.0) then
          write(*,*) 'Zero priv at:', i, j, k
          stop
        end if
      end do
    end do
  end do
  write(*,*) 'No zeros found in priv (inner domain)'

end program

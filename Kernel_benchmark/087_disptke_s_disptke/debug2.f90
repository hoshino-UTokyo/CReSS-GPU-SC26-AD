program debug2
  implicit none
  integer :: ni, nj, nk
  real, allocatable :: tkefrc_ref(:,:,:), tkefrc_in(:,:,:), tke(:,:,:)
  integer :: i, j, k, cnt

  ni = 899; nj = 899; nk = 128

  allocate(tkefrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))

  open(10, file='data/tkefrc_ref.bin', access='stream', form='unformatted')
  read(10) tkefrc_ref
  close(10)

  open(10, file='data/tkefrc_in.bin', access='stream', form='unformatted')
  read(10) tkefrc_in
  close(10)

  open(10, file='data/tke.bin', access='stream', form='unformatted')
  read(10) tke
  close(10)

  write(*,*) 'Sample values at (100,100,k):'
  write(*,*) '  k   tkefrc_in     tkefrc_ref     diff           tke'
  do k = 10, 20
    write(*,'(I4,4ES14.6)') k, tkefrc_in(100,100,k), tkefrc_ref(100,100,k), &
         tkefrc_ref(100,100,k)-tkefrc_in(100,100,k), tke(100,100,k)
  end do

  ! Count non-zero tke
  cnt = 0
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        if (tke(i,j,k) > 0.0) cnt = cnt + 1
      end do
    end do
  end do
  write(*,*) 'Non-zero tke count:', cnt

  ! Count differences
  cnt = 0
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        if (abs(tkefrc_ref(i,j,k) - tkefrc_in(i,j,k)) > 1.0e-10) cnt = cnt + 1
      end do
    end do
  end do
  write(*,*) 'Different tkefrc count:', cnt

end program

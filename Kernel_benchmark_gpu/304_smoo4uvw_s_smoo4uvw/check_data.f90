program check_data
  implicit none
  ! ni=899, so 0:ni+1 = 0:900 = 901 elements
  integer, parameter :: ni = 899, nj = 899, nk = 128
  real :: u(0:ni+1, 0:nj+1, 1:nk)
  real :: ufrc_ref(0:ni+1, 0:nj+1, 1:nk)
  integer :: ios
  
  open(unit=10, file='data/u.bin', access='stream', form='unformatted', status='old', iostat=ios)
  if (ios /= 0) stop 'Cannot open u.bin'
  read(10) u
  close(10)
  
  open(unit=10, file='data/ufrc_ref.bin', access='stream', form='unformatted', status='old', iostat=ios)
  if (ios /= 0) stop 'Cannot open ufrc_ref.bin'
  read(10) ufrc_ref
  close(10)
  
  write(*,*) 'u(100,100,10) = ', u(100,100,10)
  write(*,*) 'u(200,200,20) = ', u(200,200,20)
  write(*,*) 'ufrc_ref(100,100,10) = ', ufrc_ref(100,100,10)
  write(*,*) 'ufrc_ref(200,200,20) = ', ufrc_ref(200,200,20)
end program

program check_ref
  implicit none
  integer :: ni, nj, nk, i, j, k
  real, allocatable :: nsq8w_ref(:,:,:), qall(:,:,:), ptv(:,:,:), pt(:,:,:), a(:,:,:)
  real, allocatable :: ptbr(:,:,:), jcb8w(:,:,:)
  real :: thresq, g, dziv, gdzv, gdzv05
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  thresq = 1.0e-12
  g = 9.8
  dziv = 0.01
  gdzv = g * dziv
  gdzv05 = 0.5 * g * dziv
  
  allocate(nsq8w_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qall(0:ni+1, 0:nj+1, 1:nk))
  
  open(unit=10, file=trim(data_dir)//'/nsq8w_ref.bin', form='unformatted', access='stream')
  read(10) nsq8w_ref
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/qall.bin', form='unformatted', access='stream')
  read(10) qall
  close(10)
  
  i = 833
  j = 331
  k = 114
  
  write(*,'(A)') 'Reference data at (833, 331, 114):'
  write(*,'(A,ES15.7)') 'nsq8w_ref     = ', nsq8w_ref(i,j,k)
  write(*,'(A,ES15.7)') 'qall(i,j,k)   = ', qall(i,j,k)
  write(*,'(A,ES15.7)') 'qall(i,j,k-1) = ', qall(i,j,k-1)
  write(*,'(A,ES15.7)') 'thresq        = ', thresq
  write(*,'(A,L1)')     'qall > thresq = ', qall(i,j,k) > thresq
  
  ! Check the other kernel computed this with if or else branch
  ! If qall > thresq, uses pt, a formula
  ! Else uses ptv formula
  
end program check_ref

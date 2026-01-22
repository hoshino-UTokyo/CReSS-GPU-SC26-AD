program check_t_single
  implicit none
  integer :: ni, nj, nk, i, j, k
  real, allocatable :: ptbr(:,:,:), ptp(:,:,:), pbr(:,:,:), pp(:,:,:)
  real :: pt_val, t_val  ! single precision
  real :: p0iv, rddvcp, tlow, rd, cp, p0
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  
  ! Single precision values (computed at runtime like GPU)
  rd = 2.870000000000E+02
  cp = 1.004000000000E+03
  p0 = 1.000000000000E+05
  tlow = 2.331600036621E+02
  
  rddvcp = rd / cp
  p0iv = 1.0e0 / p0
  
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(pbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(pp(0:ni+1, 0:nj+1, 1:nk))
  
  open(unit=10, file=trim(data_dir)//'/ptbr.bin', form='unformatted', access='stream')
  read(10) ptbr
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/ptp.bin', form='unformatted', access='stream')
  read(10) ptp
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/pbr.bin', form='unformatted', access='stream')
  read(10) pbr
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/pp.bin', form='unformatted', access='stream')
  read(10) pp
  close(10)
  
  i = 833
  j = 331
  k = 114
  
  pt_val = ptbr(i,j,k) + ptp(i,j,k)
  t_val = pt_val * exp(rddvcp * log(p0iv * (pbr(i,j,k) + pp(i,j,k))))
  
  write(*,'(A)') 'Single precision T calculation at (833, 331, 114):'
  write(*,'(A,ES15.7)') 'rddvcp (computed)     = ', rddvcp
  write(*,'(A,ES15.7)') 'p0iv (computed)       = ', p0iv
  write(*,'(A,ES15.7)') 'pt_val                = ', pt_val
  write(*,'(A,ES15.7)') 't_val                 = ', t_val
  write(*,'(A,ES15.7)') 'tlow                  = ', tlow
  write(*,'(A,ES15.7)') 't_val - tlow          = ', t_val - tlow
  write(*,'(A,L1)')      't_val <= tlow         = ', t_val <= tlow
  
  ! Show bit pattern
  write(*,'(A,Z8.8)') 't_val bits            = ', transfer(t_val, 0)
  write(*,'(A,Z8.8)') 'tlow bits             = ', transfer(tlow, 0)
  
end program check_t_single

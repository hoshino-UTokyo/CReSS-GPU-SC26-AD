program check_t_calc
  implicit none
  integer :: ni, nj, nk, i, j, k
  real, allocatable :: ptbr(:,:,:), ptp(:,:,:), pbr(:,:,:), pp(:,:,:)
  real :: pt_val, t_val
  real :: p0iv, rddvcp, tlow
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  
  ! Exact values from params.txt
  rddvcp = 2.858565747738E-01
  p0iv = 9.999999747379E-06
  tlow = 2.331600036621E+02
  
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
  
  write(*,'(A)') 'T calculation at (833, 331, 114):'
  write(*,'(A,ES20.12)') 'pt_val                = ', pt_val
  write(*,'(A,ES20.12)') 'pbr(i,j,k)            = ', pbr(i,j,k)
  write(*,'(A,ES20.12)') 'pp(i,j,k)             = ', pp(i,j,k)
  write(*,'(A,ES20.12)') 'pbr+pp                = ', pbr(i,j,k) + pp(i,j,k)
  write(*,'(A,ES20.12)') 'p0iv*(pbr+pp)         = ', p0iv * (pbr(i,j,k) + pp(i,j,k))
  write(*,'(A,ES20.12)') 'log(...)              = ', log(p0iv * (pbr(i,j,k) + pp(i,j,k)))
  write(*,'(A,ES20.12)') 'rddvcp*log(...)       = ', rddvcp * log(p0iv * (pbr(i,j,k) + pp(i,j,k)))
  write(*,'(A,ES20.12)') 'exp(...)              = ', exp(rddvcp * log(p0iv * (pbr(i,j,k) + pp(i,j,k))))
  write(*,'(A,ES20.12)') 't_val                 = ', t_val
  write(*,'(A,ES20.12)') 'tlow                  = ', tlow
  write(*,'(A,ES20.12)') 't_val - tlow          = ', t_val - tlow
  write(*,'(A,L1)')      't_val <= tlow         = ', t_val <= tlow
  
end program check_t_calc

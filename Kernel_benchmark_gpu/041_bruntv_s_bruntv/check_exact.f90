program check_exact
  implicit none
  integer :: ni, nj, nk, i, j, k
  real, allocatable :: ptbr(:,:,:), ptp(:,:,:), qv(:,:,:)
  real, allocatable :: pbr(:,:,:), pp(:,:,:)
  real :: pt_orig, pt_k, pt_km1, t_val, a_val, lhcpt
  ! Use EXACT values from params.txt
  real :: rd, cp, lv0, lf0, t0, tlow
  real :: cw, ci, cwmci, epsva, p0, p0iv, rddvcp
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  
  ! Exact values from params.txt
  rd = 2.870000000000E+02
  cp = 1.004000000000E+03
  lv0 = 2.500780000000E+06
  lf0 = 3.340000000000E+05
  t0 = 2.731600036621E+02
  tlow = 2.331600036621E+02
  cw = 4.218000000000E+03
  ci = 2.106000000000E+03
  epsva = 6.219999790192E-01
  p0 = 1.000000000000E+05
  cwmci = 2.112000000000E+03
  p0iv = 9.999999747379E-06
  rddvcp = 2.858565747738E-01
  
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(pbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(pp(0:ni+1, 0:nj+1, 1:nk))
  
  open(unit=10, file=trim(data_dir)//'/ptbr.bin', form='unformatted', access='stream')
  read(10) ptbr
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/ptp.bin', form='unformatted', access='stream')
  read(10) ptp
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/qv.bin', form='unformatted', access='stream')
  read(10) qv
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
  
  ! Compute for k=114
  pt_orig = ptbr(i,j,k) + ptp(i,j,k)
  t_val = pt_orig * exp(rddvcp * log(p0iv * (pbr(i,j,k) + pp(i,j,k))))
  a_val = lv0 * exp((0.167e0 + 3.67e-4 * t_val) * log(t0 / t_val))
  if (t_val <= tlow) a_val = a_val + (lf0 + cwmci * (t_val - t0))
  lhcpt = a_val / (cp * t_val)
  pt_k = pt_orig * exp(lhcpt * qv(i,j,k))
  
  write(*,'(A)') 'Recomputed with EXACT params at k=114:'
  write(*,'(A,ES15.7)') 'pt_orig       = ', pt_orig
  write(*,'(A,ES15.7)') 't_val         = ', t_val
  write(*,'(A,ES15.7)') 'a_val(initial)= ', lv0 * exp((0.167e0 + 3.67e-4 * t_val) * log(t0 / t_val))
  write(*,'(A,ES15.7)') 'a_val(final)  = ', a_val
  write(*,'(A,L1)')     't <= tlow     = ', t_val <= tlow
  write(*,'(A,ES15.7)') 'lhcpt         = ', lhcpt
  write(*,'(A,ES15.7)') 'qv            = ', qv(i,j,k)
  write(*,'(A,ES15.7)') 'lhcpt*qv      = ', lhcpt * qv(i,j,k)
  write(*,'(A,ES15.7)') 'exp(lhcpt*qv) = ', exp(lhcpt * qv(i,j,k))
  write(*,'(A,ES15.7)') 'pt_final      = ', pt_k
  
  ! Also check for k=113
  k = 113
  pt_orig = ptbr(i,j,k) + ptp(i,j,k)
  t_val = pt_orig * exp(rddvcp * log(p0iv * (pbr(i,j,k) + pp(i,j,k))))
  a_val = lv0 * exp((0.167e0 + 3.67e-4 * t_val) * log(t0 / t_val))
  if (t_val <= tlow) a_val = a_val + (lf0 + cwmci * (t_val - t0))
  lhcpt = a_val / (cp * t_val)
  pt_km1 = pt_orig * exp(lhcpt * qv(i,j,k))
  
  write(*,'(A)')
  write(*,'(A)') 'Recomputed with EXACT params at k=113:'
  write(*,'(A,ES15.7)') 'pt_final(k-1) = ', pt_km1
  
end program check_exact

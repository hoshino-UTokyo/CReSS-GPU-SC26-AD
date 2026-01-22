program check_cpu
  implicit none
  integer :: ni, nj, nk, i, j, k
  real, allocatable :: pt(:,:,:), a(:,:,:), ptv(:,:,:)
  real, allocatable :: ptbr(:,:,:), ptp(:,:,:), qv(:,:,:), qall(:,:,:)
  real, allocatable :: pbr(:,:,:), pp(:,:,:), jcb8w(:,:,:)
  real :: thresq, g, dziv, gdzv, gdzv05, rd, cp, lv0, lf0, t0, tlow
  real :: cw, ci, cwmci, epsav, epsva, p0, p0iv, rddvcp
  real :: t_val, a_val, lhcpt, nsq8w_val
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  
  ! Physical constants from params.txt
  g = 9.8
  dziv = 0.01
  rd = 287.0
  cp = 1004.0
  lv0 = 2500780.0
  lf0 = 334000.0
  t0 = 273.16
  tlow = 233.16
  cw = 4218.0
  ci = 2106.0
  epsav = 1.6077
  epsva = 0.622
  p0 = 100000.0
  thresq = 1.0e-12
  
  gdzv = g * dziv
  gdzv05 = 0.5 * g * dziv
  rddvcp = rd / cp
  cwmci = cw - ci
  p0iv = 1.0 / p0
  
  allocate(pt(0:ni+1, 0:nj+1, 1:nk))
  allocate(a(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptv(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qall(0:ni+1, 0:nj+1, 1:nk))
  allocate(pbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(pp(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8w(0:ni+1, 0:nj+1, 1:nk))
  
  open(unit=10, file=trim(data_dir)//'/ptbr.bin', form='unformatted', access='stream')
  read(10) ptbr
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/ptp.bin', form='unformatted', access='stream')
  read(10) ptp
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/qv.bin', form='unformatted', access='stream')
  read(10) qv
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/qall.bin', form='unformatted', access='stream')
  read(10) qall
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/pbr.bin', form='unformatted', access='stream')
  read(10) pbr
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/pp.bin', form='unformatted', access='stream')
  read(10) pp
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/jcb8w.bin', form='unformatted', access='stream')
  read(10) jcb8w
  close(10)
  
  i = 833
  j = 331
  k = 114
  
  ! Compute values following CPU algorithm
  pt(i,j,k) = ptbr(i,j,k) + ptp(i,j,k)
  pt(i,j,k-1) = ptbr(i,j,k-1) + ptp(i,j,k-1)
  
  ptv(i,j,k) = pt(i,j,k) * (1.0 + epsav * qv(i,j,k)) / (1.0 + qv(i,j,k))
  ptv(i,j,k-1) = pt(i,j,k-1) * (1.0 + epsav * qv(i,j,k-1)) / (1.0 + qv(i,j,k-1))
  
  ! Compute t and a for k
  t_val = pt(i,j,k) * exp(rddvcp * log(p0iv * (pbr(i,j,k) + pp(i,j,k))))
  a_val = lv0 * exp((0.167 + 3.67e-4 * t_val) * log(t0 / t_val))
  if (t_val <= tlow) a_val = a_val + (lf0 + cwmci * (t_val - t0))
  lhcpt = a_val / (cp * t_val)
  pt(i,j,k) = pt(i,j,k) * exp(lhcpt * qv(i,j,k))
  a(i,j,k) = a_val * qv(i,j,k) / (rd * t_val)
  a(i,j,k) = (1.0 + a(i,j,k)) / (1.0 + epsva * lhcpt * a(i,j,k))
  
  ! Compute for k-1
  t_val = pt(i,j,k-1) * exp(rddvcp * log(p0iv * (pbr(i,j,k-1) + pp(i,j,k-1))))
  a_val = lv0 * exp((0.167 + 3.67e-4 * t_val) * log(t0 / t_val))
  if (t_val <= tlow) a_val = a_val + (lf0 + cwmci * (t_val - t0))
  lhcpt = a_val / (cp * t_val)
  pt(i,j,k-1) = pt(i,j,k-1) * exp(lhcpt * qv(i,j,k-1))
  a(i,j,k-1) = a_val * qv(i,j,k-1) / (rd * t_val)
  a(i,j,k-1) = (1.0 + a(i,j,k-1)) / (1.0 + epsva * lhcpt * a(i,j,k-1))
  
  write(*,'(A)') 'Recomputed values at (833, 331, 114):'
  write(*,'(A,ES15.7)') 'pt(k)         = ', pt(i,j,k)
  write(*,'(A,ES15.7)') 'pt(k-1)       = ', pt(i,j,k-1)
  write(*,'(A,ES15.7)') 'ptv(k)        = ', ptv(i,j,k)
  write(*,'(A,ES15.7)') 'ptv(k-1)      = ', ptv(i,j,k-1)
  write(*,'(A,ES15.7)') 'a(k)          = ', a(i,j,k)
  write(*,'(A,ES15.7)') 'a(k-1)        = ', a(i,j,k-1)
  
  ! Compute nsq8w
  if (qall(i,j,k) > thresq) then
    nsq8w_val = gdzv05 * ((qall(i,j,k-1) - qall(i,j,k)) &
      + (pt(i,j,k) - pt(i,j,k-1)) * (a(i,j,k-1) + a(i,j,k)) &
      / (ptbr(i,j,k-1) + ptbr(i,j,k))) / jcb8w(i,j,k)
    write(*,'(A)') 'Using if branch (qall > thresq)'
  else
    nsq8w_val = gdzv * (ptv(i,j,k) - ptv(i,j,k-1)) &
      / (jcb8w(i,j,k) * (ptbr(i,j,k-1) + ptbr(i,j,k)))
    write(*,'(A)') 'Using else branch'
  end if
  
  write(*,'(A,ES15.7)') 'nsq8w_computed= ', nsq8w_val
  
end program check_cpu

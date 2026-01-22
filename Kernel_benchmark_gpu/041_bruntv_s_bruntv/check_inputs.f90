program check_inputs
  implicit none
  integer :: ni, nj, nk, i, j, k
  real, allocatable :: ptbr(:,:,:), ptp(:,:,:), qv(:,:,:)
  real :: pt_orig, ptv_val
  real :: epsav
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  epsav = 1.6077
  
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  
  open(unit=10, file=trim(data_dir)//'/ptbr.bin', form='unformatted', access='stream')
  read(10) ptbr
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/ptp.bin', form='unformatted', access='stream')
  read(10) ptp
  close(10)
  
  open(unit=10, file=trim(data_dir)//'/qv.bin', form='unformatted', access='stream')
  read(10) qv
  close(10)
  
  i = 833
  j = 331
  k = 114
  
  pt_orig = ptbr(i,j,k) + ptp(i,j,k)
  ptv_val = pt_orig * (1.0 + epsav * qv(i,j,k)) / (1.0 + qv(i,j,k))
  
  write(*,'(A)') 'Input values at (833, 331, 114):'
  write(*,'(A,ES15.7)') 'ptbr(k)       = ', ptbr(i,j,k)
  write(*,'(A,ES15.7)') 'ptp(k)        = ', ptp(i,j,k)
  write(*,'(A,ES15.7)') 'pt_orig       = ', pt_orig
  write(*,'(A,ES15.7)') 'qv(k)         = ', qv(i,j,k)
  write(*,'(A,ES15.7)') 'ptv_computed  = ', ptv_val
  
  write(*,'(A)')
  k = 113
  pt_orig = ptbr(i,j,k) + ptp(i,j,k)
  ptv_val = pt_orig * (1.0 + epsav * qv(i,j,k)) / (1.0 + qv(i,j,k))
  write(*,'(A)') 'Input values at (833, 331, 113):'
  write(*,'(A,ES15.7)') 'ptbr(k-1)     = ', ptbr(i,j,k)
  write(*,'(A,ES15.7)') 'ptp(k-1)      = ', ptp(i,j,k)
  write(*,'(A,ES15.7)') 'pt_orig(k-1)  = ', pt_orig
  write(*,'(A,ES15.7)') 'qv(k-1)       = ', qv(i,j,k)
  write(*,'(A,ES15.7)') 'ptv_comp(k-1) = ', ptv_val
  
end program check_inputs

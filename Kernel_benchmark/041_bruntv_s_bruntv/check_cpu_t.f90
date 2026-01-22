program check_cpu_t
  implicit none
  integer :: ni, nj, nk
  real, allocatable :: t_ref(:,:)
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  
  allocate(t_ref(0:ni+1, 0:nj+1))
  
  open(unit=10, file=trim(data_dir)//'/t_ref.bin', form='unformatted', access='stream')
  read(10) t_ref
  close(10)
  
  write(*,'(A,ES15.7)') 't_ref(833,331) = ', t_ref(833,331)
  write(*,'(A,ES15.7)') 'tlow           = ', 2.3316000E+02
  write(*,'(A,L1)') 't <= tlow      = ', t_ref(833,331) <= 2.3316000E+02
  
end program check_cpu_t

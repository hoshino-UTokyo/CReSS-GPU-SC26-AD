program debug_test2
  implicit none
  integer :: ni, nj, nk, ios
  real, allocatable :: nsq8w_ref(:,:,:)
  integer :: i, j, k, extreme_cnt, boundary_extreme
  real :: threshold
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  threshold = 1.0e30
  
  allocate(nsq8w_ref(0:ni+1, 0:nj+1, 1:nk))
  
  open(unit=10, file=trim(data_dir)//'/nsq8w_ref.bin', form='unformatted', access='stream', status='old')
  read(10) nsq8w_ref
  close(10)
  
  ! Count extreme values
  extreme_cnt = 0
  boundary_extreme = 0
  
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        if (abs(nsq8w_ref(i,j,k)) > threshold) then
          extreme_cnt = extreme_cnt + 1
          if (i == 0 .or. i == ni+1 .or. j == 0 .or. j == nj+1 .or. k == 1 .or. k == nk) then
            boundary_extreme = boundary_extreme + 1
          end if
        end if
      end do
    end do
  end do
  
  write(*,'(A,I12)') 'Total extreme values (>1e30): ', extreme_cnt
  write(*,'(A,I12)') 'Extreme values at boundaries: ', boundary_extreme
  
  ! Print values at inner points only
  write(*,'(A)') ''
  write(*,'(A)') 'Sample inner values (avoiding boundaries):'
  write(*,'(A,ES15.7)') 'nsq8w_ref(10,10,10) = ', nsq8w_ref(10,10,10)
  write(*,'(A,ES15.7)') 'nsq8w_ref(100,100,50) = ', nsq8w_ref(100,100,50)
  write(*,'(A,ES15.7)') 'nsq8w_ref(450,450,64) = ', nsq8w_ref(450,450,64)
  write(*,'(A,ES15.7)') 'nsq8w_ref(800,800,100) = ', nsq8w_ref(800,800,100)
  
  ! Check range excluding boundaries
  write(*,'(A)') ''
  write(*,'(A)') 'Range for inner points (1:ni-1, 1:nj-1, 2:nk-1):'
  write(*,'(A,ES15.7)') 'Min: ', minval(nsq8w_ref(1:ni-1,1:nj-1,2:nk-1))
  write(*,'(A,ES15.7)') 'Max: ', maxval(nsq8w_ref(1:ni-1,1:nj-1,2:nk-1))
  
end program debug_test2

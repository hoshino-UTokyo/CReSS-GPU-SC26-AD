program debug_test
  implicit none
  integer :: ni, nj, nk, ios
  real, allocatable :: nsq8w_ref(:,:,:), nsq8w(:,:,:)
  real :: ref_val, val
  integer :: i, j, k, err_cnt
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  
  allocate(nsq8w_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(nsq8w(0:ni+1, 0:nj+1, 1:nk))
  
  ! Read reference
  open(unit=10, file=trim(data_dir)//'/nsq8w_ref.bin', form='unformatted', access='stream', status='old')
  read(10) nsq8w_ref
  close(10)
  
  ! Print some reference values
  write(*,'(A)') 'Reference nsq8w values:'
  write(*,'(A,ES15.7)') 'nsq8w_ref(1,1,2) = ', nsq8w_ref(1,1,2)
  write(*,'(A,ES15.7)') 'nsq8w_ref(100,100,50) = ', nsq8w_ref(100,100,50)
  write(*,'(A,ES15.7)') 'nsq8w_ref(450,450,64) = ', nsq8w_ref(450,450,64)
  
  ! Check for non-zero values
  err_cnt = 0
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        if (abs(nsq8w_ref(i,j,k)) > 1.0e-10) err_cnt = err_cnt + 1
      end do
    end do
  end do
  write(*,'(A,I12)') 'Non-zero reference values: ', err_cnt
  
  ! Check range
  write(*,'(A,ES15.7)') 'Min ref: ', minval(nsq8w_ref)
  write(*,'(A,ES15.7)') 'Max ref: ', maxval(nsq8w_ref)
  
end program debug_test

program debug_test3
  implicit none
  integer :: ni, nj, nk, ios
  real, allocatable :: nsq8w_ref(:,:,:)
  integer :: i, j, k, cnt_nan, cnt_inf
  real :: val
  character(len=256) :: data_dir
  
  data_dir = './data'
  ni = 899
  nj = 899
  nk = 128
  
  allocate(nsq8w_ref(0:ni+1, 0:nj+1, 1:nk))
  
  open(unit=10, file=trim(data_dir)//'/nsq8w_ref.bin', form='unformatted', access='stream', status='old')
  read(10) nsq8w_ref
  close(10)
  
  ! Check for NaN and Inf
  cnt_nan = 0
  cnt_inf = 0
  
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        val = nsq8w_ref(i,j,k)
        if (val /= val) cnt_nan = cnt_nan + 1  ! NaN check
        if (abs(val) > 3.4e38) cnt_inf = cnt_inf + 1
      end do
    end do
  end do
  
  write(*,'(A,I12)') 'NaN count: ', cnt_nan
  write(*,'(A,I12)') 'Very large (>3.4e38) count: ', cnt_inf
  
  ! Find first extreme value location
  write(*,'(A)') ''
  write(*,'(A)') 'First few extreme values:'
  cnt_nan = 0
  outer: do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        if (abs(nsq8w_ref(i,j,k)) > 1.0e30) then
          write(*,'(A,I4,A,I4,A,I4,A,ES15.7)') 'i=',i,', j=',j,', k=',k,': ',nsq8w_ref(i,j,k)
          cnt_nan = cnt_nan + 1
          if (cnt_nan >= 10) exit outer
        end if
      end do
    end do
  end do outer
  
end program debug_test3

program generate_data_vbcv
  implicit none
  integer, parameter :: ni = 100, nj = 100, nk = 64
  real, allocatable :: vf(:,:,:), vf_ref(:,:,:)
  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  allocate(vf(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf_ref(0:ni+1, 0:nj+1, 1:nk))

  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        vf(i,j,k) = real(i + j*100 + k*10000) * 0.001
      end do
    end do
  end do

  open(unit=10, file=trim(data_dir)//'/vf_input.bin', status='replace', access='stream', form='unformatted')
  write(10) vf
  close(10)

  vf_ref = vf
  do j = 1, nj
    do i = 1, ni-1
      vf_ref(i,j,1) = vf_ref(i,j,2)
      vf_ref(i,j,nk-1) = vf_ref(i,j,nk-2)
    end do
  end do

  open(unit=10, file=trim(data_dir)//'/vf_ref.bin', status='replace', access='stream', form='unformatted')
  write(10) vf_ref
  close(10)

  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  close(10)

  deallocate(vf, vf_ref)
  write(*,'(A)') 'Data generation complete for vbcv!'
end program generate_data_vbcv

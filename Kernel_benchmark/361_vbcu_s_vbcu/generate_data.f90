!***********************************************************************
! Data Generator for vbcu Kernel Benchmark
!***********************************************************************
program generate_data_vbcu
  implicit none

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  real, allocatable :: uf(:,:,:)
  real, allocatable :: uf_ref(:,:,:)
  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for vbcu kernel'
  write(*,'(A)') '=================================================='

  allocate(uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Initialize with varying values
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        uf(i,j,k) = real(i + j*100 + k*10000) * 0.001
      end do
    end do
  end do

  ! Save input
  open(unit=10, file=trim(data_dir)//'/uf_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) uf
  close(10)

  ! Compute reference (apply boundary conditions)
  uf_ref = uf
  do j = 1, nj-1
    do i = 1, ni
      uf_ref(i,j,1) = uf_ref(i,j,2)
      uf_ref(i,j,nk-1) = uf_ref(i,j,nk-2)
    end do
  end do

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/uf_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) uf_ref
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  close(10)

  deallocate(uf, uf_ref)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_vbcu

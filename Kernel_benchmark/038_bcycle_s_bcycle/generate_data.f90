program generate_data_bcycle
  implicit none
  integer, parameter :: ni = 100, nj = 100, kmax = 64
  ! Cyclic boundary conditions
  integer, parameter :: wbc = -1, ebc = -1, sbc = -1, nbc = -1
  integer, parameter :: nisub = 1, njsub = 1
  integer, parameter :: iwsnd = 2, iwrcv = 0, iesnd = ni-1, iercv = ni+1
  integer, parameter :: jssnd = 2, jsrcv = 0, jnsnd = nj-1, jnrcv = nj+1

  real, allocatable :: var(:,:,:), var_ref(:,:,:)
  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  allocate(var(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var_ref(0:ni+1, 0:nj+1, 1:kmax))

  ! Initialize
  do k = 1, kmax
    do j = 0, nj+1
      do i = 0, ni+1
        var(i,j,k) = real(i + j*100 + k*10000) * 0.001
      end do
    end do
  end do

  open(unit=10, file=trim(data_dir)//'/var_input.bin', status='replace', access='stream', form='unformatted')
  write(10) var
  close(10)

  ! Compute reference
  var_ref = var
  if (nisub == 1 .and. wbc == -1 .and. ebc == -1) then
    do k = 1, kmax
      do j = 0, nj+1
        var_ref(iwrcv,j,k) = var_ref(iesnd,j,k)
        var_ref(iercv,j,k) = var_ref(iwsnd,j,k)
      end do
    end do
  end if
  if (njsub == 1 .and. sbc == -1 .and. nbc == -1) then
    do k = 1, kmax
      do i = 0, ni+1
        var_ref(i,jsrcv,k) = var_ref(i,jnsnd,k)
        var_ref(i,jnrcv,k) = var_ref(i,jssnd,k)
      end do
    end do
  end if

  open(unit=10, file=trim(data_dir)//'/var_ref.bin', status='replace', access='stream', form='unformatted')
  write(10) var_ref
  close(10)

  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') kmax
  write(10, '(I12)') wbc
  write(10, '(I12)') ebc
  write(10, '(I12)') sbc
  write(10, '(I12)') nbc
  write(10, '(I12)') nisub
  write(10, '(I12)') njsub
  write(10, '(I12)') iwsnd
  write(10, '(I12)') iwrcv
  write(10, '(I12)') iesnd
  write(10, '(I12)') iercv
  write(10, '(I12)') jssnd
  write(10, '(I12)') jsrcv
  write(10, '(I12)') jnsnd
  write(10, '(I12)') jnrcv
  close(10)

  deallocate(var, var_ref)
  write(*,'(A)') 'Data generation complete for bcycle!'
end program generate_data_bcycle

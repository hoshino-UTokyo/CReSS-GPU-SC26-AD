!***********************************************************************
! Data Generator for vbcwc Kernel Benchmark
!***********************************************************************
program generate_data_vbcwc
  implicit none

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  ! Boundary condition options
  integer, parameter :: bbc = 2  ! Bottom BC (2=rigid lid)
  integer, parameter :: tbc = 2  ! Top BC (2=rigid lid)

  ! Arrays
  real, allocatable :: wc(:,:,:)
  real, allocatable :: wc_ref(:,:,:)

  integer :: i, j, k, nkm1, nkm2
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for vbcwc kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' bbc: ', bbc
  write(*,'(A,I6)') ' tbc: ', tbc

  ! Allocate
  allocate(wc(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Initialize wc with some values
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        wc(i,j,k) = sin(real(i)*0.1) * cos(real(j)*0.1) * real(k)/real(nk)
      end do
    end do
  end do

  ! Write input array
  write(*,'(A)') ' Writing input arrays...'

  open(unit=10, file=trim(data_dir)//'/wc_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) wc
  close(10)

  ! Compute reference output (same logic as kernel)
  write(*,'(A)') ' Computing reference output...'
  wc_ref = wc

  nkm1 = nk - 1
  nkm2 = nk - 2

  ! Bottom boundary conditions
  if (bbc == 2) then
    do j = 1, nj-1
      do i = 1, ni-1
        wc_ref(i,j,1) = -wc_ref(i,j,3)
        wc_ref(i,j,2) = 0.e0
      end do
    end do
  else if (bbc == 3) then
    do j = 1, nj-1
      do i = 1, ni-1
        wc_ref(i,j,1) = wc_ref(i,j,2)
      end do
    end do
  end if

  ! Top boundary conditions
  if (tbc == 2) then
    do j = 1, nj-1
      do i = 1, ni-1
        wc_ref(i,j,nk) = -wc_ref(i,j,nkm2)
        wc_ref(i,j,nkm1) = 0.e0
      end do
    end do
  else if (tbc >= 3) then
    do j = 1, nj-1
      do i = 1, ni-1
        wc_ref(i,j,nk) = wc_ref(i,j,nkm1)
      end do
    end do
  end if

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/wc_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) wc_ref
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(I12)') bbc
  write(10, '(I12)') tbc
  close(10)

  deallocate(wc, wc_ref)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_vbcwc

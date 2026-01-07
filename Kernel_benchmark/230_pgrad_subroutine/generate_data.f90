!***********************************************************************
! Data Generator for pgrad Kernel Benchmark
!***********************************************************************
program generate_data_pgrad
  implicit none

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  ! Options (simple case: no terrain, no map scale factor)
  integer, parameter :: trnopt = 0
  integer, parameter :: mpopt = 0
  integer, parameter :: mfcopt = 0
  integer, parameter :: divopt = 0

  ! Grid parameters
  real, parameter :: dx = 1000.0
  real, parameter :: dy = 1000.0
  real, parameter :: dz = 500.0
  real, parameter :: dxiv = 1.0/dx
  real, parameter :: dyiv = 1.0/dy
  real, parameter :: dziv = 1.0/dz
  real, parameter :: dts = 0.1
  real, parameter :: divndc = 0.001
  real :: dziv25, divch, divcv

  ! Arrays
  real, allocatable :: j31(:,:,:), j32(:,:,:), jcb(:,:,:)
  real, allocatable :: mf8u(:,:), mf8v(:,:)
  real, allocatable :: pp(:,:,:)
  real, allocatable :: upg(:,:,:), vpg(:,:,:), wpg(:,:,:)
  real, allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:)

  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for pgrad kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk

  dziv25 = 0.25e0 * dziv
  divch = divndc / dts * dx * dy
  divcv = divndc / dts * dz * dz

  ! Allocate
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf8u(0:ni+1, 0:nj+1))
  allocate(mf8v(0:ni+1, 0:nj+1))
  allocate(pp(0:ni+1, 0:nj+1, 1:nk))
  allocate(upg(0:ni+1, 0:nj+1, 1:nk))
  allocate(vpg(0:ni+1, 0:nj+1, 1:nk))
  allocate(wpg(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))

  ! Initialize Jacobian (flat terrain: j31=j32=0, jcb=1)
  j31 = 0.0
  j32 = 0.0
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        jcb(i,j,k) = 1.0
      end do
    end do
  end do

  ! Initialize map scale factors
  mf8u = 1.0
  mf8v = 1.0

  ! Initialize pressure perturbation
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        pp(i,j,k) = 100.0 * sin(real(i)*0.1) * cos(real(j)*0.1) * (1.0 - real(k)/real(nk))
      end do
    end do
  end do

  ! Write input arrays
  write(*,'(A)') ' Writing input arrays...'

  open(unit=10, file=trim(data_dir)//'/j31.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) j31
  close(10)

  open(unit=10, file=trim(data_dir)//'/j32.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) j32
  close(10)

  open(unit=10, file=trim(data_dir)//'/jcb.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) jcb
  close(10)

  open(unit=10, file=trim(data_dir)//'/mf8u.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) mf8u
  close(10)

  open(unit=10, file=trim(data_dir)//'/mf8v.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) mf8v
  close(10)

  open(unit=10, file=trim(data_dir)//'/pp.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) pp
  close(10)

  ! Compute reference output
  write(*,'(A)') ' Computing reference output...'
  upg = 0.0
  vpg = 0.0
  wpg = 0.0
  tmp1 = 0.0
  tmp2 = 0.0
  tmp3 = 0.0

  ! divopt=0 case
  do k=1,nk-1
    do j=1,nj-1
      do i=1,ni-1
        tmp3(i,j,k)=pp(i,j,k)
      end do
    end do
  end do

  ! Vertical gradient
  do k=2,nk-1
    do j=2,nj-2
      do i=2,ni-2
        wpg(i,j,k)=(tmp3(i,j,k-1)-tmp3(i,j,k))*dziv
      end do
    end do
  end do

  ! trnopt=0 case: multiply by jcb
  do k=2,nk-2
    do j=1,nj-1
      do i=1,ni-1
        tmp3(i,j,k)=jcb(i,j,k)*tmp3(i,j,k)
      end do
    end do
  end do

  ! mfcopt=0 case: horizontal gradients
  do k=2,nk-2
    do j=2,nj-2
      do i=2,ni-1
        upg(i,j,k)=(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv
      end do
    end do
    do j=2,nj-1
      do i=2,ni-2
        vpg(i,j,k)=(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv
      end do
    end do
  end do

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/upg_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) upg
  close(10)

  open(unit=10, file=trim(data_dir)//'/vpg_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) vpg
  close(10)

  open(unit=10, file=trim(data_dir)//'/wpg_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) wpg
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(I12)') trnopt
  write(10, '(I12)') mpopt
  write(10, '(I12)') mfcopt
  write(10, '(I12)') divopt
  write(10, '(ES20.12)') dx
  write(10, '(ES20.12)') dy
  write(10, '(ES20.12)') dz
  write(10, '(ES20.12)') dxiv
  write(10, '(ES20.12)') dyiv
  write(10, '(ES20.12)') dziv
  write(10, '(ES20.12)') dts
  close(10)

  deallocate(j31, j32, jcb, mf8u, mf8v, pp, upg, vpg, wpg)
  deallocate(tmp1, tmp2, tmp3)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_pgrad

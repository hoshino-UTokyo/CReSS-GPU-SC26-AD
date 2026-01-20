!***********************************************************************
! Data Generator for disptke Kernel Benchmark
!***********************************************************************
program generate_data_disptke
  implicit none

  real, parameter :: oned3 = 1.0 / 3.0   ! One third
  real, parameter :: eps = 1.0e-20       ! Small number

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  real, parameter :: dx = 1000.0   ! Grid distance in x
  real, parameter :: dy = 1000.0   ! Grid distance in y
  real, parameter :: dz = 100.0    ! Grid distance in z
  integer, parameter :: mpopt = 0
  integer, parameter :: mfcopt = 0
  integer, parameter :: isoopt = 1  ! Isotropic case

  real :: dz05, ds308, ln

  ! Arrays
  real, allocatable :: jcb(:,:,:)
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: rst(:,:,:)
  real, allocatable :: priv(:,:,:)
  real, allocatable :: tke(:,:,:)
  real, allocatable :: tkefrc(:,:,:)
  real, allocatable :: tkefrc_ref(:,:,:)

  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  dz05 = 0.5 * dz
  ds308 = 0.125 * dx * dy * dz

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for disptke kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dx: ', dx
  write(*,'(A,ES12.4)') ' dy: ', dy
  write(*,'(A,ES12.4)') ' dz: ', dz
  write(*,'(A,I6)') ' isoopt: ', isoopt

  ! Allocate
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(priv(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Initialize arrays
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        jcb(i,j,k) = 1.0 + 0.05 * sin(real(i)*0.1) * cos(real(j)*0.1)
        rst(i,j,k) = 1.2 * (1.0 - 0.2 * real(k)/real(nk)) * jcb(i,j,k)
        priv(i,j,k) = 1.0 + 0.5 * (1.0 - real(k)/real(nk))  ! > 1
        tke(i,j,k) = 1.0 + 0.5 * sin(real(i)*0.05) * cos(real(j)*0.05)  ! > 0
        tkefrc(i,j,k) = 0.1 * sin(real(i+j)*0.1) * (real(k)/real(nk))
      end do
    end do
  end do

  ! Initialize rmf
  do j = 0, nj+1
    do i = 0, ni+1
      rmf(i,j,1) = 1.0
      rmf(i,j,2) = 1.0
      rmf(i,j,3) = 1.0
      rmf(i,j,4) = 1.0
    end do
  end do

  ! Write input arrays
  write(*,'(A)') ' Writing input arrays...'

  open(unit=10, file=trim(data_dir)//'/jcb.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) jcb
  close(10)

  open(unit=10, file=trim(data_dir)//'/rmf.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rmf
  close(10)

  open(unit=10, file=trim(data_dir)//'/rst.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rst
  close(10)

  open(unit=10, file=trim(data_dir)//'/priv.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) priv
  close(10)

  open(unit=10, file=trim(data_dir)//'/tke.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) tke
  close(10)

  open(unit=10, file=trim(data_dir)//'/tkefrc_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) tkefrc
  close(10)

  ! Compute reference output
  write(*,'(A)') ' Computing reference output...'
  tkefrc_ref = tkefrc

  ! Isotropic case with mfcopt=0
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        ln = (priv(i,j,k) - 1.0) * exp(oned3 * log(ds308 * jcb(i,j,k))) + eps
        tkefrc_ref(i,j,k) = tkefrc_ref(i,j,k) - (0.37 * priv(i,j,k) - 0.18) &
             * rst(i,j,k) * tke(i,j,k) * sqrt(tke(i,j,k)) / ln
      end do
    end do
  end do

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/tkefrc_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) tkefrc_ref
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(ES20.12)') dx
  write(10, '(ES20.12)') dy
  write(10, '(ES20.12)') dz
  write(10, '(I12)') mpopt
  write(10, '(I12)') mfcopt
  write(10, '(I12)') isoopt
  close(10)

  deallocate(jcb, rmf, rst, priv, tke, tkefrc, tkefrc_ref)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_disptke

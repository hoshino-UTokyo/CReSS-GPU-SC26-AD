!***********************************************************************
! Data Generator for stepuv Kernel Benchmark
!***********************************************************************
program generate_data_stepuv
  implicit none

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  real, parameter :: dts = 0.1  ! Small time step

  ! Arrays
  real, allocatable :: rst8u(:,:,:), rst8v(:,:,:)
  real, allocatable :: ufrc(:,:,:), vfrc(:,:,:)
  real, allocatable :: usml(:,:,:), vsml(:,:,:)
  real, allocatable :: uf(:,:,:), vf(:,:,:)
  real, allocatable :: uf_ref(:,:,:), vf_ref(:,:,:)

  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for stepuv kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk

  ! Allocate
  allocate(rst8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(usml(0:ni+1, 0:nj+1, 1:nk))
  allocate(vsml(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Initialize density (always positive)
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        rst8u(i,j,k) = 1.2 * (1.0 - 0.3 * real(k)/real(nk))
        rst8v(i,j,k) = 1.2 * (1.0 - 0.3 * real(k)/real(nk))
      end do
    end do
  end do

  ! Initialize forcing terms
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        ufrc(i,j,k) = 0.1 * sin(real(i)*0.1) * cos(real(j)*0.1) * (1.0 - real(k)/real(nk))
        vfrc(i,j,k) = 0.1 * cos(real(i)*0.1) * sin(real(j)*0.1) * (1.0 - real(k)/real(nk))
        usml(i,j,k) = 0.01 * sin(real(i+j)*0.05)
        vsml(i,j,k) = 0.01 * cos(real(i+j)*0.05)
      end do
    end do
  end do

  ! Initialize velocity
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        uf(i,j,k) = 10.0 * sin(real(i)*0.05) * (real(k)/real(nk))
        vf(i,j,k) = 5.0 * cos(real(j)*0.05) * (real(k)/real(nk))
      end do
    end do
  end do

  ! Write input arrays
  write(*,'(A)') ' Writing input arrays...'

  open(unit=10, file=trim(data_dir)//'/rst8u.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rst8u
  close(10)

  open(unit=10, file=trim(data_dir)//'/rst8v.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rst8v
  close(10)

  open(unit=10, file=trim(data_dir)//'/ufrc.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) ufrc
  close(10)

  open(unit=10, file=trim(data_dir)//'/vfrc.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) vfrc
  close(10)

  open(unit=10, file=trim(data_dir)//'/usml.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) usml
  close(10)

  open(unit=10, file=trim(data_dir)//'/vsml.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) vsml
  close(10)

  open(unit=10, file=trim(data_dir)//'/uf_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) uf
  close(10)

  open(unit=10, file=trim(data_dir)//'/vf_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) vf
  close(10)

  ! Compute reference output
  write(*,'(A)') ' Computing reference output...'
  uf_ref = uf
  vf_ref = vf

  do k=2,nk-2
    do j=2,nj-2
      do i=2,ni-1
        uf_ref(i,j,k)=uf_ref(i,j,k)+dts*(ufrc(i,j,k)+usml(i,j,k))/rst8u(i,j,k)
      end do
    end do
    do j=2,nj-1
      do i=2,ni-2
        vf_ref(i,j,k)=vf_ref(i,j,k)+dts*(vfrc(i,j,k)+vsml(i,j,k))/rst8v(i,j,k)
      end do
    end do
  end do

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/uf_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) uf_ref
  close(10)

  open(unit=10, file=trim(data_dir)//'/vf_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) vf_ref
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(ES20.12)') dts
  close(10)

  deallocate(rst8u, rst8v, ufrc, vfrc, usml, vsml, uf, vf, uf_ref, vf_ref)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_stepuv

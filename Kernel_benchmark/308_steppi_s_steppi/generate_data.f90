!***********************************************************************
! Data Generator for steppi Kernel Benchmark
!***********************************************************************
program generate_data_steppi
  implicit none

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  real, parameter :: dts = 0.1  ! Small time step

  ! Arrays
  real, allocatable :: jcb(:,:,:)
  real, allocatable :: fp(:,:,:)
  real, allocatable :: ppf(:,:,:)
  real, allocatable :: ppf_ref(:,:,:)

  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for steppi kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dts: ', dts

  ! Allocate
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(fp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ppf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ppf_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Initialize jcb (Jacobian - always positive)
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        jcb(i,j,k) = 1.0 + 0.1 * sin(real(i)*0.1) * cos(real(j)*0.1)
      end do
    end do
  end do

  ! Initialize fp (forcing term)
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        fp(i,j,k) = 100.0 * sin(real(i)*0.05) * cos(real(j)*0.05) * (1.0 - real(k)/real(nk))
      end do
    end do
  end do

  ! Initialize ppf (pressure perturbation at future - initial value)
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        ppf(i,j,k) = 50.0 * cos(real(i)*0.03) * sin(real(j)*0.03) * (real(k)/real(nk))
      end do
    end do
  end do

  ! Write input arrays
  write(*,'(A)') ' Writing input arrays...'

  open(unit=10, file=trim(data_dir)//'/jcb.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) jcb
  close(10)

  open(unit=10, file=trim(data_dir)//'/fp.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) fp
  close(10)

  open(unit=10, file=trim(data_dir)//'/ppf_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) ppf
  close(10)

  ! Compute reference output
  write(*,'(A)') ' Computing reference output...'
  ppf_ref = ppf

  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        ppf_ref(i,j,k) = ppf_ref(i,j,k) + dts * fp(i,j,k) / jcb(i,j,k)
      end do
    end do
  end do

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/ppf_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) ppf_ref
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(ES20.12)') dts
  close(10)

  deallocate(jcb, fp, ppf, ppf_ref)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_steppi

!***********************************************************************
! Data Generator for steps Kernel Benchmark (simplified dry case)
!***********************************************************************
program generate_data_steps
  implicit none

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  real, parameter :: dtb = 1.0  ! Time step
  integer, parameter :: advopt = 4  ! Advection option

  ! Arrays
  real, allocatable :: rst(:,:,:)
  real, allocatable :: ptpp(:,:,:)
  real, allocatable :: ptfrc(:,:,:)
  real, allocatable :: dtdrst(:,:,:)
  real, allocatable :: ptpf_ref(:,:,:)

  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for steps kernel (simplified)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dtb: ', dtb
  write(*,'(A,I6)') ' advopt: ', advopt

  ! Allocate
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(dtdrst(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpf_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Initialize rst (base state density x Jacobian - always positive)
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        rst(i,j,k) = 1.2 * (1.0 - 0.3 * real(k)/real(nk)) + 0.1
      end do
    end do
  end do

  ! Initialize ptpp (potential temperature perturbation at past)
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        ptpp(i,j,k) = 2.0 * sin(real(i)*0.05) * cos(real(j)*0.05) * (1.0 - real(k)/real(nk))
      end do
    end do
  end do

  ! Initialize ptfrc (potential temperature forcing)
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        ptfrc(i,j,k) = 0.01 * sin(real(i+j)*0.1) * (real(k)/real(nk))
      end do
    end do
  end do

  ! Write input arrays
  write(*,'(A)') ' Writing input arrays...'

  open(unit=10, file=trim(data_dir)//'/rst.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rst
  close(10)

  open(unit=10, file=trim(data_dir)//'/ptpp.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) ptpp
  close(10)

  open(unit=10, file=trim(data_dir)//'/ptfrc.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) ptfrc
  close(10)

  ! Compute reference output
  write(*,'(A)') ' Computing reference output...'
  dtdrst = 0.0
  ptpf_ref = 0.0

  ! Compute dtdrst (advopt > 3 case)
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        dtdrst(i,j,k) = dtb / rst(i,j,k)
      end do
    end do
  end do

  ! Compute ptpf
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        ptpf_ref(i,j,k) = ptpp(i,j,k) + ptfrc(i,j,k) * dtdrst(i,j,k)
      end do
    end do
  end do

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/ptpf_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) ptpf_ref
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(ES20.12)') dtb
  write(10, '(I12)') advopt
  close(10)

  deallocate(rst, ptpp, ptfrc, dtdrst, ptpf_ref)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_steps

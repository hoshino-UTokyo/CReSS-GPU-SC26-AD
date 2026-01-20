!***********************************************************************
! Data Generator for stepwi Kernel Benchmark
!***********************************************************************
program generate_data_stepwi
  implicit none

  real, parameter :: g = 9.80665  ! Gravitational acceleration

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  real, parameter :: dts = 0.1       ! Small time step
  real, parameter :: dziv = 0.01     ! Inverse of dz
  real, parameter :: weicoe = 0.5    ! Weighting coefficient
  integer, parameter :: buyopt = 0   ! Buoyancy option

  ! Arrays
  real, allocatable :: rst8w(:,:,:), rbr(:,:,:), rcsq(:,:,:)
  real, allocatable :: jcb(:,:,:), rst(:,:,:)
  real, allocatable :: fw(:,:,:), wf(:,:,:), wc(:,:,:)
  real, allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:)

  real :: g05, sbsqzi, sbsqg5, a, mm, nn
  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for stepwi kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dts:    ', dts
  write(*,'(A,ES12.4)') ' dziv:   ', dziv
  write(*,'(A,ES12.4)') ' weicoe: ', weicoe
  write(*,'(A,I6)') ' buyopt: ', buyopt

  ! Allocate
  allocate(rst8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rcsq(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(fw(0:ni+1, 0:nj+1, 1:nk))
  allocate(wf(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))

  ! Initialize arrays
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        rst8w(i,j,k) = 1.2 * (1.0 - 0.2 * real(k)/real(nk)) + 0.1
        rbr(i,j,k) = 1.2 * (1.0 - 0.2 * real(k)/real(nk)) + 0.1
        rcsq(i,j,k) = 400.0 * rbr(i,j,k)  ! ~sound speed squared
        jcb(i,j,k) = 1.0 + 0.05 * sin(real(i)*0.1) * cos(real(j)*0.1)
        rst(i,j,k) = rbr(i,j,k) * jcb(i,j,k)
        fw(i,j,k) = 0.1 * sin(real(i)*0.05) * cos(real(j)*0.05) * sin(real(k)*0.1)
        wf(i,j,k) = 0.5 * cos(real(i)*0.03) * sin(real(j)*0.03) * cos(real(k)*0.1)
      end do
    end do
  end do

  ! Write input arrays
  write(*,'(A)') ' Writing input arrays...'

  open(unit=10, file=trim(data_dir)//'/rst8w.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rst8w
  close(10)

  open(unit=10, file=trim(data_dir)//'/rbr.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rbr
  close(10)

  open(unit=10, file=trim(data_dir)//'/rcsq.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rcsq
  close(10)

  open(unit=10, file=trim(data_dir)//'/jcb.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) jcb
  close(10)

  open(unit=10, file=trim(data_dir)//'/rst.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rst
  close(10)

  open(unit=10, file=trim(data_dir)//'/fw_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) fw
  close(10)

  open(unit=10, file=trim(data_dir)//'/wf_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) wf
  close(10)

  ! Compute reference output
  write(*,'(A)') ' Computing reference output...'
  wc = 0.0
  tmp1 = 0.0
  tmp2 = 0.0
  tmp3 = 0.0

  g05 = 0.5 * g
  sbsqzi = dts * dts * weicoe * weicoe * dziv
  sbsqg5 = 0.5 * g * dts * dts * weicoe * weicoe

  ! Step 1: Update wf
  do k = 3, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        wf(i,j,k) = wf(i,j,k) + dts * fw(i,j,k) / rst8w(i,j,k)
      end do
    end do
  end do

  ! Step 2: Compute intermediate values
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        tmp1(i,j,k) = g05 * rbr(i,j,k) / rcsq(i,j,k)
        tmp2(i,j,k) = dziv / jcb(i,j,k)
      end do
    end do
  end do

  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        fw(i,j,k) = tmp1(i,j,k) + tmp2(i,j,k)
        wc(i,j,k) = tmp1(i,j,k) - tmp2(i,j,k)
      end do
    end do
  end do

  ! Step 3: Compute tridiagonal matrix coefficients (buyopt=0)
  do k = 3, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        a = sbsqzi / rst8w(i,j,k)
        mm = a * rcsq(i,j,k-1)
        nn = a * rcsq(i,j,k)
        tmp1(i,j,k) = -mm * fw(i,j,k-1)
        tmp2(i,j,k) = 1.0 + (nn * fw(i,j,k) - mm * wc(i,j,k-1))
        tmp3(i,j,k) = nn * wc(i,j,k)
      end do
    end do
  end do

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/tmp1_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) tmp1
  close(10)

  open(unit=10, file=trim(data_dir)//'/tmp2_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) tmp2
  close(10)

  open(unit=10, file=trim(data_dir)//'/tmp3_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) tmp3
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(ES20.12)') dts
  write(10, '(ES20.12)') dziv
  write(10, '(ES20.12)') weicoe
  write(10, '(I12)') buyopt
  close(10)

  deallocate(rst8w, rbr, rcsq, jcb, rst, fw, wf, wc, tmp1, tmp2, tmp3)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_stepwi

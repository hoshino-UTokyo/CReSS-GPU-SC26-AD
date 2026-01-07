!***********************************************************************
! Data Generator for phy2cnt Kernel Benchmark
!***********************************************************************
program generate_data_phy2cnt
  implicit none

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  ! Physics options (trnopt=1 for terrain case)
  integer, parameter :: sthopt = 1
  integer, parameter :: trnopt = 1
  integer, parameter :: mpopt = 0
  integer, parameter :: mfcopt = 0

  ! Arrays
  real, allocatable :: j31(:,:,:), j32(:,:,:), jcb8w(:,:,:)
  real, allocatable :: mf(:,:)
  real, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
  real, allocatable :: wc(:,:,:)
  real, allocatable :: mf25(:,:), j31u2(:,:,:), j32v2(:,:,:)

  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for phy2cnt kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6,A,I6,A,I6)') ' Options: sthopt=', sthopt, ', trnopt=', trnopt, &
       ', mpopt=', mpopt, ', mfcopt=', mfcopt

  ! Allocate
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(w(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf25(0:ni+1, 0:nj+1))
  allocate(j31u2(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32v2(0:ni+1, 0:nj+1, 1:nk))

  ! Initialize Jacobian arrays
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        j31(i,j,k) = 0.01 * sin(real(i)*0.1) * (1.0 - real(k)/real(nk))
        j32(i,j,k) = 0.01 * sin(real(j)*0.1) * (1.0 - real(k)/real(nk))
        jcb8w(i,j,k) = 1.0 + 0.001 * real(k)
      end do
    end do
  end do

  ! Initialize map scale factor
  do j = 0, nj+1
    do i = 0, ni+1
      mf(i,j) = 1.0 + 0.0001 * real(i+j)
    end do
  end do

  ! Initialize velocity fields
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        u(i,j,k) = 10.0 * sin(real(i)*0.05) * cos(real(j)*0.05) * (real(k)/real(nk))
        v(i,j,k) = 5.0 * cos(real(i)*0.05) * sin(real(j)*0.05) * (real(k)/real(nk))
        w(i,j,k) = 0.5 * sin(real(i+j)*0.1) * (real(k)/real(nk))
      end do
    end do
  end do

  ! Save input arrays
  write(*,'(A)') ' Writing input arrays...'

  open(unit=10, file=trim(data_dir)//'/j31.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) j31
  close(10)

  open(unit=10, file=trim(data_dir)//'/j32.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) j32
  close(10)

  open(unit=10, file=trim(data_dir)//'/jcb8w.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) jcb8w
  close(10)

  open(unit=10, file=trim(data_dir)//'/mf.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) mf
  close(10)

  open(unit=10, file=trim(data_dir)//'/u.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) u
  close(10)

  open(unit=10, file=trim(data_dir)//'/v.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) v
  close(10)

  open(unit=10, file=trim(data_dir)//'/w.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) w
  close(10)

  ! Compute reference output
  write(*,'(A)') ' Computing reference output...'
  wc = 0.0
  j31u2 = 0.0
  j32v2 = 0.0
  mf25 = 0.0

  if(trnopt.eq.0) then
    if(sthopt.eq.0) then
      do k=2,nk-1
        do j=1,nj-1
          do i=1,ni-1
            wc(i,j,k)=w(i,j,k)
          end do
        end do
      end do
    else if(sthopt.ge.1) then
      do k=2,nk-1
        do j=1,nj-1
          do i=1,ni-1
            wc(i,j,k)=w(i,j,k)/jcb8w(i,j,k)
          end do
        end do
      end do
    end if
  else
    do k=2,nk-1
      do j=1,nj-1
        do i=1,ni
          j31u2(i,j,k)=(u(i,j,k-1)+u(i,j,k))*j31(i,j,k)
        end do
      end do
      do j=1,nj
        do i=1,ni-1
          j32v2(i,j,k)=(v(i,j,k-1)+v(i,j,k))*j32(i,j,k)
        end do
      end do
    end do

    if(mfcopt.eq.0) then
      do k=2,nk-1
        do j=1,nj-1
          do i=1,ni-1
            wc(i,j,k)=(.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k)) &
                 +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
          end do
        end do
      end do
    else
      if(mpopt.eq.0.or.mpopt.eq.10) then
        do k=2,nk-1
          do j=1,nj-1
            do i=1,ni-1
              wc(i,j,k)=(.25e0*(mf(i,j)*(j31u2(i,j,k)+j31u2(i+1,j,k)) &
                   +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
            end do
          end do
        end do
      else if(mpopt.eq.5) then
        do k=2,nk-1
          do j=1,nj-1
            do i=1,ni-1
              wc(i,j,k)=(.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k)) &
                   +mf(i,j)*(j32v2(i,j,k)+j32v2(i,j+1,k))) &
                   +w(i,j,k))/jcb8w(i,j,k)
            end do
          end do
        end do
      else
        do j=1,nj-1
          do i=1,ni-1
            mf25(i,j)=.25e0*mf(i,j)
          end do
        end do
        do k=2,nk-1
          do j=1,nj-1
            do i=1,ni-1
              wc(i,j,k)=(mf25(i,j)*((j31u2(i,j,k)+j31u2(i+1,j,k)) &
                   +(j32v2(i,j,k)+j32v2(i,j+1,k)))+w(i,j,k))/jcb8w(i,j,k)
            end do
          end do
        end do
      end if
    end if
  end if

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/wc_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) wc
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(I12)') sthopt
  write(10, '(I12)') trnopt
  write(10, '(I12)') mpopt
  write(10, '(I12)') mfcopt
  close(10)

  deallocate(j31, j32, jcb8w, mf, u, v, w, wc)
  deallocate(mf25, j31u2, j32v2)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_phy2cnt

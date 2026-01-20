!***********************************************************************
! Data Generator for vbcw Kernel Benchmark
!***********************************************************************
program generate_data_vbcw
  implicit none

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64

  ! Physics options (test case: bbc=2, tbc=2, mpopt=0, mfcopt=0)
  integer, parameter :: bbc = 2
  integer, parameter :: tbc = 2
  integer, parameter :: mpopt = 0
  integer, parameter :: mfcopt = 0

  ! Arrays
  real, allocatable :: j31(:,:,:), j32(:,:,:), jcb8w(:,:,:)
  real, allocatable :: mf(:,:)
  real, allocatable :: uf(:,:,:), vf(:,:,:), wc(:,:,:)
  real, allocatable :: wf(:,:,:), wf_ref(:,:,:)
  real, allocatable :: mf25(:,:), j31u2(:,:), j32v2(:,:)

  integer :: i, j, k
  integer :: nkm1, nkm2, nkm3
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for vbcw kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6,A,I6,A,I6)') ' Options: bbc=', bbc, ', tbc=', tbc, &
       ', mpopt=', mpopt, ', mfcopt=', mfcopt

  ! Allocate arrays
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc(0:ni+1, 0:nj+1, 1:nk))
  allocate(wf(0:ni+1, 0:nj+1, 1:nk))
  allocate(wf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf25(0:ni+1, 0:nj+1))
  allocate(j31u2(0:ni+1, 0:nj+1))
  allocate(j32v2(0:ni+1, 0:nj+1))

  ! Initialize Jacobian arrays (terrain-following coords)
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
        uf(i,j,k) = 10.0 * sin(real(i)*0.05) * cos(real(j)*0.05) * (real(k)/real(nk))
        vf(i,j,k) = 5.0 * cos(real(i)*0.05) * sin(real(j)*0.05) * (real(k)/real(nk))
        wc(i,j,k) = 0.1 * sin(real(i+j)*0.1) * (real(k)/real(nk))
        wf(i,j,k) = real(i + j*100 + k*10000) * 0.00001
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

  open(unit=10, file=trim(data_dir)//'/uf.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) uf
  close(10)

  open(unit=10, file=trim(data_dir)//'/vf.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) vf
  close(10)

  open(unit=10, file=trim(data_dir)//'/wc.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) wc
  close(10)

  open(unit=10, file=trim(data_dir)//'/wf_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) wf
  close(10)

  ! Compute reference output
  write(*,'(A)') ' Computing reference output...'
  wf_ref = wf
  mf25 = 0.0
  j31u2 = 0.0
  j32v2 = 0.0

  nkm1 = nk - 1
  nkm2 = nk - 2
  nkm3 = nk - 3

  ! Common variable (only if certain conditions)
  if((bbc.eq.2.or.tbc.eq.2).and.(mfcopt.eq.1 &
       .and.(mpopt.ne.0.and.mpopt.ne.5.and.mpopt.ne.10))) then
    do j=1,nj-1
      do i=1,ni-1
        mf25(i,j)=.25e0*mf(i,j)
      end do
    end do
  end if

  ! Bottom boundary conditions (bbc=2)
  if(bbc.eq.2) then
    do j=1,nj-1
      do i=1,ni
        j31u2(i,j)=(uf(i,j,2)+uf(i,j,3))*j31(i,j,3)
      end do
    end do

    do j=1,nj
      do i=1,ni-1
        j32v2(i,j)=(vf(i,j,2)+vf(i,j,3))*j32(i,j,3)
      end do
    end do

    if(mfcopt.eq.0) then
      do j=1,nj-1
        do i=1,ni-1
          wf_ref(i,j,1)=.25e0 &
               *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
        end do
      end do
    else
      if(mpopt.eq.0.or.mpopt.eq.10) then
        do j=1,nj-1
          do i=1,ni-1
            wf_ref(i,j,1)=.25e0*(mf(i,j)*(j31u2(i,j)+j31u2(i+1,j)) &
                 +(j32v2(i,j)+j32v2(i,j+1)))
          end do
        end do
      else if(mpopt.eq.5) then
        do j=1,nj-1
          do i=1,ni-1
            wf_ref(i,j,1)=.25e0*((j31u2(i,j)+j31u2(i+1,j)) &
                 +mf(i,j)*(j32v2(i,j)+j32v2(i,j+1)))
          end do
        end do
      else
        do j=1,nj-1
          do i=1,ni-1
            wf_ref(i,j,1)=mf25(i,j) &
                 *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
          end do
        end do
      end if
    end if

    do j=1,nj-1
      do i=1,ni-1
        wf_ref(i,j,1)=-jcb8w(i,j,3)*wc(i,j,3)-wf_ref(i,j,1)
      end do
    end do

    do j=1,nj-1
      do i=1,ni
        j31u2(i,j)=(uf(i,j,1)+uf(i,j,2))*j31(i,j,2)
      end do
    end do

    do j=1,nj
      do i=1,ni-1
        j32v2(i,j)=(vf(i,j,1)+vf(i,j,2))*j32(i,j,2)
      end do
    end do

    if(mfcopt.eq.0) then
      do j=1,nj-1
        do i=1,ni-1
          wf_ref(i,j,2)=-.25e0 &
               *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
        end do
      end do
    else
      if(mpopt.eq.0.or.mpopt.eq.10) then
        do j=1,nj-1
          do i=1,ni-1
            wf_ref(i,j,2)=-.25e0*(mf(i,j)*(j31u2(i,j)+j31u2(i+1,j)) &
                 +(j32v2(i,j)+j32v2(i,j+1)))
          end do
        end do
      else if(mpopt.eq.5) then
        do j=1,nj-1
          do i=1,ni-1
            wf_ref(i,j,2)=-.25e0*((j31u2(i,j)+j31u2(i+1,j)) &
                 +mf(i,j)*(j32v2(i,j)+j32v2(i,j+1)))
          end do
        end do
      else
        do j=1,nj-1
          do i=1,ni-1
            wf_ref(i,j,2)=-mf25(i,j) &
                 *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
          end do
        end do
      end if
    end if

  else if(bbc.eq.3) then
    do j=1,nj-1
      do i=1,ni-1
        wf_ref(i,j,1)=wf_ref(i,j,2)
      end do
    end do
  end if

  ! Top boundary conditions (tbc=2)
  if(tbc.eq.2) then
    do j=1,nj-1
      do i=1,ni
        j31u2(i,j)=(uf(i,j,nkm3)+uf(i,j,nkm2))*j31(i,j,nkm2)
      end do
    end do

    do j=1,nj
      do i=1,ni-1
        j32v2(i,j)=(vf(i,j,nkm3)+vf(i,j,nkm2))*j32(i,j,nkm2)
      end do
    end do

    if(mfcopt.eq.0) then
      do j=1,nj-1
        do i=1,ni-1
          wf_ref(i,j,nk)=.25e0 &
               *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
        end do
      end do
    else
      if(mpopt.eq.0.or.mpopt.eq.10) then
        do j=1,nj-1
          do i=1,ni-1
            wf_ref(i,j,nk)=.25e0*(mf(i,j)*(j31u2(i,j)+j31u2(i+1,j)) &
                 +(j32v2(i,j)+j32v2(i,j+1)))
          end do
        end do
      else if(mpopt.eq.5) then
        do j=1,nj-1
          do i=1,ni-1
            wf_ref(i,j,nk)=.25e0*((j31u2(i,j)+j31u2(i+1,j)) &
                 +mf(i,j)*(j32v2(i,j)+j32v2(i,j+1)))
          end do
        end do
      else
        do j=1,nj-1
          do i=1,ni-1
            wf_ref(i,j,nk)=mf25(i,j) &
                 *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
          end do
        end do
      end if
    end if

    do j=1,nj-1
      do i=1,ni-1
        wf_ref(i,j,nk)=-jcb8w(i,j,nkm2)*wc(i,j,nkm2)-wf_ref(i,j,nk)
      end do
    end do

    do j=1,nj-1
      do i=1,ni-1
        wf_ref(i,j,nkm1)=0.e0
      end do
    end do

  else if(tbc.eq.3) then
    do j=1,nj-1
      do i=1,ni-1
        wf_ref(i,j,nk)=wf_ref(i,j,nkm1)
      end do
    end do
  end if

  ! Write reference output
  open(unit=10, file=trim(data_dir)//'/wf_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) wf_ref
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(I12)') bbc
  write(10, '(I12)') tbc
  write(10, '(I12)') mpopt
  write(10, '(I12)') mfcopt
  close(10)

  deallocate(j31, j32, jcb8w, mf, uf, vf, wc, wf, wf_ref)
  deallocate(mf25, j31u2, j32v2)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_vbcw

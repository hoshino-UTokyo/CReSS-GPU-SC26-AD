!***********************************************************************
! Data Generator for bc4news Kernel Benchmark
!***********************************************************************
program generate_data_bc4news
  implicit none

  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: kmax = 64

  ! Boundary indices
  integer, parameter :: isw = 0
  integer, parameter :: ise = 101  ! ni+1
  integer, parameter :: jss = 0
  integer, parameter :: jsn = 101  ! nj+1

  ! Boundary condition options (non-1 to enable corner processing)
  integer, parameter :: wbc = 2
  integer, parameter :: ebc = 2
  integer, parameter :: sbc = 2
  integer, parameter :: nbc = 2

  ! Domain decomposition (single domain - all corners active)
  integer, parameter :: ebsw = 1
  integer, parameter :: ebse = 1
  integer, parameter :: ebnw = 1
  integer, parameter :: ebne = 1
  integer, parameter :: isub = 0
  integer, parameter :: jsub = 0
  integer, parameter :: nisub = 1
  integer, parameter :: njsub = 1

  ! Derived indices
  integer :: iswp1, isem1, jssp1, jsnm1

  ! Arrays
  real, allocatable :: var(:,:,:), var_ref(:,:,:)
  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for bc4news kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', kmax=', kmax

  iswp1 = isw + 1
  isem1 = ise - 1
  jssp1 = jss + 1
  jsnm1 = jsn - 1

  allocate(var(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var_ref(0:ni+1, 0:nj+1, 1:kmax))

  ! Initialize array with varying values
  do k = 1, kmax
    do j = 0, nj+1
      do i = 0, ni+1
        var(i,j,k) = real(i + j*100 + k*10000) * 0.001
      end do
    end do
  end do

  ! Save input
  open(unit=10, file=trim(data_dir)//'/var_input.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) var
  close(10)

  ! Compute reference (apply corner boundary conditions)
  var_ref = var

  if(abs(wbc).ne.1.or.abs(ebc).ne.1.or.abs(sbc).ne.1.or.abs(nbc).ne.1) then

    ! SW corner
    if(ebsw.eq.1.and.isub.eq.0.and.jsub.eq.0) then
      do k=1,kmax
        var_ref(isw,jss,k)=.5e0*(var_ref(iswp1,jss,k)+var_ref(isw,jssp1,k))
      end do
    end if

    ! SE corner
    if(ebse.eq.1.and.isub.eq.nisub-1.and.jsub.eq.0) then
      do k=1,kmax
        var_ref(ise,jss,k)=.5e0*(var_ref(isem1,jss,k)+var_ref(ise,jssp1,k))
      end do
    end if

    ! NW corner
    if(ebnw.eq.1.and.isub.eq.0.and.jsub.eq.njsub-1) then
      do k=1,kmax
        var_ref(isw,jsn,k)=.5e0*(var_ref(iswp1,jsn,k)+var_ref(isw,jsnm1,k))
      end do
    end if

    ! NE corner
    if(ebne.eq.1.and.isub.eq.nisub-1.and.jsub.eq.njsub-1) then
      do k=1,kmax
        var_ref(ise,jsn,k)=.5e0*(var_ref(isem1,jsn,k)+var_ref(ise,jsnm1,k))
      end do
    end if

  end if

  ! Write reference
  open(unit=10, file=trim(data_dir)//'/var_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) var_ref
  close(10)

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') kmax
  write(10, '(I12)') isw
  write(10, '(I12)') ise
  write(10, '(I12)') jss
  write(10, '(I12)') jsn
  write(10, '(I12)') wbc
  write(10, '(I12)') ebc
  write(10, '(I12)') sbc
  write(10, '(I12)') nbc
  write(10, '(I12)') ebsw
  write(10, '(I12)') ebse
  write(10, '(I12)') ebnw
  write(10, '(I12)') ebne
  write(10, '(I12)') isub
  write(10, '(I12)') jsub
  write(10, '(I12)') nisub
  write(10, '(I12)') njsub
  close(10)

  deallocate(var, var_ref)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_bc4news

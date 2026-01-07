!***********************************************************************
! Data Generator for copy3d Kernel Benchmark
!***********************************************************************
program generate_data_copy3d
  implicit none

  ! Grid dimensions
  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64
  integer, parameter :: imin = 0, imax = ni+1
  integer, parameter :: jmin = 0, jmax = nj+1
  integer, parameter :: kmin = 1, kmax = nk

  ! Arrays
  real, allocatable :: invar(:,:,:)
  real, allocatable :: outvar(:,:,:)

  integer :: i, j, k
  character(len=256) :: data_dir

  data_dir = './data'
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for copy3d kernel'
  write(*,'(A)') '=================================================='

  allocate(invar(imin:imax, jmin:jmax, kmin:kmax))
  allocate(outvar(imin:imax, jmin:jmax, kmin:kmax))

  ! Initialize input with varying values
  do k = kmin, kmax
    do j = jmin, jmax
      do i = imin, imax
        invar(i,j,k) = real(i + j*100 + k*10000) * 0.001
      end do
    end do
  end do

  ! Compute reference output (simple copy)
  outvar = invar

  ! Write parameters
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(I12)') imin
  write(10, '(I12)') imax
  write(10, '(I12)') jmin
  write(10, '(I12)') jmax
  write(10, '(I12)') kmin
  write(10, '(I12)') kmax
  close(10)

  ! Write arrays
  open(unit=10, file=trim(data_dir)//'/invar.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) invar
  close(10)

  open(unit=10, file=trim(data_dir)//'/outvar_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) outvar
  close(10)

  deallocate(invar, outvar)

  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='

end program generate_data_copy3d

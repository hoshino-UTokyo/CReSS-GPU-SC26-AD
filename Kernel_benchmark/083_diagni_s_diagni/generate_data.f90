!***********************************************************************
! Data Generator for diagni Kernel Benchmark
!***********************************************************************
!
! Generates synthetic test data for the diagni kernel benchmark.
! Can be used for testing when actual simulation data is not available.
!
! Usage: ./generate_data
!
!***********************************************************************
program generate_data_diagni
  implicit none

  ! Grid dimensions (adjust as needed)
  integer, parameter :: ni = 100
  integer, parameter :: nj = 100
  integer, parameter :: nk = 64
  integer, parameter :: nqi = 4  ! Number of ice categories
  integer, parameter :: nni = 4  ! Number of concentration categories
  integer, parameter :: haiopt = 1  ! Include hail (0 or 1)

  ! Physical constants (from m_comphy and m_commath)
  real, parameter :: cc = 3.141592e0
  real, parameter :: ns0 = 1.8e6
  real, parameter :: ng0 = 1.1e6
  real, parameter :: nh0 = 1.1e6
  real, parameter :: rhos = 1.e2
  real, parameter :: rhog = 4.e2
  real, parameter :: rhoh = 4.e2
  real, parameter :: ms0 = 1.76e-10
  real, parameter :: mg0 = 7.e-10
  real, parameter :: mh0 = 7.e-10
  real, parameter :: mimax = 1.4e-10
  real, parameter :: msmax = 6.e-4
  real, parameter :: mgmax = 8.e-3
  real, parameter :: mhmax = 8.e-3

  ! Arrays
  real, allocatable :: rbr(:,:,:)
  real, allocatable :: qice(:,:,:,:)
  real, allocatable :: nidia(:,:,:,:)

  ! Local variables
  integer :: i, j, k, n
  real :: rbv
  real :: miiv, msmiv2, ms0iv2, mgmiv2, mg0iv2, mhmiv2, mh0iv2
  real :: cdiaqs, cdiaqg, cdiaqh

  character(len=256) :: data_dir

  data_dir = './data'

  ! Create data directory
  call execute_command_line('mkdir -p '//trim(data_dir), wait=.true.)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Generating test data for diagni kernel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' nqi=', nqi, ', nni=', nni
  write(*,'(A,I6)') ' haiopt=', haiopt

  ! Allocate arrays
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qice(0:ni+1, 0:nj+1, 1:nk, 1:nqi))
  allocate(nidia(0:ni+1, 0:nj+1, 1:nk, 1:nni))

  ! Initialize input data with realistic values
  write(*,'(A)') ' Initializing arrays...'

  ! Base state density (typical atmospheric values: ~1.2 kg/m3 at surface)
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        ! Exponential decrease with height
        rbr(i,j,k) = 1.2 * exp(-real(k-1) / 30.0)
      end do
    end do
  end do

  ! Ice mixing ratios (typical values: 1e-6 to 1e-3 kg/kg)
  do n = 1, nqi
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          ! Create some spatial variation
          qice(i,j,k,n) = 1.0e-5 * (1.0 + 0.1 * sin(real(i)*0.1) * cos(real(j)*0.1))
          ! Zero out below freezing level
          if (k < 10) qice(i,j,k,n) = 0.0
        end do
      end do
    end do
  end do

  ! Write parameters
  write(*,'(A)') ' Writing params.txt...'
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I12)') ni
  write(10, '(I12)') nj
  write(10, '(I12)') nk
  write(10, '(I12)') nqi
  write(10, '(I12)') nni
  write(10, '(I12)') haiopt
  close(10)

  ! Write input arrays
  write(*,'(A)') ' Writing rbr.bin...'
  open(unit=10, file=trim(data_dir)//'/rbr.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) rbr
  close(10)

  write(*,'(A)') ' Writing qice.bin...'
  open(unit=10, file=trim(data_dir)//'/qice.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) qice
  close(10)

  ! Compute reference output
  write(*,'(A)') ' Computing reference output...'

  ! Set derived constants
  miiv = 4.e0 / (3.e0 * mimax)
  msmiv2 = 1.e-2 / msmax
  ms0iv2 = 1.e2 / ms0
  mgmiv2 = 1.e-2 / mgmax
  mg0iv2 = 1.e2 / mg0
  mhmiv2 = 1.e-2 / mhmax
  mh0iv2 = 1.e2 / mh0
  cdiaqs = ns0*ns0*ns0 / (cc*rhos)
  cdiaqg = ng0*ng0*ng0 / (cc*rhog)
  cdiaqh = nh0*nh0*nh0 / (cc*rhoh)

  nidia = 0.0

  if (haiopt == 0) then

    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          rbv = 1.e0 / rbr(i,j,k)

          ! Cloud ice
          nidia(i,j,k,1) = miiv * qice(i,j,k,1)

          ! Snow
          nidia(i,j,k,2) = sqrt(sqrt(cdiaqs*rbr(i,j,k)*qice(i,j,k,2))) * rbv
          nidia(i,j,k,2) = min(max(nidia(i,j,k,2), msmiv2*qice(i,j,k,2)), ms0iv2*qice(i,j,k,2))

          ! Graupel
          nidia(i,j,k,3) = sqrt(sqrt(cdiaqg*rbr(i,j,k)*qice(i,j,k,3))) * rbv
          nidia(i,j,k,3) = min(max(nidia(i,j,k,3), mgmiv2*qice(i,j,k,3)), mg0iv2*qice(i,j,k,3))
        end do
      end do
    end do

  else

    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          rbv = 1.e0 / rbr(i,j,k)

          ! Cloud ice
          nidia(i,j,k,1) = miiv * qice(i,j,k,1)

          ! Snow
          nidia(i,j,k,2) = sqrt(sqrt(cdiaqs*rbr(i,j,k)*qice(i,j,k,2))) * rbv
          nidia(i,j,k,2) = min(max(nidia(i,j,k,2), msmiv2*qice(i,j,k,2)), ms0iv2*qice(i,j,k,2))

          ! Graupel
          nidia(i,j,k,3) = sqrt(sqrt(cdiaqg*rbr(i,j,k)*qice(i,j,k,3))) * rbv
          nidia(i,j,k,3) = min(max(nidia(i,j,k,3), mgmiv2*qice(i,j,k,3)), mg0iv2*qice(i,j,k,3))

          ! Hail
          nidia(i,j,k,4) = sqrt(sqrt(cdiaqh*rbr(i,j,k)*qice(i,j,k,4))) * rbv
          nidia(i,j,k,4) = min(max(nidia(i,j,k,4), mhmiv2*qice(i,j,k,4)), mh0iv2*qice(i,j,k,4))
        end do
      end do
    end do

  end if

  ! Write reference output
  write(*,'(A)') ' Writing nidia_ref.bin...'
  open(unit=10, file=trim(data_dir)//'/nidia_ref.bin', status='replace', &
       access='stream', form='unformatted')
  write(10) nidia
  close(10)

  ! Cleanup
  deallocate(rbr, qice, nidia)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') ' Files created in: '//trim(data_dir)
  write(*,'(A)') '=================================================='

end program generate_data_diagni

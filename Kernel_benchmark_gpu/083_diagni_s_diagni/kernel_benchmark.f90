!***********************************************************************
! Kernel Benchmark: diagni (s_diagni) - GPU version
!***********************************************************************
!
! Source: Src/diagni.f90
! Description: Calculate diagnostic concentrations for all ice hydrometeor
!              categories (cloud ice, snow, graupel, hail) from mixing ratios.
!
!***********************************************************************
program kernel_benchmark_diagni
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk, nqi, nni

  ! Parameters
  integer :: haiopt

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

  ! Input arrays
  real, allocatable :: rbr(:,:,:)
  real, allocatable :: qice(:,:,:,:)

  ! Output array
  real, allocatable :: nidia(:,:,:,:)

  ! Reference output for validation
  real, allocatable :: nidia_ref(:,:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
  integer :: error_count
  logical :: validation_passed

  ! Loop variables
  integer :: iter, i, j, k, n

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', &
       ni, nj, nk, nqi, nni, haiopt)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: diagni (GPU)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' nqi=', nqi, ', nni=', nni
  write(*,'(A,I6)') ' haiopt=', haiopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qice(0:ni+1, 0:nj+1, 1:nk, 1:nqi))
  allocate(nidia(0:ni+1, 0:nj+1, 1:nk, 1:nni))
  allocate(nidia_ref(0:ni+1, 0:nj+1, 1:nk, 1:nni))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_4d(trim(data_dir)//'/qice.bin', qice, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_4d(trim(data_dir)//'/nidia_ref.bin', nidia_ref, 0, ni+1, 0, nj+1, 1, nk, 1, nni)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_diagni(haiopt, ni, nj, nk, nqi, nni, rbr, qice, nidia)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    !$acc wait
    t_start = omp_get_wtime()

    call kernel_diagni(haiopt, ni, nj, nk, nqi, nni, rbr, qice, nidia)

    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0

  do n = 1, nni
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          rel_error = abs(nidia(i,j,k,n) - nidia_ref(i,j,k,n))
          if (abs(nidia_ref(i,j,k,n)) > 1.0e-20) then
            rel_error = rel_error / abs(nidia_ref(i,j,k,n))
          end if
          if (rel_error > max_error) max_error = rel_error
          if (rel_error > tolerance) error_count = error_count + 1
        end do
      end do
    end do
  end do

  validation_passed = (error_count == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:    ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:        ', minval(times) * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:        ', maxval(times) * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:      ', t_total * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(rbr, qice, nidia, nidia_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: diagni (GPU version)
  ! Calculate diagnostic concentrations for all ice hydrometeor categories
  !=====================================================================
  subroutine kernel_diagni(haiopt, ni, nj, nk, nqi, nni, rbr, qice, nidia)
    implicit none

    integer, intent(in) :: haiopt
    integer, intent(in) :: ni, nj, nk, nqi, nni

    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qice(0:ni+1, 0:nj+1, 1:nk, 1:nqi)

    real, intent(out) :: nidia(0:ni+1, 0:nj+1, 1:nk, 1:nni)

    ! Local variables
    integer :: i, j, k
    real :: rbv

    ! Derived constants
    real :: miiv, msmiv2, ms0iv2, mgmiv2, mg0iv2, mhmiv2, mh0iv2
    real :: cdiaqs, cdiaqg, cdiaqh

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

    if (haiopt == 0) then

      !$acc kernels
      !$acc loop independent collapse(3) private(rbv)
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
      !$acc end kernels

    else

      !$acc kernels
      !$acc loop independent collapse(3) private(rbv)
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
      !$acc end kernels

    end if

  end subroutine kernel_diagni

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    character(len=256) :: config_file
    integer :: ios
    logical :: exists

    ! Default values
    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    config_file = 'benchmark.conf'
    inquire(file=config_file, exist=exists)

    if (exists) then
      open(unit=10, file=config_file, status='old', iostat=ios)
      if (ios == 0) then
        read(10, '(A)', iostat=ios) data_dir
        read(10, *, iostat=ios) num_iter
        read(10, *, iostat=ios) warmup_iter
        read(10, *, iostat=ios) tol
        close(10)
      end if
    end if

  end subroutine read_config

  !=====================================================================
  ! Parameter reader
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk, nqi, nni, haiopt)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, nqi, nni, haiopt

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open parameter file: ', trim(filename)
      stop 1
    end if

    do while (.true.)
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))
        select case (trim(key))
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('nqi')
            read(val, *) nqi
          case ('nni')
            read(val, *) nni
          case ('haiopt')
            read(val, *) haiopt
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary array readers
  !=====================================================================
  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_3d

  subroutine read_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_4d

end program kernel_benchmark_diagni

!***********************************************************************
! Kernel Benchmark: diagnw (s_diagnw) - GPU version
!***********************************************************************
!
! Source: Src/diagnw.f90
! Description: Compute diagnostic concentrations of cloud water and rain
!              water based on base state density and water hydrometeor
!              mixing ratios.
!
!***********************************************************************
program kernel_benchmark_diagnw
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk, nqw, nnw

  ! Physical constants from m_comphy
  real, parameter :: rhow = 1.0e3       ! Density of water
  real, parameter :: nr0 = 8.0e6        ! Parameter of rain water size distribution
  real, parameter :: mr0 = 2.7e-10      ! Mass of minimum rain water
  real, parameter :: mrmax = 2.68e-4    ! Mass of maximum rain water
  real, parameter :: nclcst = 1.0e8     ! Constant concentrations of cloud

  ! Mathematical constant from m_commath
  real, parameter :: cc = 3.141592      ! Pi

  ! Input arrays
  real, allocatable :: rbr(:,:,:)
  real, allocatable :: qwtr(:,:,:,:)

  ! Output array
  real, allocatable :: nwdia(:,:,:,:)

  ! Reference output for validation
  real, allocatable :: nwdia_ref(:,:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, nqw, nnw)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: diagnw (GPU)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' Water categories: nqw=', nqw, ', nnw=', nnw
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qwtr(0:ni+1, 0:nj+1, 1:nk, 1:nqw))
  allocate(nwdia(0:ni+1, 0:nj+1, 1:nk, 1:nnw))
  allocate(nwdia_ref(0:ni+1, 0:nj+1, 1:nk, 1:nnw))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_4d(trim(data_dir)//'/qwtr.bin', qwtr, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_4d(trim(data_dir)//'/nwdia_ref.bin', nwdia_ref, 0, ni+1, 0, nj+1, 1, nk, 1, nnw)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    nwdia = 0.0
    call kernel_diagnw(ni, nj, nk, nqw, nnw, cc, rhow, nr0, mr0, mrmax, nclcst, &
                       rbr, qwtr, nwdia)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    nwdia = 0.0

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_diagnw(ni, nj, nk, nqw, nnw, cc, rhow, nr0, mr0, mrmax, nclcst, &
                       rbr, qwtr, nwdia)
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

  do n = 1, nnw
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          rel_error = abs(nwdia(i,j,k,n) - nwdia_ref(i,j,k,n))
          if (abs(nwdia_ref(i,j,k,n)) > 1.0e-10) then
            rel_error = rel_error / abs(nwdia_ref(i,j,k,n))
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
  deallocate(rbr, qwtr, nwdia, nwdia_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: diagnw (GPU version)
  ! Compute diagnostic concentrations of water hydrometeors
  !=====================================================================
  subroutine kernel_diagnw(ni, nj, nk, nqw, nnw, cc, rhow, nr0, mr0, mrmax, nclcst, &
                           rbr, qwtr, nwdia)
    implicit none

    integer, intent(in) :: ni, nj, nk, nqw, nnw
    real, intent(in) :: cc, rhow, nr0, mr0, mrmax, nclcst
    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qwtr(0:ni+1, 0:nj+1, 1:nk, 1:nqw)
    real, intent(out) :: nwdia(0:ni+1, 0:nj+1, 1:nk, 1:nnw)

    integer :: i, j, k
    real :: rbv
    real :: mrmiv2, mr0iv2, cdiaqr

    ! Set the common used variables
    mrmiv2 = 1.0e-2 / mrmax
    mr0iv2 = 1.0e2 / mr0
    cdiaqr = nr0 * nr0 * nr0 / (cc * rhow)

    !$acc kernels
    !$acc loop independent collapse(3) private(rbv)
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          ! Calculate the inverse of base state density
          rbv = 1.0 / rbr(i,j,k)

          ! Get the diagnostic concentrations of cloud water
          nwdia(i,j,k,1) = nclcst * rbv

          ! Get the diagnostic concentrations of rain water
          nwdia(i,j,k,2) = sqrt(sqrt(cdiaqr * rbr(i,j,k) * qwtr(i,j,k,2))) * rbv
          nwdia(i,j,k,2) = min(max(nwdia(i,j,k,2), mrmiv2*qwtr(i,j,k,2)), &
                               mr0iv2*qwtr(i,j,k,2))
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_diagnw

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    inquire(file='benchmark.conf', exist=exists)
    if (exists) then
      open(unit=10, file='benchmark.conf', status='old', iostat=ios)
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
  subroutine read_parameters(filename, ni, nj, nk, nqw, nnw)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, nqw, nnw

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    ni = 1
    nj = 1
    nk = 1
    nqw = 2
    nnw = 2

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
          case ('nqw')
            read(val, *) nqw
          case ('nnw')
            read(val, *) nnw
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary 3D array reader
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

  !=====================================================================
  ! Binary 4D array reader
  !=====================================================================
  subroutine read_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, n1, n2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, n1, n2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2, n1:n2)

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

end program kernel_benchmark_diagnw

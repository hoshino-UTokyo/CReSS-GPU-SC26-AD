!***********************************************************************
! Kernel Benchmark: steps (s_steps) - Simplified dry case
!***********************************************************************
!
! Source: Src/steps.f90
! Description: Advances scalar variables (potential temperature) to next
!              time step using forcing terms. This is a simplified version
!              testing the core ptp update (gwmopt=0, fmois=dry).
!
!***********************************************************************
program kernel_benchmark_steps
  use omp_lib
  implicit none

  ! Grid dimensions
  integer :: ni, nj, nk

  ! Time step
  real :: dtb
  real :: dtb2

  ! Options
  integer :: advopt

  ! Arrays
  real, allocatable :: rst(:,:,:)
  real, allocatable :: ptpp(:,:,:)
  real, allocatable :: ptfrc(:,:,:)
  real, allocatable :: dtdrst(:,:,:)
  real, allocatable :: ptpf(:,:,:)
  real, allocatable :: ptpf_ref(:,:,:)

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
  integer :: iter, i, j, k

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, dtb, advopt)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: steps (simplified)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dtb: ', dtb
  write(*,'(A,I6)') ' advopt: ', advopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(dtdrst(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptpp.bin', ptpp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptfrc.bin', ptfrc, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ptpf_ref.bin', ptpf_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ptpf = 0.0
    dtdrst = 0.0
    call kernel_steps(advopt, dtb, ni, nj, nk, rst, ptpp, ptfrc, dtdrst, ptpf)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ptpf = 0.0
    dtdrst = 0.0

    t_start = omp_get_wtime()
    call kernel_steps(advopt, dtb, ni, nj, nk, rst, ptpp, ptfrc, dtdrst, ptpf)
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

  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(ptpf(i,j,k) - ptpf_ref(i,j,k))
        if (abs(ptpf_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(ptpf_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  validation_passed = (error_count == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', minval(times) * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', maxval(times) * 1000.0d0, ' ms'
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
  deallocate(rst, ptpp, ptfrc, dtdrst, ptpf, ptpf_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: steps (simplified for dry case with gwmopt=0)
  !=====================================================================
  subroutine kernel_steps(advopt, dtb, ni, nj, nk, rst, ptpp, ptfrc, dtdrst, ptpf)
    implicit none

    integer, intent(in) :: advopt
    real, intent(in) :: dtb
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptpp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: dtdrst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: ptpf(0:ni+1, 0:nj+1, 1:nk)

    real :: dtb2
    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! Compute dtdrst
    if (advopt <= 3) then
      dtb2 = 2.0 * dtb
      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-2
            dtdrst(i,j,k) = dtb2 / rst(i,j,k)
          end do
        end do
        !$omp end do
      end do
    else
      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-2
            dtdrst(i,j,k) = dtb / rst(i,j,k)
          end do
        end do
        !$omp end do
      end do
    end if

    ! Solve ptp to next time step (gwmopt=0 case)
    do k = 2, nk-2
      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-2
        do i = 2, ni-2
          ptpf(i,j,k) = ptpp(i,j,k) + ptfrc(i,j,k) * dtdrst(i,j,k)
        end do
      end do
      !$omp end do
    end do

    !$omp end parallel

  end subroutine kernel_steps

  !=====================================================================
  ! I/O subroutines
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    character(len=256) :: config_file
    integer :: ios
    logical :: exists

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

  subroutine read_parameters(filename, ni, nj, nk, dtb, advopt)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, advopt
    real, intent(out) :: dtb
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
          case ('dtb')
            read(val, *) dtb
          case ('advopt')
            read(val, *) advopt
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

end program kernel_benchmark_steps

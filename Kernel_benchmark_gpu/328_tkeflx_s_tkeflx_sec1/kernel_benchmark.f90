!***********************************************************************
! GPU Kernel Benchmark: tkeflx sec1 (s_tkeflx terrain preprocessing)
!***********************************************************************
!
! Source: Src/tkeflx.f90
! Description: Calculate j31*tke and j32*tke products for terrain-following
!              coordinate transformation of TKE turbulent fluxes.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_tkeflx_sec1
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Input arrays
  real, allocatable :: j31(:,:,:)
  real, allocatable :: j32(:,:,:)
  real, allocatable :: tke(:,:,:)

  ! Output arrays
  real, allocatable :: j31tke(:,:,:)
  real, allocatable :: j32tke(:,:,:)

  ! Reference outputs for validation
  real, allocatable :: j31tke_ref(:,:,:)
  real, allocatable :: j32tke_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, max_error2, rel_error
  real :: tolerance
  integer :: error_count, error_count2
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
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: tkeflx_sec1'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(j31tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(j31tke_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32tke_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tke.bin', tke, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/j31tke_ref.bin', j31tke_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32tke_ref.bin', j32tke_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    j31tke = 0.0; j32tke = 0.0
    call kernel_tkeflx_sec1(ni, nj, nk, j31, j32, tke, j31tke, j32tke)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    j31tke = 0.0; j32tke = 0.0

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_tkeflx_sec1(ni, nj, nk, j31, j32, tke, j31tke, j32tke)
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

  max_error = 0.0; error_count = 0
  do k = 2, nk-1
    do j = 1, nj-1
      do i = 2, ni-1
        rel_error = abs(j31tke(i,j,k) - j31tke_ref(i,j,k))
        if (abs(j31tke_ref(i,j,k)) > 1.0e-10) rel_error = rel_error / abs(j31tke_ref(i,j,k))
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  max_error2 = 0.0; error_count2 = 0
  do k = 2, nk-1
    do j = 2, nj-1
      do i = 1, ni-1
        rel_error = abs(j32tke(i,j,k) - j32tke_ref(i,j,k))
        if (abs(j32tke_ref(i,j,k)) > 1.0e-10) rel_error = rel_error / abs(j32tke_ref(i,j,k))
        if (rel_error > max_error2) max_error2 = rel_error
        if (rel_error > tolerance) error_count2 = error_count2 + 1
      end do
    end do
  end do

  validation_passed = (error_count == 0 .and. error_count2 == 0)

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
  write(*,'(A,ES12.4)') ' j31tke max error:   ', max_error
  write(*,'(A,I12)')    ' j31tke error count: ', error_count
  write(*,'(A,ES12.4)') ' j32tke max error:   ', max_error2
  write(*,'(A,I12)')    ' j32tke error count: ', error_count2
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  deallocate(j31, j32, tke, j31tke, j32tke, j31tke_ref, j32tke_ref, times)
  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: tkeflx_sec1 - GPU version using OpenACC
  !=====================================================================
  subroutine kernel_tkeflx_sec1(ni, nj, nk, j31, j32, tke, j31tke, j32tke)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: j31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tke(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: j31tke(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: j32tke(0:ni+1, 0:nj+1, 1:nk)
    integer :: i, j, k

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 2, ni-1
          j31tke(i,j,k) = ((tke(i-1,j,k-1)+tke(i,j,k-1)) &
                          +(tke(i-1,j,k)+tke(i,j,k))) * j31(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      !$acc loop independent
      do j = 2, nj-1
        !$acc loop independent
        do i = 1, ni-1
          j32tke(i,j,k) = ((tke(i,j-1,k-1)+tke(i,j,k-1)) &
                          +(tke(i,j-1,k)+tke(i,j,k))) * j32(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_tkeflx_sec1

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios
    logical :: exists
    data_dir = './data'; num_iter = 10; warmup_iter = 2; tol = 1.0e-5
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
  subroutine read_parameters(filename, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    character(len=256) :: line, key, val
    integer :: ios, eq_pos
    ni = 1; nj = 1; nk = 1
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
          case ('ni'); read(val, *) ni
          case ('nj'); read(val, *) nj
          case ('nk'); read(val, *) nk
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !=====================================================================
  ! Binary array reader
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

end program kernel_benchmark_tkeflx_sec1

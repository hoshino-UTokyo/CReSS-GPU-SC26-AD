!***********************************************************************
! Kernel Benchmark: defomten sec1 (velocity at w-points for deformation tensor)
!***********************************************************************
!
! Source: Src/defomten.f90
! Description: Compute s11 and s22 by summing velocity components at
!              adjacent k-levels (w-points) for the deformation tensor.
!
!***********************************************************************
program kernel_benchmark_defomten_sec1
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Boundary parameters
  integer :: jsouth, jnorth, iwest, ieast

  ! Input arrays
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)

  ! Output arrays
  real, allocatable :: s11(:,:,:)
  real, allocatable :: s22(:,:,:)

  ! Reference outputs for validation
  real, allocatable :: s11_ref(:,:,:)
  real, allocatable :: s22_ref(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', &
       ni, nj, nk, jsouth, jnorth, iwest, ieast)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: defomten_sec1'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6,A,I6,A,I6)') ' Boundaries: jsouth=', jsouth, &
       ', jnorth=', jnorth, ', iwest=', iwest, ', ieast=', ieast
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(s11(0:ni+1, 0:nj+1, 1:nk))
  allocate(s22(0:ni+1, 0:nj+1, 1:nk))
  allocate(s11_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(s22_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/s11_ref.bin', s11_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s22_ref.bin', s22_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    s11 = 0.0; s22 = 0.0
    call kernel_defomten_sec1(ni, nj, nk, jsouth, jnorth, iwest, ieast, &
         u, v, s11, s22)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    s11 = 0.0; s22 = 0.0

    t_start = omp_get_wtime()
    call kernel_defomten_sec1(ni, nj, nk, jsouth, jnorth, iwest, ieast, &
         u, v, s11, s22)
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
    do j = jsouth, nj-jnorth
      do i = 1, ni
        rel_error = abs(s11(i,j,k) - s11_ref(i,j,k))
        if (abs(s11_ref(i,j,k)) > 1.0e-10) rel_error = rel_error / abs(s11_ref(i,j,k))
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  max_error2 = 0.0; error_count2 = 0
  do k = 2, nk-1
    do j = 1, nj
      do i = iwest, ni-ieast
        rel_error = abs(s22(i,j,k) - s22_ref(i,j,k))
        if (abs(s22_ref(i,j,k)) > 1.0e-10) rel_error = rel_error / abs(s22_ref(i,j,k))
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
  write(*,'(A,ES12.4)') ' s11 max error:   ', max_error
  write(*,'(A,I12)')    ' s11 error count: ', error_count
  write(*,'(A,ES12.4)') ' s22 max error:   ', max_error2
  write(*,'(A,I12)')    ' s22 error count: ', error_count2
  write(*,'(A,ES12.4)') ' Tolerance:        ', tolerance
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  deallocate(u, v, s11, s22, s11_ref, s22_ref, times)
  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: defomten_sec1 - OpenMP version
  !=====================================================================
  subroutine kernel_defomten_sec1(ni, nj, nk, jsouth, jnorth, iwest, ieast, &
       u, v, s11, s22)
    implicit none
    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: jsouth, jnorth, iwest, ieast
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: s11(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: s22(0:ni+1, 0:nj+1, 1:nk)
    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    do k = 2, nk-1

      !$omp do schedule(runtime) private(i,j)
      do j = jsouth, nj-jnorth
        do i = 1, ni
          s11(i,j,k) = u(i,j,k-1) + u(i,j,k)
        end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj
        do i = iwest, ni-ieast
          s22(i,j,k) = v(i,j,k-1) + v(i,j,k)
        end do
      end do
      !$omp end do

    end do

    !$omp end parallel

  end subroutine kernel_defomten_sec1

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
  subroutine read_parameters(filename, ni, nj, nk, jsouth, jnorth, iwest, ieast)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, jsouth, jnorth, iwest, ieast
    character(len=256) :: line, key, val
    integer :: ios, eq_pos
    ni = 1; nj = 1; nk = 1
    jsouth = 1; jnorth = 1; iwest = 1; ieast = 1
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
          case ('ni');     read(val, *) ni
          case ('nj');     read(val, *) nj
          case ('nk');     read(val, *) nk
          case ('jsouth'); read(val, *) jsouth
          case ('jnorth'); read(val, *) jnorth
          case ('iwest');  read(val, *) iwest
          case ('ieast');  read(val, *) ieast
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

end program kernel_benchmark_defomten_sec1

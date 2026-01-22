!***********************************************************************
! GPU Kernel Benchmark: buoytke (s_buoytke)
!***********************************************************************
!
! Description: Calculate buoyancy production term for turbulent kinetic
!              energy equation using Brunt-Vaisala frequency and eddy
!              diffusivity. GPU version using OpenACC.
!
!***********************************************************************
program kernel_benchmark_buoytke
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Input arrays
  real, allocatable :: jcb8w(:,:,:)   ! Jacobian at w points
  real, allocatable :: nsq8w(:,:,:)   ! Brunt-Vaisala frequency squared at w points
  real, allocatable :: rkv8s(:,:,:)   ! Vertical eddy diffusivity at scalar points

  ! Input/output arrays
  real, allocatable :: tkefrc(:,:,:)  ! Forcing term in TKE equation
  real, allocatable :: tmp1(:,:,:)    ! Temporary array

  ! Reference outputs for validation
  real, allocatable :: tkefrc_ref(:,:,:)
  real, allocatable :: tmp1_ref(:,:,:)

  ! Backup arrays for iteration
  real, allocatable :: tkefrc_in(:,:,:)
  real, allocatable :: tmp1_in(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error, max_error_tmp1
  real :: tolerance
  integer :: error_count, error_count_tmp1
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
  write(*,'(A)') ' GPU Kernel Benchmark: buoytke'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(jcb8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(nsq8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkv8s(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/jcb8w.bin', jcb8w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/nsq8w.bin', nsq8w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rkv8s.bin', rkv8s, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tkefrc_in.bin', tkefrc_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp1_in.bin', tmp1_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/tkefrc_ref.bin', tkefrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp1_ref.bin', tmp1_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    tkefrc = tkefrc_in
    tmp1 = tmp1_in
    call kernel_buoytke(ni, nj, nk, jcb8w, nsq8w, rkv8s, tkefrc, tmp1)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    tkefrc = tkefrc_in
    tmp1 = tmp1_in

    call cpu_time(t_start)
    call kernel_buoytke(ni, nj, nk, jcb8w, nsq8w, rkv8s, tkefrc, tmp1)
    !$acc wait
    call cpu_time(t_end)
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'

  ! Validate tkefrc
  max_error = 0.0
  error_count = 0
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(tkefrc(i,j,k) - tkefrc_ref(i,j,k))
        if (abs(tkefrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(tkefrc_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate tmp1
  max_error_tmp1 = 0.0
  error_count_tmp1 = 0
  do k = 2, nk-1
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(tmp1(i,j,k) - tmp1_ref(i,j,k))
        if (abs(tmp1_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(tmp1_ref(i,j,k))
        end if
        if (rel_error > max_error_tmp1) max_error_tmp1 = rel_error
        if (rel_error > tolerance) error_count_tmp1 = error_count_tmp1 + 1
      end do
    end do
  end do

  validation_passed = (error_count == 0 .and. error_count_tmp1 == 0)

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
  write(*,'(A,ES12.4)') ' tkefrc max error:   ', max_error
  write(*,'(A,I12)') ' tkefrc error count: ', error_count
  write(*,'(A,ES12.4)') ' tmp1 max error:     ', max_error_tmp1
  write(*,'(A,I12)') ' tmp1 error count:   ', error_count_tmp1
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(jcb8w, nsq8w, rkv8s, tkefrc, tmp1)
  deallocate(tkefrc_ref, tmp1_ref, tkefrc_in, tmp1_in, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: buoytke - GPU version using OpenACC
  !=====================================================================
  subroutine kernel_buoytke(ni, nj, nk, jcb8w, nsq8w, rkv8s, tkefrc, tmp1)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: jcb8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: nsq8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rkv8s(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tkefrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    ! First pass: compute tmp1
    !$acc kernels
    do k = 2, nk-1
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          tmp1(i,j,k) = jcb8w(i,j,k) * nsq8w(i,j,k) &
                        * (rkv8s(i,j,k-1) + rkv8s(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    ! Second pass: update tkefrc using tmp1
    !$acc kernels
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          tkefrc(i,j,k) = tkefrc(i,j,k) - (tmp1(i,j,k) + tmp1(i,j,k+1))
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_buoytke

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
  subroutine read_parameters(filename, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ni = 1
    nj = 1
    nk = 1

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

end program kernel_benchmark_buoytke

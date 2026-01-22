!***********************************************************************
! GPU Kernel Benchmark: vbcu (s_vbcu)
!***********************************************************************
!
! Source: Src/vbcu.f90
! Description: Sets vertical boundary conditions for x-velocity component
!              by copying values from adjacent levels at bottom and top.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_vbcu
  use omp_lib
  implicit none

  integer :: ni, nj, nk
  real, allocatable :: uf(:,:,:)
  real, allocatable :: uf_ref(:,:,:)
  real, allocatable :: uf_input(:,:,:)

  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  real :: max_error, rel_error, tolerance
  integer :: error_count
  logical :: validation_passed
  integer :: iter, i, j, k

  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: vbcu'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  allocate(uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_input(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/uf_in.bin', uf_input, 0, ni+1, 0, nj+1, 1, nk)
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/uf_ref.bin', uf_ref, 0, ni+1, 0, nj+1, 1, nk)

  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    uf = uf_input
    call kernel_vbcu(ni, nj, nk, uf)
    !$acc wait
  end do

  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  do iter = 1, num_iterations
    uf = uf_input
    !$acc wait
    t_start = omp_get_wtime()
    call kernel_vbcu(ni, nj, nk, uf)
    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do
  t_avg = t_total / dble(num_iterations)

  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0
  do k = 1, nk
    do j = 1, nj-1
      do i = 1, ni
        rel_error = abs(uf(i,j,k) - uf_ref(i,j,k))
        if (abs(uf_ref(i,j,k)) > 1.0e-20) rel_error = rel_error / abs(uf_ref(i,j,k))
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do
  validation_passed = (error_count == 0)

  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,ES12.4)') ' Max error:    ', max_error
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  deallocate(uf, uf_ref, uf_input, times)
  if (.not. validation_passed) stop 1

contains

  subroutine kernel_vbcu(ni, nj, nk, uf)
    integer, intent(in) :: ni, nj, nk
    real, intent(inout) :: uf(0:ni+1, 0:nj+1, 1:nk)
    integer :: i, j, nkm1, nkm2

    nkm1 = nk - 1
    nkm2 = nk - 2

    ! Bottom boundary
    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni
        uf(i,j,1) = uf(i,j,2)
      end do
    end do
    !$acc end kernels

    ! Top boundary
    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni
        uf(i,j,nkm1) = uf(i,j,nkm2)
      end do
    end do
    !$acc end kernels

  end subroutine kernel_vbcu

  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    character(len=256) :: arg
    integer :: ios, nargs
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    nargs = command_argument_count()
    if (nargs >= 1) then
      call get_command_argument(1, arg)
      data_dir = trim(arg)
      return
    end if

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

  subroutine read_parameters(filename, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) stop 'Cannot open params file'
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

  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) stop 'Cannot open data file'
    read(10) arr
    close(10)
  end subroutine read_array_3d

end program kernel_benchmark_vbcu

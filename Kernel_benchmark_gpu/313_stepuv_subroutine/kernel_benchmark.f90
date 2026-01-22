!***********************************************************************
! GPU Kernel Benchmark: stepuv (s_stepuv)
!***********************************************************************
!
! Source: Src/stepuv.f90
! Description: Time integration of u and v velocity components using
!              forcing terms and acoustic mode contributions.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_stepuv
  use omp_lib
  implicit none

  ! Grid dimensions
  integer :: ni, nj, nk

  ! Time step
  real :: dts

  ! Arrays
  real, allocatable :: rst8u(:,:,:), rst8v(:,:,:)
  real, allocatable :: ufrc(:,:,:), vfrc(:,:,:)
  real, allocatable :: usml(:,:,:), vsml(:,:,:)
  real, allocatable :: uf(:,:,:), vf(:,:,:)
  real, allocatable :: uf_input(:,:,:), vf_input(:,:,:)
  real, allocatable :: uf_ref(:,:,:), vf_ref(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, dts)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: stepuv'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dts: ', dts
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rst8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(usml(0:ni+1, 0:nj+1, 1:nk))
  allocate(vsml(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_input(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf_input(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rst8u.bin', rst8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8v.bin', rst8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ufrc.bin', ufrc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc.bin', vfrc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/usml.bin', usml, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vsml.bin', vsml, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/uf_in.bin', uf_input, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vf_in.bin', vf_input, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/uf_ref.bin', uf_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vf_ref.bin', vf_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    uf = uf_input
    vf = vf_input
    call kernel_stepuv(dts, ni, nj, nk, rst8u, rst8v, ufrc, vfrc, usml, vsml, uf, vf)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    uf = uf_input
    vf = vf_input

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_stepuv(dts, ni, nj, nk, rst8u, rst8v, ufrc, vfrc, usml, vsml, uf, vf)
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

  ! Validate uf
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-1
        rel_error = abs(uf(i,j,k) - uf_ref(i,j,k))
        if (abs(uf_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(uf_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate vf
  do k = 2, nk-2
    do j = 2, nj-1
      do i = 2, ni-2
        rel_error = abs(vf(i,j,k) - vf_ref(i,j,k))
        if (abs(vf_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(vf_ref(i,j,k))
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
  deallocate(rst8u, rst8v, ufrc, vfrc, usml, vsml)
  deallocate(uf, vf, uf_input, vf_input, uf_ref, vf_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! GPU Kernel: stepuv (OpenACC)
  !=====================================================================
  subroutine kernel_stepuv(dts, ni, nj, nk, rst8u, rst8v, ufrc, vfrc, usml, vsml, uf, vf)
    implicit none

    real, intent(in) :: dts
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: rst8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ufrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: usml(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vsml(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: uf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vf(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    ! Update u velocity
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-1
          uf(i,j,k) = uf(i,j,k) + dts * (ufrc(i,j,k) + usml(i,j,k)) / rst8u(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Update v velocity
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-1
        !$acc loop independent
        do i = 2, ni-2
          vf(i,j,k) = vf(i,j,k) + dts * (vfrc(i,j,k) + vsml(i,j,k)) / rst8v(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_stepuv

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

  subroutine read_parameters(filename, ni, nj, nk, dts)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    real, intent(out) :: dts
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
          case ('dts')
            read(val, *) dts
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

end program kernel_benchmark_gpu_stepuv

!***********************************************************************
! GPU Kernel Benchmark: roughnxt (s_roughnxt)
!***********************************************************************
!
! Source: Src/roughnxt.f90
! Description: Calculates roughness lengths for momentum (z0m) and heat (z0h)
!              over ocean surfaces based on surface stress and wind speed.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_roughnxt
  use omp_lib
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup, num_iterations
  real(8) :: tolerance

  ! Dimensions
  integer :: ni, nj

  ! Physical constants from m_comphy
  real, parameter :: z0min = 1.5e-5

  ! Input arrays
  integer, allocatable :: land(:,:)
  real, allocatable :: va(:,:), cm(:,:)

  ! Input/Output arrays
  real, allocatable :: z0m(:,:), z0h(:,:)
  real, allocatable :: z0m_init(:,:), z0h_init(:,:)

  ! Reference arrays
  real, allocatable :: z0m_ref(:,:), z0h_ref(:,:)

  ! Timing variables
  real(8) :: start_time, end_time, elapsed_time
  real(8) :: total_time, min_time, max_time, avg_time
  real(8), allocatable :: times(:)

  ! Validation variables
  integer :: error_count, total_errors
  real(8) :: max_rel_error

  ! Loop variables
  integer :: iter, i, j

  ! Config file
  integer :: unit_conf, ios

  ! Read configuration
  unit_conf = 10
  open(unit_conf, file='benchmark.conf', status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print '(A)', 'Error: Cannot open benchmark.conf'
    stop 1
  end if
  read(unit_conf, '(A)') data_dir
  read(unit_conf, *) num_iterations
  read(unit_conf, *) num_warmup
  read(unit_conf, *) tolerance
  close(unit_conf)

  data_dir = trim(adjustl(data_dir))

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: roughnxt'
  write(*,'(A)') '=================================================='
  write(*,'(A,A)') ' Data directory: ', trim(data_dir)

  ! Read parameters
  call read_params(trim(data_dir) // '/params.txt')

  write(*,'(A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj
  write(*,'(A,I6)') ' Warmup iterations: ', num_warmup
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,ES12.4)') ' Tolerance: ', tolerance
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(land(0:ni+1, 0:nj+1))
  allocate(va(0:ni+1, 0:nj+1))
  allocate(cm(0:ni+1, 0:nj+1))
  allocate(z0m(0:ni+1, 0:nj+1))
  allocate(z0h(0:ni+1, 0:nj+1))
  allocate(z0m_init(0:ni+1, 0:nj+1))
  allocate(z0h_init(0:ni+1, 0:nj+1))
  allocate(z0m_ref(0:ni+1, 0:nj+1))
  allocate(z0h_ref(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Read input arrays
  write(*,'(A)') ' Loading input data...'
  call read_array_2d_int(trim(data_dir) // '/land.bin', land, ni, nj)
  call read_array_2d(trim(data_dir) // '/va.bin', va, ni, nj)
  call read_array_2d(trim(data_dir) // '/cm.bin', cm, ni, nj)
  call read_array_2d(trim(data_dir) // '/z0m_in.bin', z0m_init, ni, nj)
  call read_array_2d(trim(data_dir) // '/z0h_in.bin', z0h_init, ni, nj)
  call read_array_2d(trim(data_dir) // '/z0m_ref.bin', z0m_ref, ni, nj)
  call read_array_2d(trim(data_dir) // '/z0h_ref.bin', z0h_ref, ni, nj)
  write(*,'(A)') ' Data loaded successfully'

  ! Warmup iterations (includes GPU JIT compilation)
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, num_warmup
    call reset_arrays()
    call kernel_roughnxt()
    !$acc wait
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  total_time = 0.0d0
  min_time = huge(1.0d0)
  max_time = 0.0d0

  do iter = 1, num_iterations
    call reset_arrays()

    !$acc wait
    start_time = omp_get_wtime()
    call kernel_roughnxt()
    !$acc wait
    end_time = omp_get_wtime()

    elapsed_time = end_time - start_time
    times(iter) = elapsed_time
    total_time = total_time + elapsed_time
    min_time = min(min_time, elapsed_time)
    max_time = max(max_time, elapsed_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results
  call validate_results(total_errors, max_rel_error)

  ! Report results
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', avg_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', total_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', min_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', max_time * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_rel_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', total_errors
  if (total_errors == 0) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! Cleanup
  deallocate(land, va, cm, z0m, z0h, z0m_init, z0h_init, z0m_ref, z0h_ref, times)

  if (total_errors > 0) stop 1

contains

  !---------------------------------------------------------------------
  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    integer :: unit_num, ios_local
    character(len=256) :: line
    character(len=64) :: key
    character(len=128) :: value_str
    integer :: eq_pos

    unit_num = 20
    open(unit_num, file=filename, status='old', action='read', iostat=ios_local)
    if (ios_local /= 0) then
      print *, 'Error opening params file: ', trim(filename)
      stop 1
    end if

    do
      read(unit_num, '(A)', iostat=ios_local) line
      if (ios_local /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      eq_pos = index(line, '=')
      if (eq_pos == 0) cycle
      key = adjustl(line(1:eq_pos-1))
      value_str = adjustl(line(eq_pos+1:))

      select case (trim(key))
        case ('ni'); read(value_str, *) ni
        case ('nj'); read(value_str, *) nj
      end select
    end do
    close(unit_num)
  end subroutine read_params

  !---------------------------------------------------------------------
  subroutine read_array_2d(filename, arr, ni_loc, nj_loc)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni_loc, nj_loc
    real, intent(out) :: arr(0:ni_loc+1, 0:nj_loc+1)
    integer :: unit_num, ios_local

    unit_num = 30
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', convert='big_endian', iostat=ios_local)
    if (ios_local /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_2d

  !---------------------------------------------------------------------
  subroutine read_array_2d_int(filename, arr, ni_loc, nj_loc)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni_loc, nj_loc
    integer, intent(out) :: arr(0:ni_loc+1, 0:nj_loc+1)
    integer :: unit_num, ios_local

    unit_num = 30
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', convert='big_endian', iostat=ios_local)
    if (ios_local /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_2d_int

  !---------------------------------------------------------------------
  subroutine reset_arrays()
    z0m = z0m_init
    z0h = z0h_init
  end subroutine reset_arrays

  !---------------------------------------------------------------------
  subroutine validate_results(err_count, max_err)
    integer, intent(out) :: err_count
    real(8), intent(out) :: max_err

    integer :: i_loc, j_loc, err_z0m, err_z0h
    real(8) :: rel_err, abs_val, max_err_z0m, max_err_z0h

    err_z0m = 0
    err_z0h = 0
    max_err_z0m = 0.0d0
    max_err_z0h = 0.0d0

    ! Validate z0m
    do j_loc = 1, nj-1
      do i_loc = 1, ni-1
        abs_val = abs(dble(z0m_ref(i_loc,j_loc)))
        if (abs_val > 1.0d-30) then
          rel_err = abs(dble(z0m(i_loc,j_loc)) - dble(z0m_ref(i_loc,j_loc))) / abs_val
        else
          rel_err = abs(dble(z0m(i_loc,j_loc)) - dble(z0m_ref(i_loc,j_loc)))
        end if
        max_err_z0m = max(max_err_z0m, rel_err)
        if (rel_err > tolerance) err_z0m = err_z0m + 1
      end do
    end do

    ! Validate z0h
    do j_loc = 1, nj-1
      do i_loc = 1, ni-1
        abs_val = abs(dble(z0h_ref(i_loc,j_loc)))
        if (abs_val > 1.0d-30) then
          rel_err = abs(dble(z0h(i_loc,j_loc)) - dble(z0h_ref(i_loc,j_loc))) / abs_val
        else
          rel_err = abs(dble(z0h(i_loc,j_loc)) - dble(z0h_ref(i_loc,j_loc)))
        end if
        max_err_z0h = max(max_err_z0h, rel_err)
        if (rel_err > tolerance) err_z0h = err_z0h + 1
      end do
    end do

    write(*,'(A,I0,A,ES12.4)') ' z0m: errors=', err_z0m, ', max_rel_err=', max_err_z0m
    write(*,'(A,I0,A,ES12.4)') ' z0h: errors=', err_z0h, ', max_rel_err=', max_err_z0h

    err_count = err_z0m + err_z0h
    max_err = max(max_err_z0m, max_err_z0h)
  end subroutine validate_results

  !---------------------------------------------------------------------
  subroutine kernel_roughnxt()
    integer :: i_k, j_k
    real :: ust

    !$acc kernels
    !$acc loop independent
    do j_k = 1, nj-1
      !$acc loop independent private(ust)
      do i_k = 1, ni-1
        if (land(i_k,j_k) < 3) then
          ust = cm(i_k,j_k) * va(i_k,j_k)

          if (ust < 1.08e0) then
            z0m(i_k,j_k) = max(-34.7e-6 + 8.28e-4 * ust, z0min)
          else
            z0m(i_k,j_k) = max(-0.277e-2 + 3.39e-3 * ust, z0min)
          end if

          z0h(i_k,j_k) = z0m(i_k,j_k)
        end if
      end do
    end do
    !$acc end kernels

  end subroutine kernel_roughnxt

end program kernel_benchmark_gpu_roughnxt

program kernel_benchmark
  use omp_lib
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
  integer :: iter

  ! Config file
  integer :: unit_conf

  ! Read configuration
  unit_conf = 10
  open(unit_conf, file='benchmark.conf', status='old', action='read')
  read(unit_conf, '(A)') data_dir
  read(unit_conf, *) num_warmup
  read(unit_conf, *) num_iterations
  read(unit_conf, *) tolerance
  close(unit_conf)

  data_dir = trim(adjustl(data_dir))
  print '(A)', '=== Roughnxt Kernel Benchmark ==='
  print '(A,A)', 'Data directory: ', trim(data_dir)
  print '(A,I0)', 'Warmup iterations: ', num_warmup
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,ES10.2)', 'Tolerance: ', tolerance

  ! Read parameters
  call read_params(trim(data_dir) // '/params.txt')

  print '(A)', ''
  print '(A)', '--- Parameters ---'
  print '(A,I0,A,I0)', 'ni x nj = ', ni, ' x ', nj

  ! Allocate arrays
  allocate(land(0:ni+1, 0:nj+1))
  allocate(va(0:ni+1, 0:nj+1))
  allocate(cm(0:ni+1, 0:nj+1))
  allocate(z0m(0:ni+1, 0:nj+1))
  allocate(z0h(0:ni+1, 0:nj+1))
  allocate(z0m_ref(0:ni+1, 0:nj+1))
  allocate(z0h_ref(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Read input arrays
  print '(A)', ''
  print '(A)', '--- Loading input data ---'
  call read_array_2d_int(trim(data_dir) // '/land.bin', land, ni, nj)
  call read_array_2d(trim(data_dir) // '/va.bin', va, ni, nj)
  call read_array_2d(trim(data_dir) // '/cm.bin', cm, ni, nj)
  call read_array_2d(trim(data_dir) // '/z0m_in.bin', z0m, ni, nj)
  call read_array_2d(trim(data_dir) // '/z0h_in.bin', z0h, ni, nj)
  call read_array_2d(trim(data_dir) // '/z0m_ref.bin', z0m_ref, ni, nj)
  call read_array_2d(trim(data_dir) // '/z0h_ref.bin', z0h_ref, ni, nj)
  print '(A)', 'Data loaded successfully'

  ! Warmup iterations
  print '(A)', ''
  print '(A)', '--- Warmup ---'
  do iter = 1, num_warmup
    ! Reset input arrays
    call read_array_2d(trim(data_dir) // '/z0m_in.bin', z0m, ni, nj)
    call read_array_2d(trim(data_dir) // '/z0h_in.bin', z0h, ni, nj)
    call run_kernel()
  end do
  print '(A,I0,A)', 'Completed ', num_warmup, ' warmup iterations'

  ! Benchmark iterations
  print '(A)', ''
  print '(A)', '--- Benchmark ---'
  total_time = 0.0d0
  min_time = huge(1.0d0)
  max_time = 0.0d0

  do iter = 1, num_iterations
    ! Reset input arrays
    call read_array_2d(trim(data_dir) // '/z0m_in.bin', z0m, ni, nj)
    call read_array_2d(trim(data_dir) // '/z0h_in.bin', z0h, ni, nj)

    start_time = omp_get_wtime()
    call run_kernel()
    end_time = omp_get_wtime()

    elapsed_time = end_time - start_time
    times(iter) = elapsed_time
    total_time = total_time + elapsed_time
    min_time = min(min_time, elapsed_time)
    max_time = max(max_time, elapsed_time)
  end do

  avg_time = total_time / num_iterations

  print '(A)', ''
  print '(A)', '=== Timing Results ==='
  print '(A,ES12.4,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,ES12.4,A)', 'Min time:     ', min_time * 1000.0d0, ' ms'
  print '(A,ES12.4,A)', 'Max time:     ', max_time * 1000.0d0, ' ms'
  print '(A,ES12.4,A)', 'Total time:   ', total_time * 1000.0d0, ' ms'

  ! Validation
  print '(A)', ''
  print '(A)', '=== Validation ==='
  total_errors = 0

  call validate_array('z0m', z0m, z0m_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array('z0h', z0h, z0h_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  print '(A)', ''
  if (total_errors == 0) then
    print '(A)', '*** VALIDATION PASSED ***'
  else
    print '(A,I0,A)', '*** VALIDATION FAILED: ', total_errors, ' errors ***'
  end if

  ! Cleanup
  deallocate(land, va, cm, z0m, z0h, z0m_ref, z0h_ref, times)

contains

  subroutine run_kernel()
    integer :: i, j
    real :: ust

    !$omp parallel default(shared)

    !$omp do schedule(runtime) private(i,j,ust)
    do j = 1, nj-1
      do i = 1, ni-1
        if (land(i,j) < 3) then
          ust = cm(i,j) * va(i,j)

          if (ust < 1.08e0) then
            z0m(i,j) = max(-34.7e-6 + 8.28e-4 * ust, z0min)
          else
            z0m(i,j) = max(-0.277e-2 + 3.39e-3 * ust, z0min)
          end if

          z0h(i,j) = z0m(i,j)
        end if
      end do
    end do
    !$omp end do

    !$omp end parallel
  end subroutine run_kernel

  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    integer :: unit_num, ios
    character(len=256) :: line
    character(len=64) :: key
    character(len=128) :: value_str
    integer :: eq_pos

    unit_num = 20
    open(unit_num, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening params file: ', trim(filename)
      stop 1
    end if

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
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

  subroutine read_array_2d(filename, arr, ni, nj)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni, nj
    real, intent(out) :: arr(0:ni+1, 0:nj+1)
    integer :: unit_num, ios

    unit_num = 30
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', convert='big_endian', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_2d

  subroutine read_array_2d_int(filename, arr, ni, nj)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni, nj
    integer, intent(out) :: arr(0:ni+1, 0:nj+1)
    integer :: unit_num, ios

    unit_num = 30
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', convert='big_endian', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_2d_int

  subroutine validate_array(name, arr, ref, tol, err_count, max_err)
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(0:ni+1, 0:nj+1)
    real, intent(in) :: ref(0:ni+1, 0:nj+1)
    real(8), intent(in) :: tol
    integer, intent(out) :: err_count
    real(8), intent(out) :: max_err

    integer :: i, j
    real(8) :: rel_err, abs_val

    err_count = 0
    max_err = 0.0d0

    do j = 1, nj-1
      do i = 1, ni-1
        abs_val = abs(dble(ref(i,j)))
        if (abs_val > 1.0d-30) then
          rel_err = abs(dble(arr(i,j)) - dble(ref(i,j))) / abs_val
        else
          rel_err = abs(dble(arr(i,j)) - dble(ref(i,j)))
        end if
        max_err = max(max_err, rel_err)
        if (rel_err > tol) err_count = err_count + 1
      end do
    end do

    print '(A,A,A,I0,A,ES12.4)', 'Array ', trim(name), ': errors=', err_count, &
          ', max_rel_err=', max_err
  end subroutine validate_array

end program kernel_benchmark

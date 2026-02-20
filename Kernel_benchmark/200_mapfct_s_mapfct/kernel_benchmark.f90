program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup, num_iterations
  real(8) :: tolerance

  ! Dimensions
  integer :: ni, nj

  ! Projection parameters
  integer :: mpopt
  real :: disr, dxiv, dyiv
  real :: dxv625, dyv625

  ! Mathematical constants (from m_commath)
  real, parameter :: d2r = 3.141592e0 / 180.e0
  real, parameter :: eps = 1.0e-20

  ! cpj array (map projection parameters)
  real :: cpj(1:7)

  ! Input arrays
  real, allocatable :: x(:)
  real, allocatable :: lat(:,:)
  real, allocatable :: tmp1_in(:,:)

  ! Output arrays
  real, allocatable :: mf(:,:), mf8u(:,:), mf8v(:,:)
  real, allocatable :: rmf(:,:,:), rmf8u(:,:,:), rmf8v(:,:,:)
  real, allocatable :: tmp1(:,:)

  ! Reference arrays
  real, allocatable :: mf_ref(:,:), mf8u_ref(:,:), mf8v_ref(:,:)
  real, allocatable :: rmf_ref(:,:,:), rmf8u_ref(:,:,:), rmf8v_ref(:,:,:)
  real, allocatable :: tmp1_ref(:,:)

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
  print '(A)', '=== Mapfct Kernel Benchmark ==='
  print '(A,A)', 'Data directory: ', trim(data_dir)
  print '(A,I0)', 'Warmup iterations: ', num_warmup
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,ES10.2)', 'Tolerance: ', tolerance

  ! Read parameters
  call read_params(trim(data_dir) // '/params.txt')

  print '(A)', ''
  print '(A)', '--- Parameters ---'
  print '(A,I0,A,I0)', 'ni x nj = ', ni, ' x ', nj
  print '(A,I0)', 'mpopt = ', mpopt

  ! Allocate arrays
  allocate(x(0:ni+1))
  allocate(lat(0:ni+1, 0:nj+1))
  allocate(tmp1_in(0:ni+1, 0:nj+1))
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(mf8u(0:ni+1, 0:nj+1))
  allocate(mf8v(0:ni+1, 0:nj+1))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rmf8u(0:ni+1, 0:nj+1, 1:3))
  allocate(rmf8v(0:ni+1, 0:nj+1, 1:3))
  allocate(tmp1(0:ni+1, 0:nj+1))
  allocate(mf_ref(0:ni+1, 0:nj+1))
  allocate(mf8u_ref(0:ni+1, 0:nj+1))
  allocate(mf8v_ref(0:ni+1, 0:nj+1))
  allocate(rmf_ref(0:ni+1, 0:nj+1, 1:4))
  allocate(rmf8u_ref(0:ni+1, 0:nj+1, 1:3))
  allocate(rmf8v_ref(0:ni+1, 0:nj+1, 1:3))
  allocate(tmp1_ref(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Read input arrays
  print '(A)', ''
  print '(A)', '--- Loading input data ---'
  call read_array_1d_7(trim(data_dir) // '/cpj.bin', cpj)
  call read_array_1d(trim(data_dir) // '/x.bin', x, ni)
  call read_array_2d(trim(data_dir) // '/lat.bin', lat, ni, nj)
  call read_array_2d(trim(data_dir) // '/tmp1_in.bin', tmp1_in, ni, nj)

  ! Read reference arrays
  call read_array_2d(trim(data_dir) // '/mf_ref.bin', mf_ref, ni, nj)
  call read_array_2d(trim(data_dir) // '/mf8u_ref.bin', mf8u_ref, ni, nj)
  call read_array_2d(trim(data_dir) // '/mf8v_ref.bin', mf8v_ref, ni, nj)
  call read_array_3d_4(trim(data_dir) // '/rmf_ref.bin', rmf_ref, ni, nj)
  call read_array_3d_3(trim(data_dir) // '/rmf8u_ref.bin', rmf8u_ref, ni, nj)
  call read_array_3d_3(trim(data_dir) // '/rmf8v_ref.bin', rmf8v_ref, ni, nj)
  call read_array_2d(trim(data_dir) // '/tmp1_ref.bin', tmp1_ref, ni, nj)
  print '(A)', 'Data loaded successfully'

  ! Warmup iterations
  print '(A)', ''
  print '(A)', '--- Warmup ---'
  do iter = 1, num_warmup
    mf = 0.0
    mf8u = 0.0
    mf8v = 0.0
    rmf = 0.0
    rmf8u = 0.0
    rmf8v = 0.0
    tmp1 = tmp1_in
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
    mf = 0.0
    mf8u = 0.0
    mf8v = 0.0
    rmf = 0.0
    rmf8u = 0.0
    rmf8v = 0.0
    tmp1 = tmp1_in

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

  call validate_array_2d('mf', mf, mf_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array_2d('mf8u', mf8u, mf8u_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array_2d('mf8v', mf8v, mf8v_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_array_3d_4('rmf', rmf, rmf_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_rmf8u('rmf8u', rmf8u, rmf8u_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_rmf8v('rmf8v', rmf8v, rmf8v_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count
  call validate_tmp1('tmp1', tmp1, tmp1_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  print '(A)', ''
  if (total_errors == 0) then
    print '(A)', '*** VALIDATION PASSED ***'
  else
    print '(A,I0,A)', '*** VALIDATION FAILED: ', total_errors, ' errors ***'
  end if

  ! Cleanup
  deallocate(x, lat, tmp1_in, mf, mf8u, mf8v, rmf, rmf8u, rmf8v, tmp1)
  deallocate(mf_ref, mf8u_ref, mf8v_ref, rmf_ref, rmf8u_ref, rmf8v_ref, tmp1_ref, times)

contains

  subroutine run_kernel()
    integer :: i, j

    !$omp parallel default(shared)

    ! Calculate map scale factor at scalar points (mpopt=0: spherical)
    if (mpopt == 0 .or. mpopt == 10) then
      !$omp do schedule(runtime) private(i,j)
      do j = 0, nj
        do i = 0, ni
          mf(i,j) = 1.0e0 / (cos(lat(i,j) * d2r) + eps)
        end do
      end do
      !$omp end do
    else if (mpopt == 4) then
      !$omp do schedule(runtime) private(i,j)
      do j = 0, nj
        do i = 0, ni
          mf(i,j) = 1.0e0
        end do
      end do
      !$omp end do
    else if (mpopt == 5) then
      !$omp do schedule(runtime) private(i,j)
      do j = 0, nj
        do i = 0, ni
          mf(i,j) = 1.0e0 / (x(i) + disr)
        end do
      end do
      !$omp end do
    end if

    ! Calculate map scale factors at u and v points
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 1, ni
        mf8u(i,j) = 0.5e0 * (mf(i-1,j) + mf(i,j))
      end do
    end do
    !$omp end do

    !$omp do schedule(runtime) private(i,j)
    do j = 1, nj
      do i = 0, ni
        mf8v(i,j) = 0.5e0 * (mf(i,j-1) + mf(i,j))
      end do
    end do
    !$omp end do

    ! Calculate rmf
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 0, ni
        rmf(i,j,1) = mf(i,j) * mf(i,j)
        rmf(i,j,2) = 1.0e0 / mf(i,j)
        rmf(i,j,3) = rmf(i,j,2) * rmf(i,j,2)
        rmf(i,j,4) = sqrt(rmf(i,j,2))
      end do
    end do
    !$omp end do

    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 1, ni
        rmf8u(i,j,1) = mf8u(i,j) * mf8u(i,j)
        rmf8u(i,j,2) = 1.0e0 / mf8u(i,j)
      end do
    end do
    !$omp end do

    !$omp do schedule(runtime) private(i,j)
    do j = 1, nj
      do i = 0, ni
        rmf8v(i,j,1) = mf8v(i,j) * mf8v(i,j)
        rmf8v(i,j,2) = 1.0e0 / mf8v(i,j)
      end do
    end do
    !$omp end do

    ! Calculate differential of map scale factor x 0.0625
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 1, ni
        tmp1(i,j) = (mf(i,j) - mf(i-1,j)) * dxv625
      end do
    end do
    !$omp end do

    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 1, ni-1
        rmf8u(i,j,3) = tmp1(i,j) + tmp1(i+1,j)
      end do
    end do
    !$omp end do

    !$omp do schedule(runtime) private(i,j)
    do j = 1, nj
      do i = 0, ni
        tmp1(i,j) = (mf(i,j) - mf(i,j-1)) * dyv625
      end do
    end do
    !$omp end do

    !$omp do schedule(runtime) private(i,j)
    do j = 1, nj-1
      do i = 0, ni
        rmf8v(i,j,3) = tmp1(i,j) + tmp1(i,j+1)
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
        case ('mpopt'); read(value_str, *) mpopt
        case ('disr'); read(value_str, *) disr
        case ('dxiv'); read(value_str, *) dxiv
        case ('dyiv'); read(value_str, *) dyiv
        case ('dxv625'); read(value_str, *) dxv625
        case ('dyv625'); read(value_str, *) dyv625
      end select
    end do
    close(unit_num)
  end subroutine read_params

  subroutine read_array_1d(filename, arr, ni)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni
    real, intent(out) :: arr(0:ni+1)
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
  end subroutine read_array_1d

  subroutine read_array_1d_7(filename, arr)
    character(len=*), intent(in) :: filename
    real, intent(out) :: arr(1:7)
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
  end subroutine read_array_1d_7

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

  subroutine read_array_3d_3(filename, arr, ni, nj)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni, nj
    real, intent(out) :: arr(0:ni+1, 0:nj+1, 1:3)
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
  end subroutine read_array_3d_3

  subroutine read_array_3d_4(filename, arr, ni, nj)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni, nj
    real, intent(out) :: arr(0:ni+1, 0:nj+1, 1:4)
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
  end subroutine read_array_3d_4

  subroutine validate_array_2d(name, arr, ref, tol, err_count, max_err)
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

    do j = 0, nj
      do i = 0, ni
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
  end subroutine validate_array_2d

  subroutine validate_array_3d_3(name, arr, ref, tol, err_count, max_err)
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:3)
    real(8), intent(in) :: tol
    integer, intent(out) :: err_count
    real(8), intent(out) :: max_err

    integer :: i, j, k
    real(8) :: rel_err, abs_val

    err_count = 0
    max_err = 0.0d0

    do k = 1, 3
      do j = 0, nj
        do i = 0, ni
          abs_val = abs(dble(ref(i,j,k)))
          if (abs_val > 1.0d-30) then
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k))) / abs_val
          else
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k)))
          end if
          max_err = max(max_err, rel_err)
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

    print '(A,A,A,I0,A,ES12.4)', 'Array ', trim(name), ': errors=', err_count, &
          ', max_rel_err=', max_err
  end subroutine validate_array_3d_3

  subroutine validate_array_3d_4(name, arr, ref, tol, err_count, max_err)
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:4)
    real(8), intent(in) :: tol
    integer, intent(out) :: err_count
    real(8), intent(out) :: max_err

    integer :: i, j, k
    real(8) :: rel_err, abs_val

    err_count = 0
    max_err = 0.0d0

    do k = 1, 4
      do j = 0, nj
        do i = 0, ni
          abs_val = abs(dble(ref(i,j,k)))
          if (abs_val > 1.0d-30) then
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k))) / abs_val
          else
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k)))
          end if
          max_err = max(max_err, rel_err)
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

    print '(A,A,A,I0,A,ES12.4)', 'Array ', trim(name), ': errors=', err_count, &
          ', max_rel_err=', max_err
  end subroutine validate_array_3d_4

  subroutine validate_rmf8u(name, arr, ref, tol, err_count, max_err)
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:3)
    real(8), intent(in) :: tol
    integer, intent(out) :: err_count
    real(8), intent(out) :: max_err

    integer :: i, j, k
    real(8) :: rel_err, abs_val

    err_count = 0
    max_err = 0.0d0

    ! Components 1,2: j=0,nj, i=1,ni
    do k = 1, 2
      do j = 0, nj
        do i = 1, ni
          abs_val = abs(dble(ref(i,j,k)))
          if (abs_val > 1.0d-30) then
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k))) / abs_val
          else
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k)))
          end if
          max_err = max(max_err, rel_err)
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

    ! Component 3: j=0,nj, i=1,ni-1
    k = 3
    do j = 0, nj
      do i = 1, ni-1
        abs_val = abs(dble(ref(i,j,k)))
        if (abs_val > 1.0d-30) then
          rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k))) / abs_val
        else
          rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k)))
        end if
        max_err = max(max_err, rel_err)
        if (rel_err > tol) err_count = err_count + 1
      end do
    end do

    print '(A,A,A,I0,A,ES12.4)', 'Array ', trim(name), ': errors=', err_count, &
          ', max_rel_err=', max_err
  end subroutine validate_rmf8u

  subroutine validate_rmf8v(name, arr, ref, tol, err_count, max_err)
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:3)
    real(8), intent(in) :: tol
    integer, intent(out) :: err_count
    real(8), intent(out) :: max_err

    integer :: i, j, k
    real(8) :: rel_err, abs_val

    err_count = 0
    max_err = 0.0d0

    ! Components 1,2: j=1,nj, i=0,ni
    do k = 1, 2
      do j = 1, nj
        do i = 0, ni
          abs_val = abs(dble(ref(i,j,k)))
          if (abs_val > 1.0d-30) then
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k))) / abs_val
          else
            rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k)))
          end if
          max_err = max(max_err, rel_err)
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

    ! Component 3: j=1,nj-1, i=0,ni
    k = 3
    do j = 1, nj-1
      do i = 0, ni
        abs_val = abs(dble(ref(i,j,k)))
        if (abs_val > 1.0d-30) then
          rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k))) / abs_val
        else
          rel_err = abs(dble(arr(i,j,k)) - dble(ref(i,j,k)))
        end if
        max_err = max(max_err, rel_err)
        if (rel_err > tol) err_count = err_count + 1
      end do
    end do

    print '(A,A,A,I0,A,ES12.4)', 'Array ', trim(name), ': errors=', err_count, &
          ', max_rel_err=', max_err
  end subroutine validate_rmf8v

  subroutine validate_tmp1(name, arr, ref, tol, err_count, max_err)
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

    ! Final tmp1 computation: j=1,nj, i=0,ni
    do j = 1, nj
      do i = 0, ni
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
  end subroutine validate_tmp1

end program kernel_benchmark

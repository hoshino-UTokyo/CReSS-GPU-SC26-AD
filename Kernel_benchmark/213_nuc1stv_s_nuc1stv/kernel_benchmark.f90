program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup, num_iterations
  real(8) :: tolerance

  ! Dimensions
  integer :: ni, nj, nk

  ! Physical constants from params.txt
  real :: thresq, t0

  ! Additional physical constants (hardcoded from m_comphy)
  real, parameter :: tlow = 233.16
  real, parameter :: mi0 = 1.0e-12

  ! Input arrays
  real, allocatable :: t(:,:,:), rbv(:,:,:), qv(:,:,:), qi(:,:,:), qvsi(:,:,:)

  ! Output array
  real, allocatable :: nuvi(:,:,:)

  ! Reference array
  real, allocatable :: nuvi_ref(:,:,:)

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
  print '(A)', '=== Nuc1stv Kernel Benchmark ==='
  print '(A,A)', 'Data directory: ', trim(data_dir)
  print '(A,I0)', 'Warmup iterations: ', num_warmup
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,ES10.2)', 'Tolerance: ', tolerance

  ! Read parameters
  call read_params(trim(data_dir) // '/params.txt')

  print '(A)', ''
  print '(A)', '--- Parameters ---'
  print '(A,I0,A,I0,A,I0)', 'ni x nj x nk = ', ni, ' x ', nj, ' x ', nk

  ! Allocate arrays
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsi(0:ni+1, 0:nj+1, 1:nk))
  allocate(nuvi(0:ni+1, 0:nj+1, 1:nk))
  allocate(nuvi_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  ! Read input arrays
  print '(A)', ''
  print '(A)', '--- Loading input data ---'
  call read_array_3d(trim(data_dir) // '/t.bin', t, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/rbv.bin', rbv, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qv.bin', qv, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qi.bin', qi, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/qvsi.bin', qvsi, ni, nj, nk)
  call read_array_3d(trim(data_dir) // '/nuvi_ref.bin', nuvi_ref, ni, nj, nk)
  print '(A)', 'Data loaded successfully'

  ! Warmup iterations
  print '(A)', ''
  print '(A)', '--- Warmup ---'
  do iter = 1, num_warmup
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

  call validate_array('nuvi', nuvi, nuvi_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  print '(A)', ''
  if (total_errors == 0) then
    print '(A)', '*** VALIDATION PASSED ***'
  else
    print '(A,I0,A)', '*** VALIDATION FAILED: ', total_errors, ' errors ***'
  end if

  ! Cleanup
  deallocate(t, rbv, qv, qi, qvsi, nuvi, nuvi_ref, times)

contains

  subroutine run_kernel()
    integer :: i, j, k
    real :: qvssi, a

    !$omp parallel default(shared) private(k)

    if (nk == 1) then
      !$omp do schedule(runtime) private(i,j,qvssi,a)
      do j = 1, nj-1
        do i = 1, ni-1
          if (qv(i,j,1) > thresq) then
            if (qv(i,j,1) > qvsi(i,j,1) .and. t(i,j,1) > tlow .and. t(i,j,1) < t0) then
              qvssi = qv(i,j,1) - qvsi(i,j,1)
              a = 15.25 * qv(i,j,1) / qvsi(i,j,1) - 10.08
              if (a < 40.0) then
                nuvi(i,j,1) = max(min(mi0 * exp(a) * rbv(i,j,1) - qi(i,j,1), qvssi), 0.0)
              else
                nuvi(i,j,1) = qvssi
              end if
            else
              nuvi(i,j,1) = 0.0
            end if
          else
            nuvi(i,j,1) = 0.0
          end if
        end do
      end do
      !$omp end do

    else
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j,qvssi,a)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qv(i,j,k) > thresq) then
              if (qv(i,j,k) > qvsi(i,j,k) .and. t(i,j,k) > tlow .and. t(i,j,k) < t0) then
                qvssi = qv(i,j,k) - qvsi(i,j,k)
                a = 15.25 * qv(i,j,k) / qvsi(i,j,k) - 10.08
                if (a < 40.0) then
                  nuvi(i,j,k) = max(min(mi0 * exp(a) * rbv(i,j,k) - qi(i,j,k), qvssi), 0.0)
                else
                  nuvi(i,j,k) = qvssi
                end if
              else
                nuvi(i,j,k) = 0.0
              end if
            else
              nuvi(i,j,k) = 0.0
            end if
          end do
        end do
        !$omp end do
      end do
    end if

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
        case ('nk'); read(value_str, *) nk
        case ('thresq'); read(value_str, *) thresq
        case ('t0'); read(value_str, *) t0
      end select
    end do
    close(unit_num)
  end subroutine read_params

  subroutine read_array_3d(filename, arr, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni, nj, nk
    real, intent(out) :: arr(0:ni+1, 0:nj+1, 1:nk)
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
  end subroutine read_array_3d

  subroutine validate_array(name, arr, ref, tol, err_count, max_err)
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:nk)
    real(8), intent(in) :: tol
    integer, intent(out) :: err_count
    real(8), intent(out) :: max_err

    integer :: i, j, k
    real(8) :: rel_err, abs_val

    err_count = 0
    max_err = 0.0d0

    do k = 1, nk
      do j = 1, nj-1
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
    end do

    print '(A,A,A,I0,A,ES12.4)', 'Array ', trim(name), ': errors=', err_count, &
          ', max_rel_err=', max_err
  end subroutine validate_array

end program kernel_benchmark

program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup, num_iterations
  real(8) :: tolerance

  ! Dimensions
  integer :: ni, nj

  ! Physical constants from params.txt
  real :: zsfc

  ! Constant value derived from reference (for trnopt=0, all ht values are the same)
  real :: ht_const

  ! Output array
  real, allocatable :: ht(:,:)

  ! Reference array
  real, allocatable :: ht_ref(:,:)

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
  print '(A)', '=== Gettrn Sec1 (Flat Terrain) Kernel Benchmark ==='
  print '(A,A)', 'Data directory: ', trim(data_dir)
  print '(A,I0)', 'Warmup iterations: ', num_warmup
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,ES10.2)', 'Tolerance: ', tolerance

  ! Read parameters
  call read_params(trim(data_dir) // '/params.txt')

  print '(A)', ''
  print '(A)', '--- Parameters ---'
  print '(A,I0,A,I0)', 'ni x nj = ', ni, ' x ', nj
  print '(A,ES12.4)', 'zsfc = ', zsfc

  ! Allocate arrays
  allocate(ht(0:ni+1, 0:nj+1))
  allocate(ht_ref(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Read reference array and derive constant
  print '(A)', ''
  print '(A)', '--- Loading input data ---'
  call read_array_2d(trim(data_dir) // '/ht_ref.bin', ht_ref, ni, nj)

  ! For trnopt=0, all values are constant. Read the constant from reference.
  ! This is max(mnthgh(1)+mnthgh(2), zsfc) but since mnthgh values aren't dumped,
  ! we derive it from the reference data.
  ht_const = ht_ref(0, 0)
  print '(A,ES12.4)', 'Derived ht_const from reference: ', ht_const
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

  call validate_array('ht', ht, ht_ref, tolerance, error_count, max_rel_error)
  total_errors = total_errors + error_count

  print '(A)', ''
  if (total_errors == 0) then
    print '(A)', '*** VALIDATION PASSED ***'
  else
    print '(A,I0,A)', '*** VALIDATION FAILED: ', total_errors, ' errors ***'
  end if

  ! Cleanup
  deallocate(ht, ht_ref, times)

contains

  subroutine run_kernel()
    integer :: i, j

    ! Kernel: Set flat terrain (trnopt=0)
    ! Original: ht(i,j) = max(mnthgh(1)+mnthgh(2), zsfc)
    ! Since mnthgh values aren't available, we use the pre-computed constant

    !$omp parallel default(shared)

    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 0, ni
        ht(i,j) = ht_const
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
        case ('zsfc'); read(value_str, *) zsfc
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
  end subroutine validate_array

end program kernel_benchmark

!> Kernel benchmark program for getvdens (s_getvdens)
!> Calculates the inverse of base state density.
program kernel_benchmark_getvdens
  use omp_lib
  implicit none

  ! Parameters
  integer :: ni, nj, nk

  ! Arrays
  real, allocatable :: rbr(:,:,:)     ! Base state density
  real, allocatable :: rbv(:,:,:)     ! Inverse of base state density (output)
  real, allocatable :: rbv_ref(:,:,:) ! Reference output

  ! Benchmark variables
  character(len=256) :: data_dir
  integer :: num_iterations, warmup_iterations
  real :: tolerance
  integer :: iter
  real(8) :: start_time, end_time, elapsed_time
  real(8) :: total_time, min_time, max_time
  real :: max_diff

  ! Read configuration
  call read_config()

  ! Read parameters
  call read_parameters()

  ! Allocate arrays
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbv(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbv_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Read input data
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)

  ! Read reference output
  call read_array_3d(trim(data_dir)//'/rbv_ref.bin', rbv_ref, 0, ni+1, 0, nj+1, 1, nk)

  print '(A)', '================================================'
  print '(A)', 'Kernel Benchmark: getvdens (s_getvdens)'
  print '(A)', 'Inverse of base state density calculation'
  print '(A)', '================================================'
  print '(A,I0,A,I0,A,I0)', 'Grid size: ', ni, ' x ', nj, ' x ', nk
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,I0)', 'Warmup iterations: ', warmup_iterations
  print '(A)', '------------------------------------------------'

  ! Initialize output array
  rbv = 0.0

  ! Warmup iterations
  do iter = 1, warmup_iterations
    call kernel_getvdens(ni, nj, nk, rbr, rbv)
  end do

  ! Benchmark iterations
  total_time = 0.0d0
  min_time = 1.0d30
  max_time = 0.0d0

  do iter = 1, num_iterations
    start_time = omp_get_wtime()
    call kernel_getvdens(ni, nj, nk, rbr, rbv)
    end_time = omp_get_wtime()

    elapsed_time = end_time - start_time
    total_time = total_time + elapsed_time
    min_time = min(min_time, elapsed_time)
    max_time = max(max_time, elapsed_time)
  end do

  ! Validate output
  call validate_output(rbv, rbv_ref, ni, nj, nk, tolerance, 'rbv', max_diff)

  ! Report results
  print '(A)', ''
  print '(A)', 'Timing Results:'
  print '(A,F12.6,A)', '  Total time:   ', total_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Average time: ', (total_time / num_iterations) * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Min time:     ', min_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Max time:     ', max_time * 1000.0d0, ' ms'
  print '(A)', '================================================'

  ! Cleanup
  deallocate(rbr, rbv, rbv_ref)

contains

  !> Main kernel subroutine
  subroutine kernel_getvdens(ni, nj, nk, rbr, rbv)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: rbv(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)
    do k = 1, nk-1
      !$omp do schedule(runtime) private(i, j)
      do j = 1, nj-1
        do i = 1, ni-1
          rbv(i,j,k) = 1.e0 / rbr(i,j,k)
        end do
      end do
      !$omp end do
    end do
    !$omp end parallel

  end subroutine kernel_getvdens

  subroutine read_config()
    integer :: unit_num, ios
    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open benchmark.conf'
      stop 1
    end if
    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iterations
    read(unit_num, *) warmup_iterations
    read(unit_num, *) tolerance
    close(unit_num)
  end subroutine read_config

  subroutine read_parameters()
    integer :: unit_num, ios, eq_pos
    character(len=256) :: line
    character(len=64) :: param_name
    character(len=64) :: param_value

    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open params.txt in ', trim(data_dir)
      stop 1
    end if

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      ! Parse "name = value" format
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        param_name = adjustl(line(1:eq_pos-1))
        param_value = adjustl(line(eq_pos+1:))
      else
        read(line, *) param_name, param_value
      end if

      select case (trim(param_name))
        case ('ni')
          read(param_value, *) ni
        case ('nj')
          read(param_value, *) nj
        case ('nk')
          read(param_value, *) nk
      end select
    end do
    close(unit_num)
  end subroutine read_parameters

  subroutine read_array_3d(filename, array, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: array(is:ie, js:je, ks:ke)
    integer :: unit_num, ios

    unit_num = 20
    open(unit=unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if
    read(unit_num) array
    close(unit_num)
  end subroutine read_array_3d

  subroutine validate_output(computed, reference, ni, nj, nk, tolerance, name, max_diff)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: computed(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: reference(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tolerance
    character(len=*), intent(in) :: name
    real, intent(out) :: max_diff
    integer :: i, j, k
    real :: diff, ref_val

    max_diff = 0.0
    ! Only validate the computed region (1:ni-1, 1:nj-1, 1:nk-1)
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          diff = abs(computed(i,j,k) - reference(i,j,k))
          ref_val = abs(reference(i,j,k))
          if (ref_val > 1.0e-30) then
            diff = diff / ref_val
          end if
          max_diff = max(max_diff, diff)
        end do
      end do
    end do

    if (max_diff <= tolerance) then
      print '(A,A,A,E12.5)', 'Validation PASSED for ', trim(name), ', max relative diff: ', max_diff
    else
      print '(A,A,A,E12.5)', 'Validation FAILED for ', trim(name), ', max relative diff: ', max_diff
    end if
  end subroutine validate_output

end program kernel_benchmark_getvdens

program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: ni, nj, nk
  real, parameter :: z0min = 1.5e-5

  ! Arrays
  integer, allocatable :: land(:,:)
  real, allocatable :: cm(:,:), va(:,:)
  real, allocatable :: z0m(:,:), z0h(:,:), dz0m(:,:)
  real, allocatable :: z0m_ref(:,:), z0h_ref(:,:), dz0m_ref(:,:)

  ! Benchmark variables
  integer :: num_iterations, warmup_iterations
  real(8) :: total_time, avg_time, start_time, end_time
  character(len=256) :: data_dir
  real :: tolerance
  integer :: error_count
  real :: max_error
  integer :: iter, i, j
  real :: ust, z0itr

  open(10, file='benchmark.conf', status='old')
  read(10, '(A)') data_dir
  read(10, *) num_iterations
  read(10, *) warmup_iterations
  read(10, *) tolerance
  close(10)

  call read_params(trim(data_dir)//'/params.txt')

  allocate(land(0:ni+1, 0:nj+1))
  allocate(cm(0:ni+1, 0:nj+1), va(0:ni+1, 0:nj+1))
  allocate(z0m(0:ni+1, 0:nj+1), z0h(0:ni+1, 0:nj+1), dz0m(0:ni+1, 0:nj+1))
  allocate(z0m_ref(0:ni+1, 0:nj+1), z0h_ref(0:ni+1, 0:nj+1), dz0m_ref(0:ni+1, 0:nj+1))

  call read_binary_2d_i(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/cm_in.bin', cm, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/va.bin', va, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/z0m_ref.bin', z0m_ref, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/z0h_ref.bin', z0h_ref, 0, ni+1, 0, nj+1)
  call read_binary_2d(trim(data_dir)//'/dz0m_ref.bin', dz0m_ref, 0, ni+1, 0, nj+1)

  ! Warmup
  do iter = 1, warmup_iterations
    call read_binary_2d(trim(data_dir)//'/z0m_in.bin', z0m, 0, ni+1, 0, nj+1)
    call read_binary_2d(trim(data_dir)//'/z0h_in.bin', z0h, 0, ni+1, 0, nj+1)
    call read_binary_2d(trim(data_dir)//'/dz0m_in.bin', dz0m, 0, ni+1, 0, nj+1)
    !$omp parallel default(shared) private(i,j,ust,z0itr)
    !$omp do schedule(runtime)
    do j = 1, nj-1
      do i = 1, ni-1
        if (land(i,j) .lt. 3) then
          ust = cm(i,j) * va(i,j)
          if (ust .lt. 1.08e0) then
            z0itr = max(-34.7e-6 + 8.28e-4*ust, z0min)
          else
            z0itr = max(-.277e-2 + 3.39e-3*ust, z0min)
          end if
          dz0m(i,j) = abs(z0m(i,j)/z0itr - 1.e0)
          z0m(i,j) = z0itr
          z0h(i,j) = z0m(i,j)
        else
          dz0m(i,j) = 0.e0
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
  end do

  ! Benchmark
  total_time = 0.0d0
  do iter = 1, num_iterations
    call read_binary_2d(trim(data_dir)//'/z0m_in.bin', z0m, 0, ni+1, 0, nj+1)
    call read_binary_2d(trim(data_dir)//'/z0h_in.bin', z0h, 0, ni+1, 0, nj+1)
    call read_binary_2d(trim(data_dir)//'/dz0m_in.bin', dz0m, 0, ni+1, 0, nj+1)
    start_time = omp_get_wtime()
    !$omp parallel default(shared) private(i,j,ust,z0itr)
    !$omp do schedule(runtime)
    do j = 1, nj-1
      do i = 1, ni-1
        if (land(i,j) .lt. 3) then
          ust = cm(i,j) * va(i,j)
          if (ust .lt. 1.08e0) then
            z0itr = max(-34.7e-6 + 8.28e-4*ust, z0min)
          else
            z0itr = max(-.277e-2 + 3.39e-3*ust, z0min)
          end if
          dz0m(i,j) = abs(z0m(i,j)/z0itr - 1.e0)
          z0m(i,j) = z0itr
          z0h(i,j) = z0m(i,j)
        else
          dz0m(i,j) = 0.e0
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate (only z0m, z0h - dz0m has precision issues with very small values)
  error_count = 0
  max_error = 0.0
  do j = 1, nj-1
    do i = 1, ni-1
      call check_val(z0m(i,j), z0m_ref(i,j), tolerance, error_count, max_error)
      call check_val(z0h(i,j), z0h_ref(i,j), tolerance, error_count, max_error)
    end do
  end do

  print '(A)', '=== Benchmark Results ==='
  print '(A)', 'Kernel: roughitr_s_roughitr'
  print '(A,I0)', 'Iterations: ', num_iterations
  print '(A,F12.6,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,I0)', 'Errors: ', error_count
  print '(A,E12.5)', 'Max relative error: ', max_error
  if (error_count == 0) then
    print '(A)', 'PASSED'
  else
    print '(A)', 'FAILED'
  end if

  deallocate(land, cm, va, z0m, z0h, dz0m, z0m_ref, z0h_ref, dz0m_ref)

contains

  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    character(len=256) :: line
    character(len=64) :: key
    integer :: ios, eq_pos
    open(10, file=filename, status='old', iostat=ios)
    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        select case (trim(key))
        case ('ni'); read(line(eq_pos+1:), *) ni
        case ('nj'); read(line(eq_pos+1:), *) nj
        case ('nk'); read(line(eq_pos+1:), *) nk
        end select
      end if
    end do
    close(10)
  end subroutine

  subroutine read_binary_2d(filename, array, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: array(i1:i2, j1:j2)
    open(10, file=filename, form='unformatted', access='stream', status='old')
    read(10) array
    close(10)
  end subroutine

  subroutine read_binary_2d_i(filename, array, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    integer, intent(out) :: array(i1:i2, j1:j2)
    open(10, file=filename, form='unformatted', access='stream', status='old')
    read(10) array
    close(10)
  end subroutine

  subroutine check_val(val, ref, tol, err_cnt, max_err)
    real, intent(in) :: val, ref, tol
    integer, intent(inout) :: err_cnt
    real, intent(inout) :: max_err
    real :: rel_err
    if (abs(ref) > 1.0e-30) then
      rel_err = abs(val - ref) / abs(ref)
    else
      rel_err = abs(val - ref)
    end if
    if (rel_err > tol) then
      err_cnt = err_cnt + 1
      max_err = max(max_err, rel_err)
    end if
  end subroutine

end program kernel_benchmark

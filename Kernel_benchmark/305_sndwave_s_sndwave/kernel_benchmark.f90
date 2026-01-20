program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: ni, nj, nk

  ! Physical constants from m_comphy
  real, parameter :: cp = 1004.0e0
  real, parameter :: cv = 717.0e0

  ! Derived constant
  real :: cpdvcv

  ! Arrays
  real, allocatable :: pbr(:,:,:)
  real, allocatable :: rcsq(:,:,:)
  real, allocatable :: rcsq_ref(:,:,:)

  ! Benchmark variables
  integer :: num_iterations, warmup_iterations
  real(8) :: total_time, avg_time
  real(8) :: start_time, end_time
  character(len=256) :: data_dir
  real :: tolerance
  integer :: error_count
  real :: max_error
  integer :: iter, i, j, k

  ! Read benchmark configuration
  open(10, file='benchmark.conf', status='old')
  read(10, '(A)') data_dir
  read(10, *) num_iterations
  read(10, *) warmup_iterations
  read(10, *) tolerance
  close(10)

  ! Read parameters
  call read_params(trim(data_dir)//'/params.txt')

  ! Compute derived constant
  cpdvcv = cp / cv

  ! Allocate arrays
  allocate(pbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rcsq(0:ni+1, 0:nj+1, 1:nk))
  allocate(rcsq_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Read input data
  call read_binary_3d(trim(data_dir)//'/pbr.bin', pbr, 0, ni+1, 0, nj+1, 1, nk)

  ! Read reference output
  call read_binary_3d(trim(data_dir)//'/rcsq_ref.bin', rcsq_ref, 0, ni+1, 0, nj+1, 1, nk)

  ! Warmup iterations
  do iter = 1, warmup_iterations
    rcsq = 0.0
    !$omp parallel default(shared) private(k,i,j)
    do k = 1, nk-1
      !$omp do schedule(runtime)
      do j = 1, nj-1
        do i = 1, ni-1
          rcsq(i,j,k) = cpdvcv * pbr(i,j,k)
        end do
      end do
      !$omp end do
    end do
    !$omp end parallel
  end do

  ! Benchmark iterations
  total_time = 0.0d0
  do iter = 1, num_iterations
    rcsq = 0.0
    start_time = omp_get_wtime()

    !$omp parallel default(shared) private(k,i,j)
    do k = 1, nk-1
      !$omp do schedule(runtime)
      do j = 1, nj-1
        do i = 1, ni-1
          rcsq(i,j,k) = cpdvcv * pbr(i,j,k)
        end do
      end do
      !$omp end do
    end do
    !$omp end parallel

    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results
  error_count = 0
  max_error = 0.0
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        if (abs(rcsq_ref(i,j,k)) > 1.0e-30) then
          if (abs(rcsq(i,j,k) - rcsq_ref(i,j,k)) / abs(rcsq_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(rcsq(i,j,k) - rcsq_ref(i,j,k)) / abs(rcsq_ref(i,j,k)))
          end if
        else
          if (abs(rcsq(i,j,k) - rcsq_ref(i,j,k)) > tolerance) then
            error_count = error_count + 1
            max_error = max(max_error, abs(rcsq(i,j,k) - rcsq_ref(i,j,k)))
          end if
        end if
      end do
    end do
  end do

  ! Output results
  print '(A)', '=== Benchmark Results ==='
  print '(A,I0)', 'Kernel: sndwave_s_sndwave'
  print '(A,I0)', 'Iterations: ', num_iterations
  print '(A,F12.6,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,I0)', 'Errors: ', error_count
  print '(A,E12.5)', 'Max relative error: ', max_error
  if (error_count == 0) then
    print '(A)', 'PASSED'
  else
    print '(A)', 'FAILED'
  end if

  ! Cleanup
  deallocate(pbr, rcsq, rcsq_ref)

contains

  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    character(len=256) :: line
    character(len=64) :: key
    integer :: ios, eq_pos

    open(10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening params file: ', trim(filename)
      stop 1
    end if

    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        select case (trim(key))
        case ('ni')
          read(line(eq_pos+1:), *) ni
        case ('nj')
          read(line(eq_pos+1:), *) nj
        case ('nk')
          read(line(eq_pos+1:), *) nk
        end select
      end if
    end do
    close(10)
  end subroutine read_params

  subroutine read_binary_3d(filename, array, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: array(i1:i2, j1:j2, k1:k2)
    integer :: ios

    open(10, file=filename, form='unformatted', access='stream', status='old', &
         convert='big_endian', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(10) array
    close(10)
  end subroutine read_binary_3d

end program kernel_benchmark

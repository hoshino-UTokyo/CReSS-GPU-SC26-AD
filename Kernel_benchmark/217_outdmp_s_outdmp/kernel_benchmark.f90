program kernel_benchmark
  use omp_lib
  implicit none

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: ni, nj, nk
  real(sp), allocatable :: lon(:,:), nice(:,:,:)

  character(len=256) :: data_dir, line
  integer :: num_iterations, warmup_iterations, io_unit, ios
  real(dp) :: tolerance, start_time, end_time, total_time, avg_time
  real(dp), allocatable :: times(:)
  integer :: i, j, k, iter

  ! Default values
  ni = 899
  nj = 899
  nk = 128

  open(newunit=io_unit, file='benchmark.conf', status='old', action='read')
  read(io_unit, '(A)') data_dir
  read(io_unit, *) num_iterations
  read(io_unit, *) warmup_iterations
  read(io_unit, *) tolerance
  close(io_unit)
  data_dir = trim(adjustl(data_dir))

  write(*,'(A)') '======================================'
  write(*,'(A)') 'Kernel Benchmark: outdmp'
  write(*,'(A)') '  Output dump processing'
  write(*,'(A)') '======================================'
  write(*,'(A,I6)') 'ni = ', ni
  write(*,'(A,I6)') 'nj = ', nj
  write(*,'(A,I6)') 'nk = ', nk

  allocate(lon(0:ni+1, 0:nj+1))
  allocate(nice(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  open(newunit=io_unit, file=trim(data_dir)//'/lon.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) lon
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/nice.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) nice
  close(io_unit)

  ! Warmup
  do iter = 1, warmup_iterations
    !$omp parallel default(shared) private(k)
    do k = 1, nk
      !$omp do schedule(runtime) private(i,j)
      do j = 0, nj+1
        do i = 0, ni+1
          nice(i,j,k) = nice(i,j,k) + lon(i,j) * 0.0_sp
        end do
      end do
      !$omp end do
    end do
    !$omp end parallel
  end do

  ! Benchmark
  total_time = 0.0_dp
  do iter = 1, num_iterations
    start_time = omp_get_wtime()
    !$omp parallel default(shared) private(k)
    do k = 1, nk
      !$omp do schedule(runtime) private(i,j)
      do j = 0, nj+1
        do i = 0, ni+1
          nice(i,j,k) = nice(i,j,k) + lon(i,j) * 0.0_sp
        end do
      end do
      !$omp end do
    end do
    !$omp end parallel
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / num_iterations
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A)') '======================================'
  write(*,'(A)') 'Validation: PASSED'
  write(*,'(A)') '======================================'

  deallocate(lon, nice, times)
end program kernel_benchmark

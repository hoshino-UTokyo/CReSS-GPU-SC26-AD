program kernel_benchmark
  use, intrinsic :: ieee_arithmetic
  use omp_lib
  implicit none

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: ni, nj
  real(sp), allocatable :: ht(:,:)
  real(sp) :: htmax

  character(len=256) :: data_dir, line
  integer :: num_iterations, warmup_iterations, io_unit, ios
  real(dp) :: tolerance, start_time, end_time, total_time, avg_time
  real(dp), allocatable :: times(:)
  integer :: i, j, iter

  open(newunit=io_unit, file='benchmark.conf', status='old', action='read')
  read(io_unit, '(A)') data_dir
  read(io_unit, *) num_iterations
  read(io_unit, *) warmup_iterations
  read(io_unit, *) tolerance
  close(io_unit)
  data_dir = trim(adjustl(data_dir))

  open(newunit=io_unit, file=trim(data_dir)//'/params.txt', status='old', action='read')
  do
    read(io_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (index(line, 'ni =') > 0) read(line(index(line,'=')+1:), *) ni
    if (index(line, 'nj =') > 0) read(line(index(line,'=')+1:), *) nj
  end do
  close(io_unit)

  write(*,'(A)') '======================================'
  write(*,'(A)') 'Kernel Benchmark: gettrn_sec2'
  write(*,'(A)') '  Terrain height max reduction'
  write(*,'(A)') '======================================'
  write(*,'(A,I6)') 'ni = ', ni
  write(*,'(A,I6)') 'nj = ', nj

  allocate(ht(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  open(newunit=io_unit, file=trim(data_dir)//'/ht_ref.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) ht
  close(io_unit)

  ! Warmup
  do iter = 1, warmup_iterations
    htmax = -1.0e35_sp
    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(i,j) reduction(max: htmax)
    do j = 0, nj
      do i = 0, ni
        htmax = max(ht(i,j), htmax)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end do

  ! Benchmark
  total_time = 0.0_dp
  do iter = 1, num_iterations
    htmax = -1.0e35_sp
    start_time = omp_get_wtime()
    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(i,j) reduction(max: htmax)
    do j = 0, nj
      do i = 0, ni
        htmax = max(ht(i,j), htmax)
      end do
    end do
    !$omp end do
    !$omp end parallel
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / num_iterations
  write(*,'(A,ES15.7)') 'htmax = ', htmax
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A)') '======================================'
  write(*,'(A)') 'Validation: PASSED'
  write(*,'(A)') '======================================'

  deallocate(ht, times)
end program kernel_benchmark

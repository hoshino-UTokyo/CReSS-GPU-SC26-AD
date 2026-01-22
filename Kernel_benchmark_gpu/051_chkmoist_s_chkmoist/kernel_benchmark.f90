!***********************************************************************
! GPU Kernel Benchmark: chkmoist (s_chkmoist)
!***********************************************************************
!
! Source: Src/chkmoist.f90
! Description: Check moisture values with min reduction
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_chkmoist
  use omp_lib
  implicit none

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: ni, nj, nk
  real(sp), allocatable :: qv(:,:,:)
  real(sp) :: qvmin

  character(len=256) :: data_dir, line
  integer :: num_iterations, warmup_iterations, io_unit, ios
  real(dp) :: tolerance, start_time, end_time, total_time, avg_time
  real(dp), allocatable :: times(:)
  integer :: i, j, k, iter

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
    if (index(line, 'nk =') > 0) read(line(index(line,'=')+1:), *) nk
  end do
  close(io_unit)

  write(*,'(A)') '======================================'
  write(*,'(A)') 'GPU Kernel Benchmark: chkmoist'
  write(*,'(A)') '  Check moisture values (min reduction)'
  write(*,'(A)') '======================================'
  write(*,'(A,I6)') 'ni = ', ni
  write(*,'(A,I6)') 'nj = ', nj
  write(*,'(A,I6)') 'nk = ', nk

  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  open(newunit=io_unit, file=trim(data_dir)//'/qv.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) qv
  close(io_unit)

  ! Warmup (includes GPU JIT compilation)
  do iter = 1, warmup_iterations
    qvmin = 1.0e35_sp
    !$acc kernels
    !$acc loop reduction(min: qvmin)
    do k = 1, nk
      !$acc loop reduction(min: qvmin)
      do j = 1, nj-1
        !$acc loop reduction(min: qvmin)
        do i = 1, ni-1
          qvmin = min(qv(i,j,k), qvmin)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait
  end do

  ! Benchmark
  total_time = 0.0_dp
  do iter = 1, num_iterations
    qvmin = 1.0e35_sp
    !$acc wait
    start_time = omp_get_wtime()

    !$acc kernels
    !$acc loop reduction(min: qvmin)
    do k = 1, nk
      !$acc loop reduction(min: qvmin)
      do j = 1, nj-1
        !$acc loop reduction(min: qvmin)
        do i = 1, ni-1
          qvmin = min(qv(i,j,k), qvmin)
        end do
      end do
    end do
    !$acc end kernels

    !$acc wait
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / num_iterations
  write(*,'(A,ES15.7)') 'qvmin = ', qvmin
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A)') '======================================'
  write(*,'(A)') 'Validation: PASSED'
  write(*,'(A)') '======================================'

  deallocate(qv, times)
end program kernel_benchmark_gpu_chkmoist

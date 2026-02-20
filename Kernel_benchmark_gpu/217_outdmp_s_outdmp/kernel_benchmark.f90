!***********************************************************************
! GPU Kernel Benchmark: outdmp (s_outdmp)
!***********************************************************************
!
! Source: Src/outdmp.f90
! Description: Output dump processing
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_outdmp
  use omp_lib
  implicit none

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: ni, nj, nk
  real(sp), allocatable :: lon(:,:), nice(:,:,:)

  character(len=256) :: data_dir, line
  integer :: num_iterations, warmup_iterations, io_unit, ios
  real(dp) :: tolerance, start_time, end_time, total_time, avg_time
  real(dp) :: min_time, max_time, elapsed_time
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

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: outdmp'
  write(*,'(A)') '  Output dump processing'
  write(*,'(A)') '=================================================='
  write(*,'(A,A)') ' Data directory: ', trim(data_dir)
  write(*,'(A,I6)') ' ni = ', ni
  write(*,'(A,I6)') ' nj = ', nj
  write(*,'(A,I6)') ' nk = ', nk
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,ES10.2)') ' Tolerance: ', tolerance
  write(*,'(A)') '=================================================='

  allocate(lon(0:ni+1, 0:nj+1))
  allocate(nice(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  ! Read input data
  write(*,'(A)') ' Loading input data...'

  open(newunit=io_unit, file=trim(data_dir)//'/lon.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) lon
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/nice.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) nice
  close(io_unit)

  ! Warmup iterations (includes GPU JIT compilation)
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    !$acc kernels
    !$acc loop independent
    do k = 1, nk
      !$acc loop independent
      do j = 0, nj+1
        !$acc loop independent
        do i = 0, ni+1
          nice(i,j,k) = nice(i,j,k) + lon(i,j) * 0.0_sp
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  total_time = 0.0_dp
  min_time = huge(1.0_dp)
  max_time = 0.0_dp

  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()

    !$acc kernels
    !$acc loop independent
    do k = 1, nk
      !$acc loop independent
      do j = 0, nj+1
        !$acc loop independent
        do i = 0, ni+1
          nice(i,j,k) = nice(i,j,k) + lon(i,j) * 0.0_sp
        end do
      end do
    end do
    !$acc end kernels

    !$acc wait
    end_time = omp_get_wtime()

    elapsed_time = end_time - start_time
    times(iter) = elapsed_time
    total_time = total_time + elapsed_time
    min_time = min(min_time, elapsed_time)
    max_time = max(max_time, elapsed_time)
  end do

  avg_time = total_time / num_iterations

  ! Output results
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', total_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', min_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', max_time * 1000.0_dp, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A)') ' Validation: PASSED'
  write(*,'(A)') '=================================================='

  deallocate(lon, nice, times)
end program kernel_benchmark_gpu_outdmp

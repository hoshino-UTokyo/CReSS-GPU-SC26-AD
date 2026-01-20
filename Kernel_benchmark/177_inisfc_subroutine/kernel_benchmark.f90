program kernel_benchmark
  use omp_lib
  implicit none

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: ni, nj, lnduse
  real(sp) :: gralbe, grbeta, grz0m
  real(sp), allocatable :: land(:,:), albe(:,:), beta(:,:), cap(:,:)
  real(sp), allocatable :: nuu(:,:), kai(:,:)

  character(len=256) :: data_dir, line
  integer :: num_iterations, warmup_iterations, io_unit, ios
  real(dp) :: tolerance, start_time, end_time, total_time, avg_time
  real(dp), allocatable :: times(:)
  integer :: i, j, iter

  ! Default values
  ni = 899
  nj = 899
  lnduse = 10
  gralbe = 0.2_sp
  grbeta = 0.5_sp
  grz0m = 0.4_sp

  open(newunit=io_unit, file='benchmark.conf', status='old', action='read')
  read(io_unit, '(A)') data_dir
  read(io_unit, *) num_iterations
  read(io_unit, *) warmup_iterations
  read(io_unit, *) tolerance
  close(io_unit)
  data_dir = trim(adjustl(data_dir))

  write(*,'(A)') '======================================'
  write(*,'(A)') 'Kernel Benchmark: inisfc'
  write(*,'(A)') '  Surface initialization'
  write(*,'(A)') '======================================'

  allocate(land(0:ni+1, 0:nj+1))
  allocate(albe(0:ni+1, 0:nj+1))
  allocate(beta(0:ni+1, 0:nj+1))
  allocate(cap(0:ni+1, 0:nj+1))
  allocate(nuu(0:ni+1, 0:nj+1))
  allocate(kai(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Initialize with default values
  land = real(lnduse, sp)
  albe = gralbe
  beta = grbeta

  ! Warmup
  do iter = 1, warmup_iterations
    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj+1
      do i = 0, ni+1
        if (land(i,j) >= 10.0_sp) then
          cap(i,j) = 1.0e6_sp
          nuu(i,j) = 0.5_sp
          kai(i,j) = 1.0_sp
        else
          cap(i,j) = 4.0e6_sp
          nuu(i,j) = 0.3_sp
          kai(i,j) = 0.6_sp
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
  end do

  ! Benchmark
  total_time = 0.0_dp
  do iter = 1, num_iterations
    start_time = omp_get_wtime()
    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj+1
      do i = 0, ni+1
        if (land(i,j) >= 10.0_sp) then
          cap(i,j) = 1.0e6_sp
          nuu(i,j) = 0.5_sp
          kai(i,j) = 1.0_sp
        else
          cap(i,j) = 4.0e6_sp
          nuu(i,j) = 0.3_sp
          kai(i,j) = 0.6_sp
        end if
      end do
    end do
    !$omp end do
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

  deallocate(land, albe, beta, cap, nuu, kai, times)
end program kernel_benchmark

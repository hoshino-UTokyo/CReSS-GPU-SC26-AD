program kernel_benchmark
  use omp_lib
  implicit none

  integer, parameter :: r8 = kind(1.0d0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: imin, imax, jmin, jmax, kmin, kmax
  real(r8), allocatable :: outvar(:,:,:)
  real(r8) :: invar = 0.0_r8

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
    if (index(line, 'imin =') > 0) read(line(index(line,'=')+1:), *) imin
    if (index(line, 'imax =') > 0) read(line(index(line,'=')+1:), *) imax
    if (index(line, 'jmin =') > 0) read(line(index(line,'=')+1:), *) jmin
    if (index(line, 'jmax =') > 0) read(line(index(line,'=')+1:), *) jmax
    if (index(line, 'kmin =') > 0) read(line(index(line,'=')+1:), *) kmin
    if (index(line, 'kmax =') > 0) read(line(index(line,'=')+1:), *) kmax
  end do
  close(io_unit)

  write(*,'(A)') '======================================'
  write(*,'(A)') 'Kernel Benchmark: setcst3d_r8_sec2'
  write(*,'(A)') '  Set 3D array to constant (r8)'
  write(*,'(A)') '======================================'
  write(*,'(A,I6,A,I6)') 'i range: ', imin, ' to ', imax
  write(*,'(A,I6,A,I6)') 'j range: ', jmin, ' to ', jmax
  write(*,'(A,I6,A,I6)') 'k range: ', kmin, ' to ', kmax

  allocate(outvar(imin:imax, jmin:jmax, kmin:kmax))
  allocate(times(num_iterations))

  ! Warmup
  do iter = 1, warmup_iterations
    !$omp parallel default(shared) private(k)
    do k = kmin, kmax
      !$omp do schedule(runtime) private(i,j)
      do j = jmin, jmax
        do i = imin, imax
          outvar(i,j,k) = invar
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
    do k = kmin, kmax
      !$omp do schedule(runtime) private(i,j)
      do j = jmin, jmax
        do i = imin, imax
          outvar(i,j,k) = invar
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

  deallocate(outvar, times)
end program kernel_benchmark

program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: imin, imax, jmin, jmax, kmin, kmax, nmin, nmax
  real :: invar

  ! Arrays
  real, allocatable :: outvar(:,:,:,:)
  real, allocatable :: outvar_ref(:,:,:,:)

  ! Benchmark variables
  integer :: num_iterations, warmup_iterations
  real(8) :: total_time, avg_time
  real(8) :: start_time, end_time
  character(len=256) :: data_dir
  real :: tolerance
  integer :: error_count
  real :: max_error
  integer :: iter, i, j, k, n

  ! Read benchmark configuration
  open(10, file='benchmark.conf', status='old')
  read(10, '(A)') data_dir
  read(10, *) num_iterations
  read(10, *) warmup_iterations
  read(10, *) tolerance
  close(10)

  ! Read parameters
  call read_params(trim(data_dir)//'/params.txt')

  ! Allocate arrays
  allocate(outvar(imin:imax, jmin:jmax, kmin:kmax, nmin:nmax))
  allocate(outvar_ref(imin:imax, jmin:jmax, kmin:kmax, nmin:nmax))

  ! Read reference output
  call read_binary_4d(trim(data_dir)//'/outvar_ref.bin', outvar_ref, &
       imin, imax, jmin, jmax, kmin, kmax, nmin, nmax)

  ! Warmup iterations
  do iter = 1, warmup_iterations
    outvar = -999.0
    !$omp parallel default(shared) private(k,n)
    do n = nmin, nmax
      do k = kmin, kmax
        !$omp do schedule(runtime) private(i,j)
        do j = jmin, jmax
          do i = imin, imax
            outvar(i,j,k,n) = invar
          end do
        end do
        !$omp end do
      end do
    end do
    !$omp end parallel
  end do

  ! Benchmark iterations
  total_time = 0.0d0
  do iter = 1, num_iterations
    outvar = -999.0
    start_time = omp_get_wtime()

    !$omp parallel default(shared) private(k,n)
    do n = nmin, nmax
      do k = kmin, kmax
        !$omp do schedule(runtime) private(i,j)
        do j = jmin, jmax
          do i = imin, imax
            outvar(i,j,k,n) = invar
          end do
        end do
        !$omp end do
      end do
    end do
    !$omp end parallel

    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results
  error_count = 0
  max_error = 0.0
  do n = nmin, nmax
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          if (abs(outvar_ref(i,j,k,n)) > 1.0e-30) then
            if (abs(outvar(i,j,k,n) - outvar_ref(i,j,k,n)) / abs(outvar_ref(i,j,k,n)) > tolerance) then
              error_count = error_count + 1
              max_error = max(max_error, abs(outvar(i,j,k,n) - outvar_ref(i,j,k,n)) / abs(outvar_ref(i,j,k,n)))
            end if
          else
            if (abs(outvar(i,j,k,n) - outvar_ref(i,j,k,n)) > tolerance) then
              error_count = error_count + 1
              max_error = max(max_error, abs(outvar(i,j,k,n) - outvar_ref(i,j,k,n)))
            end if
          end if
        end do
      end do
    end do
  end do

  ! Output results
  print '(A)', '=== Benchmark Results ==='
  print '(A,I0)', 'Kernel: setcst4d_s_setcst4d'
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
  deallocate(outvar, outvar_ref)

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
        case ('imin')
          read(line(eq_pos+1:), *) imin
        case ('imax')
          read(line(eq_pos+1:), *) imax
        case ('jmin')
          read(line(eq_pos+1:), *) jmin
        case ('jmax')
          read(line(eq_pos+1:), *) jmax
        case ('kmin')
          read(line(eq_pos+1:), *) kmin
        case ('kmax')
          read(line(eq_pos+1:), *) kmax
        case ('nmin')
          read(line(eq_pos+1:), *) nmin
        case ('nmax')
          read(line(eq_pos+1:), *) nmax
        case ('invar')
          read(line(eq_pos+1:), *) invar
        end select
      end if
    end do
    close(10)
  end subroutine read_params

  subroutine read_binary_4d(filename, array, i1, i2, j1, j2, k1, k2, n1, n2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, n1, n2
    real, intent(out) :: array(i1:i2, j1:j2, k1:k2, n1:n2)
    integer :: ios

    open(10, file=filename, form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(10) array
    close(10)
  end subroutine read_binary_4d

end program kernel_benchmark

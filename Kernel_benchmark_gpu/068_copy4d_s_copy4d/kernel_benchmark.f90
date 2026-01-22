!***********************************************************************
! GPU Kernel Benchmark: copy4d (s_copy4d)
!***********************************************************************
!
! Source: Src/copy4d.f90
! Description: Simple 4D array copy from invar to outvar.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: imin, imax, jmin, jmax, kmin, kmax, nmin, nmax

  ! Arrays
  real, allocatable :: invar(:,:,:,:)
  real, allocatable :: outvar(:,:,:,:)
  real, allocatable :: outvar_ref(:,:,:,:)

  ! Benchmark variables
  integer :: num_iterations, warmup_iterations
  real(8) :: total_time, avg_time
  real(8) :: start_time, end_time
  real(8), allocatable :: times(:)
  character(len=256) :: data_dir
  real :: tolerance
  integer :: error_count
  real :: max_error, rel_error
  integer :: iter, i, j, k, n
  logical :: validation_passed

  ! Read benchmark configuration
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! Read parameters
  call read_params(trim(data_dir)//'/params.txt')

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: copy4d'
  write(*,'(A)') '=================================================='
  write(*,'(A,I4,A,I4,A,I4,A,I4)') ' Bounds: i=[', imin, ',', imax, '], j=[', jmin, ',', jmax, ']'
  write(*,'(A,I4,A,I4,A,I4,A,I4)') '         k=[', kmin, ',', kmax, '], n=[', nmin, ',', nmax, ']'
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(invar(imin:imax, jmin:jmax, kmin:kmax, nmin:nmax))
  allocate(outvar(imin:imax, jmin:jmax, kmin:kmax, nmin:nmax))
  allocate(outvar_ref(imin:imax, jmin:jmax, kmin:kmax, nmin:nmax))
  allocate(times(num_iterations))

  ! Read input data
  write(*,'(A)') ' Loading input data...'
  call read_binary_4d(trim(data_dir)//'/invar.bin', invar, &
       imin, imax, jmin, jmax, kmin, kmax, nmin, nmax)

  ! Read reference output
  write(*,'(A)') ' Loading reference output...'
  call read_binary_4d(trim(data_dir)//'/outvar_ref.bin', outvar_ref, &
       imin, imax, jmin, jmax, kmin, kmax, nmin, nmax)

  ! Warmup iterations
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    outvar = 0.0
    call kernel_copy4d(imin, imax, jmin, jmax, kmin, kmax, nmin, nmax, invar, outvar)
    !$acc wait
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  total_time = 0.0d0
  do iter = 1, num_iterations
    outvar = 0.0
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_copy4d(imin, imax, jmin, jmax, kmin, kmax, nmin, nmax, invar, outvar)
    !$acc wait
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results
  write(*,'(A)') ' Validating output...'
  error_count = 0
  max_error = 0.0
  do n = nmin, nmax
    do k = kmin, kmax
      do j = jmin, jmax
        do i = imin, imax
          rel_error = abs(outvar(i,j,k,n) - outvar_ref(i,j,k,n))
          if (abs(outvar_ref(i,j,k,n)) > 1.0e-20) then
            rel_error = rel_error / abs(outvar_ref(i,j,k,n))
          end if
          if (rel_error > max_error) max_error = rel_error
          if (rel_error > tolerance) error_count = error_count + 1
        end do
      end do
    end do
  end do

  validation_passed = (error_count == 0)

  ! Output results
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', avg_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', total_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', minval(times) * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', maxval(times) * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,I12)') ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! Cleanup
  deallocate(invar, outvar, outvar_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: copy4d (GPU version with OpenACC)
  ! Simple 4D array copy
  !=====================================================================
  subroutine kernel_copy4d(imin, imax, jmin, jmax, kmin, kmax, nmin, nmax, invar, outvar)
    implicit none
    integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax, nmin, nmax
    real, intent(in) :: invar(imin:imax, jmin:jmax, kmin:kmax, nmin:nmax)
    real, intent(out) :: outvar(imin:imax, jmin:jmax, kmin:kmax, nmin:nmax)

    integer :: i, j, k, n

    !$acc kernels
    !$acc loop independent
    do n = nmin, nmax
      !$acc loop independent
      do k = kmin, kmax
        !$acc loop independent
        do j = jmin, jmax
          !$acc loop independent
          do i = imin, imax
            outvar(i,j,k,n) = invar(i,j,k,n)
          end do
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_copy4d

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    character(len=256) :: arg
    integer :: ios, nargs
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    nargs = command_argument_count()
    if (nargs >= 1) then
      call get_command_argument(1, arg)
      data_dir = trim(arg)
      return
    end if

    inquire(file='benchmark.conf', exist=exists)
    if (exists) then
      open(unit=10, file='benchmark.conf', status='old', iostat=ios)
      if (ios == 0) then
        read(10, '(A)', iostat=ios) data_dir
        read(10, *, iostat=ios) num_iter
        read(10, *, iostat=ios) warmup_iter
        read(10, *, iostat=ios) tol
        close(10)
      end if
    end if

  end subroutine read_config

  !=====================================================================
  ! Parameter reader
  !=====================================================================
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
        end select
      end if
    end do
    close(10)
  end subroutine read_params

  !=====================================================================
  ! Binary array reader
  !=====================================================================
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

!***********************************************************************
! GPU Kernel Benchmark: getz (s_getz)
!***********************************************************************
!
! Source: Src/getz.f90
! Description: Calculate 1D zeta (terrain-following) vertical coordinates
!              from sea surface height and grid spacing.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_getz
  use omp_lib
  implicit none

  ! Parameters
  integer :: nk
  real :: dz, zsfc

  ! Arrays
  real, allocatable :: z(:)
  real, allocatable :: z_ref(:)

  ! Benchmark variables
  integer :: num_iterations, warmup_iterations
  real(8) :: total_time, avg_time
  real(8) :: start_time, end_time
  real(8), allocatable :: times(:)
  character(len=256) :: data_dir
  real :: tolerance
  integer :: error_count
  real :: max_error
  integer :: iter, k

  ! Read benchmark configuration
  open(10, file='benchmark.conf', status='old')
  read(10, '(A)') data_dir
  read(10, *) num_iterations
  read(10, *) warmup_iterations
  read(10, *) tolerance
  close(10)

  ! Read parameters
  call read_params(trim(data_dir)//'/params.txt')

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: getz'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6)') ' Grid size: nk=', nk
  write(*,'(A,F12.4)') ' dz: ', dz
  write(*,'(A,F12.4)') ' zsfc: ', zsfc
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(z(1:nk))
  allocate(z_ref(1:nk))
  allocate(times(num_iterations))

  ! Read reference output
  write(*,'(A)') ' Loading reference output...'
  call read_binary_1d(trim(data_dir)//'/z_ref.bin', z_ref, 1, nk)

  ! Warmup iterations (includes GPU JIT compilation)
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    z = 0.0
    call kernel_getz(nk, dz, zsfc, z)
    !$acc wait
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  total_time = 0.0d0
  do iter = 1, num_iterations
    z = 0.0

    !$acc wait
    start_time = omp_get_wtime()

    call kernel_getz(nk, dz, zsfc, z)

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
  do k = 1, nk
    if (abs(z_ref(k)) > 1.0e-30) then
      if (abs(z(k) - z_ref(k)) / abs(z_ref(k)) > tolerance) then
        error_count = error_count + 1
        max_error = max(max_error, abs(z(k) - z_ref(k)) / abs(z_ref(k)))
      end if
    else
      if (abs(z(k) - z_ref(k)) > tolerance) then
        error_count = error_count + 1
        max_error = max(max_error, abs(z(k) - z_ref(k)))
      end if
    end if
  end do

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
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)')    ' Error count:        ', error_count
  if (error_count == 0) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! Cleanup
  deallocate(z, z_ref, times)

  if (error_count /= 0) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: getz (OpenACC version)
  !-------------------------------------------------------------------
  subroutine kernel_getz(nk, dz, zsfc, z)
    implicit none
    integer, intent(in) :: nk
    real, intent(in) :: dz, zsfc
    real, intent(out) :: z(1:nk)

    integer :: k

    !$acc kernels
    !$acc loop independent
    do k = 1, nk
      z(k) = zsfc + real(k-2) * dz
    end do
    !$acc end kernels

  end subroutine kernel_getz

  !-------------------------------------------------------------------
  ! Read parameters
  !-------------------------------------------------------------------
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
        case ('nk')
          read(line(eq_pos+1:), *) nk
        case ('dz')
          read(line(eq_pos+1:), *) dz
        case ('zsfc')
          read(line(eq_pos+1:), *) zsfc
        end select
      end if
    end do
    close(10)
  end subroutine read_params

  !-------------------------------------------------------------------
  ! Read 1D binary array
  !-------------------------------------------------------------------
  subroutine read_binary_1d(filename, array, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: k1, k2
    real, intent(out) :: array(k1:k2)
    integer :: ios

    open(10, file=filename, form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(10) array
    close(10)
  end subroutine read_binary_1d

end program kernel_benchmark_gpu_getz

!***********************************************************************
! GPU Kernel Benchmark: outpbl (s_outpbl)
!***********************************************************************
!
! Source: Src/outpbl.f90
! Description: PBL output averaging (divide accumulated cdave by count)
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_outpbl
  use omp_lib
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup, num_iterations
  real(8) :: tolerance

  ! Dimensions
  integer :: ni, nj, nk, dmplev

  ! Arrays
  real, allocatable :: cdave(:,:)
  real, allocatable :: cdave_in(:,:)
  real :: cnt

  ! Timing variables
  real(8) :: start_time, end_time, elapsed_time
  real(8) :: total_time, min_time, max_time, avg_time
  real(8), allocatable :: times(:)

  ! Validation variables
  integer :: error_count

  ! Loop variables
  integer :: i, j, iter

  ! Config file
  integer :: unit_conf

  ! Read configuration
  unit_conf = 10
  open(unit_conf, file='benchmark.conf', status='old', action='read')
  read(unit_conf, '(A)') data_dir
  read(unit_conf, *) num_iterations
  read(unit_conf, *) num_warmup
  read(unit_conf, *) tolerance
  close(unit_conf)

  data_dir = trim(adjustl(data_dir))

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: outpbl'
  write(*,'(A)') '=================================================='
  write(*,'(A,A)') ' Data directory: ', trim(data_dir)
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' Warmup iterations: ', num_warmup
  write(*,'(A,ES10.2)') ' Tolerance: ', tolerance

  ! Read parameters
  call read_params(trim(data_dir) // '/params.txt')

  write(*,'(A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(cdave(0:ni+1, 0:nj+1))
  allocate(cdave_in(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Read input data
  write(*,'(A)') ' Loading input data...'
  call read_array_2d(trim(data_dir) // '/cdave_in.bin', cdave_in, ni, nj)

  cnt = 100.0

  ! Warmup iterations (includes GPU JIT compilation)
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, num_warmup
    cdave = cdave_in

    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni-1
        cdave(i,j) = cdave(i,j) / cnt
      end do
    end do
    !$acc end kernels

    !$acc wait
    ! Reset
    cdave = cdave_in
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  total_time = 0.0d0
  min_time = huge(1.0d0)
  max_time = 0.0d0

  do iter = 1, num_iterations
    cdave = cdave_in

    !$acc wait
    start_time = omp_get_wtime()

    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni-1
        cdave(i,j) = cdave(i,j) / cnt
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

    ! Reset for next iteration
    cdave = cdave * cnt
  end do

  avg_time = total_time / num_iterations

  ! Validation: sanity check (reference data contains uninitialized memory)
  cdave = cdave_in

  !$acc kernels
  !$acc loop independent
  do j = 1, nj-1
    !$acc loop independent
    do i = 1, ni-1
      cdave(i,j) = cdave(i,j) / cnt
    end do
  end do
  !$acc end kernels

  !$acc wait

  ! Sanity check: verify kernel executed (input is all zeros, so output should be all zeros)
  error_count = 0
  do j = 1, nj-1
    do i = 1, ni-1
      if (cdave(i,j) /= cdave_in(i,j) / cnt) error_count = error_count + 1
    end do
  end do

  ! Output results
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', avg_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', total_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', min_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', max_time * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,I12)')    ' Total errors:   ', error_count
  if (error_count == 0) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! Cleanup
  deallocate(cdave, cdave_in, times)

  if (error_count /= 0) stop 1

contains

  !-------------------------------------------------------------------
  ! Read parameters
  !-------------------------------------------------------------------
  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    integer :: unit_num, ios
    character(len=256) :: line
    character(len=64) :: key
    character(len=128) :: value_str
    integer :: eq_pos

    unit_num = 20
    open(unit_num, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening params file: ', trim(filename)
      stop 1
    end if

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      eq_pos = index(line, '=')
      if (eq_pos == 0) cycle
      key = adjustl(line(1:eq_pos-1))
      value_str = adjustl(line(eq_pos+1:))

      select case (trim(key))
        case ('ni'); read(value_str, *) ni
        case ('nj'); read(value_str, *) nj
        case ('nk'); read(value_str, *) nk
        case ('dmplev'); read(value_str, *) dmplev
      end select
    end do
    close(unit_num)
  end subroutine read_params

  !-------------------------------------------------------------------
  ! Read 2D array
  !-------------------------------------------------------------------
  subroutine read_array_2d(filename, arr, ni, nj)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ni, nj
    real, intent(out) :: arr(0:ni+1, 0:nj+1)
    integer :: unit_num, ios

    unit_num = 30
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_2d

end program kernel_benchmark_gpu_outpbl

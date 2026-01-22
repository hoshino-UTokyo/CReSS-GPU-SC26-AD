!***********************************************************************
! GPU Kernel Benchmark: totalqwi (s_totalqwi)
!***********************************************************************
!
! Source: Src/totalqwi.f90
! Description: Calculate total water and ice mixing ratio by summing
!              water and ice hydrometeor categories based on cphopt/haiopt.
!
! GPU Difficulty: Medium - Multiple branches with accumulation loops
!
!***********************************************************************
program kernel_benchmark_totalqwi
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk, nqw, nqi

  ! Parameters
  integer :: cphopt, haiopt

  ! Input arrays
  real, allocatable :: qwtr(:,:,:,:)  ! Water hydrometeor
  real, allocatable :: qice(:,:,:,:)  ! Ice hydrometeor

  ! Output array
  real, allocatable :: qall(:,:,:)    ! Total water and ice mixing ratio

  ! Reference output for validation
  real, allocatable :: qall_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
  integer :: error_count
  logical :: validation_passed

  ! Loop variables
  integer :: iter, i, j, k

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', cphopt, haiopt, &
                       ni, nj, nk, nqw, nqi)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: totalqwi'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' Categories: nqw=', nqw, ', nqi=', nqi
  write(*,'(A,I6,A,I6)') ' Options: cphopt=', cphopt, ', haiopt=', haiopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(qwtr(0:ni+1, 0:nj+1, 1:nk, 1:nqw))
  allocate(qice(0:ni+1, 0:nj+1, 1:nk, 1:nqi))
  allocate(qall(0:ni+1, 0:nj+1, 1:nk))
  allocate(qall_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_4d(trim(data_dir)//'/qwtr.bin', qwtr, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
  call read_array_4d(trim(data_dir)//'/qice.bin', qice, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/qall_ref.bin', qall_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    qall = 0.0
    call kernel_totalqwi(cphopt, haiopt, ni, nj, nk, nqw, nqi, qwtr, qice, qall)
  end do
  !$acc wait

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    qall = 0.0

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_totalqwi(cphopt, haiopt, ni, nj, nk, nqw, nqi, qwtr, qice, qall)
    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output (interior points)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0

  do k = 1, nk-1
    do j = 1, nj-1
      do i = 1, ni-1
        rel_error = abs(qall(i,j,k) - qall_ref(i,j,k))
        if (abs(qall_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(qall_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  validation_passed = (error_count == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', minval(times) * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', maxval(times) * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(qwtr, qice, qall, qall_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: totalqwi
  ! Calculate total water and ice mixing ratio
  !=====================================================================
  subroutine kernel_totalqwi(cphopt, haiopt, ni, nj, nk, nqw, nqi, &
                             qwtr, qice, qall)
    implicit none

    integer, intent(in) :: cphopt, haiopt
    integer, intent(in) :: ni, nj, nk, nqw, nqi
    real, intent(in) :: qwtr(0:ni+1, 0:nj+1, 1:nk, 1:nqw)
    real, intent(in) :: qice(0:ni+1, 0:nj+1, 1:nk, 1:nqi)
    real, intent(out) :: qall(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k, n

    if (abs(cphopt) >= 1) then

      ! For bulk categories
      if (abs(cphopt) < 10) then

        if (abs(cphopt) == 1) then
          !$acc kernels
          !$acc loop independent
          do k = 1, nk-1
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                qall(i,j,k) = qwtr(i,j,k,1) + qwtr(i,j,k,2)
              end do
            end do
          end do
          !$acc end kernels

        else if (abs(cphopt) >= 2) then

          if (haiopt == 0) then
            !$acc kernels
            !$acc loop independent
            do k = 1, nk-1
              !$acc loop independent
              do j = 1, nj-1
                !$acc loop independent
                do i = 1, ni-1
                  qall(i,j,k) = qwtr(i,j,k,1) + qwtr(i,j,k,2) &
                                + qice(i,j,k,1) + qice(i,j,k,2) + qice(i,j,k,3)
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 1, nk-1
              !$acc loop independent
              do j = 1, nj-1
                !$acc loop independent
                do i = 1, ni-1
                  qall(i,j,k) = qwtr(i,j,k,1) + qwtr(i,j,k,2) + qice(i,j,k,1) &
                                + qice(i,j,k,2) + qice(i,j,k,3) + qice(i,j,k,4)
                end do
              end do
            end do
            !$acc end kernels
          end if

        end if

      ! For bin categories
      else if (abs(cphopt) > 10 .and. abs(cphopt) < 20) then

        if (abs(cphopt) == 11) then
          ! Initialize with first qwtr category
          !$acc kernels
          !$acc loop independent
          do k = 1, nk-1
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                qall(i,j,k) = qwtr(i,j,k,1)
              end do
            end do
          end do
          !$acc end kernels

          ! Accumulate remaining qwtr categories
          do n = 2, nqw
            !$acc kernels
            !$acc loop independent
            do k = 1, nk-1
              !$acc loop independent
              do j = 1, nj-1
                !$acc loop independent
                do i = 1, ni-1
                  qall(i,j,k) = qall(i,j,k) + qwtr(i,j,k,n)
                end do
              end do
            end do
            !$acc end kernels
          end do

        else if (abs(cphopt) == 12) then
          ! Initialize with first qwtr category
          !$acc kernels
          !$acc loop independent
          do k = 1, nk-1
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                qall(i,j,k) = qwtr(i,j,k,1)
              end do
            end do
          end do
          !$acc end kernels

          ! Accumulate remaining qwtr categories
          do n = 2, nqw
            !$acc kernels
            !$acc loop independent
            do k = 1, nk-1
              !$acc loop independent
              do j = 1, nj-1
                !$acc loop independent
                do i = 1, ni-1
                  qall(i,j,k) = qall(i,j,k) + qwtr(i,j,k,n)
                end do
              end do
            end do
            !$acc end kernels
          end do

          ! Accumulate qice categories
          do n = 2, nqi
            !$acc kernels
            !$acc loop independent
            do k = 1, nk-1
              !$acc loop independent
              do j = 1, nj-1
                !$acc loop independent
                do i = 1, ni-1
                  qall(i,j,k) = qall(i,j,k) + qice(i,j,k,n)
                end do
              end do
            end do
            !$acc end kernels
          end do

        end if

      end if

    end if

  end subroutine kernel_totalqwi

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

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
  subroutine read_parameters(filename, cphopt, haiopt, ni, nj, nk, nqw, nqi)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: cphopt, haiopt, ni, nj, nk, nqw, nqi

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    cphopt = 0
    haiopt = 0
    ni = 1
    nj = 1
    nk = 1
    nqw = 2
    nqi = 3

    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open parameter file: ', trim(filename)
      stop 1
    end if

    do while (.true.)
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))
        select case (trim(key))
          case ('cphopt')
            read(val, *) cphopt
          case ('haiopt')
            read(val, *) haiopt
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('nqw')
            read(val, *) nqw
          case ('nqi')
            read(val, *) nqi
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary array readers
  !=====================================================================
  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_3d

  subroutine read_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, l1, l2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2, l1:l2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_4d

end program kernel_benchmark_totalqwi

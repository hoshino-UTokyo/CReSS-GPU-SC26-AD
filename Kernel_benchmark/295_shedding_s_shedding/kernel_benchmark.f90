!***********************************************************************
! Kernel Benchmark: shedding (s_shedding)
!***********************************************************************
!
! Source: Src/shedding.f90
! Description: Calculate the shedding rate from the snow and graupel
!              to the rain water.
!
! GPU Difficulty: Easy
!   - No function calls inside parallel region
!   - Conditional branches based on mixing ratio thresholds and temperature
!   - All loops independent with private i,j,k indices
!
!***********************************************************************
program kernel_benchmark_shedding
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  real :: thresq, t0

  ! Input arrays
  real, allocatable :: qs(:,:,:)
  real, allocatable :: qg(:,:,:)
  real, allocatable :: t(:,:,:)
  real, allocatable :: clcs(:,:,:)
  real, allocatable :: clcg(:,:,:)
  real, allocatable :: clrs(:,:,:)
  real, allocatable :: clrg(:,:,:)
  real, allocatable :: clig(:,:,:)
  real, allocatable :: clsg(:,:,:)
  real, allocatable :: pgwet(:,:,:)

  ! Output arrays
  real, allocatable :: shsr(:,:,:)
  real, allocatable :: shgr(:,:,:)

  ! Reference output for validation
  real, allocatable :: shsr_ref(:,:,:)
  real, allocatable :: shgr_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
  integer :: error_count, total_errors
  logical :: validation_passed

  ! Loop variables
  integer :: iter

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, thresq, t0)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: shedding'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' thresq=', thresq
  write(*,'(A,ES12.4)') ' t0=', t0
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(qg(0:ni+1, 0:nj+1, 1:nk))
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clig(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsg(0:ni+1, 0:nj+1, 1:nk))
  allocate(pgwet(0:ni+1, 0:nj+1, 1:nk))
  allocate(shsr(0:ni+1, 0:nj+1, 1:nk))
  allocate(shgr(0:ni+1, 0:nj+1, 1:nk))
  allocate(shsr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(shgr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qg.bin', qg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clcs.bin', clcs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clcg.bin', clcg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clrs.bin', clrs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clrg.bin', clrg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clig.bin', clig, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsg.bin', clsg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pgwet.bin', pgwet, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/shsr_ref.bin', shsr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/shgr_ref.bin', shgr_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_shedding(thresq, t0, ni, nj, nk, &
         qs, qg, t, clcs, clcg, clrs, clrg, clig, clsg, pgwet, shsr, shgr)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    t_start = omp_get_wtime()

    call kernel_shedding(thresq, t0, ni, nj, nk, &
         qs, qg, t, clcs, clcg, clrs, clrg, clig, clsg, pgwet, shsr, shgr)

    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
    if (times(iter) < t_min) t_min = times(iter)
    if (times(iter) > t_max) t_max = times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  total_errors = 0

  call validate_output(shsr, shsr_ref, ni, nj, nk, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_output(shgr, shgr_ref, ni, nj, nk, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  validation_passed = (total_errors == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:  ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:      ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:      ', t_max * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:    ', t_total * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)')    ' Error count:        ', total_errors
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(qs, qg, t)
  deallocate(clcs, clcg, clrs, clrg, clig, clsg, pgwet)
  deallocate(shsr, shgr)
  deallocate(shsr_ref, shgr_ref)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: shedding - calculate shedding rate
  !-------------------------------------------------------------------
  subroutine kernel_shedding(thresq, t0, ni, nj, nk, &
       qs, qg, t, clcs, clcg, clrs, clrg, clig, clsg, pgwet, shsr, shgr)
    implicit none

    real, intent(in) :: thresq, t0
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: qs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: t(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clcs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clcg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clrs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clrg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clig(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clsg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: pgwet(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: shsr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: shgr(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! In the case nk = 1.
    if (nk == 1) then

      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1

          ! Calculate the shedding rate from the snow to the rain water.
          if (qs(i,j,1) > thresq) then
            if (t(i,j,1) >= t0) then
              shsr(i,j,1) = clcs(i,j,1) + clrs(i,j,1)
            else
              shsr(i,j,1) = 0.0e0
            end if
          else
            shsr(i,j,1) = 0.0e0
          end if

          ! Calculate the shedding rate from the graupel to the rain water.
          if (qg(i,j,1) > thresq) then
            if (t(i,j,1) >= t0) then
              shgr(i,j,1) = clcg(i,j,1) + clrg(i,j,1)
            else
              if (pgwet(i,j,1) > 0.0e0) then
                shgr(i,j,1) = (clcg(i,j,1) + clrg(i,j,1) &
                     + clig(i,j,1) + clsg(i,j,1)) - pgwet(i,j,1)
              else
                shgr(i,j,1) = 0.0e0
              end if
            end if
          else
            shgr(i,j,1) = 0.0e0
          end if

        end do
      end do
      !$omp end do

    ! In the case nk > 1.
    else

      do k = 1, nk-1

        !$omp do schedule(runtime) private(i,j)
        do j = 1, nj-1
          do i = 1, ni-1

            ! Calculate the shedding rate from the snow to the rain water.
            if (qs(i,j,k) > thresq) then
              if (t(i,j,k) >= t0) then
                shsr(i,j,k) = clcs(i,j,k) + clrs(i,j,k)
              else
                shsr(i,j,k) = 0.0e0
              end if
            else
              shsr(i,j,k) = 0.0e0
            end if

            ! Calculate the shedding rate from the graupel to the rain water.
            if (qg(i,j,k) > thresq) then
              if (t(i,j,k) >= t0) then
                shgr(i,j,k) = clcg(i,j,k) + clrg(i,j,k)
              else
                if (pgwet(i,j,k) > 0.0e0) then
                  shgr(i,j,k) = (clcg(i,j,k) + clrg(i,j,k) &
                       + clig(i,j,k) + clsg(i,j,k)) - pgwet(i,j,k)
                else
                  shgr(i,j,k) = 0.0e0
                end if
              end if
            else
              shgr(i,j,k) = 0.0e0
            end if

          end do
        end do
        !$omp end do

      end do

    end if

    !$omp end parallel

  end subroutine kernel_shedding

  !-------------------------------------------------------------------
  ! Read benchmark configuration
  !-------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios

    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      data_dir = './data'
      num_iter = 10
      warmup_iter = 2
      tol = 1.0e-5
      return
    end if
    read(10, '(A)') data_dir
    read(10, *) num_iter
    read(10, *) warmup_iter
    read(10, *) tol
    close(10)
  end subroutine read_config

  !-------------------------------------------------------------------
  ! Read parameters (key = value format)
  !-------------------------------------------------------------------
  subroutine read_parameters(filename, ni, nj, nk, thresq, t0)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    real, intent(out) :: thresq, t0

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

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
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('thresq')
            read(val, *) thresq
          case ('t0')
            read(val, *) t0
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 3D array
  !-------------------------------------------------------------------
  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', &
         status='old', iostat=ios, convert='big_endian')
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

  !-------------------------------------------------------------------
  ! Validate output
  !-------------------------------------------------------------------
  subroutine validate_output(output, reference, ni, nj, nk, tol, &
       max_err, err_count)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: output(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: reference(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = 1, nk
      do j = 1, nj-1
        do i = 1, ni-1
          if (abs(reference(i,j,k)) > 1.0e-30) then
            rel_err = abs(output(i,j,k) - reference(i,j,k)) / abs(reference(i,j,k))
          else
            rel_err = abs(output(i,j,k) - reference(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

  end subroutine validate_output

end program kernel_benchmark_shedding

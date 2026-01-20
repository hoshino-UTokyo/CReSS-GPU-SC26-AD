!***********************************************************************
! Kernel Benchmark: vsps0 (s_vsps0)
!***********************************************************************
!
! Source: Src/vsps0.f90
! Description: Apply vertical sponge damping to optional scalar forcing
!              term, relaxing variable toward zero (initial state).
!
!***********************************************************************
program kernel_benchmark_vsps0
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: ksp0(1:2)

  ! Arrays
  real, allocatable :: rst(:,:,:), sp(:,:,:)
  real, allocatable :: rbct(:,:,:,:)
  real, allocatable :: sfrc(:,:,:), sfrc_init(:,:,:), sfrc_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error
  real :: tolerance
  integer :: error_count
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
  call read_parameters(trim(data_dir)//'/params.txt')

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: vsps0'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' ksp0(1)=', ksp0(1), ', ksp0(2)=', ksp0(2)
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(sp(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbct(1:ni, 1:nj, 1:nk, 1:2))
  allocate(sfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(sfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(sfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_1d_int(trim(data_dir)//'/ksp0.bin', ksp0, 1, 2)
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/sp.bin', sp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_4d(trim(data_dir)//'/rbct.bin', rbct, 1, ni, 1, nj, 1, nk, 1, 2)
  call read_array_3d(trim(data_dir)//'/sfrc_in.bin', sfrc_init, 0, ni+1, 0, nj+1, 1, nk)

  write(*,'(A,I6,A,I6)') ' ksp0(1)=', ksp0(1), ', ksp0(2)=', ksp0(2)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/sfrc_ref.bin', sfrc_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    sfrc = sfrc_init
    call kernel_vsps0(ksp0, ni, nj, nk, rst, sp, rbct, sfrc)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    sfrc = sfrc_init
    t_start = omp_get_wtime()

    call kernel_vsps0(ksp0, ni, nj, nk, rst, sp, rbct, sfrc)

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
  call validate_output(sfrc, sfrc_ref, ni, nj, nk, tolerance, max_error, error_count)

  validation_passed = (error_count == 0)

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
  write(*,'(A,I12)')    ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(rst, sp, rbct, sfrc, sfrc_init, sfrc_ref, times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: vsps0
  !-------------------------------------------------------------------
  subroutine kernel_vsps0(ksp0, ni, nj, nk, rst, sp, rbct, sfrc)
    implicit none

    integer, intent(in) :: ksp0(1:2)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: sp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rbct(1:ni, 1:nj, 1:nk, 1:2)
    real, intent(inout) :: sfrc(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    do k = ksp0(2)-1, nk-2
      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-2
        do i = 2, ni-2
          sfrc(i,j,k) = sfrc(i,j,k) &
               - 0.5e0 * (rbct(i,j,k,2) + rbct(i,j,k+1,2)) * rst(i,j,k) * sp(i,j,k)
        end do
      end do
      !$omp end do
    end do

    !$omp end parallel

  end subroutine kernel_vsps0

  !-------------------------------------------------------------------
  ! Read benchmark configuration
  !-------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios == 0) then
      read(10, '(A)', iostat=ios) data_dir
      read(10, *, iostat=ios) num_iter
      read(10, *, iostat=ios) warmup_iter
      read(10, *, iostat=ios) tol
      close(10)
    end if
  end subroutine read_config

  !-------------------------------------------------------------------
  ! Read parameters
  !-------------------------------------------------------------------
  subroutine read_parameters(filename)
    implicit none
    character(len=*), intent(in) :: filename

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
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 1D integer array
  !-------------------------------------------------------------------
  subroutine read_array_1d_int(filename, arr, is, ie)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie
    integer, intent(out) :: arr(is:ie)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_1d_int

  !-------------------------------------------------------------------
  ! Read 3D array
  !-------------------------------------------------------------------
  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

  !-------------------------------------------------------------------
  ! Read 4D array
  !-------------------------------------------------------------------
  subroutine read_array_4d(filename, arr, is, ie, js, je, ks, ke, ls, le)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke, ls, le
    real, intent(out) :: arr(is:ie, js:je, ks:ke, ls:le)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_4d

  !-------------------------------------------------------------------
  ! Validate output
  !-------------------------------------------------------------------
  subroutine validate_output(output, reference, ni, nj, nk, tol, max_err, err_count)
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
      do j = 0, nj+1
        do i = 0, ni+1
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

end program kernel_benchmark_vsps0

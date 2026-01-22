!***********************************************************************
! Kernel Benchmark: vspdmp (s_vspdmp) - OpenACC GPU Version
!***********************************************************************
!
! Source: Src/vspdmp.f90
! Description: Calculate relaxed vertical sponge damping coefficients.
!
!***********************************************************************
program kernel_benchmark_vspdmp
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: vspopt
  real :: vspgpv, vspbar, botgpv, botbar
  real :: cc  ! pi constant

  ! Arrays
  real, allocatable :: zph(:,:,:)
  real, allocatable :: rbct(:,:,:,:), rbct_ref(:,:,:,:)
  real, allocatable :: z1dmax(:)
  integer :: ksp0(1:2), ksp0_ref(1:2)

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
  logical :: ksp0_ref_exists

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
  write(*,'(A)') ' Kernel Benchmark: vspdmp (GPU)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' vspopt=', vspopt
  write(*,'(A,ES12.4,A,ES12.4)') ' vspgpv=', vspgpv, ', vspbar=', vspbar
  write(*,'(A,ES12.4,A,ES12.4)') ' botgpv=', botgpv, ', botbar=', botbar
  write(*,'(A,ES12.4)') ' cc (pi)=', cc
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(zph(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbct(1:ni, 1:nj, 1:nk, 1:2))
  allocate(rbct_ref(1:ni, 1:nj, 1:nk, 1:2))
  allocate(z1dmax(1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/zph.bin', zph, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  inquire(file=trim(data_dir)//'/ksp0_ref.bin', exist=ksp0_ref_exists)
  if (ksp0_ref_exists) then
    call read_array_1d_int(trim(data_dir)//'/ksp0_ref.bin', ksp0_ref, 1, 2)
    write(*,'(A,I6,A,I6)') ' Reference ksp0(1)=', ksp0_ref(1), ', ksp0(2)=', ksp0_ref(2)
  else
    write(*,'(A)') ' Note: ksp0_ref.bin not found, skipping ksp0 validation'
    ksp0_ref(1) = -1
    ksp0_ref(2) = -1
  end if
  call read_array_4d(trim(data_dir)//'/rbct_ref.bin', rbct_ref, 1, ni, 1, nj, 1, nk, 1, 2)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ksp0 = 0
    rbct = 0.0
    z1dmax = 0.0
    call kernel_vspdmp(vspopt, vspgpv, vspbar, botgpv, botbar, cc, &
                       ni, nj, nk, zph, ksp0, rbct, z1dmax)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    ksp0 = 0
    rbct = 0.0
    z1dmax = 0.0

    !$acc wait
    t_start = omp_get_wtime()

    call kernel_vspdmp(vspopt, vspgpv, vspbar, botgpv, botbar, cc, &
                       ni, nj, nk, zph, ksp0, rbct, z1dmax)

    !$acc wait
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
  write(*,'(A,I6,A,I6)') ' Computed ksp0(1)=', ksp0(1), ', ksp0(2)=', ksp0(2)

  error_count = 0
  max_error = 0.0

  ! Validate ksp0 (only if reference data exists)
  if (ksp0_ref_exists) then
    if (ksp0(1) /= ksp0_ref(1)) then
      write(*,'(A,I6,A,I6)') ' WARNING: ksp0(1) mismatch: got ', ksp0(1), ' expected ', ksp0_ref(1)
      error_count = error_count + 1
    end if
    if (ksp0(2) /= ksp0_ref(2)) then
      write(*,'(A,I6,A,I6)') ' WARNING: ksp0(2) mismatch: got ', ksp0(2), ' expected ', ksp0_ref(2)
      error_count = error_count + 1
    end if
  end if

  ! Validate rbct
  call validate_output_4d(rbct, rbct_ref, ni, nj, nk, tolerance, max_error, error_count)

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
  deallocate(zph, rbct, rbct_ref, z1dmax, times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: vspdmp (OpenACC GPU version)
  !-------------------------------------------------------------------
  subroutine kernel_vspdmp(vspopt, vspgpv, vspbar, botgpv, botbar, cc, &
                           ni, nj, nk, zph, ksp0, rbct, z1dmax)
    implicit none

    integer, intent(in) :: vspopt
    real, intent(in) :: vspgpv, vspbar, botgpv, botbar
    real, intent(in) :: cc
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: zph(0:ni+1, 0:nj+1, 1:nk)
    integer, intent(out) :: ksp0(1:2)
    real, intent(out) :: rbct(1:ni, 1:nj, 1:nk, 1:2)
    real, intent(inout) :: z1dmax(1:nk)

    integer :: i, j, k
    integer :: nkm1
    real :: cgpv05, cbar05
    real, parameter :: lim36n = -1.0e+36

    nkm1 = nk - 1
    cgpv05 = 0.5e0 * vspgpv
    cbar05 = 0.5e0 * vspbar

    ksp0(1) = 2
    ksp0(2) = 2

    ! Set the maximum z physical coordinates at each plane.
    !$acc kernels
    !$acc loop independent
    do k = 3, nk
      z1dmax(k) = lim36n
    end do
    !$acc end kernels

    ! Find maximum z at each level - this requires reduction
    ! Since z1dmax reduction across (i,j) is needed, we do this sequentially per k
    do k = 3, nk
      !$acc kernels
      !$acc loop independent reduction(max:z1dmax(k))
      do j = 1, nj-1
        do i = 1, ni-1
          z1dmax(k) = max(zph(i,j,k), z1dmax(k))
        end do
      end do
      !$acc end kernels
    end do

    ! Get the lowest damping level (sequential, run on CPU)
    !$acc wait

    if (vspopt.eq.1) then
      do_k_1: do k = 3, nk
        if (z1dmax(k).gt.botgpv) then
          ksp0(1) = k - 1
          exit do_k_1
        end if
      end do do_k_1
    end if

    do_k_2: do k = 3, nk
      if (z1dmax(k).gt.botbar) then
        ksp0(2) = k - 1
        exit do_k_2
      end if
    end do do_k_2

    ! Finally get the relaxed vertical sponge damping coefficients.
    if (vspopt.eq.1) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-1
        do j = 1, nj-1
          do i = 1, ni-1
            if (zph(i,j,k).gt.botgpv) then
              rbct(i,j,k,1) = cgpv05 * (1.e0 - cos(cc * (zph(i,j,k) - botgpv) &
                   / (zph(i,j,nkm1) - botgpv)))
            else
              rbct(i,j,k,1) = 0.e0
            end if
          end do
        end do
      end do
      !$acc end kernels
    end if

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          if (zph(i,j,k).gt.botbar) then
            rbct(i,j,k,2) = cbar05 * (1.e0 - cos(cc * (zph(i,j,k) - botbar) &
                 / (zph(i,j,nkm1) - botbar)))
          else
            rbct(i,j,k,2) = 0.e0
          end if
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_vspdmp

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

    ! Set defaults
    vspopt = 0
    vspgpv = 0.0
    vspbar = 0.0
    botgpv = 0.0
    botbar = 0.0
    cc = 3.141592653589793

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
          case ('vspopt')
            read(val, *) vspopt
          case ('vspgpv')
            read(val, *) vspgpv
          case ('vspbar')
            read(val, *) vspbar
          case ('botgpv')
            read(val, *) botgpv
          case ('botbar')
            read(val, *) botbar
          case ('cc')
            read(val, *) cc
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
  ! Validate 4D output
  !-------------------------------------------------------------------
  subroutine validate_output_4d(output, reference, ni, nj, nk, tol, max_err, err_count)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: output(1:ni, 1:nj, 1:nk, 1:2)
    real, intent(in) :: reference(1:ni, 1:nj, 1:nk, 1:2)
    real, intent(in) :: tol
    real, intent(inout) :: max_err
    integer, intent(inout) :: err_count

    real :: rel_err
    integer :: i, j, k, l

    do l = 1, 2
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni
            if (abs(reference(i,j,k,l)) > 1.0e-30) then
              rel_err = abs(output(i,j,k,l) - reference(i,j,k,l)) / abs(reference(i,j,k,l))
            else
              rel_err = abs(output(i,j,k,l) - reference(i,j,k,l))
            end if
            if (rel_err > max_err) max_err = rel_err
            if (rel_err > tol) err_count = err_count + 1
          end do
        end do
      end do
    end do
  end subroutine validate_output_4d

end program kernel_benchmark_vspdmp

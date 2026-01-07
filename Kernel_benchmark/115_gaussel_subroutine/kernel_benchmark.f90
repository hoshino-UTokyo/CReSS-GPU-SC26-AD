!***********************************************************************
! Kernel Benchmark: gaussel (s_gaussel)
!***********************************************************************
!
! Source: Src/gaussel.f90
! Description: Solves tridiagonal linear systems using Gauss elimination
!              (Thomas algorithm) or partial pivoting Gauss elimination.
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_gaussel
  use omp_lib
  implicit none

  ! Array dimensions and loop bounds
  integer :: ni, nj, kmax
  integer :: istr, iend, jstr, jend, kstr, kend
  integer :: impopt

  ! Input/output arrays (coefficient matrices and vectors)
  real, allocatable :: rr(:,:,:)
  real, allocatable :: ss(:,:,:)
  real, allocatable :: tt(:,:,:)
  real, allocatable :: ff(:,:,:)
  real, allocatable :: pv(:,:,:)

  ! Arrays to hold initial input data (for re-initialization between iterations)
  real, allocatable :: rr_in(:,:,:)
  real, allocatable :: ss_in(:,:,:)
  real, allocatable :: tt_in(:,:,:)
  real, allocatable :: ff_in(:,:,:)
  real, allocatable :: pv_in(:,:,:)

  ! Reference output for validation
  real, allocatable :: rr_ref(:,:,:)
  real, allocatable :: ss_ref(:,:,:)
  real, allocatable :: tt_ref(:,:,:)
  real, allocatable :: ff_ref(:,:,:)
  real, allocatable :: pv_ref(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', &
       impopt, istr, iend, jstr, jend, kstr, kend, ni, nj, kmax)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: gaussel'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', kmax=', kmax
  write(*,'(A,I6,A,I6,A,I6,A,I6)') ' Loop bounds: i=[', istr, ',', iend, '], j=[', jstr, ',', jend, ']'
  write(*,'(A,I6,A,I6)') '              k=[', kstr, ',', kend, ']'
  write(*,'(A,I6)') ' Implicit option (impopt): ', impopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rr(0:ni+1, 0:nj+1, 1:kmax))
  allocate(ss(0:ni+1, 0:nj+1, 1:kmax))
  allocate(tt(0:ni+1, 0:nj+1, 1:kmax))
  allocate(ff(0:ni+1, 0:nj+1, 1:kmax))
  allocate(pv(0:ni+1, 0:nj+1, 1:kmax))

  allocate(rr_in(0:ni+1, 0:nj+1, 1:kmax))
  allocate(ss_in(0:ni+1, 0:nj+1, 1:kmax))
  allocate(tt_in(0:ni+1, 0:nj+1, 1:kmax))
  allocate(ff_in(0:ni+1, 0:nj+1, 1:kmax))
  allocate(pv_in(0:ni+1, 0:nj+1, 1:kmax))

  allocate(rr_ref(0:ni+1, 0:nj+1, 1:kmax))
  allocate(ss_ref(0:ni+1, 0:nj+1, 1:kmax))
  allocate(tt_ref(0:ni+1, 0:nj+1, 1:kmax))
  allocate(ff_ref(0:ni+1, 0:nj+1, 1:kmax))
  allocate(pv_ref(0:ni+1, 0:nj+1, 1:kmax))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rr_in.bin', rr_in, 0, ni+1, 0, nj+1, 1, kmax)
  call read_array_3d(trim(data_dir)//'/ss_in.bin', ss_in, 0, ni+1, 0, nj+1, 1, kmax)
  call read_array_3d(trim(data_dir)//'/tt_in.bin', tt_in, 0, ni+1, 0, nj+1, 1, kmax)
  call read_array_3d(trim(data_dir)//'/ff_in.bin', ff_in, 0, ni+1, 0, nj+1, 1, kmax)
  call read_array_3d(trim(data_dir)//'/pv_in.bin', pv_in, 0, ni+1, 0, nj+1, 1, kmax)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/rr_ref.bin', rr_ref, 0, ni+1, 0, nj+1, 1, kmax)
  call read_array_3d(trim(data_dir)//'/ss_ref.bin', ss_ref, 0, ni+1, 0, nj+1, 1, kmax)
  call read_array_3d(trim(data_dir)//'/tt_ref.bin', tt_ref, 0, ni+1, 0, nj+1, 1, kmax)
  call read_array_3d(trim(data_dir)//'/ff_ref.bin', ff_ref, 0, ni+1, 0, nj+1, 1, kmax)
  call read_array_3d(trim(data_dir)//'/pv_ref.bin', pv_ref, 0, ni+1, 0, nj+1, 1, kmax)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ! Re-initialize arrays from input data
    rr = rr_in
    ss = ss_in
    tt = tt_in
    ff = ff_in
    pv = pv_in
    call kernel_gaussel(impopt, istr, iend, jstr, jend, kstr, kend, &
         ni, nj, kmax, rr, ss, tt, ff, pv)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ! Re-initialize arrays from input data
    rr = rr_in
    ss = ss_in
    tt = tt_in
    ff = ff_in
    pv = pv_in

    t_start = omp_get_wtime()

    call kernel_gaussel(impopt, istr, iend, jstr, jend, kstr, kend, &
         ni, nj, kmax, rr, ss, tt, ff, pv)

    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output (primary output is ff - the solved vector)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0

  do k = kstr, kend
    do j = jstr, jend
      do i = istr, iend
        rel_error = abs(ff(i,j,k) - ff_ref(i,j,k))
        if (abs(ff_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(ff_ref(i,j,k))
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
  deallocate(rr, ss, tt, ff, pv)
  deallocate(rr_in, ss_in, tt_in, ff_in, pv_in)
  deallocate(rr_ref, ss_ref, tt_ref, ff_ref, pv_ref)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: gaussel
  ! Solve tridiagonal equation with Gauss elimination
  !=====================================================================
  subroutine kernel_gaussel(impopt, istr, iend, jstr, jend, kstr, kend, &
       ni, nj, kmax, rr, ss, tt, ff, pv)
    implicit none

    integer, intent(in) :: impopt
    integer, intent(in) :: istr, iend, jstr, jend, kstr, kend
    integer, intent(in) :: ni, nj, kmax

    real, intent(inout) :: rr(0:ni+1, 0:nj+1, 1:kmax)
    real, intent(inout) :: ss(0:ni+1, 0:nj+1, 1:kmax)
    real, intent(inout) :: tt(0:ni+1, 0:nj+1, 1:kmax)
    real, intent(inout) :: ff(0:ni+1, 0:nj+1, 1:kmax)
    real, intent(inout) :: pv(0:ni+1, 0:nj+1, 1:kmax)

    integer :: i, j, k
    integer :: kem1, kem2
    integer :: kp, kpm1, kpp1
    real :: aa

    kem1 = kend - 1
    kem2 = kend - 2

    !$omp parallel default(shared) private(k)

    ! Solve the tridiagonal equation with the Gauss elimination.
    if (impopt == 1) then

      ! Perform the forward eliminations.
      if (kstr == kend) then

        !$omp do schedule(runtime) private(i,j)
        do j = jstr, jend
          do i = istr, iend
            ff(i,j,kstr) = ff(i,j,kstr) / ss(i,j,kstr)
          end do
        end do
        !$omp end do

      else

        !$omp do schedule(runtime) private(i,j,aa)
        do j = jstr, jend
          do i = istr, iend
            aa = 1.0e0 / ss(i,j,kstr)
            ss(i,j,kstr) = tt(i,j,kstr) * aa
            ff(i,j,kstr) = ff(i,j,kstr) * aa
          end do
        end do
        !$omp end do

      end if

      if (kstr + 2 <= kend) then
        do k = kstr + 1, kend - 1
          !$omp do schedule(runtime) private(i,j,aa)
          do j = jstr, jend
            do i = istr, iend
              aa = 1.0e0 / (ss(i,j,k) - rr(i,j,k) * ss(i,j,k-1))
              ss(i,j,k) = tt(i,j,k) * aa
              ff(i,j,k) = (ff(i,j,k) - rr(i,j,k) * ff(i,j,k-1)) * aa
            end do
          end do
          !$omp end do
        end do
      end if

      if (kstr + 1 <= kend) then
        !$omp do schedule(runtime) private(i,j,aa)
        do j = jstr, jend
          do i = istr, iend
            aa = 1.0e0 / (ss(i,j,kend) - rr(i,j,kend) * ss(i,j,kem1))
            ff(i,j,kend) = (ff(i,j,kend) - rr(i,j,kend) * ff(i,j,kem1)) * aa
          end do
        end do
        !$omp end do
      end if

      ! Perform the back substitutions.
      if (kstr + 1 <= kend) then
        do k = kend - 1, kstr, -1
          !$omp do schedule(runtime) private(i,j)
          do j = jstr, jend
            do i = istr, iend
              ff(i,j,k) = ff(i,j,k) - ss(i,j,k) * ff(i,j,k+1)
            end do
          end do
          !$omp end do
        end do
      end if

    ! Solve with partial pivoting Gauss elimination
    else if (impopt == 2) then

      ! Perform the forward eliminations.
      !$omp do schedule(runtime) private(i,j,aa)
      do j = 2, nj - 2
        do i = 2, ni - 2
          aa = 1.0e0 / ss(i,j,kstr)
          ss(i,j,kstr) = tt(i,j,kstr) * aa
          ff(i,j,kstr) = ff(i,j,kstr) * aa
        end do
      end do
      !$omp end do

      !$omp do schedule(runtime)
      do k = kstr, kend
        pv(0,0,k) = real(k) + 0.1e0
      end do
      !$omp end do

      do k = kstr, kend
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj - 2
          do i = 2, ni - 2
            pv(i,j,k) = pv(0,0,k)
          end do
        end do
        !$omp end do
      end do

      do k = kstr + 1, kend - 2
        !$omp do schedule(runtime) private(i,j,kp,kpm1,aa)
        do j = 2, nj - 2
          do i = 2, ni - 2
            if (abs(ss(i,j,int(pv(i,j,k)))) < abs(rr(i,j,k+1))) then
              aa = pv(i,j,k)
              pv(i,j,k) = pv(i,j,k+1)
              pv(i,j,k+1) = aa
            end if

            kp = int(pv(i,j,k))
            kpm1 = int(pv(i,j,k-1))

            aa = 1.0e0 / (ss(i,j,kp) - rr(i,j,kp) * ss(i,j,kpm1))
            ss(i,j,kp) = tt(i,j,kp) * aa
            ff(i,j,kp) = (ff(i,j,kp) - rr(i,j,kp) * ff(i,j,kpm1)) * aa
          end do
        end do
        !$omp end do
      end do

      !$omp do schedule(runtime) private(i,j,kpm1,aa)
      do j = 2, nj - 2
        do i = 2, ni - 2
          kpm1 = int(pv(i,j,kem2))

          aa = 1.0e0 / (ss(i,j,kem1) - rr(i,j,kem1) * ss(i,j,kpm1))
          ss(i,j,kem1) = tt(i,j,kem1) * aa
          ff(i,j,kem1) = (ff(i,j,kem1) - rr(i,j,kem1) * ff(i,j,kpm1)) * aa

          aa = 1.0e0 / (ss(i,j,kend) - rr(i,j,kend) * ss(i,j,kem1))
          ff(i,j,kend) = (ff(i,j,kend) - rr(i,j,kend) * ff(i,j,kem1)) * aa
        end do
      end do
      !$omp end do

      ! Perform the back substitutions.
      do k = kend - 1, kstr, -1
        !$omp do schedule(runtime) private(i,j,kp,kpp1)
        do j = 2, nj - 2
          do i = 2, ni - 2
            kp = int(pv(i,j,k))
            kpp1 = int(pv(i,j,k+1))
            ff(i,j,kp) = ff(i,j,kp) - ss(i,j,kp) * ff(i,j,kpp1)
          end do
        end do
        !$omp end do
      end do

    end if

    !$omp end parallel

  end subroutine kernel_gaussel

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    character(len=256) :: config_file
    integer :: ios
    logical :: exists

    ! Default values
    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    config_file = 'benchmark.conf'
    inquire(file=config_file, exist=exists)

    if (exists) then
      open(unit=10, file=config_file, status='old', iostat=ios)
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
  subroutine read_parameters(filename, impopt, istr, iend, jstr, jend, &
       kstr, kend, ni, nj, kmax)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: impopt, istr, iend, jstr, jend, kstr, kend
    integer, intent(out) :: ni, nj, kmax

    integer :: ios, fpimpopt
    character(len=256) :: line
    character(len=64) :: name
    character(len=64) :: val_str

    ! Default values
    impopt = 1
    istr = 1
    iend = 1
    jstr = 1
    jend = 1
    kstr = 1
    kend = 1
    ni = 1
    nj = 1
    kmax = 1

    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open parameter file: ', trim(filename)
      stop 1
    end if

    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit

      ! Parse "name = value" format
      if (index(line, '=') > 0) then
        read(line, *) name
        val_str = line(index(line, '=')+1:)

        select case (trim(name))
        case ('fpimpopt')
          read(val_str, *) fpimpopt
          ! Convert fpimpopt to impopt (typically impopt = 1 for Thomas algorithm)
          impopt = 1
        case ('istr')
          read(val_str, *) istr
        case ('iend')
          read(val_str, *) iend
        case ('jstr')
          read(val_str, *) jstr
        case ('jend')
          read(val_str, *) jend
        case ('kstr')
          read(val_str, *) kstr
        case ('kend')
          read(val_str, *) kend
        case ('ni')
          read(val_str, *) ni
        case ('nj')
          read(val_str, *) nj
        case ('kmax')
          read(val_str, *) kmax
        end select
      end if
    end do

    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary array reader
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

end program kernel_benchmark_gaussel

!***********************************************************************
! Kernel Benchmark: adjstuv_sec2 (s_adjstuv - second parallel section)
!***********************************************************************
!
! Source: Src/adjstuv.f90
! Description: Applies adjustment value to u and v velocity components
!              at boundary faces.
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_adjstuv_sec2
  use omp_lib
  implicit none

  ! Parameters
  integer :: wbc, ebc, mpopt, mfcopt
  integer :: ni, nj, nk
  real :: adj
  integer :: ebe, ebn, ebs, ebw
  integer :: iend, istr, isub, jend, jstr, jsub
  integer :: nisub, njsub

  ! Input arrays
  real, allocatable :: rst8u(:,:,:), rst8v(:,:,:)

  ! Input/output arrays
  real, allocatable :: uf(:,:,:), vf(:,:,:)
  real, allocatable :: uf_in(:,:,:), vf_in(:,:,:)
  real, allocatable :: uf_ref(:,:,:), vf_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
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
  call read_parameters(trim(data_dir)//'/params.txt')

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: adjstuv_sec2'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' adj=', adj
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rst8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rst8u.bin', rst8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8v.bin', rst8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/uf_in.bin', uf_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vf_in.bin', vf_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/uf_ref.bin', uf_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vf_ref.bin', vf_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    uf = uf_in
    vf = vf_in
    call kernel_adjstuv_sec2(wbc, ebc, ni, nj, nk, &
         istr, iend, jstr, jend, &
         ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
         adj, rst8u, rst8v, uf, vf)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    uf = uf_in
    vf = vf_in

    t_start = omp_get_wtime()
    call kernel_adjstuv_sec2(wbc, ebc, ni, nj, nk, &
         istr, iend, jstr, jend, &
         ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
         adj, rst8u, rst8v, uf, vf)
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  total_errors = 0

  call validate_3d(uf, uf_ref, ni, nj, nk, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(vf, vf_ref, ni, nj, nk, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  validation_passed = (total_errors == 0)

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
  write(*,'(A,I12)') ' Error count:        ', total_errors
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(rst8u, rst8v)
  deallocate(uf, vf, uf_in, vf_in, uf_ref, vf_ref)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: adjstuv section 2
  ! Applies adjustment value to u and v at boundary faces
  !=====================================================================
  subroutine kernel_adjstuv_sec2(wbc, ebc, ni, nj, nk, &
       istr, iend, jstr, jend, &
       ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
       adj, rst8u, rst8v, uf, vf)
    implicit none

    integer, intent(in) :: wbc, ebc
    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: istr, iend, jstr, jend
    integer, intent(in) :: ebw, ebe, ebs, ebn
    integer, intent(in) :: isub, jsub, nisub, njsub
    real, intent(in) :: adj
    real, intent(in) :: rst8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: uf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vf(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$omp parallel default(shared)

    if (ebw == 1 .and. isub == 0 .and. abs(wbc) /= 1) then

      !$omp do schedule(runtime) private(j,k)
      do k = 1, nk-1
        do j = jstr, jend
          uf(1,j,k) = uf(1,j,k) + adj / rst8u(1,j,k)
        end do
      end do
      !$omp end do

    end if

    if (ebe == 1 .and. isub == nisub-1 .and. abs(ebc) /= 1) then

      !$omp do schedule(runtime) private(j,k)
      do k = 1, nk-1
        do j = jstr, jend
          uf(ni,j,k) = uf(ni,j,k) - adj / rst8u(ni,j,k)
        end do
      end do
      !$omp end do

    end if

    if (ebs == 1 .and. jsub == 0) then

      !$omp do schedule(runtime) private(i,k)
      do k = 1, nk-1
        do i = istr, iend
          vf(i,1,k) = vf(i,1,k) + adj / rst8v(i,1,k)
        end do
      end do
      !$omp end do

    end if

    if (ebn == 1 .and. jsub == njsub-1) then

      !$omp do schedule(runtime) private(i,k)
      do k = 1, nk-1
        do i = istr, iend
          vf(i,nj,k) = vf(i,nj,k) - adj / rst8v(i,nj,k)
        end do
      end do
      !$omp end do

    end if

    !$omp end parallel

  end subroutine kernel_adjstuv_sec2

  !=====================================================================
  ! Validation subroutine
  !=====================================================================
  subroutine validate_3d(arr, ref, ni, nj, nk, tol, max_err, err_count)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:nk)
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
          rel_err = abs(arr(i,j,k) - ref(i,j,k))
          if (abs(ref(i,j,k)) > 1.0e-10) then
            rel_err = rel_err / abs(ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

  end subroutine validate_3d

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
  subroutine read_parameters(filename)
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
          case ('wbc');    read(val, *) wbc
          case ('ebc');    read(val, *) ebc
          case ('mpopt');  read(val, *) mpopt
          case ('mfcopt'); read(val, *) mfcopt
          case ('ni');     read(val, *) ni
          case ('nj');     read(val, *) nj
          case ('nk');     read(val, *) nk
          case ('adj');    read(val, *) adj
          case ('ebe');    read(val, *) ebe
          case ('ebn');    read(val, *) ebn
          case ('ebs');    read(val, *) ebs
          case ('ebw');    read(val, *) ebw
          case ('iend');   read(val, *) iend
          case ('istr');   read(val, *) istr
          case ('isub');   read(val, *) isub
          case ('jend');   read(val, *) jend
          case ('jstr');   read(val, *) jstr
          case ('jsub');   read(val, *) jsub
          case ('nisub');  read(val, *) nisub
          case ('njsub');  read(val, *) njsub
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

end program kernel_benchmark_adjstuv_sec2

!***********************************************************************
! Kernel Benchmark: bcbase (s_bcbase)
!***********************************************************************
!
! Source: Src/bcbase.f90
! Description: Set bottom and top boundary conditions for base state
!              variables (ubr, vbr, ptbr, qvbr, ptvbr, pibr, pbr, rbr)
!              using extrapolation.
!
!***********************************************************************
program kernel_benchmark_bcbase
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: smtopt
  real :: g_val, cp_val, rd_val, p0_val

  ! Input arrays
  real, allocatable :: zph8s(:,:,:)

  ! Input/output arrays
  real, allocatable :: ubr(:,:,:), vbr(:,:,:), pbr(:,:,:), ptbr(:,:,:)
  real, allocatable :: qvbr(:,:,:), rbr(:,:,:), pibr(:,:,:), ptvbr(:,:,:)

  ! Reference output for validation
  real, allocatable :: ubr_ref(:,:,:), vbr_ref(:,:,:), pbr_ref(:,:,:), ptbr_ref(:,:,:)
  real, allocatable :: qvbr_ref(:,:,:), rbr_ref(:,:,:), pibr_ref(:,:,:), ptvbr_ref(:,:,:)

  ! Backup arrays for iteration
  real, allocatable :: ubr_in(:,:,:), vbr_in(:,:,:), pbr_in(:,:,:), ptbr_in(:,:,:)
  real, allocatable :: qvbr_in(:,:,:), rbr_in(:,:,:), pibr_in(:,:,:), ptvbr_in(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error, arr_error
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
  call read_parameters(trim(data_dir)//'/params.txt', smtopt, ni, nj, nk, &
                       g_val, cp_val, rd_val, p0_val)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: bcbase'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' smtopt=', smtopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(zph8s(0:ni+1, 0:nj+1, 1:nk))
  allocate(ubr(0:ni+1, 0:nj+1, 1:nk), vbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(pbr(0:ni+1, 0:nj+1, 1:nk), ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvbr(0:ni+1, 0:nj+1, 1:nk), rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(pibr(0:ni+1, 0:nj+1, 1:nk), ptvbr(0:ni+1, 0:nj+1, 1:nk))

  allocate(ubr_ref(0:ni+1, 0:nj+1, 1:nk), vbr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(pbr_ref(0:ni+1, 0:nj+1, 1:nk), ptbr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvbr_ref(0:ni+1, 0:nj+1, 1:nk), rbr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(pibr_ref(0:ni+1, 0:nj+1, 1:nk), ptvbr_ref(0:ni+1, 0:nj+1, 1:nk))

  allocate(ubr_in(0:ni+1, 0:nj+1, 1:nk), vbr_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(pbr_in(0:ni+1, 0:nj+1, 1:nk), ptbr_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvbr_in(0:ni+1, 0:nj+1, 1:nk), rbr_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(pibr_in(0:ni+1, 0:nj+1, 1:nk), ptvbr_in(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/zph8s.bin', zph8s, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ubr_in.bin', ubr_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vbr_in.bin', vbr_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pbr_in.bin', pbr_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptbr_in.bin', ptbr_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvbr_in.bin', qvbr_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rbr_in.bin', rbr_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pibr_in.bin', pibr_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptvbr_in.bin', ptvbr_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ubr_ref.bin', ubr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vbr_ref.bin', vbr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pbr_ref.bin', pbr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptbr_ref.bin', ptbr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvbr_ref.bin', qvbr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rbr_ref.bin', rbr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pibr_ref.bin', pibr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptvbr_ref.bin', ptvbr_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ubr = ubr_in; vbr = vbr_in; pbr = pbr_in; ptbr = ptbr_in
    qvbr = qvbr_in; rbr = rbr_in; pibr = pibr_in; ptvbr = ptvbr_in
    call kernel_bcbase(ni, nj, nk, g_val, cp_val, rd_val, p0_val, &
                       zph8s, ubr, vbr, pbr, ptbr, qvbr, rbr, pibr, ptvbr)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ubr = ubr_in; vbr = vbr_in; pbr = pbr_in; ptbr = ptbr_in
    qvbr = qvbr_in; rbr = rbr_in; pibr = pibr_in; ptvbr = ptvbr_in

    t_start = omp_get_wtime()
    call kernel_bcbase(ni, nj, nk, g_val, cp_val, rd_val, p0_val, &
                       zph8s, ubr, vbr, pbr, ptbr, qvbr, rbr, pibr, ptvbr)
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
  error_count = 0

  call validate_array(ubr, ubr_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error
  call validate_array(vbr, vbr_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error
  call validate_array(pbr, pbr_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error
  call validate_array(ptbr, ptbr_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error
  call validate_array(qvbr, qvbr_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error
  call validate_array(rbr, rbr_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error
  call validate_array(pibr, pibr_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error
  call validate_array(ptvbr, ptvbr_ref, ni, nj, nk, tolerance, arr_error, error_count)
  if (arr_error > max_error) max_error = arr_error

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
  deallocate(zph8s)
  deallocate(ubr, vbr, pbr, ptbr, qvbr, rbr, pibr, ptvbr)
  deallocate(ubr_ref, vbr_ref, pbr_ref, ptbr_ref, qvbr_ref, rbr_ref, pibr_ref, ptvbr_ref)
  deallocate(ubr_in, vbr_in, pbr_in, ptbr_in, qvbr_in, rbr_in, pibr_in, ptvbr_in)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: bcbase
  ! Set bottom and top boundary conditions for base state variables
  !=====================================================================
  subroutine kernel_bcbase(ni, nj, nk, g, cp, rd, p0, &
                           zph8s, ubr, vbr, pbr, ptbr, qvbr, rbr, pibr, ptvbr)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: g, cp, rd, p0
    real, intent(in) :: zph8s(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: ubr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: pbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: ptbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: qvbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: pibr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: ptvbr(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j
    integer :: nkm1, nkm2
    real :: gdvcp2, cpdvrd

    nkm1 = nk - 1
    nkm2 = nk - 2
    gdvcp2 = 2.0 * g / cp
    cpdvrd = cp / rd

    !$omp parallel default(shared)

    ! Set the bottom and the top boundary conditions for the base state velocity.
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 1, ni
        ubr(i,j,1) = ubr(i,j,2)
        ubr(i,j,nkm1) = ubr(i,j,nkm2)
      end do
    end do
    !$omp end do

    !$omp do schedule(runtime) private(i,j)
    do j = 1, nj
      do i = 0, ni
        vbr(i,j,1) = vbr(i,j,2)
        vbr(i,j,nkm1) = vbr(i,j,nkm2)
      end do
    end do
    !$omp end do

    ! Set the bottom and the top boundary conditions for the base state
    ! potential temperature and water vapor mixing ratio.
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 0, ni
        ptbr(i,j,1) = ptbr(i,j,2)
        ptbr(i,j,nkm1) = ptbr(i,j,nkm2)
        qvbr(i,j,1) = qvbr(i,j,2)
        qvbr(i,j,nkm1) = qvbr(i,j,nkm2)
      end do
    end do
    !$omp end do

    ! Set the bottom and the top boundary conditions for the base state
    ! virtual potential temperature.
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 0, ni
        ptvbr(i,j,1) = ptvbr(i,j,2)
        ptvbr(i,j,nkm1) = ptvbr(i,j,nkm2)
      end do
    end do
    !$omp end do

    ! Set the bottom and the top boundary conditions for the base state
    ! Exner function.
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 0, ni
        pibr(i,j,1) = pibr(i,j,2) &
          + gdvcp2 * (zph8s(i,j,2) - zph8s(i,j,1)) / (ptvbr(i,j,1) + ptvbr(i,j,2))
        pibr(i,j,nkm1) = pibr(i,j,nkm2) &
          - gdvcp2 * (zph8s(i,j,nkm1) - zph8s(i,j,nkm2)) / (ptvbr(i,j,nkm2) + ptvbr(i,j,nkm1))
      end do
    end do
    !$omp end do

    ! Set the bottom and the top boundary conditions for the base state
    ! pressure and the base state density.
    !$omp do schedule(runtime) private(i,j)
    do j = 0, nj
      do i = 0, ni
        pbr(i,j,1) = p0 * exp(cpdvrd * log(pibr(i,j,1)))
        pbr(i,j,nkm1) = p0 * exp(cpdvrd * log(pibr(i,j,nkm1)))
        rbr(i,j,1) = pbr(i,j,1) / (rd * ptvbr(i,j,1) * pibr(i,j,1))
        rbr(i,j,nkm1) = pbr(i,j,nkm1) / (rd * ptvbr(i,j,nkm1) * pibr(i,j,nkm1))
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine kernel_bcbase

  !=====================================================================
  ! Validate array
  !=====================================================================
  subroutine validate_array(arr, arr_ref, ni, nj, nk, tol, max_err, err_count)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: arr_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(inout) :: err_count

    integer :: i, j, k
    real :: rel_error

    max_err = 0.0
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          rel_error = abs(arr(i,j,k) - arr_ref(i,j,k))
          if (abs(arr_ref(i,j,k)) > 1.0e-10) then
            rel_error = rel_error / abs(arr_ref(i,j,k))
          end if
          if (rel_error > max_err) max_err = rel_error
          if (rel_error > tol) err_count = err_count + 1
        end do
      end do
    end do
  end subroutine validate_array

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
  subroutine read_parameters(filename, smtopt, ni, nj, nk, g, cp, rd, p0)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: smtopt, ni, nj, nk
    real, intent(out) :: g, cp, rd, p0

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    smtopt = 0
    ni = 1
    nj = 1
    nk = 1
    g = 9.8
    cp = 1004.0
    rd = 287.0
    p0 = 100000.0

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
          case ('smtopt')
            read(val, *) smtopt
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('nkm1')
            read(val, *) nk
            nk = nk + 1  ! nk = nkm1 + 1
          case ('g')
            read(val, *) g
          case ('cp')
            read(val, *) cp
          case ('rd')
            read(val, *) rd
          case ('p0')
            read(val, *) p0
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

end program kernel_benchmark_bcbase

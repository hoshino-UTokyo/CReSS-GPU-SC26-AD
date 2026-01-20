!***********************************************************************
! Kernel Benchmark: forcesfc (s_forcesfc)
!***********************************************************************
!
! Source: Src/forcesfc.f90
! Description: Get the surface flux to bottom boundary for potential
!              temperature, water vapor, and velocity components.
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_forcesfc
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  character(len=5) :: fmois
  real :: epsav

  ! Input arrays (3D)
  real, allocatable :: j31(:,:,:)
  real, allocatable :: j32(:,:,:)
  real, allocatable :: ptbr(:,:,:)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)
  real, allocatable :: w(:,:,:)
  real, allocatable :: ptp(:,:,:)
  real, allocatable :: qv(:,:,:)
  real, allocatable :: ptv(:,:,:)

  ! Input arrays (2D)
  real, allocatable :: qvsfc(:,:)
  real, allocatable :: ce(:,:)
  real, allocatable :: ct(:,:)
  real, allocatable :: cq(:,:)

  ! Output arrays
  real, allocatable :: ufrc(:,:,:)
  real, allocatable :: vfrc(:,:,:)
  real, allocatable :: ptfrc(:,:,:)
  real, allocatable :: qvfrc(:,:,:)

  ! Reference output for validation
  real, allocatable :: ufrc_ref(:,:,:)
  real, allocatable :: vfrc_ref(:,:,:)
  real, allocatable :: ptfrc_ref(:,:,:)
  real, allocatable :: qvfrc_ref(:,:,:)

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
  integer :: iter, i, j

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', &
       fmois, ni, nj, nk, epsav)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: forcesfc'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,A)') ' fmois=', trim(fmois)
  write(*,'(A,ES12.4)') ' epsav=', epsav
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(w(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptv(0:ni+1, 0:nj+1, 1:nk))

  allocate(qvsfc(0:ni+1, 0:nj+1))
  allocate(ce(0:ni+1, 0:nj+1))
  allocate(ct(0:ni+1, 0:nj+1))
  allocate(cq(0:ni+1, 0:nj+1))

  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvfrc(0:ni+1, 0:nj+1, 1:nk))

  allocate(ufrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvfrc_ref(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptv.bin', ptv, 0, ni+1, 0, nj+1, 1, nk)

  call read_array_2d(trim(data_dir)//'/qvsfc.bin', qvsfc, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/ce.bin', ce, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/ct.bin', ct, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/cq.bin', cq, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ufrc_ref.bin', ufrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_ref.bin', vfrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptfrc_ref.bin', ptfrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvfrc_ref.bin', qvfrc_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_forcesfc(fmois, ni, nj, nk, epsav, j31, j32, ptbr, u, v, w, &
         ptp, qv, ptv, qvsfc, ce, ct, cq, ufrc, vfrc, ptfrc, qvfrc)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    t_start = omp_get_wtime()

    call kernel_forcesfc(fmois, ni, nj, nk, epsav, j31, j32, ptbr, u, v, w, &
         ptp, qv, ptv, qvsfc, ce, ct, cq, ufrc, vfrc, ptfrc, qvfrc)

    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output (only k=1 level is computed)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0

  ! Validate ptfrc at k=1
  do j = 1, nj-1
    do i = 1, ni-1
      rel_error = abs(ptfrc(i,j,1) - ptfrc_ref(i,j,1))
      if (abs(ptfrc_ref(i,j,1)) > 1.0e-20) then
        rel_error = rel_error / abs(ptfrc_ref(i,j,1))
      end if
      if (rel_error > max_error) max_error = rel_error
      if (rel_error > tolerance) error_count = error_count + 1
    end do
  end do

  ! Validate qvfrc at k=1
  do j = 1, nj-1
    do i = 1, ni-1
      rel_error = abs(qvfrc(i,j,1) - qvfrc_ref(i,j,1))
      if (abs(qvfrc_ref(i,j,1)) > 1.0e-20) then
        rel_error = rel_error / abs(qvfrc_ref(i,j,1))
      end if
      if (rel_error > max_error) max_error = rel_error
      if (rel_error > tolerance) error_count = error_count + 1
    end do
  end do

  ! Validate ufrc at k=1
  do j = 1, nj-1
    do i = 2, ni-1
      rel_error = abs(ufrc(i,j,1) - ufrc_ref(i,j,1))
      if (abs(ufrc_ref(i,j,1)) > 1.0e-20) then
        rel_error = rel_error / abs(ufrc_ref(i,j,1))
      end if
      if (rel_error > max_error) max_error = rel_error
      if (rel_error > tolerance) error_count = error_count + 1
    end do
  end do

  ! Validate vfrc at k=1
  do j = 2, nj-1
    do i = 1, ni-1
      rel_error = abs(vfrc(i,j,1) - vfrc_ref(i,j,1))
      if (abs(vfrc_ref(i,j,1)) > 1.0e-20) then
        rel_error = rel_error / abs(vfrc_ref(i,j,1))
      end if
      if (rel_error > max_error) max_error = rel_error
      if (rel_error > tolerance) error_count = error_count + 1
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
  deallocate(j31, j32, ptbr, u, v, w, ptp, qv, ptv)
  deallocate(qvsfc, ce, ct, cq)
  deallocate(ufrc, vfrc, ptfrc, qvfrc)
  deallocate(ufrc_ref, vfrc_ref, ptfrc_ref, qvfrc_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: forcesfc
  ! Get surface flux to bottom boundary
  !=====================================================================
  subroutine kernel_forcesfc(fmois, ni, nj, nk, epsav, j31, j32, ptbr, &
       u, v, w, ptp, qv, ptv, qvsfc, ce, ct, cq, ufrc, vfrc, ptfrc, qvfrc)
    implicit none

    character(len=5), intent(in) :: fmois
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: epsav

    real, intent(in) :: j31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptv(0:ni+1, 0:nj+1, 1:nk)

    real, intent(in) :: qvsfc(0:ni+1, 0:nj+1)
    real, intent(in) :: ce(0:ni+1, 0:nj+1)
    real, intent(in) :: ct(0:ni+1, 0:nj+1)
    real, intent(in) :: cq(0:ni+1, 0:nj+1)

    real, intent(out) :: ufrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: vfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: ptfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: qvfrc(0:ni+1, 0:nj+1, 1:nk)

    ! Local variables
    integer :: i, j
    real :: j318u, j328v
    real :: xcomp, ycomp, zcomp

    !$omp parallel default(shared)

    ! Get the surface flux for the potential temperature
    if (fmois(1:3) == 'dry') then

      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          ptfrc(i,j,1) = ct(i,j)*(ptv(i,j,2) - ptv(i,j,1))
        end do
      end do
      !$omp end do

    else if (fmois(1:5) == 'moist') then

      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          ptfrc(i,j,1) = ct(i,j)*((ptbr(i,j,2) + ptp(i,j,2)) &
               - ptv(i,j,1)*(1.e0 + qvsfc(i,j))/(1.e0 + epsav*qvsfc(i,j)))
        end do
      end do
      !$omp end do

    end if

    ! Get the surface flux for water vapor mixing ratio
    if (fmois(1:3) == 'dry') then

      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          qvfrc(i,j,1) = 0.e0
        end do
      end do
      !$omp end do

    else if (fmois(1:5) == 'moist') then

      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          qvfrc(i,j,1) = cq(i,j)*(qv(i,j,2) - qvsfc(i,j))
        end do
      end do
      !$omp end do

    end if

    ! Get the surface flux for the x components of velocity
    !$omp do schedule(runtime) private(i,j,j318u,xcomp,zcomp)
    do j = 1, nj-1
      do i = 2, ni-1
        j318u = j31(i,j,2) + j31(i,j,3)
        xcomp = 1.e0/sqrt(4.e0 + j318u*j318u)
        zcomp = .125e0*j318u*xcomp
        ufrc(i,j,1) = (ce(i-1,j) + ce(i,j))*(u(i,j,2)*xcomp &
             + ((w(i-1,j,2) + w(i,j,3)) + (w(i-1,j,3) + w(i,j,2)))*zcomp)
      end do
    end do
    !$omp end do

    ! Get the surface flux for the y components of velocity
    !$omp do schedule(runtime) private(i,j,j328v,ycomp,zcomp)
    do j = 2, nj-1
      do i = 1, ni-1
        j328v = j32(i,j,2) + j32(i,j,3)
        ycomp = 1.e0/sqrt(4.e0 + j328v*j328v)
        zcomp = .125e0*j328v*ycomp
        vfrc(i,j,1) = (ce(i,j-1) + ce(i,j))*(v(i,j,2)*ycomp &
             + ((w(i,j-1,2) + w(i,j,3)) + (w(i,j-1,3) + w(i,j,2)))*zcomp)
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine kernel_forcesfc

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
  subroutine read_parameters(filename, fmois, ni, nj, nk, epsav)
    character(len=*), intent(in) :: filename
    character(len=5), intent(out) :: fmois
    integer, intent(out) :: ni, nj, nk
    real, intent(out) :: epsav

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    fmois = 'dry  '

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
          case ('fmois')
            read(val, '(A)') fmois
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('epsav')
            read(val, *) epsav
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary array readers
  !=====================================================================
  subroutine read_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: arr(i1:i2, j1:j2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_2d

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

end program kernel_benchmark_forcesfc

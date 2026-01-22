!***********************************************************************
! Kernel Benchmark: disptke (s_disptke) - GPU Version
!***********************************************************************
!
! Source: Src/disptke.f90
! Description: Calculate TKE dissipation term using turbulent length scale,
!              with different formulations for isotropic/anisotropic cases.
!
!***********************************************************************
program kernel_benchmark_disptke
  use omp_lib
  implicit none

  ! Constants from commath
  real, parameter :: oned3 = 1.0 / 3.0   ! One third
  real, parameter :: eps = 1.0e-20       ! Small number

  ! Grid dimensions
  integer :: ni, nj, nk

  ! Parameters
  real :: dx, dy, dz
  integer :: mpopt, mfcopt, isoopt

  ! Derived
  real :: dz05, ds308

  ! Arrays
  real, allocatable :: jcb(:,:,:)
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: rst(:,:,:)
  real, allocatable :: priv(:,:,:)
  real, allocatable :: tke(:,:,:)
  real, allocatable :: tkefrc(:,:,:)
  real, allocatable :: tkefrc_input(:,:,:)
  real, allocatable :: tkefrc_ref(:,:,:)

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
  integer :: max_i, max_j, max_k

  ! Loop variables
  integer :: iter, i, j, k

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, dx, dy, dz, mpopt, mfcopt, isoopt, ds308, dz05)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: disptke (GPU)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dx: ', dx
  write(*,'(A,ES12.4)') ' dy: ', dy
  write(*,'(A,ES12.4)') ' dz: ', dz
  write(*,'(A,I6)') ' isoopt: ', isoopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(priv(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_input(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d4(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/priv.bin', priv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tke.bin', tke, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tkefrc_in.bin', tkefrc_input, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/tkefrc_ref.bin', tkefrc_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    tkefrc = tkefrc_input
    call kernel_disptke(isoopt, mfcopt, mpopt, dz05, ds308, ni, nj, nk, jcb, rmf, rst, priv, tke, tkefrc)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    tkefrc = tkefrc_input

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_disptke(isoopt, mfcopt, mpopt, dz05, ds308, ni, nj, nk, jcb, rmf, rst, priv, tke, tkefrc)
    !$acc wait
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
  max_i = 2
  max_j = 2
  max_k = 2

  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(tkefrc(i,j,k) - tkefrc_ref(i,j,k))
        if (abs(tkefrc_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(tkefrc_ref(i,j,k))
        end if
        if (rel_error > max_error) then
          max_error = rel_error
          max_i = i
          max_j = j
          max_k = k
        end if
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do
  write(*,'(A,I0,A,I0,A,I0)') ' Max error at (i,j,k)=(', max_i, ',', max_j, ',', max_k, ')'
  write(*,'(A,ES15.8)') ' GPU value: ', tkefrc(max_i, max_j, max_k)
  write(*,'(A,ES15.8)') ' Ref value: ', tkefrc_ref(max_i, max_j, max_k)
  write(*,'(A,ES15.8)') ' Input val: ', tkefrc_input(max_i, max_j, max_k)
  write(*,'(A,ES15.8)') ' priv:      ', priv(max_i, max_j, max_k)
  write(*,'(A,ES15.8)') ' tke:       ', tke(max_i, max_j, max_k)
  write(*,'(A,ES15.8)') ' rst:       ', rst(max_i, max_j, max_k)
  write(*,'(A,ES15.8)') ' jcb:       ', jcb(max_i, max_j, max_k)
  write(*,'(A,ES15.8)') ' rmf(i,j,2):', rmf(max_i, max_j, 2)
  ! CPU-side debug calculation
  block
    real :: ln_cpu, term_cpu
    ln_cpu = (priv(max_i,max_j,max_k) - 1.0) * exp(oned3 * log(ds308 * rmf(max_i,max_j,2) * jcb(max_i,max_j,max_k))) + eps
    term_cpu = (0.37 * priv(max_i,max_j,max_k) - 0.18) * rst(max_i,max_j,max_k) * tke(max_i,max_j,max_k) * sqrt(tke(max_i,max_j,max_k)) / ln_cpu
    write(*,'(A,ES15.8)') ' CPU ln:    ', ln_cpu
    write(*,'(A,ES15.8)') ' CPU term:  ', term_cpu
    write(*,'(A,ES15.8)') ' CPU result:', tkefrc_input(max_i,max_j,max_k) - term_cpu
  end block

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
  deallocate(jcb, rmf, rst, priv, tke, tkefrc, tkefrc_input, tkefrc_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: disptke (GPU version)
  !=====================================================================
  subroutine kernel_disptke(isoopt, mfcopt, mpopt, dz05, ds308, ni, nj, nk, jcb, rmf, rst, priv, tke, tkefrc)
    implicit none

    integer, intent(in) :: isoopt, mfcopt, mpopt
    real, intent(in) :: dz05, ds308
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: jcb(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rmf(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: priv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tke(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tkefrc(0:ni+1, 0:nj+1, 1:nk)

    real :: ln
    integer :: i, j, k

    ! Isotropic case
    if (isoopt == 1) then

      if (mfcopt == 0) then

        !$acc kernels
        !$acc loop independent collapse(3) private(ln)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-2
              ln = (priv(i,j,k) - 1.0) * exp(oned3 * log(ds308 * jcb(i,j,k))) + eps
              tkefrc(i,j,k) = tkefrc(i,j,k) - (0.37 * priv(i,j,k) - 0.18) &
                   * rst(i,j,k) * tke(i,j,k) * sqrt(tke(i,j,k)) / ln
            end do
          end do
        end do
        !$acc end kernels

      else

        if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then

          !$acc kernels
          !$acc loop independent collapse(3) private(ln)
          do k = 2, nk-2
            do j = 2, nj-2
              do i = 2, ni-2
                ln = (priv(i,j,k) - 1.0) * exp(oned3 * log(ds308 * rmf(i,j,2) * jcb(i,j,k))) + eps
                tkefrc(i,j,k) = tkefrc(i,j,k) - (0.37 * priv(i,j,k) - 0.18) &
                     * rst(i,j,k) * tke(i,j,k) * sqrt(tke(i,j,k)) / ln
              end do
            end do
          end do
          !$acc end kernels

        else

          !$acc kernels
          !$acc loop independent collapse(3) private(ln)
          do k = 2, nk-2
            do j = 2, nj-2
              do i = 2, ni-2
                ln = (priv(i,j,k) - 1.0) * exp(oned3 * log(ds308 * rmf(i,j,3) * jcb(i,j,k))) + eps
                tkefrc(i,j,k) = tkefrc(i,j,k) - (0.37 * priv(i,j,k) - 0.18) &
                     * rst(i,j,k) * tke(i,j,k) * sqrt(tke(i,j,k)) / ln
              end do
            end do
          end do
          !$acc end kernels

        end if

      end if

    ! Anisotropic case
    else if (isoopt == 2) then

      !$acc kernels
      !$acc loop independent collapse(3) private(ln)
      do k = 2, nk-2
        do j = 2, nj-2
          do i = 2, ni-2
            ln = (priv(i,j,k) - 1.0) * jcb(i,j,k) * dz05 + eps
            tkefrc(i,j,k) = tkefrc(i,j,k) - (0.37 * priv(i,j,k) - 0.18) &
                 * rst(i,j,k) * tke(i,j,k) * sqrt(tke(i,j,k)) / ln
          end do
        end do
      end do
      !$acc end kernels

    end if

  end subroutine kernel_disptke

  !=====================================================================
  ! I/O subroutines
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

  subroutine read_parameters(filename, ni, nj, nk, dx, dy, dz, mpopt, mfcopt, isoopt, ds308_out, dz05_out)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, mpopt, mfcopt, isoopt
    real, intent(out) :: dx, dy, dz, ds308_out, dz05_out
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
          case ('dx')
            read(val, *) dx
          case ('dy')
            read(val, *) dy
          case ('dz')
            read(val, *) dz
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
          case ('isoopt')
            read(val, *) isoopt
          case ('ds308')
            read(val, *) ds308_out
          case ('dz05')
            read(val, *) dz05_out
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

  subroutine read_array_2d4(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: arr(i1:i2, j1:j2, 1:4)
    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_2d4

end program kernel_benchmark_disptke

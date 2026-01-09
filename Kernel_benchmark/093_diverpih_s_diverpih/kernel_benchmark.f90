!***********************************************************************
! Kernel Benchmark: diverpih (s_diverpih)
!***********************************************************************
!
! Source: Src/diverpih.f90
! Description: Calculate horizontal divergence for pressure equation
!              with horizontally explicit and vertically implicit method
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_diverpih
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: trnopt, mpopt, mfcopt
  real :: dziv

  ! Input arrays
  real, allocatable :: j31(:,:,:)
  real, allocatable :: j32(:,:,:)
  real, allocatable :: jcb8u(:,:,:)
  real, allocatable :: jcb8v(:,:,:)
  real, allocatable :: mf(:,:)
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: rmf8u(:,:,:)
  real, allocatable :: rmf8v(:,:,:)
  real, allocatable :: rcsq(:,:,:)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)

  ! Output array
  real, allocatable :: pdiv(:,:,:)
  real, allocatable :: pdiv_init(:,:,:)

  ! Work arrays
  real, allocatable :: tmp1(:,:,:)
  real, allocatable :: tmp2(:,:,:)
  real, allocatable :: tmp3(:,:,:)

  ! Reference output for validation
  real, allocatable :: pdiv_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
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
       trnopt, mpopt, mfcopt, dziv, ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: diverpih'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6,A,I6)') ' Options: trnopt=', trnopt, ', mpopt=', mpopt, ', mfcopt=', mfcopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rmf8u(0:ni+1, 0:nj+1, 1:3))
  allocate(rmf8v(0:ni+1, 0:nj+1, 1:3))
  allocate(rcsq(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(pdiv(0:ni+1, 0:nj+1, 1:nk))
  allocate(pdiv_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
  allocate(pdiv_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8u.bin', jcb8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8v.bin', jcb8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/mf.bin', mf, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call read_array_3d(trim(data_dir)//'/rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call read_array_3d(trim(data_dir)//'/rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  call read_array_3d(trim(data_dir)//'/rcsq.bin', rcsq, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pdiv_in.bin', pdiv_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/pdiv_ref.bin', pdiv_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    pdiv = pdiv_init
    call kernel_diverpih(trnopt, mpopt, mfcopt, dziv, ni, nj, nk, &
         j31, j32, jcb8u, jcb8v, mf, rmf, rmf8u, rmf8v, &
         rcsq, u, v, pdiv, tmp1, tmp2, tmp3)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    pdiv = pdiv_init
    t_start = omp_get_wtime()

    call kernel_diverpih(trnopt, mpopt, mfcopt, dziv, ni, nj, nk, &
         j31, j32, jcb8u, jcb8v, mf, rmf, rmf8u, rmf8v, &
         rcsq, u, v, pdiv, tmp1, tmp2, tmp3)

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
  call validate_output(pdiv, pdiv_ref, ni, nj, nk, tolerance, &
       max_error, error_count, validation_passed)

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
  deallocate(j31, j32, jcb8u, jcb8v)
  deallocate(mf, rmf, rmf8u, rmf8v)
  deallocate(rcsq, u, v)
  deallocate(pdiv, pdiv_init, pdiv_ref)
  deallocate(tmp1, tmp2, tmp3)
  deallocate(times)

contains

  !-------------------------------------------------------------------
  ! Kernel: diverpih
  !-------------------------------------------------------------------
  subroutine kernel_diverpih(trnopt, mpopt, mfcopt, dziv, ni, nj, nk, &
       j31, j32, jcb8u, jcb8v, mf, rmf, rmf8u, rmf8v, &
       rcsq, u, v, pdiv, tmp1, tmp2, tmp3)
    implicit none

    integer, intent(in) :: trnopt, mpopt, mfcopt
    real, intent(in) :: dziv
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: mf(0:ni+1,0:nj+1)
    real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
    real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
    real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
    real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: pdiv(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)

    integer :: nkm1
    integer :: i, j, k

    nkm1 = nk - 1

    !$omp parallel default(shared) private(k)

    ! For the flat terrain case
    if (trnopt == 0) then
      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-2
            pdiv(i,j,k) = rcsq(i,j,k) * pdiv(i,j,k)
          end do
        end do
        !$omp end do
      end do

    ! For the curved grid case
    else
      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-1
            tmp1(i,j,k) = (u(i,j,k) + u(i,j,k-1)) * j31(i,j,k)
          end do
        end do
        !$omp end do

        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-1
          do i = 2, ni-2
            tmp2(i,j,k) = (v(i,j,k) + v(i,j,k-1)) * j32(i,j,k)
          end do
        end do
        !$omp end do
      end do

      if (mfcopt == 0) then
        do k = 2, nk-2
          !$omp do schedule(runtime) private(i,j)
          do j = 2, nj-2
            do i = 2, ni-2
              tmp3(i,j,k) = 0.25e0 * ((tmp1(i,j,k) + tmp1(i+1,j,k)) &
                   + (tmp2(i,j,k) + tmp2(i,j+1,k)))
            end do
          end do
          !$omp end do
        end do

      else
        if (mpopt == 0 .or. mpopt == 10) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                tmp3(i,j,k) = 0.25e0 * (mf(i,j) * (tmp1(i,j,k) + tmp1(i+1,j,k)) &
                     + (tmp2(i,j,k) + tmp2(i,j+1,k)))
              end do
            end do
            !$omp end do
          end do

        else if (mpopt == 5) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                tmp3(i,j,k) = 0.25e0 * ((tmp1(i,j,k) + tmp1(i+1,j,k)) &
                     + mf(i,j) * (tmp2(i,j,k) + tmp2(i,j+1,k)))
              end do
            end do
            !$omp end do
          end do

        else
          !$omp do schedule(runtime) private(i,j)
          do j = 2, nj-2
            do i = 2, ni-2
              pdiv(i,j,nk) = 0.25e0 * mf(i,j)
            end do
          end do
          !$omp end do

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                tmp3(i,j,k) = pdiv(i,j,nk) * ((tmp1(i,j,k) + tmp1(i+1,j,k)) &
                     + (tmp2(i,j,k) + tmp2(i,j+1,k)))
              end do
            end do
            !$omp end do
          end do
        end if
      end if

      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-2
        do i = 2, ni-2
          tmp3(i,j,nkm1) = 0.0e0
        end do
      end do
      !$omp end do

      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-2
            pdiv(i,j,k) = rcsq(i,j,k) &
                 * (pdiv(i,j,k) + (tmp3(i,j,k) - tmp3(i,j,k+1)) * dziv)
          end do
        end do
        !$omp end do
      end do
    end if

    !$omp end parallel

  end subroutine kernel_diverpih

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
  subroutine read_parameters(filename, trnopt, mpopt, mfcopt, dziv, ni, nj, nk)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: trnopt, mpopt, mfcopt, ni, nj, nk
    real, intent(out) :: dziv

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
          case ('trnopt')
            read(val, *) trnopt
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
          case ('dziv')
            read(val, *) dziv
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
  ! Read 2D array
  !-------------------------------------------------------------------
  subroutine read_array_2d(filename, arr, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: arr(is:ie, js:je)

    open(unit=10, file=filename, access='stream', form='unformatted', status='old')
    read(10) arr
    close(10)
  end subroutine read_array_2d

  !-------------------------------------------------------------------
  ! Read 3D array
  !-------------------------------------------------------------------
  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)

    open(unit=10, file=filename, access='stream', form='unformatted', status='old')
    read(10) arr
    close(10)
  end subroutine read_array_3d

  !-------------------------------------------------------------------
  ! Validate output
  !-------------------------------------------------------------------
  subroutine validate_output(output, reference, ni, nj, nk, tol, &
       max_err, err_count, passed)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: output(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: reference(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count
    logical, intent(out) :: passed

    real :: rel_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = 2, nk-2
      do j = 2, nj-2
        do i = 2, ni-2
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

    passed = (err_count == 0)
  end subroutine validate_output

end program kernel_benchmark_diverpih

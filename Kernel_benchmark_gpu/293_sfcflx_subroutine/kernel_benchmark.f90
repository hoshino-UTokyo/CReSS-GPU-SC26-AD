!-----------------------------------------------------------------------
! Kernel Benchmark Program: sfcflx (GPU Version)
! Source: Src/sfcflx.f90
! Subroutine: s_sfcflx
! Description: Calculate exchange coefficients for surface momentum,
!              heat and moisture fluxes from bulk coefficients
!-----------------------------------------------------------------------
program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters (dimensions)
  integer :: ni, nj, nk

  ! Physical constants
  real :: rddwkp

  ! Arrays
  real, allocatable :: za(:,:)
  real, allocatable :: rbr(:,:,:)
  integer, allocatable :: land(:,:)
  real, allocatable :: kai(:,:)
  real, allocatable :: z0m(:,:), z0h(:,:)
  real, allocatable :: va(:,:)
  real, allocatable :: rch(:,:)
  real, allocatable :: cm(:,:), ch(:,:)
  real, allocatable :: ce(:,:), ct(:,:), cq(:,:)
  real, allocatable :: cm_ref(:,:), ch_ref(:,:)
  real, allocatable :: ce_ref(:,:), ct_ref(:,:), cq_ref(:,:)

  ! Benchmark variables
  character(len=256) :: data_dir
  integer :: num_iterations, warmup_iterations
  real :: tolerance
  double precision :: start_time, end_time, total_time, avg_time
  integer :: iter, errors
  real :: max_error

  ! Read benchmark configuration
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! Read parameters from dump
  call read_parameters(trim(data_dir)//'/params.txt')

  ! Allocate arrays
  allocate(za(0:ni+1, 0:nj+1))
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(land(0:ni+1, 0:nj+1))
  allocate(kai(0:ni+1, 0:nj+1))
  allocate(z0m(0:ni+1, 0:nj+1))
  allocate(z0h(0:ni+1, 0:nj+1))
  allocate(va(0:ni+1, 0:nj+1))
  allocate(rch(0:ni+1, 0:nj+1))
  allocate(cm(0:ni+1, 0:nj+1))
  allocate(ch(0:ni+1, 0:nj+1))
  allocate(ce(0:ni+1, 0:nj+1))
  allocate(ct(0:ni+1, 0:nj+1))
  allocate(cq(0:ni+1, 0:nj+1))
  allocate(cm_ref(0:ni+1, 0:nj+1))
  allocate(ch_ref(0:ni+1, 0:nj+1))
  allocate(ce_ref(0:ni+1, 0:nj+1))
  allocate(ct_ref(0:ni+1, 0:nj+1))
  allocate(cq_ref(0:ni+1, 0:nj+1))

  ! Read input arrays
  call read_array_2d(trim(data_dir)//'/za.bin', za, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d_int(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/kai.bin', kai, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/z0m.bin', z0m, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/z0h.bin', z0h, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/va.bin', va, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/rch.bin', rch, 0, ni+1, 0, nj+1)

  ! Read cm and ch (output from bulksfc, input to sfcflx final calculation)
  ! These are stored as reference but we use them as input
  call read_array_2d(trim(data_dir)//'/cm_ref.bin', cm, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/ch_ref.bin', ch, 0, ni+1, 0, nj+1)

  ! Read reference output for validation (ce, ct, cq)
  call read_array_2d(trim(data_dir)//'/ce_ref.bin', ce_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/ct_ref.bin', ct_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/cq_ref.bin', cq_ref, 0, ni+1, 0, nj+1)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: sfcflx (GPU Version)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' rddwkp: ', rddwkp
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  ! Warmup iterations
  do iter = 1, warmup_iterations
    call kernel_sfcflx(ni, nj, nk, za, rbr, land, kai, z0m, z0h, va, rch, &
                       cm, ch, ce, ct, cq, rddwkp)
  end do
  !$acc wait

  ! Benchmark iterations
  total_time = 0.0d0
  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_sfcflx(ni, nj, nk, za, rbr, land, kai, z0m, z0h, va, rch, &
                       cm, ch, ce, ct, cq, rddwkp)
    !$acc wait
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results
  call validate_results(ce, ce_ref, ct, ct_ref, cq, cq_ref, &
                        ni, nj, tolerance, errors, max_error)

  ! Output results
  print '(A)',        '========================================'
  print '(A)',        'Kernel: sfcflx (s_sfcflx)'
  print '(A)',        '========================================'
  print '(A,I0,A,I0)', 'Grid size: ', ni, ' x ', nj
  print '(A,I0)',     'Iterations: ', num_iterations
  print '(A,F12.6,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', 'Total time: ', total_time, ' s'
  print '(A,I0)',     'Validation errors: ', errors
  print '(A,ES12.5)', 'Max error: ', max_error
  print '(A,ES12.5)', 'Tolerance: ', tolerance
  if (errors == 0) then
    print '(A)',      'Validation: PASSED'
  else
    print '(A)',      'Validation: FAILED'
    stop 1
  end if
  print '(A)',        '========================================'

  ! Cleanup
  deallocate(za, rbr, land, kai, z0m, z0h, va, rch)
  deallocate(cm, ch, ce, ct, cq)
  deallocate(cm_ref, ch_ref, ce_ref, ct_ref, cq_ref)

contains

  !---------------------------------------------------------------------
  ! Kernel subroutine: sfcflx (GPU Version)
  ! Note: This only benchmarks the final exchange coefficient calculation
  ! cm and ch are read as input (pre-computed by bulksfc)
  !---------------------------------------------------------------------
  subroutine kernel_sfcflx(ni, nj, nk, za, rbr, land, kai, z0m, z0h, va, rch, &
                           cm, ch, ce, ct, cq, rddwkp)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: za(0:ni+1, 0:nj+1)
    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    integer, intent(in) :: land(0:ni+1, 0:nj+1)
    real, intent(in) :: kai(0:ni+1, 0:nj+1)
    real, intent(in) :: z0m(0:ni+1, 0:nj+1)
    real, intent(in) :: z0h(0:ni+1, 0:nj+1)
    real, intent(in) :: va(0:ni+1, 0:nj+1)
    real, intent(in) :: rch(0:ni+1, 0:nj+1)
    real, intent(in) :: cm(0:ni+1, 0:nj+1)
    real, intent(in) :: ch(0:ni+1, 0:nj+1)
    real, intent(out) :: ce(0:ni+1, 0:nj+1)
    real, intent(out) :: ct(0:ni+1, 0:nj+1)
    real, intent(out) :: cq(0:ni+1, 0:nj+1)
    real, intent(in) :: rddwkp

    integer :: i, j
    real :: a

    ! Calculate exchange coefficients from bulk coefficients
    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni-1
        a = rbr(i,j,2) * cm(i,j) * va(i,j)

        ce(i,j) = a * cm(i,j)
        ct(i,j) = a * ch(i,j)

        if (land(i,j) .lt. 0) then
          cq(i,j) = ct(i,j) / (1.e0 + rddwkp * ch(i,j))
        else
          cq(i,j) = ct(i,j)
        end if
      end do
    end do
    !$acc end kernels

  end subroutine kernel_sfcflx

  !---------------------------------------------------------------------
  ! Read configuration file
  !---------------------------------------------------------------------
  subroutine read_config(data_dir, num_iterations, warmup_iterations, tolerance)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iterations, warmup_iterations
    real, intent(out) :: tolerance

    open(10, file='benchmark.conf', status='old')
    read(10, '(A)') data_dir
    read(10, *) num_iterations
    read(10, *) warmup_iterations
    read(10, *) tolerance
    close(10)
  end subroutine read_config

  !---------------------------------------------------------------------
  ! Read parameters from dump
  !---------------------------------------------------------------------
  subroutine read_parameters(filename)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=256) :: line
    character(len=64) :: name
    integer :: ios

    open(10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: Cannot open ', trim(filename)
      stop 1
    end if

    do while (.true.)
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle

      if (index(line, 'ni =') > 0) then
        read(line(index(line,'=')+1:), *) ni
      else if (index(line, 'nj =') > 0) then
        read(line(index(line,'=')+1:), *) nj
      else if (index(line, 'nk =') > 0) then
        read(line(index(line,'=')+1:), *) nk
      else if (index(line, 'rddwkp =') > 0) then
        read(line(index(line,'=')+1:), *) rddwkp
      end if
    end do

    close(10)
  end subroutine read_parameters

  !---------------------------------------------------------------------
  ! Read 2D real array
  !---------------------------------------------------------------------
  subroutine read_array_2d(filename, arr, i1, i2, j1, j2)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: arr(i1:i2, j1:j2)
    integer :: ios

    open(10, file=filename, access='stream', form='unformatted', status='old', &
         convert='big_endian', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: Cannot open ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_2d

  !---------------------------------------------------------------------
  ! Read 2D integer array
  !---------------------------------------------------------------------
  subroutine read_array_2d_int(filename, arr, i1, i2, j1, j2)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    integer, intent(out) :: arr(i1:i2, j1:j2)
    integer :: ios

    open(10, file=filename, access='stream', form='unformatted', status='old', &
         convert='big_endian', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: Cannot open ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_2d_int

  !---------------------------------------------------------------------
  ! Read 3D real array
  !---------------------------------------------------------------------
  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: ios

    open(10, file=filename, access='stream', form='unformatted', status='old', &
         convert='big_endian', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: Cannot open ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

  !---------------------------------------------------------------------
  ! Validate results
  !---------------------------------------------------------------------
  subroutine validate_results(ce, ce_ref, ct, ct_ref, cq, cq_ref, &
                              ni, nj, tolerance, errors, max_error)
    implicit none
    real, intent(in) :: ce(0:ni+1, 0:nj+1)
    real, intent(in) :: ce_ref(0:ni+1, 0:nj+1)
    real, intent(in) :: ct(0:ni+1, 0:nj+1)
    real, intent(in) :: ct_ref(0:ni+1, 0:nj+1)
    real, intent(in) :: cq(0:ni+1, 0:nj+1)
    real, intent(in) :: cq_ref(0:ni+1, 0:nj+1)
    integer, intent(in) :: ni, nj
    real, intent(in) :: tolerance
    integer, intent(out) :: errors
    real, intent(out) :: max_error

    integer :: i, j
    real :: rel_err

    errors = 0
    max_error = 0.0

    do j = 1, nj-1
      do i = 1, ni-1
        ! Check ce
        if (abs(ce_ref(i,j)) > 1.0e-30) then
          rel_err = abs(ce(i,j) - ce_ref(i,j)) / abs(ce_ref(i,j))
        else
          rel_err = abs(ce(i,j) - ce_ref(i,j))
        end if
        if (rel_err > max_error) max_error = rel_err
        if (rel_err > tolerance) errors = errors + 1

        ! Check ct
        if (abs(ct_ref(i,j)) > 1.0e-30) then
          rel_err = abs(ct(i,j) - ct_ref(i,j)) / abs(ct_ref(i,j))
        else
          rel_err = abs(ct(i,j) - ct_ref(i,j))
        end if
        if (rel_err > max_error) max_error = rel_err
        if (rel_err > tolerance) errors = errors + 1

        ! Check cq
        if (abs(cq_ref(i,j)) > 1.0e-30) then
          rel_err = abs(cq(i,j) - cq_ref(i,j)) / abs(cq_ref(i,j))
        else
          rel_err = abs(cq(i,j) - cq_ref(i,j))
        end if
        if (rel_err > max_error) max_error = rel_err
        if (rel_err > tolerance) errors = errors + 1
      end do
    end do
  end subroutine validate_results

end program kernel_benchmark

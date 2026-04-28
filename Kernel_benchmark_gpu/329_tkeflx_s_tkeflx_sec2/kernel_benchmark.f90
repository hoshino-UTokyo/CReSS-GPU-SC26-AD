!***********************************************************************
! GPU Kernel Benchmark: tkeflx sec2 (s_tkeflx main TKE turbulent flux)
!***********************************************************************
!
! Source: Src/tkeflx.f90
! Description: Calculate the turbulent fluxes for the turbulent kinetic
!              energy (x, y, z components). Includes conditional branches
!              based on mfcopt, mpopt, and trnopt.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_tkeflx_sec2
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Option flags
  integer :: trnopt, mfcopt, mpopt

  ! Scalar parameters
  real :: dxv5, dyv5, dzv5, dzv125

  ! Input arrays
  real, allocatable :: jcb(:,:,:)
  real, allocatable :: tke(:,:,:)
  real, allocatable :: rkh(:,:,:)
  real, allocatable :: rkv(:,:,:)
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: j31tke(:,:,:)
  real, allocatable :: j32tke(:,:,:)

  ! Output arrays
  real, allocatable :: h1(:,:,:)
  real, allocatable :: h2(:,:,:)
  real, allocatable :: h3(:,:,:)
  real, allocatable :: jcbtke(:,:,:)

  ! Reference outputs for validation
  real, allocatable :: h1_ref(:,:,:)
  real, allocatable :: h2_ref(:,:,:)
  real, allocatable :: h3_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error_h1, max_error_h2, max_error_h3, rel_error
  real :: tolerance
  integer :: error_count_h1, error_count_h2, error_count_h3
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
       ni, nj, nk, trnopt, mfcopt, mpopt, dxv5, dyv5, dzv5, dzv125)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: tkeflx_sec2'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' trnopt: ', trnopt
  write(*,'(A,I6)') ' mfcopt: ', mfcopt
  write(*,'(A,I6)') ' mpopt:  ', mpopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkh(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkv(0:ni+1, 0:nj+1, 1:nk))
  allocate(h1(0:ni+1, 0:nj+1, 1:nk))
  allocate(h2(0:ni+1, 0:nj+1, 1:nk))
  allocate(h3(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcbtke(0:ni+1, 0:nj+1, 1:nk))
  allocate(h1_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(h2_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(h3_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  ! Conditional allocations
  if (mfcopt == 1) then
    allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  end if
  if (trnopt >= 1) then
    allocate(j31tke(0:ni+1, 0:nj+1, 1:nk))
    allocate(j32tke(0:ni+1, 0:nj+1, 1:nk))
  end if

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tke.bin', tke, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rkh.bin', rkh, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rkv.bin', rkv, 0, ni+1, 0, nj+1, 1, nk)

  if (mfcopt == 1) then
    call read_array_3d(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  end if
  if (trnopt >= 1) then
    call read_array_3d(trim(data_dir)//'/j31tke.bin', j31tke, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/j32tke.bin', j32tke, 0, ni+1, 0, nj+1, 1, nk)
  end if

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/h1_ref.bin', h1_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/h2_ref.bin', h2_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/h3_ref.bin', h3_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    h1 = 0.0; h2 = 0.0; h3 = 0.0; jcbtke = 0.0
    call kernel_tkeflx_sec2(ni, nj, nk, trnopt, mfcopt, mpopt, &
         dxv5, dyv5, dzv5, dzv125, &
         jcb, tke, rkh, rkv, rmf, j31tke, j32tke, &
         h1, h2, h3, jcbtke)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    h1 = 0.0; h2 = 0.0; h3 = 0.0; jcbtke = 0.0

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_tkeflx_sec2(ni, nj, nk, trnopt, mfcopt, mpopt, &
         dxv5, dyv5, dzv5, dzv125, &
         jcb, tke, rkh, rkv, rmf, j31tke, j32tke, &
         h1, h2, h3, jcbtke)
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

  ! Validate h1
  max_error_h1 = 0.0; error_count_h1 = 0
  do k = 1, nk-1
    do j = 2, nj-2
      do i = 2, ni-1
        rel_error = abs(h1(i,j,k) - h1_ref(i,j,k))
        if (abs(h1_ref(i,j,k)) > 1.0e-10) rel_error = rel_error / abs(h1_ref(i,j,k))
        if (rel_error > max_error_h1) max_error_h1 = rel_error
        if (rel_error > tolerance) error_count_h1 = error_count_h1 + 1
      end do
    end do
  end do

  ! Validate h2
  max_error_h2 = 0.0; error_count_h2 = 0
  do k = 1, nk-1
    do j = 2, nj-1
      do i = 2, ni-2
        rel_error = abs(h2(i,j,k) - h2_ref(i,j,k))
        if (abs(h2_ref(i,j,k)) > 1.0e-10) rel_error = rel_error / abs(h2_ref(i,j,k))
        if (rel_error > max_error_h2) max_error_h2 = rel_error
        if (rel_error > tolerance) error_count_h2 = error_count_h2 + 1
      end do
    end do
  end do

  ! Validate h3
  max_error_h3 = 0.0; error_count_h3 = 0
  do k = 2, nk-1
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(h3(i,j,k) - h3_ref(i,j,k))
        if (abs(h3_ref(i,j,k)) > 1.0e-10) rel_error = rel_error / abs(h3_ref(i,j,k))
        if (rel_error > max_error_h3) max_error_h3 = rel_error
        if (rel_error > tolerance) error_count_h3 = error_count_h3 + 1
      end do
    end do
  end do

  validation_passed = (error_count_h1 == 0 .and. error_count_h2 == 0 &
                       .and. error_count_h3 == 0)

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
  write(*,'(A,ES12.4)') ' h1 max error:   ', max_error_h1
  write(*,'(A,I12)')    ' h1 error count: ', error_count_h1
  write(*,'(A,ES12.4)') ' h2 max error:   ', max_error_h2
  write(*,'(A,I12)')    ' h2 error count: ', error_count_h2
  write(*,'(A,ES12.4)') ' h3 max error:   ', max_error_h3
  write(*,'(A,I12)')    ' h3 error count: ', error_count_h3
  write(*,'(A,ES12.4)') ' Tolerance:      ', tolerance
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(jcb, tke, rkh, rkv, h1, h2, h3, jcbtke)
  deallocate(h1_ref, h2_ref, h3_ref, times)
  if (allocated(rmf)) deallocate(rmf)
  if (allocated(j31tke)) deallocate(j31tke)
  if (allocated(j32tke)) deallocate(j32tke)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: tkeflx_sec2 - GPU version using OpenACC
  !=====================================================================
  subroutine kernel_tkeflx_sec2(ni, nj, nk, trnopt, mfcopt, mpopt, &
       dxv5, dyv5, dzv5, dzv125, &
       jcb, tke, rkh, rkv, rmf, j31tke, j32tke, &
       h1, h2, h3, jcbtke)
    implicit none

    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: trnopt, mfcopt, mpopt
    real, intent(in) :: dxv5, dyv5, dzv5, dzv125
    real, intent(in) :: jcb(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tke(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rkh(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rkv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rmf(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: j31tke(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32tke(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: h1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: h2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: h3(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: jcbtke(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    ! Compute jcbtke = tke * jcb
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          jcbtke(i,j,k) = tke(i,j,k) * jcb(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Calculate the x components of the turbulent fluxes
    if (mfcopt == 1 .and. mpopt == 5) then

      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 2, nj-2
          !$acc loop independent
          do i = 1, ni-1
            h3(i,j,k) = rmf(i,j,2) * rkh(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      if (trnopt == 0) then

        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-1
              h1(i,j,k) = (h3(i-1,j,k) + h3(i,j,k)) &
                * (jcbtke(i,j,k) - jcbtke(i-1,j,k)) * dxv5
            end do
          end do
        end do
        !$acc end kernels

      else if (trnopt >= 1) then

        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-1
              h1(i,j,k) = (h3(i-1,j,k) + h3(i,j,k)) &
                * ((jcbtke(i,j,k) - jcbtke(i-1,j,k)) * dxv5 &
                + (j31tke(i,j,k+1) - j31tke(i,j,k)) * dzv125)
            end do
          end do
        end do
        !$acc end kernels

      end if

    else

      if (trnopt == 0) then

        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-1
              h1(i,j,k) = (rkh(i-1,j,k) + rkh(i,j,k)) &
                * (jcbtke(i,j,k) - jcbtke(i-1,j,k)) * dxv5
            end do
          end do
        end do
        !$acc end kernels

      else if (trnopt >= 1) then

        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-1
              h1(i,j,k) = (rkh(i-1,j,k) + rkh(i,j,k)) &
                * ((jcbtke(i,j,k) - jcbtke(i-1,j,k)) * dxv5 &
                + (j31tke(i,j,k+1) - j31tke(i,j,k)) * dzv125)
            end do
          end do
        end do
        !$acc end kernels

      end if

    end if

    ! Calculate the y components of the turbulent fluxes
    if (mfcopt == 1 .and. (mpopt == 0 .or. mpopt == 10)) then

      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 2, ni-2
            h3(i,j,k) = rmf(i,j,2) * rkh(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      if (trnopt == 0) then

        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 2, nj-1
            !$acc loop independent
            do i = 2, ni-2
              h2(i,j,k) = (h3(i,j-1,k) + h3(i,j,k)) &
                * (jcbtke(i,j,k) - jcbtke(i,j-1,k)) * dyv5
            end do
          end do
        end do
        !$acc end kernels

      else if (trnopt >= 1) then

        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 2, nj-1
            !$acc loop independent
            do i = 2, ni-2
              h2(i,j,k) = (h3(i,j-1,k) + h3(i,j,k)) &
                * ((jcbtke(i,j,k) - jcbtke(i,j-1,k)) * dyv5 &
                + (j32tke(i,j,k+1) - j32tke(i,j,k)) * dzv125)
            end do
          end do
        end do
        !$acc end kernels

      end if

    else

      if (trnopt == 0) then

        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 2, nj-1
            !$acc loop independent
            do i = 2, ni-2
              h2(i,j,k) = (rkh(i,j-1,k) + rkh(i,j,k)) &
                * (jcbtke(i,j,k) - jcbtke(i,j-1,k)) * dyv5
            end do
          end do
        end do
        !$acc end kernels

      else if (trnopt >= 1) then

        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 2, nj-1
            !$acc loop independent
            do i = 2, ni-2
              h2(i,j,k) = (rkh(i,j-1,k) + rkh(i,j,k)) &
                * ((jcbtke(i,j,k) - jcbtke(i,j-1,k)) * dyv5 &
                + (j32tke(i,j,k+1) - j32tke(i,j,k)) * dzv125)
            end do
          end do
        end do
        !$acc end kernels

      end if

    end if

    ! Get the z components of the turbulent fluxes
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          h3(i,j,k) = (rkv(i,j,k-1) + rkv(i,j,k)) &
            * (tke(i,j,k) - tke(i,j,k-1)) * dzv5
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_tkeflx_sec2

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios
    logical :: exists
    data_dir = './data'; num_iter = 10; warmup_iter = 2; tol = 1.0e-5
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
  ! Parameter reader (extended for scalar options)
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk, &
       trnopt, mfcopt, mpopt, dxv5, dyv5, dzv5, dzv125)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    integer, intent(out) :: trnopt, mfcopt, mpopt
    real, intent(out) :: dxv5, dyv5, dzv5, dzv125
    character(len=256) :: line, key, val
    integer :: ios, eq_pos
    ni = 1; nj = 1; nk = 1
    trnopt = 0; mfcopt = 0; mpopt = 0
    dxv5 = 0.0; dyv5 = 0.0; dzv5 = 0.0; dzv125 = 0.0
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
          case ('ni');     read(val, *) ni
          case ('nj');     read(val, *) nj
          case ('nk');     read(val, *) nk
          case ('trnopt'); read(val, *) trnopt
          case ('mfcopt'); read(val, *) mfcopt
          case ('mpopt');  read(val, *) mpopt
          case ('dxv5');   read(val, *) dxv5
          case ('dyv5');   read(val, *) dyv5
          case ('dzv5');   read(val, *) dzv5
          case ('dzv125'); read(val, *) dzv125
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

end program kernel_benchmark_tkeflx_sec2

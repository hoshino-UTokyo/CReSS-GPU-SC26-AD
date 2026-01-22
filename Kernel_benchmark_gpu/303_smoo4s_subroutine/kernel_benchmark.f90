!***********************************************************************
! Kernel Benchmark: smoo4s (s_smoo4s) - GPU OpenACC Version
!***********************************************************************
!
! Source: Src/smoo4s.f90
! Description: Perform 4th order numerical smoothing for scalar variable
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_smoo4s
  use omp_lib
  use openacc
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: smtopt
  integer :: iwest, ieast, jsouth, jnorth
  real :: smhcoe, smvcoe

  ! Input arrays
  real, allocatable :: rbr(:,:,:)
  real, allocatable :: s(:,:,:)

  ! Input/Output array
  real, allocatable :: sfrc(:,:,:)
  real, allocatable :: sfrc_init(:,:,:)

  ! Work arrays
  real, allocatable :: rbrs(:,:,:)
  real, allocatable :: rbrs2(:,:,:)
  real, allocatable :: tmp1(:,:,:)
  real, allocatable :: tmp2(:,:,:)
  real, allocatable :: tmp3(:,:,:)

  ! Reference output for validation
  real, allocatable :: sfrc_ref(:,:,:)

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
       ni, nj, nk, smtopt, iwest, ieast, jsouth, jnorth, smhcoe, smvcoe)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: smoo4s (GPU OpenACC)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' smtopt = ', smtopt
  write(*,'(A,I6,A,I6)') ' iwest=', iwest, ', ieast=', ieast
  write(*,'(A,I6,A,I6)') ' jsouth=', jsouth, ', jnorth=', jnorth
  write(*,'(A,ES12.4)') ' smhcoe = ', smhcoe
  write(*,'(A,ES12.4)') ' smvcoe = ', smvcoe
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' GPU device: ', acc_get_device_num(acc_device_default)
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(s(0:ni+1, 0:nj+1, 1:nk))
  allocate(sfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(sfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbrs(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbrs2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
  allocate(sfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s.bin', s, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/sfrc_in.bin', sfrc_init, 0, ni+1, 0, nj+1, 1, nk)

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
    rbrs = 0.0
    rbrs2 = 0.0
    tmp1 = 0.0
    tmp2 = 0.0
    tmp3 = 0.0
    call kernel_smoo4s(smtopt, iwest, ieast, jsouth, jnorth, &
         smhcoe, smvcoe, ni, nj, nk, rbr, s, sfrc, &
         rbrs, rbrs2, tmp1, tmp2, tmp3)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    sfrc = sfrc_init
    rbrs = 0.0
    rbrs2 = 0.0
    tmp1 = 0.0
    tmp2 = 0.0
    tmp3 = 0.0

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_smoo4s(smtopt, iwest, ieast, jsouth, jnorth, &
         smhcoe, smvcoe, ni, nj, nk, rbr, s, sfrc, &
         rbrs, rbrs2, tmp1, tmp2, tmp3)
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

  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(sfrc(i,j,k) - sfrc_ref(i,j,k))
        if (abs(sfrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(sfrc_ref(i,j,k))
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
  deallocate(rbr, s, sfrc, sfrc_init)
  deallocate(rbrs, rbrs2, tmp1, tmp2, tmp3)
  deallocate(sfrc_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: smoo4s
  ! Perform 4th order numerical smoothing for scalar variable
  !=====================================================================
  subroutine kernel_smoo4s(smtopt, iwest, ieast, jsouth, jnorth, &
       smhcoe, smvcoe, ni, nj, nk, rbr, s, sfrc, &
       rbrs, rbrs2, tmp1, tmp2, tmp3)
    implicit none

    integer, intent(in) :: smtopt
    integer, intent(in) :: iwest, ieast, jsouth, jnorth
    real, intent(in) :: smhcoe, smvcoe
    integer, intent(in) :: ni, nj, nk

    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: s(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: sfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: rbrs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: rbrs2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp3(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    integer :: nkm1, nkm2

    nkm1 = nk - 1
    nkm2 = nk - 2

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-1
      do j = jsouth, nj-jnorth
        do i = iwest, ni-ieast
          rbrs(i,j,k) = rbr(i,j,k) * s(i,j,k)
          rbrs2(i,j,k) = 2.0 * rbrs(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 1+iwest, ni-1-ieast
          tmp1(i,j,k) = (rbrs(i+1,j,k) + rbrs(i-1,j,k)) - rbrs2(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 1+jsouth, nj-1-jnorth
        do i = 2, ni-2
          tmp2(i,j,k) = (rbrs(i,j+1,k) + rbrs(i,j-1,k)) - rbrs2(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 2, ni-2
          tmp3(i,j,k) = (rbrs(i,j,k+1) + rbrs(i,j,k-1)) - rbrs2(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 2, ni-2
          sfrc(i,j,k) = sfrc(i,j,k) &
               + smhcoe * (tmp1(i,j,k) + tmp2(i,j,k)) + smvcoe * tmp3(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    if (mod(smtopt, 10) == 2) then

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-2-jnorth
          do i = 2+iwest, ni-2-ieast
            sfrc(i,j,k) = sfrc(i,j,k) &
                 + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
                 + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k))))
          end do
        end do
      end do
      !$acc end kernels

    else

      !$acc kernels
      !$acc loop independent collapse(2)
      do j = 2+jsouth, nj-2-jnorth
        do i = 2+iwest, ni-2-ieast
          tmp3(i,j,1) = tmp3(i,j,2)
          tmp3(i,j,nkm1) = tmp3(i,j,nkm2)
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-2-jnorth
          do i = 2+iwest, ni-2-ieast
            sfrc(i,j,k) = sfrc(i,j,k) &
                 + (smvcoe * (tmp3(i,j,k) - (tmp3(i,j,k+1) + tmp3(i,j,k-1))) &
                 + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
                 + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k)))))
          end do
        end do
      end do
      !$acc end kernels

    end if

  end subroutine kernel_smoo4s

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
  subroutine read_parameters(filename, ni, nj, nk, smtopt, &
       iwest, ieast, jsouth, jnorth, smhcoe, smvcoe)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, smtopt
    integer, intent(out) :: iwest, ieast, jsouth, jnorth
    real, intent(out) :: smhcoe, smvcoe

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
          case ('smtopt')
            read(val, *) smtopt
          case ('iwest')
            read(val, *) iwest
          case ('ieast')
            read(val, *) ieast
          case ('jsouth')
            read(val, *) jsouth
          case ('jnorth')
            read(val, *) jnorth
          case ('smhcoe')
            read(val, *) smhcoe
          case ('smvcoe')
            read(val, *) smvcoe
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

end program kernel_benchmark_smoo4s

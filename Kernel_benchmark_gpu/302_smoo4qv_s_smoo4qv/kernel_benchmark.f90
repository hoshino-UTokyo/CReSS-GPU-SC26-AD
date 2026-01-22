!***********************************************************************
! GPU Kernel Benchmark: smoo4qv (s_smoo4qv)
!***********************************************************************
!
! Source: Src/smoo4qv.f90
! Description: 4th order smoothing for water vapor mixing ratio
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_smoo4qv
  use omp_lib
  implicit none

  ! Grid dimensions
  integer :: ni, nj, nk, nkm1, nkm2

  ! Control parameters
  integer :: smtopt, iwest, ieast, jsouth, jnorth
  real :: smhcoe, smvcoe

  ! Input arrays
  real, allocatable :: qv(:,:,:)      ! Water vapor mixing ratio
  real, allocatable :: qvbr(:,:,:)    ! Base state qv
  real, allocatable :: rbr(:,:,:)     ! Base state density

  ! Input/output arrays
  real, allocatable :: qvfrc(:,:,:)   ! Forcing term
  real, allocatable :: rbrqv(:,:,:)   ! rbr x (qv - qvbr)
  real, allocatable :: rbrqv2(:,:,:)  ! 2.0 x rbrqv
  real, allocatable :: tmp1(:,:,:)    ! Work array
  real, allocatable :: tmp2(:,:,:)    ! Work array
  real, allocatable :: tmp3(:,:,:)    ! Work array

  ! Reference arrays
  real, allocatable :: qvfrc_ref(:,:,:)
  real, allocatable :: rbrqv_ref(:,:,:)
  real, allocatable :: rbrqv2_ref(:,:,:)
  real, allocatable :: tmp1_ref(:,:,:)
  real, allocatable :: tmp2_ref(:,:,:)
  real, allocatable :: tmp3_ref(:,:,:)

  ! Initial state arrays (for reset)
  real, allocatable :: qvfrc_init(:,:,:)
  real, allocatable :: rbrqv_init(:,:,:)
  real, allocatable :: rbrqv2_init(:,:,:)
  real, allocatable :: tmp1_init(:,:,:)
  real, allocatable :: tmp2_init(:,:,:)
  real, allocatable :: tmp3_init(:,:,:)

  ! Benchmark parameters
  character(len=256) :: data_dir
  integer :: warmup_iterations, benchmark_iterations
  real :: tolerance

  ! Timing variables
  real(8), allocatable :: times(:)
  real(8) :: start_time, end_time
  real(8) :: avg_time, min_time, max_time, total_time

  ! Validation variables
  real :: max_rel_error
  integer :: error_count, total_elements

  ! Loop variables
  integer :: i, j, k, iter

  ! Read benchmark configuration
  call read_config(data_dir, warmup_iterations, benchmark_iterations, tolerance)

  ! Read parameters
  call read_params(data_dir, ni, nj, nk, smtopt, iwest, ieast, jsouth, jnorth, &
                   smhcoe, smvcoe)

  ! Derived values
  nkm1 = nk - 1
  nkm2 = nk - 2

  ! Print benchmark info
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: smoo4qv'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I3,A,I3,A,I3,A,I3)') ' iwest=', iwest, ', ieast=', ieast, &
        ', jsouth=', jsouth, ', jnorth=', jnorth
  write(*,'(A,I3)') ' smtopt=', smtopt
  write(*,'(A,E12.4,A,E12.4)') ' smhcoe=', smhcoe, ', smvcoe=', smvcoe
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', benchmark_iterations
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(qv(0:ni+1,0:nj+1,1:nk))
  allocate(qvbr(0:ni+1,0:nj+1,1:nk))
  allocate(rbr(0:ni+1,0:nj+1,1:nk))
  allocate(qvfrc(0:ni+1,0:nj+1,1:nk))
  allocate(rbrqv(0:ni+1,0:nj+1,1:nk))
  allocate(rbrqv2(0:ni+1,0:nj+1,1:nk))
  allocate(tmp1(0:ni+1,0:nj+1,1:nk))
  allocate(tmp2(0:ni+1,0:nj+1,1:nk))
  allocate(tmp3(0:ni+1,0:nj+1,1:nk))
  allocate(qvfrc_ref(0:ni+1,0:nj+1,1:nk))
  allocate(rbrqv_ref(0:ni+1,0:nj+1,1:nk))
  allocate(rbrqv2_ref(0:ni+1,0:nj+1,1:nk))
  allocate(tmp1_ref(0:ni+1,0:nj+1,1:nk))
  allocate(tmp2_ref(0:ni+1,0:nj+1,1:nk))
  allocate(tmp3_ref(0:ni+1,0:nj+1,1:nk))
  allocate(qvfrc_init(0:ni+1,0:nj+1,1:nk))
  allocate(rbrqv_init(0:ni+1,0:nj+1,1:nk))
  allocate(rbrqv2_init(0:ni+1,0:nj+1,1:nk))
  allocate(tmp1_init(0:ni+1,0:nj+1,1:nk))
  allocate(tmp2_init(0:ni+1,0:nj+1,1:nk))
  allocate(tmp3_init(0:ni+1,0:nj+1,1:nk))
  allocate(times(benchmark_iterations))

  ! Load input data
  write(*,'(A)') ' Loading input data...'
  call read_3d_array(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/qvbr.bin', qvbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/qvfrc_in.bin', qvfrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/rbrqv_in.bin', rbrqv_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/rbrqv2_in.bin', rbrqv2_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/tmp1_in.bin', tmp1_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/tmp2_in.bin', tmp2_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/tmp3_in.bin', tmp3_init, 0, ni+1, 0, nj+1, 1, nk)

  ! Load reference output
  write(*,'(A)') ' Loading reference output...'
  call read_3d_array(trim(data_dir)//'/qvfrc_ref.bin', qvfrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/rbrqv_ref.bin', rbrqv_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/rbrqv2_ref.bin', rbrqv2_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/tmp1_ref.bin', tmp1_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/tmp2_ref.bin', tmp2_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/tmp3_ref.bin', tmp3_ref, 0, ni+1, 0, nj+1, 1, nk)

  ! Warmup iterations
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call reset_arrays()
    call kernel_smoo4qv()
    !$acc wait
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  do iter = 1, benchmark_iterations
    call reset_arrays()
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_smoo4qv()
    !$acc wait
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
  end do

  ! Validate results
  write(*,'(A)') ' Validating output...'
  max_rel_error = 0.0
  error_count = 0
  total_elements = 0

  call validate_3d_array(qvfrc, qvfrc_ref, 0, ni+1, 0, nj+1, 1, nk, &
                         tolerance, max_rel_error, error_count, total_elements)
  call validate_3d_array(rbrqv, rbrqv_ref, 0, ni+1, 0, nj+1, 1, nk, &
                         tolerance, max_rel_error, error_count, total_elements)
  call validate_3d_array(rbrqv2, rbrqv2_ref, 0, ni+1, 0, nj+1, 1, nk, &
                         tolerance, max_rel_error, error_count, total_elements)
  call validate_3d_array(tmp1, tmp1_ref, 0, ni+1, 0, nj+1, 1, nk, &
                         tolerance, max_rel_error, error_count, total_elements)
  call validate_3d_array(tmp2, tmp2_ref, 0, ni+1, 0, nj+1, 1, nk, &
                         tolerance, max_rel_error, error_count, total_elements)
  call validate_3d_array(tmp3, tmp3_ref, 0, ni+1, 0, nj+1, 1, nk, &
                         tolerance, max_rel_error, error_count, total_elements)

  ! Calculate timing statistics
  total_time = sum(times) * 1000.0d0  ! Convert to ms
  avg_time = total_time / benchmark_iterations
  min_time = minval(times) * 1000.0d0
  max_time = maxval(times) * 1000.0d0

  ! Print results
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:    ', avg_time, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:        ', min_time, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:        ', max_time, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:      ', total_time, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,E12.4)') ' Max relative error: ', max_rel_error
  write(*,'(A,E12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', error_count
  if (error_count == 0) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! Cleanup
  deallocate(qv, qvbr, rbr, qvfrc, rbrqv, rbrqv2, tmp1, tmp2, tmp3)
  deallocate(qvfrc_ref, rbrqv_ref, rbrqv2_ref, tmp1_ref, tmp2_ref, tmp3_ref)
  deallocate(qvfrc_init, rbrqv_init, rbrqv2_init, tmp1_init, tmp2_init, tmp3_init)
  deallocate(times)

  if (error_count > 0) stop 1

contains

  subroutine reset_arrays()
    qvfrc = qvfrc_init
    rbrqv = rbrqv_init
    rbrqv2 = rbrqv2_init
    tmp1 = tmp1_init
    tmp2 = tmp2_init
    tmp3 = tmp3_init
  end subroutine reset_arrays

  subroutine kernel_smoo4qv()
    integer :: i, j, k
    integer :: nkm1_local, nkm2_local

    nkm1_local = nk - 1
    nkm2_local = nk - 2

    ! Compute rbrqv and rbrqv2
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = jsouth, nj-jnorth
        !$acc loop independent
        do i = iwest, ni-ieast
          rbrqv(i,j,k) = rbr(i,j,k) * (qv(i,j,k) - qvbr(i,j,k))
          rbrqv2(i,j,k) = 2.e0 * rbrqv(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    ! Compute tmp1, tmp2, tmp3
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 1+iwest, ni-1-ieast
          tmp1(i,j,k) = (rbrqv(i+1,j,k) + rbrqv(i-1,j,k)) - rbrqv2(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 1+jsouth, nj-1-jnorth
        !$acc loop independent
        do i = 2, ni-2
          tmp2(i,j,k) = (rbrqv(i,j+1,k) + rbrqv(i,j-1,k)) - rbrqv2(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          tmp3(i,j,k) = (rbrqv(i,j,k+1) + rbrqv(i,j,k-1)) - rbrqv2(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    ! Apply 2nd order smoothing to qvfrc
    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          qvfrc(i,j,k) = qvfrc(i,j,k) &
            + smhcoe * (tmp1(i,j,k) + tmp2(i,j,k)) + smvcoe * tmp3(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    ! Apply 4th order smoothing
    if (mod(smtopt,10) == 2) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-2
        !$acc loop independent
        do j = 2+jsouth, nj-2-jnorth
          !$acc loop independent
          do i = 2+iwest, ni-2-ieast
            qvfrc(i,j,k) = qvfrc(i,j,k) &
              + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
              + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k))))
          end do
        end do
      end do
      !$acc end kernels
    else
      !$acc kernels
      !$acc loop independent
      do j = 2+jsouth, nj-2-jnorth
        !$acc loop independent
        do i = 2+iwest, ni-2-ieast
          tmp3(i,j,1) = tmp3(i,j,2)
          tmp3(i,j,nkm1_local) = tmp3(i,j,nkm2_local)
        end do
      end do
      !$acc end kernels
      !$acc wait

      !$acc kernels
      !$acc loop independent
      do k = 2, nk-2
        !$acc loop independent
        do j = 2+jsouth, nj-2-jnorth
          !$acc loop independent
          do i = 2+iwest, ni-2-ieast
            qvfrc(i,j,k) = qvfrc(i,j,k) &
              + (smvcoe * (tmp3(i,j,k) - (tmp3(i,j,k+1) + tmp3(i,j,k-1))) &
              + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
              + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k)))))
          end do
        end do
      end do
      !$acc end kernels
    end if

  end subroutine kernel_smoo4qv

  subroutine read_config(data_dir, warmup_iters, bench_iters, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: warmup_iters, bench_iters
    real, intent(out) :: tol
    integer :: unit_num

    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read')
    read(unit_num, '(A)') data_dir
    read(unit_num, *) bench_iters
    read(unit_num, *) warmup_iters
    read(unit_num, *) tol
    close(unit_num)
  end subroutine read_config

  subroutine read_params(data_dir, ni, nj, nk, smtopt, iwest, ieast, jsouth, jnorth, &
                         smhcoe, smvcoe)
    character(len=*), intent(in) :: data_dir
    integer, intent(out) :: ni, nj, nk, smtopt, iwest, ieast, jsouth, jnorth
    real, intent(out) :: smhcoe, smvcoe
    character(len=256) :: line, key
    integer :: unit_num, ios, eq_pos
    real :: val

    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', action='read')
    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        read(line(eq_pos+1:), *) val
        select case(trim(key))
          case('ni'); ni = int(val)
          case('nj'); nj = int(val)
          case('nk'); nk = int(val)
          case('smtopt'); smtopt = int(val)
          case('iwest'); iwest = int(val)
          case('ieast'); ieast = int(val)
          case('jsouth'); jsouth = int(val)
          case('jnorth'); jnorth = int(val)
          case('smhcoe'); smhcoe = val
          case('smvcoe'); smvcoe = val
        end select
      end if
    end do
    close(unit_num)
  end subroutine read_params

  subroutine read_3d_array(filename, arr, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: unit_num

    unit_num = 12
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) arr
    close(unit_num)
  end subroutine read_3d_array

  subroutine validate_3d_array(computed, reference, is, ie, js, je, ks, ke, &
                               tol, max_err, err_count, total_count)
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(in) :: computed(is:ie, js:je, ks:ke)
    real, intent(in) :: reference(is:ie, js:je, ks:ke)
    real, intent(in) :: tol
    real, intent(inout) :: max_err
    integer, intent(inout) :: err_count, total_count
    real :: rel_err, abs_ref, abs_err
    real :: abs_tol
    integer :: i, j, k

    ! Absolute tolerance for very small values (single precision noise floor)
    ! Values smaller than 1e-8 are essentially at the noise floor for single precision
    abs_tol = 1.0e-8

    do k = ks, ke
      do j = js, je
        do i = is, ie
          total_count = total_count + 1
          abs_ref = abs(reference(i,j,k))
          abs_err = abs(computed(i,j,k) - reference(i,j,k))

          ! Use combined tolerance: pass if either abs or rel error is small
          if (abs_ref > abs_tol) then
            rel_err = abs_err / abs_ref
          else
            ! For very small values, use absolute error comparison
            if (abs_err < abs_tol) then
              rel_err = 0.0
            else
              rel_err = abs_err
            end if
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do
  end subroutine validate_3d_array

end program kernel_benchmark_gpu_smoo4qv

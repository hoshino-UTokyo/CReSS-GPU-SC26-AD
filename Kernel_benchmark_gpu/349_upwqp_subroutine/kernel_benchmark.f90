!> Kernel benchmark program for upwqp (s_upwqp) (GPU Version)
!> Calculates sedimentation flux and precipitation for optional
!> precipitation mixing ratio using upwind scheme.
program kernel_benchmark_upwqp
  use openacc
  implicit none

  ! Parameters from dump
  integer :: advopt
  real :: dziv
  integer :: ni, nj, nk
  real :: dtp

  ! Arrays
  real, allocatable :: rbr(:,:,:)      ! Base state density
  real, allocatable :: rst(:,:,:)      ! Base state density x Jacobian
  real, allocatable :: uq(:,:,:)       ! Terminal velocity
  real, allocatable :: qpf(:,:,:)      ! Precipitation mixing ratio
  real, allocatable :: qpf_save(:,:,:) ! Save for reset
  real, allocatable :: precip(:,:,:)   ! Precipitation and accumulation
  real, allocatable :: precip_save(:,:,:)
  real, allocatable :: qpflx(:,:,:)    ! Fallout flux
  real, allocatable :: qpflx_save(:,:,:)
  real, allocatable :: qpf_ref(:,:,:)    ! Reference output
  real, allocatable :: precip_ref(:,:,:)
  real, allocatable :: qpflx_ref(:,:,:)

  ! Benchmark variables
  character(len=256) :: data_dir
  integer :: num_iterations, warmup_iterations
  real :: tolerance
  integer :: iter
  real(8) :: start_time, end_time, elapsed_time
  real(8) :: total_time, min_time, max_time
  real :: max_diff
  integer :: count_start, count_end, count_rate, count_max

  ! Physical constant
  real, parameter :: rhow = 1.0e3  ! Water density

  ! Read configuration
  call read_config()

  ! Read parameters
  call read_parameters()

  ! Allocate arrays
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(uq(0:ni+1, 0:nj+1, 1:nk))
  allocate(qpf(0:ni+1, 0:nj+1, 1:nk))
  allocate(qpf_save(0:ni+1, 0:nj+1, 1:nk))
  allocate(precip(0:ni+1, 0:nj+1, 1:2))
  allocate(precip_save(0:ni+1, 0:nj+1, 1:2))
  allocate(qpflx(0:ni+1, 0:nj+1, 1:nk))
  allocate(qpflx_save(0:ni+1, 0:nj+1, 1:nk))
  allocate(qpf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(precip_ref(0:ni+1, 0:nj+1, 1:2))
  allocate(qpflx_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Read input data
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/uq.bin', uq, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qpf_in.bin', qpf_save, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d_2(trim(data_dir)//'/precip_in.bin', precip_save, 0, ni+1, 0, nj+1, 1, 2)
  call read_array_3d(trim(data_dir)//'/qpflx_in.bin', qpflx_save, 0, ni+1, 0, nj+1, 1, nk)

  ! Read reference output
  call read_array_3d(trim(data_dir)//'/qpf_ref.bin', qpf_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d_2(trim(data_dir)//'/precip_ref.bin', precip_ref, 0, ni+1, 0, nj+1, 1, 2)
  call read_array_3d(trim(data_dir)//'/qpflx_ref.bin', qpflx_ref, 0, ni+1, 0, nj+1, 1, nk)

  print '(A)', '================================================'
  print '(A)', 'Kernel Benchmark: upwqp (s_upwqp) (GPU)'
  print '(A)', 'Sedimentation flux and precipitation calculation'
  print '(A)', '================================================'
  print '(A,I0,A,I0,A,I0)', 'Grid size: ', ni, ' x ', nj, ' x ', nk
  print '(A,I0)', 'advopt: ', advopt
  print '(A,E12.5)', 'dziv: ', dziv
  print '(A,E12.5)', 'dtp: ', dtp
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,I0)', 'Warmup iterations: ', warmup_iterations
  print '(A,I0)', 'GPU device: ', acc_get_device_num(acc_device_default)
  print '(A)', '------------------------------------------------'

  ! Warmup iterations
  do iter = 1, warmup_iterations
    qpf = qpf_save
    precip = precip_save
    qpflx = qpflx_save
    call kernel_upwqp(advopt, dziv, dtp, ni, nj, nk, rbr, rst, uq, &
                      qpf, precip, qpflx, rhow)
  end do

  ! Benchmark iterations
  total_time = 0.0d0
  min_time = 1.0d30
  max_time = 0.0d0

  call system_clock(count_rate=count_rate, count_max=count_max)

  do iter = 1, num_iterations
    qpf = qpf_save
    precip = precip_save
    qpflx = qpflx_save

    !$acc wait
    call system_clock(count_start)
    call kernel_upwqp(advopt, dziv, dtp, ni, nj, nk, rbr, rst, uq, &
                      qpf, precip, qpflx, rhow)
    !$acc wait
    call system_clock(count_end)

    start_time = dble(count_start) / dble(count_rate)
    end_time = dble(count_end) / dble(count_rate)
    elapsed_time = end_time - start_time
    total_time = total_time + elapsed_time
    min_time = min(min_time, elapsed_time)
    max_time = max(max_time, elapsed_time)
  end do

  ! Validate output
  call validate_output(qpf, qpf_ref, ni, nj, nk, tolerance, 'qpf', max_diff)
  call validate_output_2(precip, precip_ref, ni, nj, tolerance, 'precip', max_diff)
  call validate_output(qpflx, qpflx_ref, ni, nj, nk, tolerance, 'qpflx', max_diff)

  ! Report results
  print '(A)', ''
  print '(A)', 'Timing Results:'
  print '(A,F12.6,A)', '  Total time:   ', total_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Average time: ', (total_time / num_iterations) * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Min time:     ', min_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Max time:     ', max_time * 1000.0d0, ' ms'
  print '(A)', '================================================'

  ! Cleanup
  deallocate(rbr, rst, uq, qpf, qpf_save, precip, precip_save)
  deallocate(qpflx, qpflx_save, qpf_ref, precip_ref, qpflx_ref)

contains

  !> Main kernel subroutine (OpenACC version)
  subroutine kernel_upwqp(advopt, dziv, dtp, ni, nj, nk, rbr, rst, uq, &
                          qpf, precip, qpflx, rhow)
    implicit none
    integer, intent(in) :: advopt
    real, intent(in) :: dziv, dtp
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: uq(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: qpf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: precip(0:ni+1, 0:nj+1, 1:2)
    real, intent(inout) :: qpflx(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rhow

    integer :: i, j, k
    integer :: nkm1, nkm2
    real :: dtp05, dzvdt, rwiv05

    nkm1 = nk - 1
    nkm2 = nk - 2
    dtp05 = 0.5e0 * dtp
    dzvdt = dziv * dtp
    rwiv05 = 0.5e0 / rhow

    ! Calculate flux
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          qpflx(i,j,k) = rbr(i,j,k) * uq(i,j,k) * qpf(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Update mixing ratio
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-2
      do j = 1, nj-1
        do i = 1, ni-1
          qpf(i,j,k) = max(qpf(i,j,k) + &
                       (qpflx(i,j,k+1) - qpflx(i,j,k)) / rst(i,j,k) * dzvdt, 0.e0)
        end do
      end do
    end do
    !$acc end kernels

    if (advopt <= 3) then
      !$acc kernels
      !$acc loop independent collapse(2)
      do j = 1, nj-1
        do i = 1, ni-1
          qpf(i,j,nkm1) = qpf(i,j,nkm2)
          precip(i,j,1) = (qpflx(i,j,1) + qpflx(i,j,2)) * rwiv05
          precip(i,j,2) = precip(i,j,2) + precip(i,j,1) * dtp05
        end do
      end do
      !$acc end kernels
    else
      !$acc kernels
      !$acc loop independent collapse(2)
      do j = 1, nj-1
        do i = 1, ni-1
          qpf(i,j,nkm1) = qpf(i,j,nkm2)
          precip(i,j,1) = (qpflx(i,j,1) + qpflx(i,j,2)) * rwiv05
          precip(i,j,2) = precip(i,j,2) + precip(i,j,1) * dtp
        end do
      end do
      !$acc end kernels
    end if

  end subroutine kernel_upwqp

  subroutine read_config()
    integer :: unit_num, ios
    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open benchmark.conf'
      stop 1
    end if
    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iterations
    read(unit_num, *) warmup_iterations
    read(unit_num, *) tolerance
    close(unit_num)
  end subroutine read_config

  subroutine read_parameters()
    integer :: unit_num, ios, eq_pos
    character(len=256) :: line
    character(len=64) :: param_name
    character(len=64) :: param_value

    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open params.txt in ', trim(data_dir)
      stop 1
    end if

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      ! Parse "name = value" format
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        param_name = adjustl(line(1:eq_pos-1))
        param_value = adjustl(line(eq_pos+1:))
      else
        read(line, *) param_name, param_value
      end if

      select case (trim(param_name))
        case ('advopt')
          read(param_value, *) advopt
        case ('dziv')
          read(param_value, *) dziv
        case ('ni')
          read(param_value, *) ni
        case ('nj')
          read(param_value, *) nj
        case ('nk')
          read(param_value, *) nk
        case ('dtp')
          read(param_value, *) dtp
      end select
    end do
    close(unit_num)
  end subroutine read_parameters

  subroutine read_array_3d(filename, array, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: array(is:ie, js:je, ks:ke)
    integer :: unit_num, ios

    unit_num = 20
    open(unit=unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if
    read(unit_num) array
    close(unit_num)
  end subroutine read_array_3d

  subroutine read_array_3d_2(filename, array, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: array(is:ie, js:je, ks:ke)
    integer :: unit_num, ios

    unit_num = 20
    open(unit=unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if
    read(unit_num) array
    close(unit_num)
  end subroutine read_array_3d_2

  subroutine validate_output(computed, reference, ni, nj, nk, tolerance, name, max_diff)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: computed(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: reference(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tolerance
    character(len=*), intent(in) :: name
    real, intent(out) :: max_diff
    integer :: i, j, k
    real :: diff, ref_val

    max_diff = 0.0
    ! Only validate the computed region (1:ni-1, 1:nj-1, 1:nk-1)
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          diff = abs(computed(i,j,k) - reference(i,j,k))
          ref_val = abs(reference(i,j,k))
          if (ref_val > 1.0e-30) then
            diff = diff / ref_val
          end if
          max_diff = max(max_diff, diff)
        end do
      end do
    end do

    if (max_diff <= tolerance) then
      print '(A,A,A,E12.5)', 'Validation PASSED for ', trim(name), ', max relative diff: ', max_diff
    else
      print '(A,A,A,E12.5)', 'Validation FAILED for ', trim(name), ', max relative diff: ', max_diff
    end if
  end subroutine validate_output

  subroutine validate_output_2(computed, reference, ni, nj, tolerance, name, max_diff)
    integer, intent(in) :: ni, nj
    real, intent(in) :: computed(0:ni+1, 0:nj+1, 1:2)
    real, intent(in) :: reference(0:ni+1, 0:nj+1, 1:2)
    real, intent(in) :: tolerance
    character(len=*), intent(in) :: name
    real, intent(out) :: max_diff
    integer :: i, j, k
    real :: diff, ref_val

    max_diff = 0.0
    ! Only validate the computed region (1:ni-1, 1:nj-1)
    do k = 1, 2
      do j = 1, nj-1
        do i = 1, ni-1
          diff = abs(computed(i,j,k) - reference(i,j,k))
          ref_val = abs(reference(i,j,k))
          if (ref_val > 1.0e-30) then
            diff = diff / ref_val
          end if
          max_diff = max(max_diff, diff)
        end do
      end do
    end do

    if (max_diff <= tolerance) then
      print '(A,A,A,E12.5)', 'Validation PASSED for ', trim(name), ', max relative diff: ', max_diff
    else
      print '(A,A,A,E12.5)', 'Validation FAILED for ', trim(name), ', max relative diff: ', max_diff
    end if
  end subroutine validate_output_2

end program kernel_benchmark_upwqp

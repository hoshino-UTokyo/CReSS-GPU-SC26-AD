!> Kernel benchmark program for freezing (s_freezing) - GPU version
!> Calculates the freezing rate from rain water to graupel.
program kernel_benchmark_freezing
  use omp_lib
  implicit none

  ! Parameters from dump
  integer :: cphopt
  integer :: ni, nj, nk
  real :: dtb, thresq
  real :: t0, t0cel, tlow, rhow, cc

  ! Arrays
  real, allocatable :: qr(:,:,:)       ! Rain water mixing ratio
  real, allocatable :: ncr(:,:,:)      ! Concentrations of rain water
  real, allocatable :: tcel(:,:,:)     ! Ambient air temperature (Celsius)
  real, allocatable :: diaqr(:,:,:)    ! Mean diameter of rain water
  real, allocatable :: frrg(:,:,:)     ! Freezing rate (output)
  real, allocatable :: frrgn(:,:,:)    ! Freezing rate for concentrations (output)
  real, allocatable :: frrg_ref(:,:,:) ! Reference output
  real, allocatable :: frrgn_ref(:,:,:)

  ! Benchmark variables
  character(len=256) :: data_dir
  integer :: num_iterations, warmup_iterations
  real :: tolerance
  integer :: iter
  real(8) :: start_time, end_time, elapsed_time
  real(8) :: total_time, min_time, max_time
  real :: max_diff

  ! Read configuration
  call read_config()

  ! Read parameters
  call read_parameters()

  ! Allocate arrays
  allocate(qr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncr(0:ni+1, 0:nj+1, 1:nk))
  allocate(tcel(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqr(0:ni+1, 0:nj+1, 1:nk))
  allocate(frrg(0:ni+1, 0:nj+1, 1:nk))
  allocate(frrgn(0:ni+1, 0:nj+1, 1:nk))
  allocate(frrg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(frrgn_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Read input data
  call read_array_3d(trim(data_dir)//'/qr.bin', qr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncr.bin', ncr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tcel.bin', tcel, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/diaqr.bin', diaqr, 0, ni+1, 0, nj+1, 1, nk)

  ! Read reference output
  call read_array_3d(trim(data_dir)//'/frrg_ref.bin', frrg_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/frrgn_ref.bin', frrgn_ref, 0, ni+1, 0, nj+1, 1, nk)

  print '(A)', '================================================'
  print '(A)', 'Kernel Benchmark: freezing (s_freezing) - GPU'
  print '(A)', 'Rain water to graupel freezing rate calculation'
  print '(A)', '================================================'
  print '(A,I0,A,I0,A,I0)', 'Grid size: ', ni, ' x ', nj, ' x ', nk
  print '(A,I0)', 'cphopt: ', cphopt
  print '(A,E12.5)', 'dtb: ', dtb
  print '(A,E12.5)', 'thresq: ', thresq
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,I0)', 'Warmup iterations: ', warmup_iterations
  print '(A)', '------------------------------------------------'

  ! Initialize output arrays
  frrg = 0.0
  frrgn = 0.0

  ! Warmup iterations
  do iter = 1, warmup_iterations
    call kernel_freezing(cphopt, dtb, thresq, ni, nj, nk, qr, ncr, tcel, &
                         diaqr, frrg, frrgn, t0, t0cel, tlow, rhow, cc)
  end do

  ! Benchmark iterations
  total_time = 0.0d0
  min_time = 1.0d30
  max_time = 0.0d0

  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_freezing(cphopt, dtb, thresq, ni, nj, nk, qr, ncr, tcel, &
                         diaqr, frrg, frrgn, t0, t0cel, tlow, rhow, cc)
    !$acc wait
    end_time = omp_get_wtime()

    elapsed_time = end_time - start_time
    total_time = total_time + elapsed_time
    min_time = min(min_time, elapsed_time)
    max_time = max(max_time, elapsed_time)
  end do

  ! Validate output
  call validate_output(frrg, frrg_ref, ni, nj, nk, tolerance, 'frrg', max_diff)
  call validate_output(frrgn, frrgn_ref, ni, nj, nk, tolerance, 'frrgn', max_diff)

  ! Report results
  print '(A)', ''
  print '(A)', 'Timing Results:'
  print '(A,F12.6,A)', '  Total time:   ', total_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Average time: ', (total_time / num_iterations) * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Min time:     ', min_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Max time:     ', max_time * 1000.0d0, ' ms'
  print '(A)', '================================================'

  ! Cleanup
  deallocate(qr, ncr, tcel, diaqr, frrg, frrgn, frrg_ref, frrgn_ref)

contains

  !> Main kernel subroutine - GPU version
  subroutine kernel_freezing(cphopt, dtb, thresq, ni, nj, nk, qr, ncr, tcel, &
                             diaqr, frrg, frrgn, t0, t0cel, tlow, rhow, cc)
    implicit none
    integer, intent(in) :: cphopt, ni, nj, nk
    real, intent(in) :: dtb, thresq, t0, t0cel, tlow, rhow, cc
    real, intent(in) :: qr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ncr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tcel(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: diaqr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: frrg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: frrgn(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    real :: tclow, cfrrg, cfrrgn
    real :: diaqr3, a

    tclow = tlow - t0
    cfrrg = 20.e0 * 100.e0 * cc * cc * rhow * dtb
    cfrrgn = 100.e0 * cc * dtb

    if (nk == 1) then
      if (abs(cphopt) == 2) then
        !$acc kernels
        !$acc loop independent collapse(2) private(diaqr3)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qr(i,j,1) > thresq) then
              if (tcel(i,j,1) < t0cel) then
                if (tcel(i,j,1) <= tclow) then
                  frrg(i,j,1) = qr(i,j,1)
                else
                  diaqr3 = diaqr(i,j,1) * diaqr(i,j,1) * diaqr(i,j,1)
                  frrg(i,j,1) = min(cfrrg * diaqr3 * diaqr3 * &
                    (exp(-0.66e0 * tcel(i,j,1)) - 1.e0) * ncr(i,j,1), qr(i,j,1))
                end if
              else
                frrg(i,j,1) = 0.e0
              end if
            else
              frrg(i,j,1) = 0.e0
            end if
          end do
        end do
        !$acc end kernels
      else if (abs(cphopt) >= 3) then
        !$acc kernels
        !$acc loop independent collapse(2) private(diaqr3, a)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qr(i,j,1) > thresq) then
              if (tcel(i,j,1) < t0cel) then
                if (tcel(i,j,1) <= tclow) then
                  frrg(i,j,1) = qr(i,j,1)
                  frrgn(i,j,1) = ncr(i,j,1)
                else
                  a = (exp(-0.66e0 * tcel(i,j,1)) - 1.e0) * ncr(i,j,1)
                  diaqr3 = diaqr(i,j,1) * diaqr(i,j,1) * diaqr(i,j,1)
                  frrg(i,j,1) = min(cfrrg * diaqr3 * diaqr3 * a, qr(i,j,1))
                  frrgn(i,j,1) = min(cfrrgn * diaqr3 * a, ncr(i,j,1))
                end if
              else
                frrg(i,j,1) = 0.e0
                frrgn(i,j,1) = 0.e0
              end if
            else
              frrg(i,j,1) = 0.e0
              frrgn(i,j,1) = 0.e0
            end if
          end do
        end do
        !$acc end kernels
      end if
    else
      if (abs(cphopt) == 2) then
        !$acc kernels
        !$acc loop independent collapse(3) private(diaqr3)
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni-1
              if (qr(i,j,k) > thresq) then
                if (tcel(i,j,k) < t0cel) then
                  if (tcel(i,j,k) <= tclow) then
                    frrg(i,j,k) = qr(i,j,k)
                  else
                    diaqr3 = diaqr(i,j,k) * diaqr(i,j,k) * diaqr(i,j,k)
                    frrg(i,j,k) = min(qr(i,j,k), cfrrg * diaqr3 * diaqr3 * &
                      (exp(-0.66e0 * tcel(i,j,k)) - 1.e0) * ncr(i,j,k))
                  end if
                else
                  frrg(i,j,k) = 0.e0
                end if
              else
                frrg(i,j,k) = 0.e0
              end if
            end do
          end do
        end do
        !$acc end kernels
      else if (abs(cphopt) >= 3) then
        !$acc kernels
        !$acc loop independent collapse(3) private(diaqr3, a)
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni-1
              if (qr(i,j,k) > thresq) then
                if (tcel(i,j,k) < t0cel) then
                  if (tcel(i,j,k) <= tclow) then
                    frrg(i,j,k) = qr(i,j,k)
                    frrgn(i,j,k) = ncr(i,j,k)
                  else
                    a = (exp(-0.66e0 * tcel(i,j,k)) - 1.e0) * ncr(i,j,k)
                    diaqr3 = diaqr(i,j,k) * diaqr(i,j,k) * diaqr(i,j,k)
                    frrg(i,j,k) = min(cfrrg * diaqr3 * diaqr3 * a, qr(i,j,k))
                    frrgn(i,j,k) = min(cfrrgn * diaqr3 * a, ncr(i,j,k))
                  end if
                else
                  frrg(i,j,k) = 0.e0
                  frrgn(i,j,k) = 0.e0
                end if
              else
                frrg(i,j,k) = 0.e0
                frrgn(i,j,k) = 0.e0
              end if
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_freezing

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
        case ('cphopt')
          read(param_value, *) cphopt
        case ('ni')
          read(param_value, *) ni
        case ('nj')
          read(param_value, *) nj
        case ('nk')
          read(param_value, *) nk
        case ('dtb')
          read(param_value, *) dtb
        case ('thresq')
          read(param_value, *) thresq
        case ('t0')
          read(param_value, *) t0
        case ('t0cel')
          read(param_value, *) t0cel
        case ('tlow')
          read(param_value, *) tlow
        case ('rhow')
          read(param_value, *) rhow
        case ('cc')
          read(param_value, *) cc
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

end program kernel_benchmark_freezing

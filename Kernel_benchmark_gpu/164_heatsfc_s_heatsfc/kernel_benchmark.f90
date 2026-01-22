!-----------------------------------------------------------------------
! GPU Kernel Benchmark Program: heatsfc
! Source: Src/heatsfc.f90
! Subroutine: s_heatsfc
! Description: Calculate sensible and latent heat on the surface
! GPU Difficulty: Medium - 2D loop with conditionals and private variables
!-----------------------------------------------------------------------
program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: ni, nj, nk, nund
  real :: cp, t0, lv0, lf0, cwmci
  character(len=5) :: fmois

  ! Physical constants (from comphy.f90)
  real, parameter :: tlow = 233.16e0

  ! Arrays
  real, allocatable :: t(:,:,:)
  real, allocatable :: qv(:,:,:)
  real, allocatable :: qvsfc(:,:)
  real, allocatable :: ct(:,:)
  real, allocatable :: cq(:,:)
  integer, allocatable :: land(:,:)
  real, allocatable :: kai(:,:)
  real, allocatable :: tund(:,:,:)
  real, allocatable :: tice(:,:)
  real, allocatable :: hs(:,:), hs_ref(:,:)
  real, allocatable :: le(:,:), le_ref(:,:)

  ! Benchmark variables
  character(len=256) :: data_dir
  integer :: num_iterations, warmup_iterations
  real :: tolerance
  double precision :: start_time, end_time, total_time, avg_time
  integer :: iter, errors
  real :: max_error_hs, max_error_le

  ! Read benchmark configuration
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! Read parameters from dump
  call read_parameters(trim(data_dir)//'/params.txt')

  ! Allocate arrays
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsfc(0:ni+1, 0:nj+1))
  allocate(ct(0:ni+1, 0:nj+1))
  allocate(cq(0:ni+1, 0:nj+1))
  allocate(land(0:ni+1, 0:nj+1))
  allocate(kai(0:ni+1, 0:nj+1))
  allocate(tund(0:ni+1, 0:nj+1, 1:nund))
  allocate(tice(0:ni+1, 0:nj+1))
  allocate(hs(0:ni+1, 0:nj+1))
  allocate(le(0:ni+1, 0:nj+1))
  allocate(hs_ref(0:ni+1, 0:nj+1))
  allocate(le_ref(0:ni+1, 0:nj+1))

  ! Read input arrays
  call read_array_3d(trim(data_dir)//'/t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/qvsfc.bin', qvsfc, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/ct.bin', ct, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/cq.bin', cq, 0, ni+1, 0, nj+1)
  call read_array_2d_int(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/kai.bin', kai, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/tund.bin', tund, 0, ni+1, 0, nj+1, 1, nund)
  call read_array_2d(trim(data_dir)//'/tice.bin', tice, 0, ni+1, 0, nj+1)

  ! Read reference output for validation
  call read_array_2d(trim(data_dir)//'/hs_ref.bin', hs_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/le_ref.bin', le_ref, 0, ni+1, 0, nj+1)

  ! Warmup iterations
  do iter = 1, warmup_iterations
    call kernel_heatsfc(fmois, ni, nj, nk, nund, t, qv, qvsfc, ct, cq, &
                        land, kai, tund, tice, hs, le, cp, t0, lv0, lf0, cwmci)
  end do
  !$acc wait

  ! Benchmark iterations
  total_time = 0.0d0
  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_heatsfc(fmois, ni, nj, nk, nund, t, qv, qvsfc, ct, cq, &
                        land, kai, tund, tice, hs, le, cp, t0, lv0, lf0, cwmci)
    !$acc wait
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results
  call validate_results(hs, hs_ref, le, le_ref, ni, nj, tolerance, errors, max_error_hs, max_error_le)

  ! Output results
  print '(A)',        '========================================'
  print '(A)',        'GPU Kernel: heatsfc (s_heatsfc)'
  print '(A)',        '========================================'
  print '(A,I0,A,I0,A,I0)', 'Grid size: ', ni, ' x ', nj, ' x ', nk
  print '(A,I0)',     'Soil layers: ', nund
  print '(A,A)',      'Moisture mode: ', trim(fmois)
  print '(A,I0)',     'Iterations: ', num_iterations
  print '(A,F12.6,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', 'Total time: ', total_time, ' s'
  print '(A,I0)',     'Validation errors: ', errors
  print '(A,ES12.5)', 'Max error (hs): ', max_error_hs
  print '(A,ES12.5)', 'Max error (le): ', max_error_le
  print '(A,ES12.5)', 'Tolerance: ', tolerance
  if (errors == 0) then
    print '(A)',      'VALIDATION: PASSED'
  else
    print '(A)',      'VALIDATION: FAILED'
  end if
  print '(A)',        '========================================'

  ! Cleanup
  deallocate(t, qv, qvsfc, ct, cq, land, kai, tund, tice)
  deallocate(hs, le, hs_ref, le_ref)

contains

  !---------------------------------------------------------------------
  ! Kernel subroutine: heatsfc
  !---------------------------------------------------------------------
  subroutine kernel_heatsfc(fmois, ni, nj, nk, nund, t, qv, qvsfc, ct, cq, &
                            land, kai, tund, tice, hs, le, cp, t0, lv0, lf0, cwmci)
    implicit none

    character(len=5), intent(in) :: fmois
    integer, intent(in) :: ni, nj, nk, nund
    real, intent(in) :: t(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qvsfc(0:ni+1, 0:nj+1)
    real, intent(in) :: ct(0:ni+1, 0:nj+1)
    real, intent(in) :: cq(0:ni+1, 0:nj+1)
    integer, intent(in) :: land(0:ni+1, 0:nj+1)
    real, intent(in) :: kai(0:ni+1, 0:nj+1)
    real, intent(in) :: tund(0:ni+1, 0:nj+1, 1:nund)
    real, intent(in) :: tice(0:ni+1, 0:nj+1)
    real, intent(in) :: cp, t0, lv0, lf0, cwmci
    real, intent(out) :: hs(0:ni+1, 0:nj+1)
    real, intent(out) :: le(0:ni+1, 0:nj+1)

    integer :: i, j
    real :: tmsfc, lva, lsa

    if (fmois(1:3) .eq. 'dry') then

      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent private(tmsfc)
        do i = 1, ni-1
          if (land(i,j) .eq. 1) then
            tmsfc = kai(i,j) * tice(i,j) + (1.e0 - kai(i,j)) * tund(i,j,1)
            hs(i,j) = cp * ct(i,j) * (tmsfc - t(i,j,2))
          else
            hs(i,j) = cp * ct(i,j) * (tund(i,j,1) - t(i,j,2))
          end if
          le(i,j) = 0.e0
        end do
      end do
      !$acc end kernels

    else if (fmois(1:5) .eq. 'moist') then

      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent private(tmsfc, lva, lsa)
        do i = 1, ni-1
          if (land(i,j) .lt. 0) then
            lva = lv0 * exp((.167e0 + 3.67e-4*t(i,j,2)) * log(t0/t(i,j,2)))
            hs(i,j) = cp * ct(i,j) * (tund(i,j,1) - t(i,j,2))
            le(i,j) = lva * cq(i,j) * (qvsfc(i,j) - qv(i,j,2))
          else if (land(i,j) .eq. 1) then
            tmsfc = kai(i,j) * tice(i,j) + (1.e0 - kai(i,j)) * tund(i,j,1)
            lva = lv0 * exp((.167e0 + 3.67e-4*t(i,j,2)) * log(t0/t(i,j,2)))
            hs(i,j) = cp * ct(i,j) * (tmsfc - t(i,j,2))
            le(i,j) = lva * cq(i,j) * (qvsfc(i,j) - qv(i,j,2))
          else
            if (t(i,j,2) .gt. tlow) then
              lva = lv0 * exp((.167e0 + 3.67e-4*t(i,j,2)) * log(t0/t(i,j,2)))
              hs(i,j) = cp * ct(i,j) * (tund(i,j,1) - t(i,j,2))
              le(i,j) = lva * cq(i,j) * (qvsfc(i,j) - qv(i,j,2))
            else
              lsa = lv0 * exp((.167e0 + 3.67e-4*t(i,j,2)) * log(t0/t(i,j,2))) &
                  + (lf0 + cwmci * (t(i,j,2) - t0))
              hs(i,j) = cp * ct(i,j) * (tund(i,j,1) - t(i,j,2))
              le(i,j) = lsa * cq(i,j) * (qvsfc(i,j) - qv(i,j,2))
            end if
          end if
        end do
      end do
      !$acc end kernels

    end if

  end subroutine kernel_heatsfc

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: unit_num, ios

    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open benchmark.conf'
      stop 1
    end if

    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iter
    read(unit_num, *) warmup_iter
    read(unit_num, *) tol
    close(unit_num)
  end subroutine read_config

  !---------------------------------------------------------------------
  ! Read parameters from dump file
  !---------------------------------------------------------------------
  subroutine read_parameters(filename)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=256) :: line, key, value_str
    integer :: unit_num, ios, eq_pos

    unit_num = 11
    open(unit=unit_num, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        value_str = adjustl(line(eq_pos+1:))

        select case (trim(key))
          case ('ni')
            read(value_str, *) ni
          case ('nj')
            read(value_str, *) nj
          case ('nk')
            read(value_str, *) nk
          case ('nund')
            read(value_str, *) nund
          case ('cp')
            read(value_str, *) cp
          case ('t0')
            read(value_str, *) t0
          case ('lv0')
            read(value_str, *) lv0
          case ('lf0')
            read(value_str, *) lf0
          case ('cwmci')
            read(value_str, *) cwmci
          case ('fmois')
            read(value_str, '(A)') fmois
        end select
      end if
    end do

    close(unit_num)
  end subroutine read_parameters

  !---------------------------------------------------------------------
  ! Read 2D real array from binary file
  !---------------------------------------------------------------------
  subroutine read_array_2d(filename, array, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: array(is:ie, js:je)
    integer :: unit_num, ios

    unit_num = 12
    open(unit=unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(unit_num) array
    close(unit_num)
  end subroutine read_array_2d

  !---------------------------------------------------------------------
  ! Read 2D integer array from binary file
  !---------------------------------------------------------------------
  subroutine read_array_2d_int(filename, array, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    integer, intent(out) :: array(is:ie, js:je)
    integer :: unit_num, ios

    unit_num = 12
    open(unit=unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(unit_num) array
    close(unit_num)
  end subroutine read_array_2d_int

  !---------------------------------------------------------------------
  ! Read 3D real array from binary file
  !---------------------------------------------------------------------
  subroutine read_array_3d(filename, array, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: array(is:ie, js:je, ks:ke)
    integer :: unit_num, ios

    unit_num = 12
    open(unit=unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(unit_num) array
    close(unit_num)
  end subroutine read_array_3d

  !---------------------------------------------------------------------
  ! Validate results against reference
  !---------------------------------------------------------------------
  subroutine validate_results(hs, hs_ref, le, le_ref, ni, nj, tolerance, errors, max_err_hs, max_err_le)
    implicit none
    integer, intent(in) :: ni, nj
    real, intent(in) :: hs(0:ni+1, 0:nj+1), hs_ref(0:ni+1, 0:nj+1)
    real, intent(in) :: le(0:ni+1, 0:nj+1), le_ref(0:ni+1, 0:nj+1)
    real, intent(in) :: tolerance
    integer, intent(out) :: errors
    real, intent(out) :: max_err_hs, max_err_le

    integer :: i, j
    real :: rel_err, denom

    errors = 0
    max_err_hs = 0.0
    max_err_le = 0.0

    do j = 1, nj-1
      do i = 1, ni-1
        ! Check hs
        denom = max(abs(hs_ref(i,j)), 1.0e-20)
        rel_err = abs(hs(i,j) - hs_ref(i,j)) / denom
        max_err_hs = max(max_err_hs, rel_err)
        if (rel_err > tolerance) then
          errors = errors + 1
        end if

        ! Check le
        denom = max(abs(le_ref(i,j)), 1.0e-20)
        rel_err = abs(le(i,j) - le_ref(i,j)) / denom
        max_err_le = max(max_err_le, rel_err)
        if (rel_err > tolerance) then
          errors = errors + 1
        end if
      end do
    end do
  end subroutine validate_results

end program kernel_benchmark

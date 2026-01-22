!> GPU Kernel benchmark program for distrpg (s_distrpg)
!> Calculates the distribution ratio at which the collisions between
!> rain water and snow, and resets the collection rate.
!> GPU Difficulty: Medium - Multiple branches and conditional computations
program kernel_benchmark_distrpg
  use omp_lib
  implicit none

  ! Parameters from dump
  integer :: cphopt
  integer :: ni, nj, nk
  real :: thresq
  real :: t0

  ! Physical constants (from m_comphy)
  real, parameter :: rhos = 1.0e2   ! Snow density [kg/m^3]
  real, parameter :: rhow = 1.0e3   ! Water density [kg/m^3]

  ! Arrays
  real, allocatable :: t(:,:,:)       ! Air temperature
  real, allocatable :: qs(:,:,:)      ! Snow mixing ratio
  real, allocatable :: diaqr(:,:,:)   ! Mean diameter of rain water
  real, allocatable :: diaqs(:,:,:)   ! Mean diameter of snow
  real, allocatable :: clrs(:,:,:)    ! Collection rate: rain to snow
  real, allocatable :: clsr(:,:,:)    ! Collection rate: snow to rain
  real, allocatable :: clrsn(:,:,:)   ! Collection rate for concentrations
  real, allocatable :: clsrn(:,:,:)   ! Collection rate for concentrations
  real, allocatable :: clrsg(:,:,:)   ! Production rate of graupel
  real, allocatable :: clrs_save(:,:,:)
  real, allocatable :: clsr_save(:,:,:)
  real, allocatable :: clrsn_save(:,:,:)
  real, allocatable :: clsrn_save(:,:,:)
  real, allocatable :: clrsg_ref(:,:,:)
  real, allocatable :: clrs_ref(:,:,:)
  real, allocatable :: clsr_ref(:,:,:)
  real, allocatable :: clrsn_ref(:,:,:)
  real, allocatable :: clsrn_ref(:,:,:)

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
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqr(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsr(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrsn(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsrn(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrsg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrs_save(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsr_save(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrsn_save(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsrn_save(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrsg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrs_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrsn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsrn_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Read input data
  call read_array_3d(trim(data_dir)//'/t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/diaqr.bin', diaqr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/diaqs.bin', diaqs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clrs_in.bin', clrs_save, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsr_in.bin', clsr_save, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clrsn_in.bin', clrsn_save, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsrn_in.bin', clsrn_save, 0, ni+1, 0, nj+1, 1, nk)

  ! Read reference output
  call read_array_3d(trim(data_dir)//'/clrsg_ref.bin', clrsg_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clrs_ref.bin', clrs_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsr_ref.bin', clsr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clrsn_ref.bin', clrsn_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsrn_ref.bin', clsrn_ref, 0, ni+1, 0, nj+1, 1, nk)

  print '(A)', '================================================'
  print '(A)', 'GPU Kernel Benchmark: distrpg (s_distrpg)'
  print '(A)', 'Rain-snow collision distribution rate calculation'
  print '(A)', '================================================'
  print '(A,I0,A,I0,A,I0)', 'Grid size: ', ni, ' x ', nj, ' x ', nk
  print '(A,I0)', 'cphopt: ', cphopt
  print '(A,E12.5)', 'thresq: ', thresq
  print '(A,E12.5)', 't0: ', t0
  print '(A,I0)', 'Benchmark iterations: ', num_iterations
  print '(A,I0)', 'Warmup iterations: ', warmup_iterations
  print '(A)', '------------------------------------------------'

  ! Warmup iterations
  do iter = 1, warmup_iterations
    clrs = clrs_save
    clsr = clsr_save
    clrsn = clrsn_save
    clsrn = clsrn_save
    call kernel_distrpg(cphopt, thresq, ni, nj, nk, qs, t, diaqr, diaqs, &
                        clrs, clsr, clrsn, clsrn, clrsg, t0, rhos, rhow)
  end do
  !$acc wait

  ! Benchmark iterations
  total_time = 0.0d0
  min_time = 1.0d30
  max_time = 0.0d0

  do iter = 1, num_iterations
    clrs = clrs_save
    clsr = clsr_save
    clrsn = clrsn_save
    clsrn = clsrn_save

    !$acc wait
    start_time = omp_get_wtime()
    call kernel_distrpg(cphopt, thresq, ni, nj, nk, qs, t, diaqr, diaqs, &
                        clrs, clsr, clrsn, clsrn, clrsg, t0, rhos, rhow)
    !$acc wait
    end_time = omp_get_wtime()

    elapsed_time = end_time - start_time
    total_time = total_time + elapsed_time
    min_time = min(min_time, elapsed_time)
    max_time = max(max_time, elapsed_time)
  end do

  ! Validate output
  call validate_output(clrsg, clrsg_ref, ni, nj, nk, tolerance, 'clrsg', max_diff)
  call validate_output(clrs, clrs_ref, ni, nj, nk, tolerance, 'clrs', max_diff)
  call validate_output(clsr, clsr_ref, ni, nj, nk, tolerance, 'clsr', max_diff)
  call validate_output(clrsn, clrsn_ref, ni, nj, nk, tolerance, 'clrsn', max_diff)
  call validate_output(clsrn, clsrn_ref, ni, nj, nk, tolerance, 'clsrn', max_diff)

  ! Report results
  print '(A)', ''
  print '(A)', 'Timing Results:'
  print '(A,F12.6,A)', '  Total time:   ', total_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Average time: ', (total_time / num_iterations) * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Min time:     ', min_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', '  Max time:     ', max_time * 1000.0d0, ' ms'
  print '(A)', '================================================'

  ! Cleanup
  deallocate(t, qs, diaqr, diaqs)
  deallocate(clrs, clsr, clrsn, clsrn, clrsg)
  deallocate(clrs_save, clsr_save, clrsn_save, clsrn_save)
  deallocate(clrsg_ref, clrs_ref, clsr_ref, clrsn_ref, clsrn_ref)

contains

  !> Main kernel subroutine
  subroutine kernel_distrpg(cphopt, thresq, ni, nj, nk, qs, t, diaqr, diaqs, &
                            clrs, clsr, clrsn, clsrn, clrsg, t0, rhos, rhow)
    implicit none
    integer, intent(in) :: cphopt
    real, intent(in) :: thresq, t0, rhos, rhow
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: qs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: t(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: diaqr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: diaqs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: clrs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: clsr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: clrsn(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: clsrn(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: clrsg(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    real :: rhos2, rhow2
    real :: alpha, alpha1, a, b

    rhos2 = rhos * rhos
    rhow2 = rhow * rhow

    if (nk == 1) then
      if (abs(cphopt) == 2) then
        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(alpha, alpha1, a, b)
          do i = 1, ni-1
            if (qs(i,j,1) > thresq) then
              if (t(i,j,1) < t0) then
                a = diaqs(i,j,1) * diaqs(i,j,1)
                b = diaqr(i,j,1) * diaqr(i,j,1)
                a = rhos2 * a * a * a
                b = rhow2 * b * b * b
                alpha = a / (a + b)
                alpha1 = 1.e0 - alpha
                clrsg(i,j,1) = alpha1 * clrs(i,j,1)
                clrs(i,j,1) = alpha * clrs(i,j,1)
                clsr(i,j,1) = alpha1 * clsr(i,j,1)
              else
                clrsg(i,j,1) = clrs(i,j,1)
                clrs(i,j,1) = 0.e0
              end if
            else
              clrsg(i,j,1) = 0.e0
            end if
          end do
        end do
        !$acc end kernels
      else if (abs(cphopt) >= 3) then
        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(alpha, alpha1, a, b)
          do i = 1, ni-1
            if (qs(i,j,1) > thresq) then
              if (t(i,j,1) < t0) then
                a = diaqs(i,j,1) * diaqs(i,j,1)
                b = diaqr(i,j,1) * diaqr(i,j,1)
                a = rhos2 * a * a * a
                b = rhow2 * b * b * b
                alpha = a / (a + b)
                alpha1 = 1.e0 - alpha
                clrsg(i,j,1) = alpha1 * clrs(i,j,1)
                clrs(i,j,1) = alpha * clrs(i,j,1)
                clsr(i,j,1) = alpha1 * clsr(i,j,1)
                clrsn(i,j,1) = alpha1 * clrsn(i,j,1)
                clsrn(i,j,1) = alpha1 * clsrn(i,j,1)
              else
                clrsg(i,j,1) = clrs(i,j,1)
                clrs(i,j,1) = 0.e0
              end if
            else
              clrsg(i,j,1) = 0.e0
            end if
          end do
        end do
        !$acc end kernels
      end if
    else
      if (abs(cphopt) == 2) then
        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent private(alpha, alpha1, a, b)
            do i = 1, ni-1
              if (qs(i,j,k) > thresq) then
                if (t(i,j,k) < t0) then
                  a = diaqs(i,j,k) * diaqs(i,j,k)
                  b = diaqr(i,j,k) * diaqr(i,j,k)
                  a = rhos2 * a * a * a
                  b = rhow2 * b * b * b
                  alpha = a / (a + b)
                  alpha1 = 1.e0 - alpha
                  clrsg(i,j,k) = alpha1 * clrs(i,j,k)
                  clrs(i,j,k) = alpha * clrs(i,j,k)
                  clsr(i,j,k) = alpha1 * clsr(i,j,k)
                else
                  clrsg(i,j,k) = clrs(i,j,k)
                  clrs(i,j,k) = 0.e0
                end if
              else
                clrsg(i,j,k) = 0.e0
              end if
            end do
          end do
        end do
        !$acc end kernels
      else if (abs(cphopt) >= 3) then
        !$acc kernels
        !$acc loop independent
        do k = 1, nk-1
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent private(alpha, alpha1, a, b)
            do i = 1, ni-1
              if (qs(i,j,k) > thresq) then
                if (t(i,j,k) < t0) then
                  a = diaqs(i,j,k) * diaqs(i,j,k)
                  b = diaqr(i,j,k) * diaqr(i,j,k)
                  a = rhos2 * a * a * a
                  b = rhow2 * b * b * b
                  alpha = a / (a + b)
                  alpha1 = 1.e0 - alpha
                  clrsg(i,j,k) = alpha1 * clrs(i,j,k)
                  clrs(i,j,k) = alpha * clrs(i,j,k)
                  clsr(i,j,k) = alpha1 * clsr(i,j,k)
                  clrsn(i,j,k) = alpha1 * clrsn(i,j,k)
                  clsrn(i,j,k) = alpha1 * clsrn(i,j,k)
                else
                  clrsg(i,j,k) = clrs(i,j,k)
                  clrs(i,j,k) = 0.e0
                end if
              else
                clrsg(i,j,k) = 0.e0
              end if
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_distrpg

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
        case ('thresq')
          read(param_value, *) thresq
        case ('t0')
          read(param_value, *) t0
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

end program kernel_benchmark_distrpg

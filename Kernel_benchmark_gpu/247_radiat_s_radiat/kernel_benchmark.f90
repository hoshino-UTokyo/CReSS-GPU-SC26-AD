!***********************************************************************
! GPU Kernel Benchmark: radiat (s_radiat)
!***********************************************************************
!
! Source: Src/radiat.f90
! Description: Calculates short and long wave radiation including
!              zenith angle, z coordinates, and radiation fluxes.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_radiat
  use omp_lib
  implicit none

  ! Configuration
  character(len=256) :: data_dir
  integer :: num_iterations, warmup_iterations
  real :: tolerance

  ! Array dimensions
  integer :: ni, nj, nk, nkm1, nund

  ! Scalar parameters
  character(len=5) :: fmois
  integer :: cphopt
  real :: ln1013, esgm, esgm51, eqt, rchr, rcmn, sinphs, cosphs

  ! Physical constants
  real, parameter :: cc = 3.141592e0
  real, parameter :: d2r = 3.141592e0/180.e0
  real, parameter :: oned6 = 1.e0/6.e0
  real, parameter :: oned15 = 1.e0/15.e0
  real, parameter :: oned60 = 1.e0/60.e0
  real, parameter :: eps = 1.e-20
  real, parameter :: epsva = 0.622e0
  real, parameter :: sigma = 5.67e-8
  real, parameter :: sun0 = 1367.e0
  real, parameter :: el = 0.95e0
  real, parameter :: icalbe = 0.5e0
  real, parameter :: zarad = 100.e0

  ! Input arrays
  real, allocatable :: zph(:,:,:)
  real, allocatable :: lat(:,:), lon(:,:)
  real, allocatable :: p(:,:,:), t(:,:,:), qv(:,:,:)
  integer, allocatable :: land(:,:)
  real, allocatable :: albe(:,:), kai(:,:)
  real, allocatable :: tund(:,:,:), tice(:,:)
  real, allocatable :: cdl(:,:), cdm(:,:), cdh(:,:)
  real, allocatable :: fall(:,:)

  ! Inout arrays
  real, allocatable :: zph8s(:,:,:), zph8s_in(:,:,:)
  real, allocatable :: zref(:,:), zref_in(:,:)
  real, allocatable :: coseta(:,:), coseta_in(:,:)

  ! Output arrays
  real, allocatable :: rgd(:,:), rsd(:,:), rld(:,:), rlu(:,:)

  ! Reference arrays
  real, allocatable :: rgd_ref(:,:), rsd_ref(:,:), rld_ref(:,:), rlu_ref(:,:)
  real, allocatable :: zph8s_ref(:,:,:), zref_ref(:,:), coseta_ref(:,:)

  ! Timing and validation
  real(8) :: start_time, end_time, elapsed_time, total_time
  real(8) :: t_min, t_max
  real :: max_error_rgd, max_error_rsd, max_error_rld, max_error_rlu
  real :: max_error_zph8s, max_error_zref, max_error_coseta
  integer :: iter
  logical :: validation_passed

  ! Read configuration
  call read_config()

  ! Read parameters
  call read_parameters()

  ! Allocate arrays
  call allocate_arrays()

  ! Read input data
  call read_input_data()

  ! Read reference data
  call read_reference_data()

  ! Print header
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: radiat (s_radiat)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,A)') ' fmois: ', trim(fmois)
  write(*,'(A,I2)') ' cphopt: ', cphopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  ! Warmup
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call reset_inout_arrays()
    call s_radiat_kernel()
    !$acc wait
  end do

  ! Benchmark
  write(*,'(A)') ' Running benchmark iterations...'
  total_time = 0.0d0
  t_min = huge(t_min)
  t_max = 0.0d0

  do iter = 1, num_iterations
    call reset_inout_arrays()
    !$acc wait
    start_time = omp_get_wtime()
    call s_radiat_kernel()
    !$acc wait
    end_time = omp_get_wtime()
    elapsed_time = end_time - start_time
    total_time = total_time + elapsed_time
    t_min = min(t_min, elapsed_time)
    t_max = max(t_max, elapsed_time)
  end do

  elapsed_time = total_time / dble(num_iterations)

  ! Validate
  call validate_results(max_error_rgd, max_error_rsd, max_error_rld, max_error_rlu, &
                        max_error_zph8s, max_error_zref, max_error_coseta)

  ! Check validation
  validation_passed = (max_error_rgd < tolerance .and. max_error_rsd < tolerance .and. &
      max_error_rld < tolerance .and. max_error_rlu < tolerance .and. &
      max_error_zph8s < tolerance .and. max_error_zref < tolerance .and. &
      max_error_coseta < tolerance)

  ! Print results
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', elapsed_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', total_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', t_max * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.5)') ' Max error (rgd):    ', max_error_rgd
  write(*,'(A,ES12.5)') ' Max error (rsd):    ', max_error_rsd
  write(*,'(A,ES12.5)') ' Max error (rld):    ', max_error_rld
  write(*,'(A,ES12.5)') ' Max error (rlu):    ', max_error_rlu
  write(*,'(A,ES12.5)') ' Max error (zph8s):  ', max_error_zph8s
  write(*,'(A,ES12.5)') ' Max error (zref):   ', max_error_zref
  write(*,'(A,ES12.5)') ' Max error (coseta): ', max_error_coseta
  write(*,'(A,ES12.5)') ' Tolerance:          ', tolerance

  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! Cleanup
  call deallocate_arrays()

  if (.not. validation_passed) stop 1

contains

  subroutine read_config()
    integer :: unit_num, ios
    character(len=256) :: line

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
    integer :: unit_num, ios
    character(len=256) :: line, key
    character(len=256) :: value_str
    integer :: eq_pos

    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open params.txt'
      stop 1
    end if

    fmois = 'moist'  ! Default

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit

      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      eq_pos = index(line, '=')
      if (eq_pos == 0) cycle

      key = adjustl(trim(line(1:eq_pos-1)))
      value_str = adjustl(trim(line(eq_pos+1:)))

      select case (trim(key))
        case ('fmois')
          fmois = trim(value_str)
        case ('cphopt')
          read(value_str, *) cphopt
        case ('ni')
          read(value_str, *) ni
        case ('nj')
          read(value_str, *) nj
        case ('nk')
          read(value_str, *) nk
        case ('nkm1')
          read(value_str, *) nkm1
        case ('nund')
          read(value_str, *) nund
        case ('ln1013')
          read(value_str, *) ln1013
        case ('esgm')
          read(value_str, *) esgm
        case ('esgm51')
          read(value_str, *) esgm51
        case ('eqt')
          read(value_str, *) eqt
        case ('rchr')
          read(value_str, *) rchr
        case ('rcmn')
          read(value_str, *) rcmn
        case ('sinphs')
          read(value_str, *) sinphs
        case ('cosphs')
          read(value_str, *) cosphs
      end select
    end do

    close(unit_num)
  end subroutine read_parameters

  subroutine allocate_arrays()
    allocate(zph(0:ni+1, 0:nj+1, 1:nk))
    allocate(lat(0:ni+1, 0:nj+1))
    allocate(lon(0:ni+1, 0:nj+1))
    allocate(p(0:ni+1, 0:nj+1, 1:nk))
    allocate(t(0:ni+1, 0:nj+1, 1:nk))
    allocate(qv(0:ni+1, 0:nj+1, 1:nk))
    allocate(land(0:ni+1, 0:nj+1))
    allocate(albe(0:ni+1, 0:nj+1))
    allocate(kai(0:ni+1, 0:nj+1))
    allocate(tund(0:ni+1, 0:nj+1, 1:nund))
    allocate(tice(0:ni+1, 0:nj+1))
    allocate(cdl(0:ni+1, 0:nj+1))
    allocate(cdm(0:ni+1, 0:nj+1))
    allocate(cdh(0:ni+1, 0:nj+1))
    allocate(fall(0:ni+1, 0:nj+1))

    allocate(zph8s(0:ni+1, 0:nj+1, 1:nk))
    allocate(zph8s_in(0:ni+1, 0:nj+1, 1:nk))
    allocate(zref(0:ni+1, 0:nj+1))
    allocate(zref_in(0:ni+1, 0:nj+1))
    allocate(coseta(0:ni+1, 0:nj+1))
    allocate(coseta_in(0:ni+1, 0:nj+1))

    allocate(rgd(0:ni+1, 0:nj+1))
    allocate(rsd(0:ni+1, 0:nj+1))
    allocate(rld(0:ni+1, 0:nj+1))
    allocate(rlu(0:ni+1, 0:nj+1))

    allocate(rgd_ref(0:ni+1, 0:nj+1))
    allocate(rsd_ref(0:ni+1, 0:nj+1))
    allocate(rld_ref(0:ni+1, 0:nj+1))
    allocate(rlu_ref(0:ni+1, 0:nj+1))
    allocate(zph8s_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(zref_ref(0:ni+1, 0:nj+1))
    allocate(coseta_ref(0:ni+1, 0:nj+1))
  end subroutine allocate_arrays

  subroutine read_input_data()
    call read_3d_real(trim(data_dir)//'/zph.bin', zph, 0, ni+1, 0, nj+1, 1, nk)
    call read_2d_real(trim(data_dir)//'/lat.bin', lat, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/lon.bin', lon, 0, ni+1, 0, nj+1)
    call read_3d_real(trim(data_dir)//'/p.bin', p, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
    call read_2d_int(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/albe.bin', albe, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/kai.bin', kai, 0, ni+1, 0, nj+1)
    call read_3d_real(trim(data_dir)//'/tund.bin', tund, 0, ni+1, 0, nj+1, 1, nund)
    call read_2d_real(trim(data_dir)//'/tice.bin', tice, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/cdl.bin', cdl, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/cdm.bin', cdm, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/cdh.bin', cdh, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/fall.bin', fall, 0, ni+1, 0, nj+1)

    call read_3d_real(trim(data_dir)//'/zph8s_in.bin', zph8s_in, 0, ni+1, 0, nj+1, 1, nk)
    call read_2d_real(trim(data_dir)//'/zref_in.bin', zref_in, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/coseta_in.bin', coseta_in, 0, ni+1, 0, nj+1)
  end subroutine read_input_data

  subroutine read_reference_data()
    call read_2d_real(trim(data_dir)//'/rgd_ref.bin', rgd_ref, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/rsd_ref.bin', rsd_ref, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/rld_ref.bin', rld_ref, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/rlu_ref.bin', rlu_ref, 0, ni+1, 0, nj+1)
    call read_3d_real(trim(data_dir)//'/zph8s_ref.bin', zph8s_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_2d_real(trim(data_dir)//'/zref_ref.bin', zref_ref, 0, ni+1, 0, nj+1)
    call read_2d_real(trim(data_dir)//'/coseta_ref.bin', coseta_ref, 0, ni+1, 0, nj+1)
  end subroutine read_reference_data

  subroutine reset_inout_arrays()
    zph8s(:,:,:) = zph8s_in(:,:,:)
    zref(:,:) = zref_in(:,:)
    coseta(:,:) = coseta_in(:,:)
    rgd(:,:) = 0.0
    rsd(:,:) = 0.0
    rld(:,:) = 0.0
    rlu(:,:) = 0.0
  end subroutine reset_inout_arrays

  subroutine read_2d_real(filename, array, is, ie, js, je)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: array(is:ie, js:je)
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
  end subroutine read_2d_real

  subroutine read_2d_int(filename, array, is, ie, js, je)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    integer, intent(out) :: array(is:ie, js:je)
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
  end subroutine read_2d_int

  subroutine read_3d_real(filename, array, is, ie, js, je, ks, ke)
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
  end subroutine read_3d_real

  subroutine s_radiat_kernel()
    integer :: i, j, k
    real :: tlc, pa, ta, ea, cdall, absrp, dk, a, b

    ! Calculate the zenith angle (2D loop)
    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent private(tlc)
      do i = 1, ni-1
        tlc = rchr + rcmn + oned15 * lon(i,j)
        coseta(i,j) = sinphs * sin(lat(i,j) * d2r) &
                    + cosphs * cos(lat(i,j) * d2r) * cos(eqt + 15.e0 * (tlc - 12.e0) * d2r)
      end do
    end do
    !$acc end kernels

    ! Get z physical coordinates at scalar points (3D loop with outer k)
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          zph8s(i,j,k) = 0.5e0 * (zph(i,j,k) + zph(i,j,k+1))
        end do
      end do
    end do
    !$acc end kernels

    ! Reference z coordinates (2D loop)
    !$acc kernels
    !$acc loop independent
    do j = 1, nj-1
      !$acc loop independent
      do i = 1, ni-1
        zref(i,j) = min(zarad + zph(i,j,2), zph8s(i,j,nkm1))
      end do
    end do
    !$acc end kernels

    ! Dry air case
    if (fmois(1:3) == 'dry') then
      !$acc kernels
      !$acc loop independent
      do k = 1, nk-2
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(dk, pa, ta, absrp, a)
          do i = 1, ni-1
            if (zref(i,j) > zph8s(i,j,k) .and. zref(i,j) <= zph8s(i,j,k+1)) then
              dk = (zref(i,j) - zph8s(i,j,k)) / (zph8s(i,j,k+1) - zph8s(i,j,k))
              pa = (1.e0 - dk) * p(i,j,k) + dk * p(i,j,k+1)
              ta = (1.e0 - dk) * t(i,j,k) + dk * t(i,j,k+1)

              if (coseta(i,j) > 0.e0) then
                if (land(i,j) < 0) then
                  absrp = 1.e0 - (9.e0 * (1.e0 - coseta(i,j)) * albe(i,j) + albe(i,j))
                else if (land(i,j) == 1) then
                  absrp = 1.e0 - (kai(i,j) * icalbe + (1.e0 - kai(i,j)) &
                        * (9.e0 * (1.e0 - coseta(i,j)) * albe(i,j) + albe(i,j)))
                else
                  absrp = max(1.e0 - (0.5e0 * (1.e0 - coseta(i,j)) * albe(i,j) + albe(i,j)), 0.e0)
                end if

                rgd(i,j) = sun0 * (0.554e0 + 0.43e0 * exp(ln1013 / (coseta(i,j) + eps))) * coseta(i,j)
                rsd(i,j) = absrp * rgd(i,j)
              else
                rgd(i,j) = 0.e0
                rsd(i,j) = 0.e0
              end if

              a = ta * ta
              rld(i,j) = esgm51 * a * a

              if (land(i,j) == 1) then
                a = kai(i,j) * tice(i,j) + (1.e0 - kai(i,j)) * tund(i,j,1)
                a = a * a
              else
                a = tund(i,j,1) * tund(i,j,1)
              end if
              rlu(i,j) = esgm * a * a
            end if
          end do
        end do
      end do
      !$acc end kernels

    ! Moist air case
    else if (fmois(1:5) == 'moist') then

      ! No cloud physics
      if (abs(cphopt) == 0) then
        !$acc kernels
        !$acc loop independent
        do k = 1, nk-2
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent private(dk, pa, ta, ea, cdall, absrp, a, b)
            do i = 1, ni-1
              if (zref(i,j) > zph8s(i,j,k) .and. zref(i,j) <= zph8s(i,j,k+1)) then
                dk = (zref(i,j) - zph8s(i,j,k)) / (zph8s(i,j,k+1) - zph8s(i,j,k))
                pa = (1.e0 - dk) * p(i,j,k) + dk * p(i,j,k+1)
                ta = (1.e0 - dk) * t(i,j,k) + dk * t(i,j,k+1)
                ea = (1.e0 - dk) * qv(i,j,k) + dk * qv(i,j,k+1)
                ea = pa * ea / (epsva + ea)
                cdall = cdl(i,j) + cdm(i,j) + cdh(i,j)

                if (coseta(i,j) > 0.e0) then
                  if (land(i,j) < 0) then
                    absrp = 1.e0 - ((9.e0 - 3.e0 * cdall) * (1.e0 - coseta(i,j)) * albe(i,j) + albe(i,j))
                  else if (land(i,j) == 1) then
                    absrp = 1.e0 - (kai(i,j) * icalbe + (1.e0 - kai(i,j)) * ((9.e0 - 3.e0 * cdall) &
                          * (1.e0 - coseta(i,j)) * albe(i,j) + albe(i,j)))
                  else
                    absrp = max(1.e0 - ((0.5e0 - oned6 * cdall) * (1.e0 - coseta(i,j)) * albe(i,j) + albe(i,j)), 0.e0)
                  end if

                  b = 0.43e0 + 0.00016e0 * ea
                  if (ea > 3000.e0) then
                    a = 0.e0
                  else if (ea > 100.e0 .and. ea <= 3000.e0) then
                    a = 1.12e0 - b - 0.06e0 * log10(ea)
                  else
                    a = 0.554e0
                  end if

                  rgd(i,j) = sun0 * (a + b * exp(ln1013 / (coseta(i,j) + eps))) &
                           * (1.e0 - 0.7e0 * cdl(i,j)) * (1.e0 - 0.6e0 * cdm(i,j)) &
                           * (1.e0 - 0.3e0 * cdh(i,j)) * coseta(i,j)
                  rsd(i,j) = absrp * rgd(i,j)
                else
                  rgd(i,j) = 0.e0
                  rsd(i,j) = 0.e0
                end if

                a = cdl(i,j) + 0.85e0 * cdm(i,j) + 0.5e0 * cdh(i,j)
                b = ta * ta
                rld(i,j) = esgm * b * b * (1.e0 + (0.66e-2 * sqrt(ea) - 0.49e0) * (1.e0 - (0.75e0 - 0.5e-4 * ea) * a))

                if (land(i,j) == 1) then
                  b = kai(i,j) * tice(i,j) + (1.e0 - kai(i,j)) * tund(i,j,1)
                  b = b * b
                else
                  b = tund(i,j,1) * tund(i,j,1)
                end if
                rlu(i,j) = esgm * b * b
              end if
            end do
          end do
        end do
        !$acc end kernels

      ! Cloud physics case
      else
        !$acc kernels
        !$acc loop independent
        do k = 1, nk-2
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent private(dk, pa, ta, ea, cdall, absrp, a, b)
            do i = 1, ni-1
              if (zref(i,j) > zph8s(i,j,k) .and. zref(i,j) <= zph8s(i,j,k+1)) then
                dk = (zref(i,j) - zph8s(i,j,k)) / (zph8s(i,j,k+1) - zph8s(i,j,k))
                pa = (1.e0 - dk) * p(i,j,k) + dk * p(i,j,k+1)
                ta = (1.e0 - dk) * t(i,j,k) + dk * t(i,j,k+1)
                ea = (1.e0 - dk) * qv(i,j,k) + dk * qv(i,j,k+1)
                ea = pa * ea / (epsva + ea)
                cdall = cdl(i,j) + cdm(i,j) + cdh(i,j)

                if (coseta(i,j) > 0.e0) then
                  if (land(i,j) < 0) then
                    absrp = 1.e0 - ((9.e0 - 3.e0 * cdall) * (1.e0 - coseta(i,j)) * albe(i,j) + albe(i,j))
                  else if (land(i,j) == 1) then
                    absrp = 1.e0 - (kai(i,j) * icalbe + (1.e0 - kai(i,j)) * ((9.e0 - 3.e0 * cdall) &
                          * (1.e0 - coseta(i,j)) * albe(i,j) + albe(i,j)))
                  else
                    absrp = max(1.e0 - ((0.5e0 - oned6 * cdall) * (1.e0 - coseta(i,j)) * albe(i,j) + albe(i,j)), 0.e0)
                  end if

                  b = 0.43e0 + 0.00016e0 * ea
                  if (ea > 3000.e0) then
                    a = 0.e0
                  else if (ea > 100.e0 .and. ea <= 3000.e0) then
                    a = 1.12e0 - b - 0.06e0 * log10(ea)
                  else
                    a = 0.554e0
                  end if

                  rgd(i,j) = sun0 * (a + b * exp(ln1013 / (coseta(i,j) + eps))) &
                           * (1.e0 - 0.7e0 * cdl(i,j)) * (1.e0 - 0.6e0 * cdm(i,j)) &
                           * (1.e0 - 0.3e0 * cdh(i,j)) * coseta(i,j)
                  rsd(i,j) = absrp * rgd(i,j)
                else
                  rgd(i,j) = 0.e0
                  rsd(i,j) = 0.e0
                end if

                if (fall(i,j) > 0.e0) then
                  a = cdl(i,j) + 0.85e0 * cdm(i,j) + 0.5e0 * cdh(i,j) + 0.1e0 * cdall
                else
                  a = cdl(i,j) + 0.85e0 * cdm(i,j) + 0.5e0 * cdh(i,j)
                end if

                b = ta * ta
                rld(i,j) = esgm * b * b * (1.e0 + (0.66e-2 * sqrt(ea) - 0.49e0) * (1.e0 - (0.75e0 - 0.5e-4 * ea) * a))

                if (land(i,j) == 1) then
                  b = kai(i,j) * tice(i,j) + (1.e0 - kai(i,j)) * tund(i,j,1)
                  b = b * b
                else
                  b = tund(i,j,1) * tund(i,j,1)
                end if
                rlu(i,j) = esgm * b * b
              end if
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine s_radiat_kernel

  subroutine validate_results(max_err_rgd, max_err_rsd, max_err_rld, max_err_rlu, &
                              max_err_zph8s, max_err_zref, max_err_coseta)
    real, intent(out) :: max_err_rgd, max_err_rsd, max_err_rld, max_err_rlu
    real, intent(out) :: max_err_zph8s, max_err_zref, max_err_coseta
    integer :: i, j, k
    real :: err, ref_val

    max_err_rgd = 0.0
    max_err_rsd = 0.0
    max_err_rld = 0.0
    max_err_rlu = 0.0
    max_err_zph8s = 0.0
    max_err_zref = 0.0
    max_err_coseta = 0.0

    do j = 1, nj-1
      do i = 1, ni-1
        ref_val = abs(rgd_ref(i,j))
        if (ref_val > 1.0e-30) then
          err = abs(rgd(i,j) - rgd_ref(i,j)) / ref_val
          max_err_rgd = max(max_err_rgd, err)
        end if

        ref_val = abs(rsd_ref(i,j))
        if (ref_val > 1.0e-30) then
          err = abs(rsd(i,j) - rsd_ref(i,j)) / ref_val
          max_err_rsd = max(max_err_rsd, err)
        end if

        ref_val = abs(rld_ref(i,j))
        if (ref_val > 1.0e-30) then
          err = abs(rld(i,j) - rld_ref(i,j)) / ref_val
          max_err_rld = max(max_err_rld, err)
        end if

        ref_val = abs(rlu_ref(i,j))
        if (ref_val > 1.0e-30) then
          err = abs(rlu(i,j) - rlu_ref(i,j)) / ref_val
          max_err_rlu = max(max_err_rlu, err)
        end if

        ref_val = abs(zref_ref(i,j))
        if (ref_val > 1.0e-30) then
          err = abs(zref(i,j) - zref_ref(i,j)) / ref_val
          max_err_zref = max(max_err_zref, err)
        end if

        ref_val = abs(coseta_ref(i,j))
        if (ref_val > 1.0e-30) then
          err = abs(coseta(i,j) - coseta_ref(i,j)) / ref_val
          max_err_coseta = max(max_err_coseta, err)
        end if
      end do
    end do

    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          ref_val = abs(zph8s_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(zph8s(i,j,k) - zph8s_ref(i,j,k)) / ref_val
            max_err_zph8s = max(max_err_zph8s, err)
          end if
        end do
      end do
    end do
  end subroutine validate_results

  subroutine deallocate_arrays()
    deallocate(zph, lat, lon, p, t, qv)
    deallocate(land, albe, kai, tund, tice)
    deallocate(cdl, cdm, cdh, fall)
    deallocate(zph8s, zph8s_in, zref, zref_in, coseta, coseta_in)
    deallocate(rgd, rsd, rld, rlu)
    deallocate(rgd_ref, rsd_ref, rld_ref, rlu_ref)
    deallocate(zph8s_ref, zref_ref, coseta_ref)
  end subroutine deallocate_arrays

end program kernel_benchmark_gpu_radiat

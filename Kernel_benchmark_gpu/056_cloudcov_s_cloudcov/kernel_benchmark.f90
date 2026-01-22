!***********************************************************************
! GPU Kernel Benchmark: cloudcov (s_cloudcov)
!***********************************************************************
!
! Source: Src/cloudcov.f90
! Description: Calculate low/mid/high cloud cover from relative humidity
!              or hydrometeor mixing ratios with multiple interpolation levels
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_cloudcov
  use omp_lib
  implicit none

  ! Parameters
  integer :: cphopt, ni, nj, nk
  real :: dz, es0iv2
  character(len=5) :: fmois
  character(len=3) :: fproc

  ! Physical constants (from comphy.f90)
  real, parameter :: tlow = 233.16e0
  real, parameter :: t0 = 273.16e0
  real, parameter :: epsva = 0.622e0

  ! Mathematical constants (from commath.f90)
  real, parameter :: tend3 = 10.e0/3.e0
  real, parameter :: tend6 = 10.e0/6.e0
  real, parameter :: tend7 = 10.e0/7.e0

  ! Lookup tables
  real :: rcdl(0:100), rcdm(0:100), rcdh(0:100)

  ! Arrays
  real, allocatable :: zph(:,:,:), rst(:,:,:), p(:,:,:), t(:,:,:)
  real, allocatable :: qv(:,:,:), qall(:,:,:)
  real, allocatable :: cdl(:,:), cdm(:,:), cdh(:,:)
  real, allocatable :: cdl_ref(:,:), cdm_ref(:,:), cdh_ref(:,:)
  real, allocatable :: zph8s(:,:,:), qsum(:,:,:)
  real, allocatable :: rh24(:,:), rh32(:,:), rh48(:,:), rh72(:,:)
  real, allocatable :: qsuml(:,:), qsumm(:,:), qsumh(:,:)

  ! Initial values for resetting
  real, allocatable :: zph8s_init(:,:,:), qsum_init(:,:,:)
  real, allocatable :: rh24_init(:,:), rh32_init(:,:), rh48_init(:,:), rh72_init(:,:)
  real, allocatable :: qsuml_init(:,:), qsumm_init(:,:), qsumh_init(:,:)

  ! Benchmark variables
  character(len=256) :: data_dir
  integer :: num_iterations, warmup_iterations
  real :: tolerance
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)
  integer :: iter, errors
  real :: max_error_cdl, max_error_cdm, max_error_cdh

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters from dump
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt')

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: cloudcov'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,A)') ' fmois: ', trim(fmois)
  write(*,'(A,A)') ' fproc: ', trim(fproc)
  write(*,'(A,I6)') ' cphopt: ', cphopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(zph(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(p(0:ni+1, 0:nj+1, 1:nk))
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qall(0:ni+1, 0:nj+1, 1:nk))
  allocate(cdl(0:ni+1, 0:nj+1))
  allocate(cdm(0:ni+1, 0:nj+1))
  allocate(cdh(0:ni+1, 0:nj+1))
  allocate(cdl_ref(0:ni+1, 0:nj+1))
  allocate(cdm_ref(0:ni+1, 0:nj+1))
  allocate(cdh_ref(0:ni+1, 0:nj+1))
  allocate(zph8s(0:ni+1, 0:nj+1, 1:nk))
  allocate(qsum(0:ni+1, 0:nj+1, 1:nk))
  allocate(rh24(0:ni+1, 0:nj+1))
  allocate(rh32(0:ni+1, 0:nj+1))
  allocate(rh48(0:ni+1, 0:nj+1))
  allocate(rh72(0:ni+1, 0:nj+1))
  allocate(qsuml(0:ni+1, 0:nj+1))
  allocate(qsumm(0:ni+1, 0:nj+1))
  allocate(qsumh(0:ni+1, 0:nj+1))

  ! Initial value arrays
  allocate(zph8s_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(qsum_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(rh24_init(0:ni+1, 0:nj+1))
  allocate(rh32_init(0:ni+1, 0:nj+1))
  allocate(rh48_init(0:ni+1, 0:nj+1))
  allocate(rh72_init(0:ni+1, 0:nj+1))
  allocate(qsuml_init(0:ni+1, 0:nj+1))
  allocate(qsumm_init(0:ni+1, 0:nj+1))
  allocate(qsumh_init(0:ni+1, 0:nj+1))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input arrays
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/zph.bin', zph, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/p.bin', p, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qall.bin', qall, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/zph8s_in.bin', zph8s_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qsum_in.bin', qsum_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/rh24_in.bin', rh24_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/rh32_in.bin', rh32_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/rh48_in.bin', rh48_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/rh72_in.bin', rh72_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/qsuml_in.bin', qsuml_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/qsumm_in.bin', qsumm_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/qsumh_in.bin', qsumh_init, 0, ni+1, 0, nj+1)

  ! Read lookup tables
  call read_array_1d(trim(data_dir)//'/rcdl.bin', rcdl, 0, 100)
  call read_array_1d(trim(data_dir)//'/rcdm.bin', rcdm, 0, 100)
  call read_array_1d(trim(data_dir)//'/rcdh.bin', rcdh, 0, 100)

  !---------------------------------------------------------------------
  ! Read reference output for validation
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_2d(trim(data_dir)//'/cdl_ref.bin', cdl_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/cdm_ref.bin', cdm_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/cdh_ref.bin', cdh_ref, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call reset_inout_arrays()
    call kernel_cloudcov(fmois, fproc, cphopt, ni, nj, nk, dz, es0iv2, &
                         zph, rst, p, t, qv, qall, cdl, cdm, cdh, &
                         zph8s, rh24, rh32, rh48, rh72, qsum, qsuml, qsumm, qsumh, &
                         rcdl, rcdm, rcdh)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  do iter = 1, num_iterations
    call reset_inout_arrays()

    !$acc wait
    t_start = omp_get_wtime()

    call kernel_cloudcov(fmois, fproc, cphopt, ni, nj, nk, dz, es0iv2, &
                         zph, rst, p, t, qv, qall, cdl, cdm, cdh, &
                         zph8s, rh24, rh32, rh48, rh72, qsum, qsuml, qsumm, qsumh, &
                         rcdl, rcdm, rcdh)

    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate results
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  call validate_results(cdl, cdl_ref, cdm, cdm_ref, cdh, cdh_ref, ni, nj, tolerance, &
                        errors, max_error_cdl, max_error_cdm, max_error_cdh)

  !---------------------------------------------------------------------
  ! Output results
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
  write(*,'(A,ES12.4)') ' Max error (cdl): ', max_error_cdl
  write(*,'(A,ES12.4)') ' Max error (cdm): ', max_error_cdm
  write(*,'(A,ES12.4)') ' Max error (cdh): ', max_error_cdh
  write(*,'(A,ES12.4)') ' Tolerance:       ', tolerance
  write(*,'(A,I12)')    ' Error count:     ', errors
  if (errors == 0) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(zph, rst, p, t, qv, qall)
  deallocate(cdl, cdm, cdh, cdl_ref, cdm_ref, cdh_ref)
  deallocate(zph8s, qsum, rh24, rh32, rh48, rh72, qsuml, qsumm, qsumh)
  deallocate(zph8s_init, qsum_init)
  deallocate(rh24_init, rh32_init, rh48_init, rh72_init)
  deallocate(qsuml_init, qsumm_init, qsumh_init)
  deallocate(times)

  if (errors /= 0) stop 1

contains

  !---------------------------------------------------------------------
  ! Reset inout arrays to initial values
  !---------------------------------------------------------------------
  subroutine reset_inout_arrays()
    zph8s = zph8s_init
    qsum = qsum_init
    rh24 = rh24_init
    rh32 = rh32_init
    rh48 = rh48_init
    rh72 = rh72_init
    qsuml = qsuml_init
    qsumm = qsumm_init
    qsumh = qsumh_init
  end subroutine reset_inout_arrays

  !---------------------------------------------------------------------
  ! Kernel subroutine: cloudcov (OpenACC version)
  !---------------------------------------------------------------------
  subroutine kernel_cloudcov(fmois, fproc, cphopt, ni, nj, nk, dz, es0iv2, &
                             zph, rst, p, t, qv, qall, cdl, cdm, cdh, &
                             zph8s, rh24, rh32, rh48, rh72, qsum, qsuml, qsumm, qsumh, &
                             rcdl, rcdm, rcdh)
    implicit none

    character(len=5), intent(in) :: fmois
    character(len=3), intent(in) :: fproc
    integer, intent(in) :: cphopt, ni, nj, nk
    real, intent(in) :: dz, es0iv2
    real, intent(in) :: zph(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: p(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: t(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qall(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rcdl(0:100), rcdm(0:100), rcdh(0:100)
    real, intent(out) :: cdl(0:ni+1, 0:nj+1)
    real, intent(out) :: cdm(0:ni+1, 0:nj+1)
    real, intent(out) :: cdh(0:ni+1, 0:nj+1)
    real, intent(inout) :: zph8s(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: rh24(0:ni+1, 0:nj+1)
    real, intent(inout) :: rh32(0:ni+1, 0:nj+1)
    real, intent(inout) :: rh48(0:ni+1, 0:nj+1)
    real, intent(inout) :: rh72(0:ni+1, 0:nj+1)
    real, intent(inout) :: qsum(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: qsuml(0:ni+1, 0:nj+1)
    real, intent(inout) :: qsumm(0:ni+1, 0:nj+1)
    real, intent(inout) :: qsumh(0:ni+1, 0:nj+1)

    integer :: i, j, k, irh
    real :: rhsfc, rha, rhb, dk

    if (fmois(1:3) .eq. 'dry') then

      !$acc kernels
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          cdl(i,j) = 0.e0
          cdm(i,j) = 0.e0
          cdh(i,j) = 0.e0
        end do
      end do
      !$acc end kernels

    else if (fmois(1:5) .eq. 'moist') then

      ! Get z physical coordinates at scalar points
      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            zph8s(i,j,k) = .5e0 * (zph(i,j,k) + zph(i,j,k+1))
          end do
        end do
      end do
      !$acc end kernels

      if (fproc(1:2) .eq. 'rh' .or. abs(cphopt) .eq. 0) then

        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(rhsfc)
          do i = 1, ni-1
            if (t(i,j,2) .gt. tlow) then
              rhsfc = es0iv2 * qv(i,j,2) * p(i,j,2) / (epsva + qv(i,j,2)) &
                    * exp(17.269e0 * (t(i,j,2) - t0) / (35.86e0 - t(i,j,2)))
            else
              rhsfc = es0iv2 * qv(i,j,2) * p(i,j,2) / (epsva + qv(i,j,2)) &
                    * exp(21.875e0 * (t(i,j,2) - t0) / (7.66e0 - t(i,j,2)))
            end if
            rh24(i,j) = rhsfc
            rh32(i,j) = rhsfc
            rh48(i,j) = rhsfc
            rh72(i,j) = rhsfc
          end do
        end do
        !$acc end kernels

        do k = 2, nk-1
          !$acc kernels
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent private(rha, rhb, dk)
            do i = 1, ni-1
              ! Interpolate at 2400m
              if (zph8s(i,j,k-1) .lt. 2400.e0 .and. zph8s(i,j,k) .ge. 2400.e0) then
                if (t(i,j,k) .gt. tlow) then
                  rha = es0iv2 * qv(i,j,k) * p(i,j,k) / (epsva + qv(i,j,k)) &
                      * exp(17.269e0 * (t(i,j,k) - t0) / (35.86e0 - t(i,j,k)))
                else
                  rha = es0iv2 * qv(i,j,k) * p(i,j,k) / (epsva + qv(i,j,k)) &
                      * exp(21.875e0 * (t(i,j,k) - t0) / (7.66e0 - t(i,j,k)))
                end if
                if (t(i,j,k-1) .gt. tlow) then
                  rhb = es0iv2 * qv(i,j,k-1) * p(i,j,k-1) / (epsva + qv(i,j,k-1)) &
                      * exp(17.269e0 * (t(i,j,k-1) - t0) / (35.86e0 - t(i,j,k-1)))
                else
                  rhb = es0iv2 * qv(i,j,k-1) * p(i,j,k-1) / (epsva + qv(i,j,k-1)) &
                      * exp(21.875e0 * (t(i,j,k-1) - t0) / (7.66e0 - t(i,j,k-1)))
                end if
                dk = (zph8s(i,j,k) - 2400.e0) / (zph8s(i,j,k) - zph8s(i,j,k-1))
                rh24(i,j) = rha * (1.e0 - dk) + rhb * dk
              end if

              ! Interpolate at 3200m
              if (zph8s(i,j,k-1) .lt. 3200.e0 .and. zph8s(i,j,k) .ge. 3200.e0) then
                if (t(i,j,k) .gt. tlow) then
                  rha = es0iv2 * qv(i,j,k) * p(i,j,k) / (epsva + qv(i,j,k)) &
                      * exp(17.269e0 * (t(i,j,k) - t0) / (35.86e0 - t(i,j,k)))
                else
                  rha = es0iv2 * qv(i,j,k) * p(i,j,k) / (epsva + qv(i,j,k)) &
                      * exp(21.875e0 * (t(i,j,k) - t0) / (7.66e0 - t(i,j,k)))
                end if
                if (t(i,j,k-1) .gt. tlow) then
                  rhb = es0iv2 * qv(i,j,k-1) * p(i,j,k-1) / (epsva + qv(i,j,k-1)) &
                      * exp(17.269e0 * (t(i,j,k-1) - t0) / (35.86e0 - t(i,j,k-1)))
                else
                  rhb = es0iv2 * qv(i,j,k-1) * p(i,j,k-1) / (epsva + qv(i,j,k-1)) &
                      * exp(21.875e0 * (t(i,j,k-1) - t0) / (7.66e0 - t(i,j,k-1)))
                end if
                dk = (zph8s(i,j,k) - 3200.e0) / (zph8s(i,j,k) - zph8s(i,j,k-1))
                rh32(i,j) = rha * (1.e0 - dk) + rhb * dk
              end if

              ! Interpolate at 4800m
              if (zph8s(i,j,k-1) .lt. 4800.e0 .and. zph8s(i,j,k) .ge. 4800.e0) then
                if (t(i,j,k) .gt. tlow) then
                  rha = es0iv2 * qv(i,j,k) * p(i,j,k) / (epsva + qv(i,j,k)) &
                      * exp(17.269e0 * (t(i,j,k) - t0) / (35.86e0 - t(i,j,k)))
                else
                  rha = es0iv2 * qv(i,j,k) * p(i,j,k) / (epsva + qv(i,j,k)) &
                      * exp(21.875e0 * (t(i,j,k) - t0) / (7.66e0 - t(i,j,k)))
                end if
                if (t(i,j,k-1) .gt. tlow) then
                  rhb = es0iv2 * qv(i,j,k-1) * p(i,j,k-1) / (epsva + qv(i,j,k-1)) &
                      * exp(17.269e0 * (t(i,j,k-1) - t0) / (35.86e0 - t(i,j,k-1)))
                else
                  rhb = es0iv2 * qv(i,j,k-1) * p(i,j,k-1) / (epsva + qv(i,j,k-1)) &
                      * exp(21.875e0 * (t(i,j,k-1) - t0) / (7.66e0 - t(i,j,k-1)))
                end if
                dk = (zph8s(i,j,k) - 4800.e0) / (zph8s(i,j,k) - zph8s(i,j,k-1))
                rh48(i,j) = rha * (1.e0 - dk) + rhb * dk
              end if

              ! Interpolate at 7200m
              if (zph8s(i,j,k-1) .lt. 7200.e0 .and. zph8s(i,j,k) .ge. 7200.e0) then
                if (t(i,j,k) .gt. tlow) then
                  rha = es0iv2 * qv(i,j,k) * p(i,j,k) / (epsva + qv(i,j,k)) &
                      * exp(17.269e0 * (t(i,j,k) - t0) / (35.86e0 - t(i,j,k)))
                else
                  rha = es0iv2 * qv(i,j,k) * p(i,j,k) / (epsva + qv(i,j,k)) &
                      * exp(21.875e0 * (t(i,j,k) - t0) / (7.66e0 - t(i,j,k)))
                end if
                if (t(i,j,k-1) .gt. tlow) then
                  rhb = es0iv2 * qv(i,j,k-1) * p(i,j,k-1) / (epsva + qv(i,j,k-1)) &
                      * exp(17.269e0 * (t(i,j,k-1) - t0) / (35.86e0 - t(i,j,k-1)))
                else
                  rhb = es0iv2 * qv(i,j,k-1) * p(i,j,k-1) / (epsva + qv(i,j,k-1)) &
                      * exp(21.875e0 * (t(i,j,k-1) - t0) / (7.66e0 - t(i,j,k-1)))
                end if
                dk = (zph8s(i,j,k) - 7200.e0) / (zph8s(i,j,k) - zph8s(i,j,k-1))
                rh72(i,j) = rha * (1.e0 - dk) + rhb * dk
              end if
            end do
          end do
          !$acc end kernels
        end do

        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent private(irh, dk)
          do i = 1, ni-1
            irh = min(int(rh24(i,j)), 100)
            dk = rh24(i,j) - aint(rh24(i,j))
            cdl(i,j) = ((1.e0 - dk) * rcdl(irh) + dk * rcdl(irh+1)) &
                     * min(max(3200.e0 - zph(i,j,2), 0.e0), 3200.e0) / 3200.e0

            irh = min(int(rh32(i,j)), 100)
            dk = rh32(i,j) - aint(rh32(i,j))
            cdm(i,j) = ((1.e0 - dk) * rcdm(irh) + dk * rcdm(irh+1)) &
                     * min(max(4800.e0 - zph(i,j,2), 0.e0), 1600.e0) / 1600.e0

            irh = min(int(rh48(i,j)), 100)
            dk = rh48(i,j) - aint(rh48(i,j))
            cdh(i,j) = ((1.e0 - dk) * rcdm(irh) + dk * rcdm(irh+1)) &
                     * min(max(7200.e0 - zph(i,j,2), 0.e0), 2400.e0) / 2400.e0

            cdm(i,j) = .5e0 * (cdm(i,j) + cdh(i,j))

            irh = min(int(rh72(i,j)), 100)
            dk = rh72(i,j) - aint(rh72(i,j))
            cdh(i,j) = (1.e0 - dk) * rcdh(irh) + dk * rcdh(irh+1)
          end do
        end do
        !$acc end kernels

      else if (fproc(1:3) .eq. 'mix' .and. abs(cphopt) .ne. 0) then

        !$acc kernels
        !$acc loop independent
        do k = 2, nk-1
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni-1
              qsum(i,j,k) = rst(i,j,k) * qall(i,j,k) * dz
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            qsuml(i,j) = 0.e0
            qsumm(i,j) = 0.e0
            qsumh(i,j) = 0.e0
          end do
        end do
        !$acc end kernels

        ! Accumulation loop - k dependency requires sequential k processing
        ! Each (i,j) point is independent, but k must be processed in order
        do k = 2, nk-1
          !$acc kernels
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni-1
              if (zph8s(i,j,k) .lt. 2800.e0) then
                qsuml(i,j) = qsuml(i,j) + qsum(i,j,k)
              else if (zph8s(i,j,k) .ge. 2800.e0 .and. zph8s(i,j,k) .lt. 6000.e0) then
                qsumm(i,j) = qsumm(i,j) + qsum(i,j,k)
              else
                qsumh(i,j) = qsumh(i,j) + qsum(i,j,k)
              end if
            end do
          end do
          !$acc end kernels
        end do

        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            cdl(i,j) = min(tend7 * (1.e0 - exp(-15.e0 * qsuml(i,j))), 1.e0)
            cdm(i,j) = min(tend6 * (1.e0 - exp(-15.e0 * qsumm(i,j))), 1.e0)
            cdh(i,j) = min(tend3 * (1.e0 - exp(-15.e0 * qsumh(i,j))), 1.e0)
          end do
        end do
        !$acc end kernels

      end if

    end if

  end subroutine kernel_cloudcov

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

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

  !---------------------------------------------------------------------
  ! Read parameters from dump file
  !---------------------------------------------------------------------
  subroutine read_parameters(filename)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=256) :: line, key, value_str
    integer :: ios, eq_pos

    open(unit=11, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    do
      read(11, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        value_str = adjustl(line(eq_pos+1:))

        select case (trim(key))
          case ('cphopt')
            read(value_str, *) cphopt
          case ('ni')
            read(value_str, *) ni
          case ('nj')
            read(value_str, *) nj
          case ('nk')
            read(value_str, *) nk
          case ('dz')
            read(value_str, *) dz
          case ('es0iv2')
            read(value_str, *) es0iv2
          case ('fmois')
            read(value_str, '(A)') fmois
          case ('fproc')
            read(value_str, '(A)') fproc
        end select
      end if
    end do

    close(11)
  end subroutine read_parameters

  !---------------------------------------------------------------------
  ! Read 1D real array from binary file
  !---------------------------------------------------------------------
  subroutine read_array_1d(filename, array, is, ie)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie
    real, intent(out) :: array(is:ie)
    integer :: ios

    open(unit=12, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(12) array
    close(12)
  end subroutine read_array_1d

  !---------------------------------------------------------------------
  ! Read 2D real array from binary file
  !---------------------------------------------------------------------
  subroutine read_array_2d(filename, array, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: array(is:ie, js:je)
    integer :: ios

    open(unit=12, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(12) array
    close(12)
  end subroutine read_array_2d

  !---------------------------------------------------------------------
  ! Read 3D real array from binary file
  !---------------------------------------------------------------------
  subroutine read_array_3d(filename, array, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: array(is:ie, js:je, ks:ke)
    integer :: ios

    open(unit=12, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(12) array
    close(12)
  end subroutine read_array_3d

  !---------------------------------------------------------------------
  ! Validate results against reference
  !---------------------------------------------------------------------
  subroutine validate_results(cdl, cdl_ref, cdm, cdm_ref, cdh, cdh_ref, ni, nj, tolerance, &
                              errors, max_err_cdl, max_err_cdm, max_err_cdh)
    implicit none
    integer, intent(in) :: ni, nj
    real, intent(in) :: cdl(0:ni+1, 0:nj+1), cdl_ref(0:ni+1, 0:nj+1)
    real, intent(in) :: cdm(0:ni+1, 0:nj+1), cdm_ref(0:ni+1, 0:nj+1)
    real, intent(in) :: cdh(0:ni+1, 0:nj+1), cdh_ref(0:ni+1, 0:nj+1)
    real, intent(in) :: tolerance
    integer, intent(out) :: errors
    real, intent(out) :: max_err_cdl, max_err_cdm, max_err_cdh

    integer :: i, j, max_i_cdm, max_j_cdm, max_i_cdl, max_j_cdl
    real :: rel_err, denom

    errors = 0
    max_err_cdl = 0.0
    max_err_cdm = 0.0
    max_err_cdh = 0.0
    max_i_cdm = 1
    max_j_cdm = 1
    max_i_cdl = 1
    max_j_cdl = 1
    do j = 1, nj-1
      do i = 1, ni-1
        denom = max(abs(cdl_ref(i,j)), 1.0e-20)
        rel_err = abs(cdl(i,j) - cdl_ref(i,j)) / denom
        if (rel_err > max_err_cdl) then
          max_err_cdl = rel_err
          max_i_cdl = i
          max_j_cdl = j
        end if
        if (rel_err > tolerance) errors = errors + 1

        denom = max(abs(cdm_ref(i,j)), 1.0e-20)
        rel_err = abs(cdm(i,j) - cdm_ref(i,j)) / denom
        if (rel_err > max_err_cdm) then
          max_err_cdm = rel_err
          max_i_cdm = i
          max_j_cdm = j
        end if
        if (rel_err > tolerance) errors = errors + 1

        denom = max(abs(cdh_ref(i,j)), 1.0e-20)
        rel_err = abs(cdh(i,j) - cdh_ref(i,j)) / denom
        max_err_cdh = max(max_err_cdh, rel_err)
        if (rel_err > tolerance) errors = errors + 1
      end do
    end do
    print '(A,I0,A,I0)', 'Max cdl error at: (', max_i_cdl, ',', max_j_cdl, ')'
    print '(A,ES15.8,A,ES15.8)', 'cdl GPU=', cdl(max_i_cdl, max_j_cdl), ' ref=', cdl_ref(max_i_cdl, max_j_cdl)
    print '(A,I0,A,I0)', 'Max cdm error at: (', max_i_cdm, ',', max_j_cdm, ')'
    print '(A,ES15.8,A,ES15.8)', 'cdm GPU=', cdm(max_i_cdm, max_j_cdm), ' ref=', cdm_ref(max_i_cdm, max_j_cdm)
  end subroutine validate_results

end program kernel_benchmark_gpu_cloudcov

!***********************************************************************
! Kernel Benchmark: exbcpt (s_exbcpt) - GPU Version
!***********************************************************************
!
! Source: Src/exbcpt.f90
! Description: Force lateral boundary value to external boundary value
!              for potential temperature
!
! This is a boundary condition kernel that operates on domain edges.
! For single-process execution, all boundary flags are set accordingly.
!
!***********************************************************************
program kernel_benchmark_exbcpt
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters from params.txt
  character(len=108) :: exbvar
  integer :: wbc, ebc, advopt, ivstp
  real :: exnews, dt, gtinc

  ! Derived parameters
  integer :: nim1, nim2, njm1, njm2
  real :: dt2, gtinc1, gtinc2, dmpdt, tpdt

  ! MPI-related parameters (for single process)
  integer :: ebw, ebe, ebs, ebn
  integer :: isub, jsub, nisub, njsub

  ! Input arrays
  real, allocatable :: ptp(:,:,:)
  real, allocatable :: ptpp(:,:,:)
  real, allocatable :: ptcpx(:,:,:)
  real, allocatable :: ptcpy(:,:,:)
  real, allocatable :: ptpgpv(:,:,:)
  real, allocatable :: ptptd(:,:,:)

  ! Output array
  real, allocatable :: ptpf(:,:,:)
  real, allocatable :: ptpf_init(:,:,:)

  ! Reference output for validation
  real, allocatable :: ptpf_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
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
  call read_parameters(trim(data_dir)//'/params.txt')

  ! Set derived parameters
  nim1 = ni - 1
  nim2 = ni - 2
  njm1 = nj - 1
  njm2 = nj - 2
  gtinc1 = gtinc
  gtinc2 = gtinc + dt
  dt2 = 2.0e0 * dt
  if (advopt <= 3) then
    dmpdt = exnews * dt2
  else
    dmpdt = exnews * dt
  end if
  tpdt = gtinc + real(ivstp - 1) * dt

  ! Set MPI parameters for single process execution
  ebw = 1
  ebe = 1
  ebs = 1
  ebn = 1
  isub = 0
  jsub = 0
  nisub = 1
  njsub = 1

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: exbcpt (GPU)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,A)') ' exbvar: ', trim(exbvar)
  write(*,'(A,I3,A,I3,A,I3)') ' wbc=', wbc, ', ebc=', ebc, ', advopt=', advopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptcpx(1:nj, 1:nk, 1:2))
  allocate(ptcpy(1:ni, 1:nk, 1:2))
  allocate(ptpgpv(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptptd(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpf_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptpp.bin', ptpp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptcpx.bin', ptcpx, 1, nj, 1, nk, 1, 2)
  call read_array_3d(trim(data_dir)//'/ptcpy.bin', ptcpy, 1, ni, 1, nk, 1, 2)
  call read_array_3d(trim(data_dir)//'/ptpgpv.bin', ptpgpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptptd.bin', ptptd, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptpf_in.bin', ptpf_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ptpf_ref.bin', ptpf_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ptpf = ptpf_init
    call kernel_exbcpt(exbvar, wbc, ebc, advopt, ni, nj, nk, &
         nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
         isub, jsub, nisub, njsub, dt, dt2, gtinc1, gtinc2, &
         dmpdt, tpdt, ptp, ptpp, ptcpx, ptcpy, ptpgpv, ptptd, ptpf)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    ptpf = ptpf_init
    !$acc wait
    t_start = omp_get_wtime()

    call kernel_exbcpt(exbvar, wbc, ebc, advopt, ni, nj, nk, &
         nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
         isub, jsub, nisub, njsub, dt, dt2, gtinc1, gtinc2, &
         dmpdt, tpdt, ptp, ptpp, ptcpx, ptcpy, ptpgpv, ptptd, ptpf)

    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
    if (times(iter) < t_min) t_min = times(iter)
    if (times(iter) > t_max) t_max = times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  call validate_output(ptpf, ptpf_ref, ni, nj, nk, tolerance, max_error, error_count)

  validation_passed = (error_count == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:  ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:      ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:      ', t_max * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:    ', t_total * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)')    ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(ptp, ptpp, ptcpx, ptcpy, ptpgpv, ptptd)
  deallocate(ptpf, ptpf_init, ptpf_ref)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: exbcpt (GPU version with OpenACC)
  !-------------------------------------------------------------------
  subroutine kernel_exbcpt(exbvar, wbc, ebc, advopt, ni, nj, nk, &
       nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
       isub, jsub, nisub, njsub, dt, dt2, gtinc1, gtinc2, &
       dmpdt, tpdt, ptp, ptpp, ptcpx, ptcpy, ptpgpv, ptptd, ptpf)
    implicit none

    character(len=*), intent(in) :: exbvar
    integer, intent(in) :: wbc, ebc, advopt
    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: nim1, nim2, njm1, njm2
    integer, intent(in) :: ebw, ebe, ebs, ebn
    integer, intent(in) :: isub, jsub, nisub, njsub
    real, intent(in) :: dt, dt2, gtinc1, gtinc2, dmpdt, tpdt
    real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ptcpx(1:nj,1:nk,1:2)
    real, intent(in) :: ptcpy(1:ni,1:nk,1:2)
    real, intent(in) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ptptd(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)

    integer :: i, j, k
    real :: ptb1, ptb2, ptb2i, ptb2j
    real :: gamma, radwe, radsn

    ! Southwest corner
    if (ebs == 1 .and. jsub == 0) then
      if (exbvar(5:5) == '-') then
        if (abs(wbc) /= 1 .and. abs(ebc) /= 1) then
          if (advopt <= 3) then
            if (ebw == 1 .and. isub == 0) then
              !$acc kernels
              !$acc loop independent private(ptb1,ptb2i,ptb2j,radwe,radsn)
              do k = 2, nk-2
                ptb1 = ptpgpv(1,1,k) + ptptd(1,1,k) * gtinc1
                ptb2i = ptpgpv(2,1,k) + ptptd(2,1,k) * gtinc2
                ptb2j = ptpgpv(1,2,k) + ptptd(1,2,k) * gtinc2
                radwe = ((ptp(2,1,k) - ptpp(1,1,k)) - (ptb2i - ptb1)) &
                     * ptcpx(1,k,1) / (1.0e0 - ptcpx(1,k,1))
                radsn = ((ptp(1,2,k) - ptpp(1,1,k)) - (ptb2j - ptb1)) &
                     * ptcpy(1,k,1) / (1.0e0 - ptcpy(1,k,1))
                ptpf(1,1,k) = ptpp(1,1,k) + ptptd(1,1,k) * dt2 &
                     - 2.0e0 * (radwe + radsn) - dmpdt * (ptpp(1,1,k) - ptb1)
              end do
              !$acc end kernels
            end if
            if (ebe == 1 .and. isub == nisub-1) then
              !$acc kernels
              !$acc loop independent private(ptb1,ptb2i,ptb2j,radwe,radsn)
              do k = 2, nk-2
                ptb1 = ptpgpv(nim1,1,k) + ptptd(nim1,1,k) * gtinc1
                ptb2i = ptpgpv(nim2,1,k) + ptptd(nim2,1,k) * gtinc2
                ptb2j = ptpgpv(nim1,2,k) + ptptd(nim1,2,k) * gtinc2
                radwe = ((ptp(nim2,1,k) - ptpp(nim1,1,k)) - (ptb2i - ptb1)) &
                     * ptcpx(1,k,2) / (1.0e0 + ptcpx(1,k,2))
                radsn = ((ptp(nim1,2,k) - ptpp(nim1,1,k)) - (ptb2j - ptb1)) &
                     * ptcpy(nim1,k,1) / (1.0e0 - ptcpy(nim1,k,1))
                ptpf(nim1,1,k) = ptpp(nim1,1,k) + ptptd(nim1,1,k) * dt2 &
                     + 2.0e0 * (radwe - radsn) - dmpdt * (ptpp(nim1,1,k) - ptb1)
              end do
              !$acc end kernels
            end if
          else
            if (ebw == 1 .and. isub == 0) then
              !$acc kernels
              !$acc loop independent private(ptb1,ptb2i,ptb2j,radwe,radsn)
              do k = 2, nk-2
                ptb1 = ptpgpv(1,1,k) + ptptd(1,1,k) * tpdt
                ptb2i = ptpgpv(2,1,k) + ptptd(2,1,k) * tpdt
                ptb2j = ptpgpv(1,2,k) + ptptd(1,2,k) * tpdt
                radwe = ((ptpp(2,1,k) - ptpp(1,1,k)) - (ptb2i - ptb1)) * ptcpx(1,k,1)
                radsn = ((ptpp(1,2,k) - ptpp(1,1,k)) - (ptb2j - ptb1)) * ptcpy(1,k,1)
                ptpf(1,1,k) = ptpp(1,1,k) + ptptd(1,1,k) * dt &
                     - (radwe + radsn) - dmpdt * (ptpp(1,1,k) - ptb1)
              end do
              !$acc end kernels
            end if
            if (ebe == 1 .and. isub == nisub-1) then
              !$acc kernels
              !$acc loop independent private(ptb1,ptb2i,ptb2j,radwe,radsn)
              do k = 2, nk-2
                ptb1 = ptpgpv(nim1,1,k) + ptptd(nim1,1,k) * tpdt
                ptb2i = ptpgpv(nim2,1,k) + ptptd(nim2,1,k) * tpdt
                ptb2j = ptpgpv(nim1,2,k) + ptptd(nim1,2,k) * tpdt
                radwe = ((ptpp(nim2,1,k) - ptpp(nim1,1,k)) - (ptb2i - ptb1)) * ptcpx(1,k,2)
                radsn = ((ptpp(nim1,2,k) - ptpp(nim1,1,k)) - (ptb2j - ptb1)) * ptcpy(nim1,k,1)
                ptpf(nim1,1,k) = ptpp(nim1,1,k) + ptptd(nim1,1,k) * dt &
                     + (radwe - radsn) - dmpdt * (ptpp(nim1,1,k) - ptb1)
              end do
              !$acc end kernels
            end if
          end if
        end if
      end if
    end if

    ! Northwest corner
    if (ebn == 1 .and. jsub == njsub-1) then
      if (exbvar(5:5) == '-') then
        if (abs(wbc) /= 1 .and. abs(ebc) /= 1) then
          if (advopt <= 3) then
            if (ebw == 1 .and. isub == 0) then
              !$acc kernels
              !$acc loop independent private(ptb1,ptb2i,ptb2j,radwe,radsn)
              do k = 2, nk-2
                ptb1 = ptpgpv(1,njm1,k) + ptptd(1,njm1,k) * gtinc1
                ptb2i = ptpgpv(2,njm1,k) + ptptd(2,njm1,k) * gtinc2
                ptb2j = ptpgpv(1,njm2,k) + ptptd(1,njm2,k) * gtinc2
                radwe = ((ptp(2,njm1,k) - ptpp(1,njm1,k)) - (ptb2i - ptb1)) &
                     * ptcpx(njm1,k,1) / (1.0e0 - ptcpx(njm1,k,1))
                radsn = ((ptp(1,njm2,k) - ptpp(1,njm1,k)) - (ptb2j - ptb1)) &
                     * ptcpy(1,k,2) / (1.0e0 + ptcpy(1,k,2))
                ptpf(1,njm1,k) = ptpp(1,njm1,k) + ptptd(1,njm1,k) * dt2 &
                     - 2.0e0 * (radwe - radsn) - dmpdt * (ptpp(1,njm1,k) - ptb1)
              end do
              !$acc end kernels
            end if
            if (ebe == 1 .and. isub == nisub-1) then
              !$acc kernels
              !$acc loop independent private(ptb1,ptb2i,ptb2j,radwe,radsn)
              do k = 2, nk-2
                ptb1 = ptpgpv(nim1,njm1,k) + ptptd(nim1,njm1,k) * gtinc1
                ptb2i = ptpgpv(nim2,njm1,k) + ptptd(nim2,njm1,k) * gtinc2
                ptb2j = ptpgpv(nim1,njm2,k) + ptptd(nim1,njm2,k) * gtinc2
                radwe = ((ptp(nim2,njm1,k) - ptpp(nim1,njm1,k)) - (ptb2i - ptb1)) &
                     * ptcpx(njm1,k,2) / (1.0e0 + ptcpx(njm1,k,2))
                radsn = ((ptp(nim1,njm2,k) - ptpp(nim1,njm1,k)) - (ptb2j - ptb1)) &
                     * ptcpy(nim1,k,2) / (1.0e0 + ptcpy(nim1,k,2))
                ptpf(nim1,njm1,k) = ptpp(nim1,njm1,k) + ptptd(nim1,njm1,k) * dt2 &
                     + 2.0e0 * (radwe + radsn) - dmpdt * (ptpp(nim1,njm1,k) - ptb1)
              end do
              !$acc end kernels
            end if
          else
            if (ebw == 1 .and. isub == 0) then
              !$acc kernels
              !$acc loop independent private(ptb1,ptb2i,ptb2j,radwe,radsn)
              do k = 2, nk-2
                ptb1 = ptpgpv(1,njm1,k) + ptptd(1,njm1,k) * tpdt
                ptb2i = ptpgpv(2,njm1,k) + ptptd(2,njm1,k) * tpdt
                ptb2j = ptpgpv(1,njm2,k) + ptptd(1,njm2,k) * tpdt
                radwe = ((ptpp(2,njm1,k) - ptpp(1,njm1,k)) - (ptb2i - ptb1)) * ptcpx(njm1,k,1)
                radsn = ((ptpp(1,njm2,k) - ptpp(1,njm1,k)) - (ptb2j - ptb1)) * ptcpy(1,k,2)
                ptpf(1,njm1,k) = ptpp(1,njm1,k) + ptptd(1,njm1,k) * dt &
                     - (radwe - radsn) - dmpdt * (ptpp(1,njm1,k) - ptb1)
              end do
              !$acc end kernels
            end if
            if (ebe == 1 .and. isub == nisub-1) then
              !$acc kernels
              !$acc loop independent private(ptb1,ptb2i,ptb2j,radwe,radsn)
              do k = 2, nk-2
                ptb1 = ptpgpv(nim1,njm1,k) + ptptd(nim1,njm1,k) * tpdt
                ptb2i = ptpgpv(nim2,njm1,k) + ptptd(nim2,njm1,k) * tpdt
                ptb2j = ptpgpv(nim1,njm2,k) + ptptd(nim1,njm2,k) * tpdt
                radwe = ((ptpp(nim2,njm1,k) - ptpp(nim1,njm1,k)) - (ptb2i - ptb1)) * ptcpx(njm1,k,2)
                radsn = ((ptpp(nim1,njm2,k) - ptpp(nim1,njm1,k)) - (ptb2j - ptb1)) * ptcpy(nim1,k,2)
                ptpf(nim1,njm1,k) = ptpp(nim1,njm1,k) + ptptd(nim1,njm1,k) * dt &
                     + (radwe + radsn) - dmpdt * (ptpp(nim1,njm1,k) - ptb1)
              end do
              !$acc end kernels
            end if
          end if
        end if
      end if
    end if

    ! West boundary
    if (ebw == 1 .and. isub == 0) then
      if (abs(wbc) /= 1) then
        if (advopt <= 3) then
          if (exbvar(5:5) == '-') then
            !$acc kernels
            !$acc loop independent collapse(2) private(ptb1,ptb2,gamma)
            do k = 2, nk-2
              do j = 2, nj-2
                gamma = 2.0e0 * ptcpx(j,k,1) / (1.0e0 - ptcpx(j,k,1))
                ptb1 = ptpgpv(1,j,k) + ptptd(1,j,k) * gtinc1
                ptb2 = ptpgpv(2,j,k) + ptptd(2,j,k) * gtinc2
                ptpf(1,j,k) = ptpp(1,j,k) + ptptd(1,j,k) * dt2 &
                     - gamma * ((ptp(2,j,k) - ptpp(1,j,k)) - (ptb2 - ptb1)) &
                     - dmpdt * (ptpp(1,j,k) - ptb1)
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent collapse(2)
            do k = 2, nk-2
              do j = 2, nj-2
                ptpf(1,j,k) = ptpp(1,j,k) + ptptd(1,j,k) * dt2
              end do
            end do
            !$acc end kernels
          end if
        else
          if (exbvar(5:5) == '-') then
            !$acc kernels
            !$acc loop independent collapse(2) private(ptb1,ptb2)
            do k = 2, nk-2
              do j = 2, nj-2
                ptb1 = ptpgpv(1,j,k) + ptptd(1,j,k) * tpdt
                ptb2 = ptpgpv(2,j,k) + ptptd(2,j,k) * tpdt
                ptpf(1,j,k) = ptpp(1,j,k) + ptptd(1,j,k) * dt &
                     - ptcpx(j,k,1) * ((ptpp(2,j,k) - ptpp(1,j,k)) - (ptb2 - ptb1)) &
                     - dmpdt * (ptpp(1,j,k) - ptb1)
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent collapse(2)
            do k = 2, nk-2
              do j = 2, nj-2
                ptpf(1,j,k) = ptpp(1,j,k) + ptptd(1,j,k) * dt
              end do
            end do
            !$acc end kernels
          end if
        end if
      end if
    end if

    ! East boundary
    if (ebe == 1 .and. isub == nisub-1) then
      if (abs(ebc) /= 1) then
        if (advopt <= 3) then
          if (exbvar(5:5) == '-') then
            !$acc kernels
            !$acc loop independent collapse(2) private(ptb1,ptb2,gamma)
            do k = 2, nk-2
              do j = 2, nj-2
                gamma = 2.0e0 * ptcpx(j,k,2) / (1.0e0 + ptcpx(j,k,2))
                ptb1 = ptpgpv(nim1,j,k) + ptptd(nim1,j,k) * gtinc1
                ptb2 = ptpgpv(nim2,j,k) + ptptd(nim2,j,k) * gtinc2
                ptpf(nim1,j,k) = ptpp(nim1,j,k) + ptptd(nim1,j,k) * dt2 &
                     + gamma * ((ptp(nim2,j,k) - ptpp(nim1,j,k)) - (ptb2 - ptb1)) &
                     - dmpdt * (ptpp(nim1,j,k) - ptb1)
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent collapse(2)
            do k = 2, nk-2
              do j = 2, nj-2
                ptpf(nim1,j,k) = ptpp(nim1,j,k) + ptptd(nim1,j,k) * dt2
              end do
            end do
            !$acc end kernels
          end if
        else
          if (exbvar(5:5) == '-') then
            !$acc kernels
            !$acc loop independent collapse(2) private(ptb1,ptb2)
            do k = 2, nk-2
              do j = 2, nj-2
                ptb1 = ptpgpv(nim1,j,k) + ptptd(nim1,j,k) * tpdt
                ptb2 = ptpgpv(nim2,j,k) + ptptd(nim2,j,k) * tpdt
                ptpf(nim1,j,k) = ptpp(nim1,j,k) + ptptd(nim1,j,k) * dt &
                     + ptcpx(j,k,2) * ((ptpp(nim2,j,k) - ptpp(nim1,j,k)) - (ptb2 - ptb1)) &
                     - dmpdt * (ptpp(nim1,j,k) - ptb1)
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent collapse(2)
            do k = 2, nk-2
              do j = 2, nj-2
                ptpf(nim1,j,k) = ptpp(nim1,j,k) + ptptd(nim1,j,k) * dt
              end do
            end do
            !$acc end kernels
          end if
        end if
      end if
    end if

    ! South boundary
    if (ebs == 1 .and. jsub == 0) then
      if (exbvar(5:5) == '-') then
        if (advopt <= 3) then
          !$acc kernels
          !$acc loop independent collapse(2) private(ptb1,ptb2,gamma)
          do k = 2, nk-2
            do i = 2, ni-2
              gamma = 2.0e0 * ptcpy(i,k,1) / (1.0e0 - ptcpy(i,k,1))
              ptb1 = ptpgpv(i,1,k) + ptptd(i,1,k) * gtinc1
              ptb2 = ptpgpv(i,2,k) + ptptd(i,2,k) * gtinc2
              ptpf(i,1,k) = ptpp(i,1,k) + ptptd(i,1,k) * dt2 &
                   - gamma * ((ptp(i,2,k) - ptpp(i,1,k)) - (ptb2 - ptb1)) &
                   - dmpdt * (ptpp(i,1,k) - ptb1)
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent collapse(2) private(ptb1,ptb2)
          do k = 2, nk-2
            do i = 2, ni-2
              ptb1 = ptpgpv(i,1,k) + ptptd(i,1,k) * tpdt
              ptb2 = ptpgpv(i,2,k) + ptptd(i,2,k) * tpdt
              ptpf(i,1,k) = ptpp(i,1,k) + ptptd(i,1,k) * dt &
                   - ptcpy(i,k,1) * ((ptpp(i,2,k) - ptpp(i,1,k)) - (ptb2 - ptb1)) &
                   - dmpdt * (ptpp(i,1,k) - ptb1)
            end do
          end do
          !$acc end kernels
        end if
      else
        if (advopt <= 3) then
          !$acc kernels
          !$acc loop independent collapse(2)
          do k = 2, nk-2
            do i = 1, ni-1
              ptpf(i,1,k) = ptpp(i,1,k) + ptptd(i,1,k) * dt2
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent collapse(2)
          do k = 2, nk-2
            do i = 1, ni-1
              ptpf(i,1,k) = ptpp(i,1,k) + ptptd(i,1,k) * dt
            end do
          end do
          !$acc end kernels
        end if
      end if
    end if

    ! North boundary
    if (ebn == 1 .and. jsub == njsub-1) then
      if (exbvar(5:5) == '-') then
        if (advopt <= 3) then
          !$acc kernels
          !$acc loop independent collapse(2) private(ptb1,ptb2,gamma)
          do k = 2, nk-2
            do i = 2, ni-2
              gamma = 2.0e0 * ptcpy(i,k,2) / (1.0e0 + ptcpy(i,k,2))
              ptb1 = ptpgpv(i,njm1,k) + ptptd(i,njm1,k) * gtinc1
              ptb2 = ptpgpv(i,njm2,k) + ptptd(i,njm2,k) * gtinc2
              ptpf(i,njm1,k) = ptpp(i,njm1,k) + ptptd(i,njm1,k) * dt2 &
                   + gamma * ((ptp(i,njm2,k) - ptpp(i,njm1,k)) - (ptb2 - ptb1)) &
                   - dmpdt * (ptpp(i,njm1,k) - ptb1)
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent collapse(2) private(ptb1,ptb2)
          do k = 2, nk-2
            do i = 2, ni-2
              ptb1 = ptpgpv(i,njm1,k) + ptptd(i,njm1,k) * tpdt
              ptb2 = ptpgpv(i,njm2,k) + ptptd(i,njm2,k) * tpdt
              ptpf(i,njm1,k) = ptpp(i,njm1,k) + ptptd(i,njm1,k) * dt &
                   + ptcpy(i,k,2) * ((ptpp(i,njm2,k) - ptpp(i,njm1,k)) - (ptb2 - ptb1)) &
                   - dmpdt * (ptpp(i,njm1,k) - ptb1)
            end do
          end do
          !$acc end kernels
        end if
      else
        if (advopt <= 3) then
          !$acc kernels
          !$acc loop independent collapse(2)
          do k = 2, nk-2
            do i = 1, ni-1
              ptpf(i,njm1,k) = ptpp(i,njm1,k) + ptptd(i,njm1,k) * dt2
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent collapse(2)
          do k = 2, nk-2
            do i = 1, ni-1
              ptpf(i,njm1,k) = ptpp(i,njm1,k) + ptptd(i,njm1,k) * dt
            end do
          end do
          !$acc end kernels
        end if
      end if
    end if

  end subroutine kernel_exbcpt

  !-------------------------------------------------------------------
  ! Read benchmark configuration
  !-------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios == 0) then
      read(10, '(A)', iostat=ios) data_dir
      read(10, *, iostat=ios) num_iter
      read(10, *, iostat=ios) warmup_iter
      read(10, *, iostat=ios) tol
      close(10)
    end if
  end subroutine read_config

  !-------------------------------------------------------------------
  ! Read parameters
  !-------------------------------------------------------------------
  subroutine read_parameters(filename)
    implicit none
    character(len=*), intent(in) :: filename

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
          case ('exbvar')
            exbvar = trim(val)
          case ('wbc')
            read(val, *) wbc
          case ('ebc')
            read(val, *) ebc
          case ('advopt')
            read(val, *) advopt
          case ('exnews')
            read(val, *) exnews
          case ('ivstp')
            read(val, *) ivstp
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('dt')
            read(val, *) dt
          case ('gtinc')
            read(val, *) gtinc
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 3D array
  !-------------------------------------------------------------------
  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

  !-------------------------------------------------------------------
  ! Validate output
  !-------------------------------------------------------------------
  subroutine validate_output(output, reference, ni, nj, nk, tol, max_err, err_count)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: output(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: reference(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    ! Validate all boundary regions
    do k = 2, nk-2
      do j = 0, nj+1
        do i = 0, ni+1
          if (abs(reference(i,j,k)) > 1.0e-30) then
            rel_err = abs(output(i,j,k) - reference(i,j,k)) / abs(reference(i,j,k))
          else
            rel_err = abs(output(i,j,k) - reference(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do
  end subroutine validate_output

end program kernel_benchmark_exbcpt

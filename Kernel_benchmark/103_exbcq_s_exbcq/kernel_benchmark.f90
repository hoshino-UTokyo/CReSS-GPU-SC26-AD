!***********************************************************************
! Kernel Benchmark: exbcq (s_exbcq)
!***********************************************************************
!
! Source: Src/exbcq.f90
! Description: Force lateral boundary value to external boundary value
!              for optional mixing ratio (water vapor, etc.)
!
!***********************************************************************
program kernel_benchmark_exbcq
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters from params.txt
  character(len=108) :: exbvar
  integer :: wbc, ebc, advopt, ape, ivstp
  real :: exnews, dt, gtinc

  ! Derived parameters
  integer :: nim1, nim2, njm1, njm2
  real :: dt2, gtinc1, gtinc2, dmpdt, tpdt

  ! MPI-related parameters (for single process)
  integer :: ebw, ebe, ebs, ebn
  integer :: isub, jsub, nisub, njsub

  ! Input arrays
  real, allocatable :: q(:,:,:)
  real, allocatable :: qp(:,:,:)
  real, allocatable :: qcpx(:,:,:)
  real, allocatable :: qcpy(:,:,:)
  real, allocatable :: qgpv(:,:,:)
  real, allocatable :: qtd(:,:,:)

  ! Output array
  real, allocatable :: qf(:,:,:)
  real, allocatable :: qf_init(:,:,:)

  ! Reference output for validation
  real, allocatable :: qf_ref(:,:,:)

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
  write(*,'(A)') ' Kernel Benchmark: exbcq'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,A)') ' exbvar: ', trim(exbvar)
  write(*,'(A,I3,A,I3,A,I3,A,I3)') ' wbc=', wbc, ', ebc=', ebc, ', advopt=', advopt, ', ape=', ape
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(q(0:ni+1, 0:nj+1, 1:nk))
  allocate(qp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qcpx(1:nj, 1:nk, 1:2))
  allocate(qcpy(1:ni, 1:nk, 1:2))
  allocate(qgpv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qtd(0:ni+1, 0:nj+1, 1:nk))
  allocate(qf(0:ni+1, 0:nj+1, 1:nk))
  allocate(qf_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(qf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/q.bin', q, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qp.bin', qp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qcpx.bin', qcpx, 1, nj, 1, nk, 1, 2)
  call read_array_3d(trim(data_dir)//'/qcpy.bin', qcpy, 1, ni, 1, nk, 1, 2)
  call read_array_3d(trim(data_dir)//'/qgpv.bin', qgpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qtd.bin', qtd, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qf_in.bin', qf_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/qf_ref.bin', qf_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    qf = qf_init
    call kernel_exbcq(exbvar, wbc, ebc, advopt, ape, ni, nj, nk, &
         nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
         isub, jsub, nisub, njsub, dt, dt2, gtinc1, gtinc2, &
         dmpdt, tpdt, q, qp, qcpx, qcpy, qgpv, qtd, qf)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    qf = qf_init
    t_start = omp_get_wtime()

    call kernel_exbcq(exbvar, wbc, ebc, advopt, ape, ni, nj, nk, &
         nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
         isub, jsub, nisub, njsub, dt, dt2, gtinc1, gtinc2, &
         dmpdt, tpdt, q, qp, qcpx, qcpy, qgpv, qtd, qf)

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
  call validate_output(qf, qf_ref, ni, nj, nk, tolerance, max_error, error_count)

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
  deallocate(q, qp, qcpx, qcpy, qgpv, qtd)
  deallocate(qf, qf_init, qf_ref)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: exbcq
  !-------------------------------------------------------------------
  subroutine kernel_exbcq(exbvar, wbc, ebc, advopt, ape, ni, nj, nk, &
       nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
       isub, jsub, nisub, njsub, dt, dt2, gtinc1, gtinc2, &
       dmpdt, tpdt, q, qp, qcpx, qcpy, qgpv, qtd, qf)
    implicit none

    character(len=*), intent(in) :: exbvar
    integer, intent(in) :: wbc, ebc, advopt, ape
    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: nim1, nim2, njm1, njm2
    integer, intent(in) :: ebw, ebe, ebs, ebn
    integer, intent(in) :: isub, jsub, nisub, njsub
    real, intent(in) :: dt, dt2, gtinc1, gtinc2, dmpdt, tpdt
    real, intent(in) :: q(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qp(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qcpx(1:nj,1:nk,1:2)
    real, intent(in) :: qcpy(1:ni,1:nk,1:2)
    real, intent(in) :: qgpv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qtd(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: qf(0:ni+1,0:nj+1,1:nk)

    integer :: i, j, k
    real :: qb1, qb2, qb2i, qb2j
    real :: gamma, radwe, radsn

    !$omp parallel default(shared)

    ! Southwest corner
    if (ebs == 1 .and. jsub == 0) then
      if (exbvar(ape:ape) == '-') then
        if (abs(wbc) /= 1 .and. abs(ebc) /= 1) then
          if (advopt <= 3) then
            if (ebw == 1 .and. isub == 0) then
              !$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)
              do k = 2, nk-2
                qb1 = max(qgpv(1,1,k) + qtd(1,1,k) * gtinc1, 0.0e0)
                qb2i = max(qgpv(2,1,k) + qtd(2,1,k) * gtinc2, 0.0e0)
                qb2j = max(qgpv(1,2,k) + qtd(1,2,k) * gtinc2, 0.0e0)
                radwe = ((q(2,1,k) - qp(1,1,k)) - (qb2i - qb1)) &
                     * qcpx(1,k,1) / (1.0e0 - qcpx(1,k,1))
                radsn = ((q(1,2,k) - qp(1,1,k)) - (qb2j - qb1)) &
                     * qcpy(1,k,1) / (1.0e0 - qcpy(1,k,1))
                qf(1,1,k) = max(qp(1,1,k) + qtd(1,1,k) * dt2 &
                     - 2.0e0 * (radwe + radsn) - dmpdt * (qp(1,1,k) - qb1), 0.0e0)
              end do
              !$omp end do
            end if
            if (ebe == 1 .and. isub == nisub-1) then
              !$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)
              do k = 2, nk-2
                qb1 = max(qgpv(nim1,1,k) + qtd(nim1,1,k) * gtinc1, 0.0e0)
                qb2i = max(qgpv(nim2,1,k) + qtd(nim2,1,k) * gtinc2, 0.0e0)
                qb2j = max(qgpv(nim1,2,k) + qtd(nim1,2,k) * gtinc2, 0.0e0)
                radwe = ((q(nim2,1,k) - qp(nim1,1,k)) - (qb2i - qb1)) &
                     * qcpx(1,k,2) / (1.0e0 + qcpx(1,k,2))
                radsn = ((q(nim1,2,k) - qp(nim1,1,k)) - (qb2j - qb1)) &
                     * qcpy(nim1,k,1) / (1.0e0 - qcpy(nim1,k,1))
                qf(nim1,1,k) = max(qp(nim1,1,k) + qtd(nim1,1,k) * dt2 &
                     + 2.0e0 * (radwe - radsn) - dmpdt * (qp(nim1,1,k) - qb1), 0.0e0)
              end do
              !$omp end do
            end if
          else
            if (ebw == 1 .and. isub == 0) then
              !$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)
              do k = 2, nk-2
                qb1 = max(qgpv(1,1,k) + qtd(1,1,k) * tpdt, 0.0e0)
                qb2i = max(qgpv(2,1,k) + qtd(2,1,k) * tpdt, 0.0e0)
                qb2j = max(qgpv(1,2,k) + qtd(1,2,k) * tpdt, 0.0e0)
                radwe = ((qp(2,1,k) - qp(1,1,k)) - (qb2i - qb1)) * qcpx(1,k,1)
                radsn = ((qp(1,2,k) - qp(1,1,k)) - (qb2j - qb1)) * qcpy(1,k,1)
                qf(1,1,k) = max(qp(1,1,k) + qtd(1,1,k) * dt &
                     - (radwe + radsn) - dmpdt * (qp(1,1,k) - qb1), 0.0e0)
              end do
              !$omp end do
            end if
            if (ebe == 1 .and. isub == nisub-1) then
              !$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)
              do k = 2, nk-2
                qb1 = max(qgpv(nim1,1,k) + qtd(nim1,1,k) * tpdt, 0.0e0)
                qb2i = max(qgpv(nim2,1,k) + qtd(nim2,1,k) * tpdt, 0.0e0)
                qb2j = max(qgpv(nim1,2,k) + qtd(nim1,2,k) * tpdt, 0.0e0)
                radwe = ((qp(nim2,1,k) - qp(nim1,1,k)) - (qb2i - qb1)) * qcpx(1,k,2)
                radsn = ((qp(nim1,2,k) - qp(nim1,1,k)) - (qb2j - qb1)) * qcpy(nim1,k,1)
                qf(nim1,1,k) = max(qp(nim1,1,k) + qtd(nim1,1,k) * dt &
                     + (radwe - radsn) - dmpdt * (qp(nim1,1,k) - qb1), 0.0e0)
              end do
              !$omp end do
            end if
          end if
        end if
      end if
    end if

    ! Northwest corner
    if (ebn == 1 .and. jsub == njsub-1) then
      if (exbvar(ape:ape) == '-') then
        if (abs(wbc) /= 1 .and. abs(ebc) /= 1) then
          if (advopt <= 3) then
            if (ebw == 1 .and. isub == 0) then
              !$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)
              do k = 2, nk-2
                qb1 = max(qgpv(1,njm1,k) + qtd(1,njm1,k) * gtinc1, 0.0e0)
                qb2i = max(qgpv(2,njm1,k) + qtd(2,njm1,k) * gtinc2, 0.0e0)
                qb2j = max(qgpv(1,njm2,k) + qtd(1,njm2,k) * gtinc2, 0.0e0)
                radwe = ((q(2,njm1,k) - qp(1,njm1,k)) - (qb2i - qb1)) &
                     * qcpx(njm1,k,1) / (1.0e0 - qcpx(njm1,k,1))
                radsn = ((q(1,njm2,k) - qp(1,njm1,k)) - (qb2j - qb1)) &
                     * qcpy(1,k,2) / (1.0e0 + qcpy(1,k,2))
                qf(1,njm1,k) = max(qp(1,njm1,k) + qtd(1,njm1,k) * dt2 &
                     - 2.0e0 * (radwe - radsn) - dmpdt * (qp(1,njm1,k) - qb1), 0.0e0)
              end do
              !$omp end do
            end if
            if (ebe == 1 .and. isub == nisub-1) then
              !$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)
              do k = 2, nk-2
                qb1 = max(qgpv(nim1,njm1,k) + qtd(nim1,njm1,k) * gtinc1, 0.0e0)
                qb2i = max(qgpv(nim2,njm1,k) + qtd(nim2,njm1,k) * gtinc2, 0.0e0)
                qb2j = max(qgpv(nim1,njm2,k) + qtd(nim1,njm2,k) * gtinc2, 0.0e0)
                radwe = ((q(nim2,njm1,k) - qp(nim1,njm1,k)) - (qb2i - qb1)) &
                     * qcpx(njm1,k,2) / (1.0e0 + qcpx(njm1,k,2))
                radsn = ((q(nim1,njm2,k) - qp(nim1,njm1,k)) - (qb2j - qb1)) &
                     * qcpy(nim1,k,2) / (1.0e0 + qcpy(nim1,k,2))
                qf(nim1,njm1,k) = max(qp(nim1,njm1,k) + qtd(nim1,njm1,k) * dt2 &
                     + 2.0e0 * (radwe + radsn) - dmpdt * (qp(nim1,njm1,k) - qb1), 0.0e0)
              end do
              !$omp end do
            end if
          else
            if (ebw == 1 .and. isub == 0) then
              !$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)
              do k = 2, nk-2
                qb1 = max(qgpv(1,njm1,k) + qtd(1,njm1,k) * tpdt, 0.0e0)
                qb2i = max(qgpv(2,njm1,k) + qtd(2,njm1,k) * tpdt, 0.0e0)
                qb2j = max(qgpv(1,njm2,k) + qtd(1,njm2,k) * tpdt, 0.0e0)
                radwe = ((qp(2,njm1,k) - qp(1,njm1,k)) - (qb2i - qb1)) * qcpx(njm1,k,1)
                radsn = ((qp(1,njm2,k) - qp(1,njm1,k)) - (qb2j - qb1)) * qcpy(1,k,2)
                qf(1,njm1,k) = max(qp(1,njm1,k) + qtd(1,njm1,k) * dt &
                     - (radwe - radsn) - dmpdt * (qp(1,njm1,k) - qb1), 0.0e0)
              end do
              !$omp end do
            end if
            if (ebe == 1 .and. isub == nisub-1) then
              !$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)
              do k = 2, nk-2
                qb1 = max(qgpv(nim1,njm1,k) + qtd(nim1,njm1,k) * tpdt, 0.0e0)
                qb2i = max(qgpv(nim2,njm1,k) + qtd(nim2,njm1,k) * tpdt, 0.0e0)
                qb2j = max(qgpv(nim1,njm2,k) + qtd(nim1,njm2,k) * tpdt, 0.0e0)
                radwe = ((qp(nim2,njm1,k) - qp(nim1,njm1,k)) - (qb2i - qb1)) * qcpx(njm1,k,2)
                radsn = ((qp(nim1,njm2,k) - qp(nim1,njm1,k)) - (qb2j - qb1)) * qcpy(nim1,k,2)
                qf(nim1,njm1,k) = max(qp(nim1,njm1,k) + qtd(nim1,njm1,k) * dt &
                     + (radwe + radsn) - dmpdt * (qp(nim1,njm1,k) - qb1), 0.0e0)
              end do
              !$omp end do
            end if
          end if
        end if
      end if
    end if

    ! West boundary
    if (ebw == 1 .and. isub == 0) then
      if (abs(wbc) /= 1) then
        if (advopt <= 3) then
          if (exbvar(ape:ape) == '-') then
            !$omp do schedule(runtime) private(j,k,qb1,qb2,gamma)
            do k = 2, nk-2
              do j = 2, nj-2
                gamma = 2.0e0 * qcpx(j,k,1) / (1.0e0 - qcpx(j,k,1))
                qb1 = max(qgpv(1,j,k) + qtd(1,j,k) * gtinc1, 0.0e0)
                qb2 = max(qgpv(2,j,k) + qtd(2,j,k) * gtinc2, 0.0e0)
                qf(1,j,k) = max(qp(1,j,k) + qtd(1,j,k) * dt2 &
                     - gamma * ((q(2,j,k) - qp(1,j,k)) - (qb2 - qb1)) &
                     - dmpdt * (qp(1,j,k) - qb1), 0.0e0)
              end do
            end do
            !$omp end do
          else
            !$omp do schedule(runtime) private(j,k)
            do k = 2, nk-2
              do j = 2, nj-2
                qf(1,j,k) = max(qp(1,j,k) + qtd(1,j,k) * dt2, 0.0e0)
              end do
            end do
            !$omp end do
          end if
        else
          if (exbvar(ape:ape) == '-') then
            !$omp do schedule(runtime) private(j,k,qb1,qb2)
            do k = 2, nk-2
              do j = 2, nj-2
                qb1 = max(qgpv(1,j,k) + qtd(1,j,k) * tpdt, 0.0e0)
                qb2 = max(qgpv(2,j,k) + qtd(2,j,k) * tpdt, 0.0e0)
                qf(1,j,k) = max(qp(1,j,k) + qtd(1,j,k) * dt &
                     - qcpx(j,k,1) * ((qp(2,j,k) - qp(1,j,k)) - (qb2 - qb1)) &
                     - dmpdt * (qp(1,j,k) - qb1), 0.0e0)
              end do
            end do
            !$omp end do
          else
            !$omp do schedule(runtime) private(j,k)
            do k = 2, nk-2
              do j = 2, nj-2
                qf(1,j,k) = max(qp(1,j,k) + qtd(1,j,k) * dt, 0.0e0)
              end do
            end do
            !$omp end do
          end if
        end if
      end if
    end if

    ! East boundary
    if (ebe == 1 .and. isub == nisub-1) then
      if (abs(ebc) /= 1) then
        if (advopt <= 3) then
          if (exbvar(ape:ape) == '-') then
            !$omp do schedule(runtime) private(j,k,qb1,qb2,gamma)
            do k = 2, nk-2
              do j = 2, nj-2
                gamma = 2.0e0 * qcpx(j,k,2) / (1.0e0 + qcpx(j,k,2))
                qb1 = max(qgpv(nim1,j,k) + qtd(nim1,j,k) * gtinc1, 0.0e0)
                qb2 = max(qgpv(nim2,j,k) + qtd(nim2,j,k) * gtinc2, 0.0e0)
                qf(nim1,j,k) = max(qp(nim1,j,k) + qtd(nim1,j,k) * dt2 &
                     + gamma * ((q(nim2,j,k) - qp(nim1,j,k)) - (qb2 - qb1)) &
                     - dmpdt * (qp(nim1,j,k) - qb1), 0.0e0)
              end do
            end do
            !$omp end do
          else
            !$omp do schedule(runtime) private(j,k)
            do k = 2, nk-2
              do j = 2, nj-2
                qf(nim1,j,k) = max(qp(nim1,j,k) + qtd(nim1,j,k) * dt2, 0.0e0)
              end do
            end do
            !$omp end do
          end if
        else
          if (exbvar(ape:ape) == '-') then
            !$omp do schedule(runtime) private(j,k,qb1,qb2)
            do k = 2, nk-2
              do j = 2, nj-2
                qb1 = max(qgpv(nim1,j,k) + qtd(nim1,j,k) * tpdt, 0.0e0)
                qb2 = max(qgpv(nim2,j,k) + qtd(nim2,j,k) * tpdt, 0.0e0)
                qf(nim1,j,k) = max(qp(nim1,j,k) + qtd(nim1,j,k) * dt &
                     + qcpx(j,k,2) * ((qp(nim2,j,k) - qp(nim1,j,k)) - (qb2 - qb1)) &
                     - dmpdt * (qp(nim1,j,k) - qb1), 0.0e0)
              end do
            end do
            !$omp end do
          else
            !$omp do schedule(runtime) private(j,k)
            do k = 2, nk-2
              do j = 2, nj-2
                qf(nim1,j,k) = max(qp(nim1,j,k) + qtd(nim1,j,k) * dt, 0.0e0)
              end do
            end do
            !$omp end do
          end if
        end if
      end if
    end if

    ! South boundary
    if (ebs == 1 .and. jsub == 0) then
      if (exbvar(ape:ape) == '-') then
        if (advopt <= 3) then
          !$omp do schedule(runtime) private(i,k,qb1,qb2,gamma)
          do k = 2, nk-2
            do i = 2, ni-2
              gamma = 2.0e0 * qcpy(i,k,1) / (1.0e0 - qcpy(i,k,1))
              qb1 = max(qgpv(i,1,k) + qtd(i,1,k) * gtinc1, 0.0e0)
              qb2 = max(qgpv(i,2,k) + qtd(i,2,k) * gtinc2, 0.0e0)
              qf(i,1,k) = max(qp(i,1,k) + qtd(i,1,k) * dt2 &
                   - gamma * ((q(i,2,k) - qp(i,1,k)) - (qb2 - qb1)) &
                   - dmpdt * (qp(i,1,k) - qb1), 0.0e0)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(i,k,qb1,qb2)
          do k = 2, nk-2
            do i = 2, ni-2
              qb1 = max(qgpv(i,1,k) + qtd(i,1,k) * tpdt, 0.0e0)
              qb2 = max(qgpv(i,2,k) + qtd(i,2,k) * tpdt, 0.0e0)
              qf(i,1,k) = max(qp(i,1,k) + qtd(i,1,k) * dt &
                   - qcpy(i,k,1) * ((qp(i,2,k) - qp(i,1,k)) - (qb2 - qb1)) &
                   - dmpdt * (qp(i,1,k) - qb1), 0.0e0)
            end do
          end do
          !$omp end do
        end if
      else
        if (advopt <= 3) then
          !$omp do schedule(runtime) private(i,k)
          do k = 2, nk-2
            do i = 1, ni-1
              qf(i,1,k) = max(qp(i,1,k) + qtd(i,1,k) * dt2, 0.0e0)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(i,k)
          do k = 2, nk-2
            do i = 1, ni-1
              qf(i,1,k) = max(qp(i,1,k) + qtd(i,1,k) * dt, 0.0e0)
            end do
          end do
          !$omp end do
        end if
      end if
    end if

    ! North boundary
    if (ebn == 1 .and. jsub == njsub-1) then
      if (exbvar(ape:ape) == '-') then
        if (advopt <= 3) then
          !$omp do schedule(runtime) private(i,k,qb1,qb2,gamma)
          do k = 2, nk-2
            do i = 2, ni-2
              gamma = 2.0e0 * qcpy(i,k,2) / (1.0e0 + qcpy(i,k,2))
              qb1 = max(qgpv(i,njm1,k) + qtd(i,njm1,k) * gtinc1, 0.0e0)
              qb2 = max(qgpv(i,njm2,k) + qtd(i,njm2,k) * gtinc2, 0.0e0)
              qf(i,njm1,k) = max(qp(i,njm1,k) + qtd(i,njm1,k) * dt2 &
                   + gamma * ((q(i,njm2,k) - qp(i,njm1,k)) - (qb2 - qb1)) &
                   - dmpdt * (qp(i,njm1,k) - qb1), 0.0e0)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(i,k,qb1,qb2)
          do k = 2, nk-2
            do i = 2, ni-2
              qb1 = max(qgpv(i,njm1,k) + qtd(i,njm1,k) * tpdt, 0.0e0)
              qb2 = max(qgpv(i,njm2,k) + qtd(i,njm2,k) * tpdt, 0.0e0)
              qf(i,njm1,k) = max(qp(i,njm1,k) + qtd(i,njm1,k) * dt &
                   + qcpy(i,k,2) * ((qp(i,njm2,k) - qp(i,njm1,k)) - (qb2 - qb1)) &
                   - dmpdt * (qp(i,njm1,k) - qb1), 0.0e0)
            end do
          end do
          !$omp end do
        end if
      else
        if (advopt <= 3) then
          !$omp do schedule(runtime) private(i,k)
          do k = 2, nk-2
            do i = 1, ni-1
              qf(i,njm1,k) = max(qp(i,njm1,k) + qtd(i,njm1,k) * dt2, 0.0e0)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(i,k)
          do k = 2, nk-2
            do i = 1, ni-1
              qf(i,njm1,k) = max(qp(i,njm1,k) + qtd(i,njm1,k) * dt, 0.0e0)
            end do
          end do
          !$omp end do
        end if
      end if
    end if

    !$omp end parallel

  end subroutine kernel_exbcq

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
          case ('ape')
            read(val, *) ape
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

end program kernel_benchmark_exbcq

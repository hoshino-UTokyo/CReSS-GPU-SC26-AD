!***********************************************************************
! GPU Kernel Benchmark: turbuvw (s_turbuvw)
!***********************************************************************
!
! Source: Src/turbuvw.f90
! Description: Calculate velocity turbulent mixing for u, v, w components
! GPU Port: OpenACC with Unified Memory
!
! NOTE: Option values are read from params.txt (dumped from simulation)
!
!***********************************************************************
program kernel_benchmark_turbuvw
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Option parameters (read from params.txt)
  integer :: trnopt, mpopt, mfcopt, advopt
  real :: dxiv, dyiv, dziv

  ! Derived parameters
  real :: dxiv05, dyiv05, dxiv25, dyiv25

  ! Input arrays
  real, allocatable :: j31(:,:,:)
  real, allocatable :: j32(:,:,:)
  real, allocatable :: jcb(:,:,:)
  real, allocatable :: jcb8u(:,:,:)
  real, allocatable :: jcb8v(:,:,:)
  real, allocatable :: mf(:,:)
  real, allocatable :: mf8u(:,:)
  real, allocatable :: mf8v(:,:)
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: rmf8u(:,:,:)
  real, allocatable :: rmf8v(:,:,:)

  ! Input/Output arrays
  real, allocatable :: t11(:,:,:)
  real, allocatable :: t22(:,:,:)
  real, allocatable :: t33(:,:,:)
  real, allocatable :: t12(:,:,:)
  real, allocatable :: t13(:,:,:)
  real, allocatable :: t23(:,:,:)
  real, allocatable :: t31(:,:,:)
  real, allocatable :: t32(:,:,:)
  real, allocatable :: ufrc(:,:,:)
  real, allocatable :: vfrc(:,:,:)
  real, allocatable :: wfrc(:,:,:)
  real, allocatable :: tmp1(:,:,:)

  ! Initial values for re-initialization
  real, allocatable :: t11_init(:,:,:), t22_init(:,:,:), t33_init(:,:,:)
  real, allocatable :: t12_init(:,:,:), t13_init(:,:,:), t23_init(:,:,:)
  real, allocatable :: t31_init(:,:,:), t32_init(:,:,:)
  real, allocatable :: ufrc_init(:,:,:), vfrc_init(:,:,:), wfrc_init(:,:,:)
  real, allocatable :: tmp1_init(:,:,:)

  ! Reference output
  real, allocatable :: ufrc_ref(:,:,:), vfrc_ref(:,:,:), wfrc_ref(:,:,:)

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
  ! Read parameters (including option values)
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', &
       trnopt, mpopt, mfcopt, advopt, dxiv, dyiv, dziv, ni, nj, nk)

  !---------------------------------------------------------------------
  ! Initialize derived parameters
  !---------------------------------------------------------------------
  dxiv05 = 0.5e0 * dxiv
  dyiv05 = 0.5e0 * dyiv
  dxiv25 = 0.25e0 * dxiv
  dyiv25 = 0.25e0 * dyiv

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: turbuvw'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I3,A,I3,A,I3,A,I3)') ' Options: trnopt=', trnopt, &
       ', mpopt=', mpopt, ', mfcopt=', mfcopt, ', advopt=', advopt
  write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') ' dxiv=', dxiv, ', dyiv=', dyiv, ', dziv=', dziv
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(mf8u(0:ni+1, 0:nj+1))
  allocate(mf8v(0:ni+1, 0:nj+1))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rmf8u(0:ni+1, 0:nj+1, 1:3))
  allocate(rmf8v(0:ni+1, 0:nj+1, 1:3))

  allocate(t11(0:ni+1, 0:nj+1, 1:nk))
  allocate(t22(0:ni+1, 0:nj+1, 1:nk))
  allocate(t33(0:ni+1, 0:nj+1, 1:nk))
  allocate(t12(0:ni+1, 0:nj+1, 1:nk))
  allocate(t13(0:ni+1, 0:nj+1, 1:nk))
  allocate(t23(0:ni+1, 0:nj+1, 1:nk))
  allocate(t31(0:ni+1, 0:nj+1, 1:nk))
  allocate(t32(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))

  allocate(t11_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t22_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t33_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t12_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t13_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t23_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t31_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(t32_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_init(0:ni+1, 0:nj+1, 1:nk))

  allocate(ufrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc_ref(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8u.bin', jcb8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8v.bin', jcb8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/mf.bin', mf, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/mf8u.bin', mf8u, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/mf8v.bin', mf8v, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call read_array_3d(trim(data_dir)//'/rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call read_array_3d(trim(data_dir)//'/rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)

  call read_array_3d(trim(data_dir)//'/t11_in.bin', t11_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t22_in.bin', t22_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t33_in.bin', t33_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t12_in.bin', t12_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t13_in.bin', t13_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t23_in.bin', t23_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t31_in.bin', t31_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t32_in.bin', t32_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ufrc_in.bin', ufrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_in.bin', vfrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wfrc_in.bin', wfrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp1_in.bin', tmp1_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ufrc_ref.bin', ufrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_ref.bin', vfrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wfrc_ref.bin', wfrc_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    t11 = t11_init; t22 = t22_init; t33 = t33_init
    t12 = t12_init; t13 = t13_init; t23 = t23_init
    t31 = t31_init; t32 = t32_init
    ufrc = ufrc_init; vfrc = vfrc_init; wfrc = wfrc_init
    tmp1 = tmp1_init

    call kernel_turbuvw(ni, nj, nk, j31, j32, jcb, jcb8u, jcb8v, &
         mf, mf8u, mf8v, rmf, rmf8u, rmf8v, &
         t11, t22, t33, t12, t13, t23, t31, t32, &
         ufrc, vfrc, wfrc, tmp1, &
         trnopt, mpopt, mfcopt, advopt, &
         dxiv, dyiv, dziv, dxiv05, dyiv05, dxiv25, dyiv25)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    t11 = t11_init; t22 = t22_init; t33 = t33_init
    t12 = t12_init; t13 = t13_init; t23 = t23_init
    t31 = t31_init; t32 = t32_init
    ufrc = ufrc_init; vfrc = vfrc_init; wfrc = wfrc_init
    tmp1 = tmp1_init

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_turbuvw(ni, nj, nk, j31, j32, jcb, jcb8u, jcb8v, &
         mf, mf8u, mf8v, rmf, rmf8u, rmf8v, &
         t11, t22, t33, t12, t13, t23, t31, t32, &
         ufrc, vfrc, wfrc, tmp1, &
         trnopt, mpopt, mfcopt, advopt, &
         dxiv, dyiv, dziv, dxiv05, dyiv05, dxiv25, dyiv25)
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

  ! Validate ufrc
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-1
        rel_error = abs(ufrc(i,j,k) - ufrc_ref(i,j,k))
        if (abs(ufrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(ufrc_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate vfrc
  do k = 2, nk-2
    do j = 2, nj-1
      do i = 2, ni-2
        rel_error = abs(vfrc(i,j,k) - vfrc_ref(i,j,k))
        if (abs(vfrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(vfrc_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate wfrc
  do k = 2, nk-1
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(wfrc(i,j,k) - wfrc_ref(i,j,k))
        if (abs(wfrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(wfrc_ref(i,j,k))
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

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: turbuvw (full implementation with all code paths)
  ! GPU Port: OpenACC with Unified Memory
  !=====================================================================
  subroutine kernel_turbuvw(ni, nj, nk, j31, j32, jcb, jcb8u, jcb8v, &
       mf, mf8u, mf8v, rmf, rmf8u, rmf8v, &
       t11, t22, t33, t12, t13, t23, t31, t32, &
       ufrc, vfrc, wfrc, tmp1, &
       trnopt, mpopt, mfcopt, advopt, &
       dxiv, dyiv, dziv, dxiv05, dyiv05, dxiv25, dyiv25)
    implicit none

    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: j31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: mf(0:ni+1, 0:nj+1)
    real, intent(in) :: mf8u(0:ni+1, 0:nj+1)
    real, intent(in) :: mf8v(0:ni+1, 0:nj+1)
    real, intent(in) :: rmf(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: rmf8u(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: rmf8v(0:ni+1, 0:nj+1, 1:3)
    real, intent(inout) :: t11(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t22(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t33(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t12(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t13(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t23(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: t32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: ufrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: wfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    integer, intent(in) :: trnopt, mpopt, mfcopt, advopt
    real, intent(in) :: dxiv, dyiv, dziv, dxiv05, dyiv05, dxiv25, dyiv25

    integer :: i, j, k

    ! Calculate the inverse of map scale factor at dot points
    if (mfcopt == 1) then
      if (mpopt == 0 .or. mpopt == 10) then
        !$acc kernels
        !$acc loop independent
        do j = 2, nj-1
          !$acc loop independent
          do i = 2, ni-1
            t13(i,j,nk) = rmf8u(i,j-1,2) + rmf8u(i,j,2)
          end do
        end do
        !$acc end kernels
      else if (mpopt == 5) then
        !$acc kernels
        !$acc loop independent
        do j = 2, nj-1
          !$acc loop independent
          do i = 2, ni-1
            t23(i,j,nk) = rmf8v(i-1,j,2) + rmf8v(i,j,2)
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do j = 2, nj-1
          !$acc loop independent
          do i = 2, ni-1
            t13(i,j,nk) = rmf8u(i,j-1,2) + rmf8u(i,j,2)
            t23(i,j,nk) = rmf8v(i-1,j,2) + rmf8v(i,j,2)
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! Calculate the u turbulent mixing
    if (trnopt == 0) then
      if (mfcopt == 0) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 1, ni-1
              t11(i,j,k) = jcb(i,j,k) * t11(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-1
            !$acc loop independent
            do i = 2, ni-1
              tmp1(i,j,k) = (jcb8u(i,j-1,k) + jcb8u(i,j,k)) * t12(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

        if (advopt <= 3) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 2, nj-2
              !$acc loop independent
              do i = 2, ni-1
                ufrc(i,j,k) = ((t11(i,j,k) - t11(i-1,j,k)) * dxiv &
                     + (tmp1(i,j+1,k) - tmp1(i,j,k)) * dyiv05) &
                     + (t13(i,j,k+1) - t13(i,j,k)) * dziv
              end do
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 2, nj-2
              !$acc loop independent
              do i = 2, ni-1
                ufrc(i,j,k) = ufrc(i,j,k) + ((t11(i,j,k) - t11(i-1,j,k)) * dxiv &
                     + (tmp1(i,j+1,k) - tmp1(i,j,k)) * dyiv05) &
                     + (t13(i,j,k+1) - t13(i,j,k)) * dziv
              end do
            end do
          end do
          !$acc end kernels
        end if
      else
        if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then
          if (mpopt == 0 .or. mpopt == 10) then
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 1, ni-1
                  t11(i,j,k) = jcb(i,j,k) * t11(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-1
                  tmp1(i,j,k) = t13(i,j,nk) * (jcb8u(i,j-1,k) + jcb8u(i,j,k)) * t12(i,j,k)
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 1, ni-1
                  t11(i,j,k) = rmf(i,j,2) * jcb(i,j,k) * t11(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-1
                  tmp1(i,j,k) = (jcb8u(i,j-1,k) + jcb8u(i,j,k)) * t12(i,j,k)
                end do
              end do
            end do
            !$acc end kernels
          end if

          if (advopt <= 3) then
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-1
                  ufrc(i,j,k) = mf8u(i,j) * ((t11(i,j,k) - t11(i-1,j,k)) * dxiv &
                       + (tmp1(i,j+1,k) - tmp1(i,j,k)) * dyiv25) &
                       + (t13(i,j,k+1) - t13(i,j,k)) * dziv
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-1
                  ufrc(i,j,k) = ufrc(i,j,k) &
                       + mf8u(i,j) * ((t11(i,j,k) - t11(i-1,j,k)) * dxiv &
                       + (tmp1(i,j+1,k) - tmp1(i,j,k)) * dyiv25) &
                       + (t13(i,j,k+1) - t13(i,j,k)) * dziv
                end do
              end do
            end do
            !$acc end kernels
          end if
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 2, nj-2
              !$acc loop independent
              do i = 1, ni-1
                t11(i,j,k) = rmf(i,j,2) * jcb(i,j,k) * t11(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 2, nj-1
              !$acc loop independent
              do i = 2, ni-1
                tmp1(i,j,k) = t13(i,j,nk) * (jcb8u(i,j-1,k) + jcb8u(i,j,k)) * t12(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

          if (advopt <= 3) then
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-1
                  ufrc(i,j,k) = rmf8u(i,j,1) * ((t11(i,j,k) - t11(i-1,j,k)) * dxiv &
                       + (tmp1(i,j+1,k) - tmp1(i,j,k)) * dyiv25) &
                       + (t13(i,j,k+1) - t13(i,j,k)) * dziv
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-1
                  ufrc(i,j,k) = ufrc(i,j,k) &
                       + rmf8u(i,j,1) * ((t11(i,j,k) - t11(i-1,j,k)) * dxiv &
                       + (tmp1(i,j+1,k) - tmp1(i,j,k)) * dyiv25) &
                       + (t13(i,j,k+1) - t13(i,j,k)) * dziv
                end do
              end do
            end do
            !$acc end kernels
          end if
        end if
      end if
    end if

    ! Calculate the v turbulent mixing
    if (trnopt == 0) then
      if (mfcopt == 0) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 2, ni-2
              t22(i,j,k) = jcb(i,j,k) * t22(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-1
            !$acc loop independent
            do i = 2, ni-1
              tmp1(i,j,k) = (jcb8v(i-1,j,k) + jcb8v(i,j,k)) * t12(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

        if (advopt <= 3) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 2, nj-1
              !$acc loop independent
              do i = 2, ni-2
                vfrc(i,j,k) = ((t22(i,j,k) - t22(i,j-1,k)) * dyiv &
                     + (tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv05) &
                     + (t23(i,j,k+1) - t23(i,j,k)) * dziv
              end do
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 2, nj-1
              !$acc loop independent
              do i = 2, ni-2
                vfrc(i,j,k) = vfrc(i,j,k) + ((t22(i,j,k) - t22(i,j-1,k)) * dyiv &
                     + (tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv05) &
                     + (t23(i,j,k+1) - t23(i,j,k)) * dziv
              end do
            end do
          end do
          !$acc end kernels
        end if
      else
        if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then
          if (mpopt == 0 .or. mpopt == 10) then
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 1, nj-1
                !$acc loop independent
                do i = 2, ni-2
                  t22(i,j,k) = rmf(i,j,2) * jcb(i,j,k) * t22(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-1
                  tmp1(i,j,k) = (jcb8v(i-1,j,k) + jcb8v(i,j,k)) * t12(i,j,k)
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 1, nj-1
                !$acc loop independent
                do i = 2, ni-2
                  t22(i,j,k) = jcb(i,j,k) * t22(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-1
                  tmp1(i,j,k) = t23(i,j,nk) * (jcb8v(i-1,j,k) + jcb8v(i,j,k)) * t12(i,j,k)
                end do
              end do
            end do
            !$acc end kernels
          end if

          if (advopt <= 3) then
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-2
                  vfrc(i,j,k) = mf8v(i,j) * ((t22(i,j,k) - t22(i,j-1,k)) * dyiv &
                       + (tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv25) &
                       + (t23(i,j,k+1) - t23(i,j,k)) * dziv
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-2
                  vfrc(i,j,k) = vfrc(i,j,k) &
                       + mf8v(i,j) * ((t22(i,j,k) - t22(i,j-1,k)) * dyiv &
                       + (tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv25) &
                       + (t23(i,j,k+1) - t23(i,j,k)) * dziv
                end do
              end do
            end do
            !$acc end kernels
          end if
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 2, ni-2
                t22(i,j,k) = rmf(i,j,2) * jcb(i,j,k) * t22(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 2, nj-1
              !$acc loop independent
              do i = 2, ni-1
                tmp1(i,j,k) = t23(i,j,nk) * (jcb8v(i-1,j,k) + jcb8v(i,j,k)) * t12(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

          if (advopt <= 3) then
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-2
                  vfrc(i,j,k) = rmf8v(i,j,1) * ((t22(i,j,k) - t22(i,j-1,k)) * dyiv &
                       + (tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv25) &
                       + (t23(i,j,k+1) - t23(i,j,k)) * dziv
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-2
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-2
                  vfrc(i,j,k) = vfrc(i,j,k) &
                       + rmf8v(i,j,1) * ((t22(i,j,k) - t22(i,j-1,k)) * dyiv &
                       + (tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv25) &
                       + (t23(i,j,k+1) - t23(i,j,k)) * dziv
                end do
              end do
            end do
            !$acc end kernels
          end if
        end if
      end if
    end if

    ! Calculate the w turbulent mixing
    if (trnopt == 0) then
      if (mfcopt == 0) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-1
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-1
              t11(i,j,k) = (jcb8u(i,j,k-1) + jcb8u(i,j,k)) * t31(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k = 2, nk-1
          !$acc loop independent
          do j = 2, nj-1
            !$acc loop independent
            do i = 2, ni-2
              t22(i,j,k) = (jcb8v(i,j,k-1) + jcb8v(i,j,k)) * t32(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

        if (advopt <= 3) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-1
            !$acc loop independent
            do j = 2, nj-2
              !$acc loop independent
              do i = 2, ni-2
                wfrc(i,j,k) = (t33(i,j,k) - t33(i,j,k-1)) * dziv &
                     + ((t11(i+1,j,k) - t11(i,j,k)) * dxiv05 &
                     + (t22(i,j+1,k) - t22(i,j,k)) * dyiv05)
              end do
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-1
            !$acc loop independent
            do j = 2, nj-2
              !$acc loop independent
              do i = 2, ni-2
                wfrc(i,j,k) = wfrc(i,j,k) + (t33(i,j,k) - t33(i,j,k-1)) * dziv &
                     + ((t11(i+1,j,k) - t11(i,j,k)) * dxiv05 &
                     + (t22(i,j+1,k) - t22(i,j,k)) * dyiv05)
              end do
            end do
          end do
          !$acc end kernels
        end if
      else
        if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then
          if (mpopt == 0 .or. mpopt == 10) then
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-1
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-1
                  t11(i,j,k) = (jcb8u(i,j,k-1) + jcb8u(i,j,k)) * t31(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

            !$acc kernels
            !$acc loop independent
            do k = 2, nk-1
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-2
                  t22(i,j,k) = rmf8v(i,j,2) * (jcb8v(i,j,k-1) + jcb8v(i,j,k)) * t32(i,j,k)
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-1
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-1
                  t11(i,j,k) = rmf8u(i,j,2) * (jcb8u(i,j,k-1) + jcb8u(i,j,k)) * t31(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

            !$acc kernels
            !$acc loop independent
            do k = 2, nk-1
              !$acc loop independent
              do j = 2, nj-1
                !$acc loop independent
                do i = 2, ni-2
                  t22(i,j,k) = (jcb8v(i,j,k-1) + jcb8v(i,j,k)) * t32(i,j,k)
                end do
              end do
            end do
            !$acc end kernels
          end if

          if (advopt <= 3) then
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-1
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-2
                  wfrc(i,j,k) = (t33(i,j,k) - t33(i,j,k-1)) * dziv &
                       + mf(i,j) * ((t11(i+1,j,k) - t11(i,j,k)) * dxiv05 &
                       + (t22(i,j+1,k) - t22(i,j,k)) * dyiv05)
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-1
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-2
                  wfrc(i,j,k) = wfrc(i,j,k) + (t33(i,j,k) - t33(i,j,k-1)) * dziv &
                       + mf(i,j) * ((t11(i+1,j,k) - t11(i,j,k)) * dxiv05 &
                       + (t22(i,j+1,k) - t22(i,j,k)) * dyiv05)
                end do
              end do
            end do
            !$acc end kernels
          end if
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-1
            !$acc loop independent
            do j = 2, nj-2
              !$acc loop independent
              do i = 2, ni-1
                t11(i,j,k) = rmf8u(i,j,2) * (jcb8u(i,j,k-1) + jcb8u(i,j,k)) * t31(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

          !$acc kernels
          !$acc loop independent
          do k = 2, nk-1
            !$acc loop independent
            do j = 2, nj-1
              !$acc loop independent
              do i = 2, ni-2
                t22(i,j,k) = rmf8v(i,j,2) * (jcb8v(i,j,k-1) + jcb8v(i,j,k)) * t32(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

          if (advopt <= 3) then
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-1
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-2
                  wfrc(i,j,k) = (t33(i,j,k) - t33(i,j,k-1)) * dziv &
                       + rmf(i,j,1) * ((t11(i+1,j,k) - t11(i,j,k)) * dxiv05 &
                       + (t22(i,j+1,k) - t22(i,j,k)) * dyiv05)
                end do
              end do
            end do
            !$acc end kernels
          else
            !$acc kernels
            !$acc loop independent
            do k = 2, nk-1
              !$acc loop independent
              do j = 2, nj-2
                !$acc loop independent
                do i = 2, ni-2
                  wfrc(i,j,k) = wfrc(i,j,k) + (t33(i,j,k) - t33(i,j,k-1)) * dziv &
                       + rmf(i,j,1) * ((t11(i+1,j,k) - t11(i,j,k)) * dxiv05 &
                       + (t22(i,j+1,k) - t22(i,j,k)) * dyiv05)
                end do
              end do
            end do
            !$acc end kernels
          end if
        end if
      end if
    end if

  end subroutine kernel_turbuvw

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

  subroutine read_parameters(filename, trnopt, mpopt, mfcopt, advopt, &
       dxiv, dyiv, dziv, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: trnopt, mpopt, mfcopt, advopt, ni, nj, nk
    real, intent(out) :: dxiv, dyiv, dziv
    character(len=256) :: line, name, val
    integer :: ios, eq_pos

    ! Initialize defaults
    trnopt = 0; mpopt = 0; mfcopt = 1; advopt = 4
    dxiv = 1.0; dyiv = 1.0; dziv = 1.0
    ni = 100; nj = 100; nk = 50

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
        name = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))
        select case (trim(name))
          case ('trnopt')
            read(val, *) trnopt
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
          case ('advopt')
            read(val, *) advopt
          case ('dxiv')
            read(val, *) dxiv
          case ('dyiv')
            read(val, *) dyiv
          case ('dziv')
            read(val, *) dziv
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  subroutine read_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: arr(i1:i2, j1:j2)
    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_2d

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

end program kernel_benchmark_turbuvw

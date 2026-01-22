!***********************************************************************
! GPU Kernel Benchmark: smoo4uvw (s_smoo4uvw)
!***********************************************************************
!
! Source: Src/smoo4uvw.f90
! Description: Perform 4th order numerical smoothing for velocity (u,v,w)
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_smoo4uvw
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: smtopt
  integer :: iwest, ieast, jsouth, jnorth
  real :: smhcoe, smvcoe

  ! Input arrays
  real, allocatable :: jcb8u(:,:,:), jcb8v(:,:,:), jcb8w(:,:,:)
  real, allocatable :: ubr(:,:,:), vbr(:,:,:)
  real, allocatable :: rst8u(:,:,:), rst8v(:,:,:), rst8w(:,:,:)
  real, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)

  ! Input/Output arrays
  real, allocatable :: ufrc(:,:,:), vfrc(:,:,:), wfrc(:,:,:)
  real, allocatable :: ufrc_init(:,:,:), vfrc_init(:,:,:), wfrc_init(:,:,:)

  ! Work arrays
  real, allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:)
  real, allocatable :: tmp4(:,:,:), tmp5(:,:,:)

  ! Reference output for validation
  real, allocatable :: ufrc_ref(:,:,:), vfrc_ref(:,:,:), wfrc_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error_u, max_error_v, max_error_w, rel_error
  real :: tolerance
  integer :: error_count_u, error_count_v, error_count_w
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
  call read_parameters(trim(data_dir)//'/params.txt', &
       ni, nj, nk, smtopt, iwest, ieast, jsouth, jnorth, smhcoe, smvcoe)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: smoo4uvw'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' smtopt = ', smtopt
  write(*,'(A,I6,A,I6)') ' iwest=', iwest, ', ieast=', ieast
  write(*,'(A,I6,A,I6)') ' jsouth=', jsouth, ', jnorth=', jnorth
  write(*,'(A,ES12.4)') ' smhcoe = ', smhcoe
  write(*,'(A,ES12.4)') ' smvcoe = ', smvcoe
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(jcb8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(ubr(0:ni+1, 0:nj+1, 1:nk))
  allocate(vbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(w(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp4(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp5(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/jcb8u.bin', jcb8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8v.bin', jcb8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8w.bin', jcb8w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ubr.bin', ubr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vbr.bin', vbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8u.bin', rst8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8v.bin', rst8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8w.bin', rst8w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ufrc_in.bin', ufrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_in.bin', vfrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wfrc_in.bin', wfrc_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ufrc_ref.bin', ufrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_ref.bin', vfrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wfrc_ref.bin', wfrc_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    !$acc wait
    ufrc = ufrc_init
    vfrc = vfrc_init
    wfrc = wfrc_init
    tmp1 = 0.0; tmp2 = 0.0; tmp3 = 0.0; tmp4 = 0.0; tmp5 = 0.0
    !$acc wait
    call kernel_smoo4uvw(smtopt, iwest, ieast, jsouth, jnorth, &
         smhcoe, smvcoe, ni, nj, nk, &
         jcb8u, jcb8v, jcb8w, ubr, vbr, rst8u, rst8v, rst8w, &
         u, v, w, ufrc, vfrc, wfrc, tmp1, tmp2, tmp3, tmp4, tmp5)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    !$acc wait
    ufrc = ufrc_init
    vfrc = vfrc_init
    wfrc = wfrc_init
    tmp1 = 0.0; tmp2 = 0.0; tmp3 = 0.0; tmp4 = 0.0; tmp5 = 0.0
    !$acc wait

    t_start = omp_get_wtime()
    call kernel_smoo4uvw(smtopt, iwest, ieast, jsouth, jnorth, &
         smhcoe, smvcoe, ni, nj, nk, &
         jcb8u, jcb8v, jcb8w, ubr, vbr, rst8u, rst8v, rst8w, &
         u, v, w, ufrc, vfrc, wfrc, tmp1, tmp2, tmp3, tmp4, tmp5)
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

  ! Validate ufrc
  max_error_u = 0.0
  error_count_u = 0
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-1
        rel_error = abs(ufrc(i,j,k) - ufrc_ref(i,j,k))
        if (abs(ufrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(ufrc_ref(i,j,k))
        end if
        if (rel_error > max_error_u) max_error_u = rel_error
        if (rel_error > tolerance) error_count_u = error_count_u + 1
      end do
    end do
  end do

  ! Validate vfrc
  max_error_v = 0.0
  error_count_v = 0
  do k = 2, nk-2
    do j = 2, nj-1
      do i = 2, ni-2
        rel_error = abs(vfrc(i,j,k) - vfrc_ref(i,j,k))
        if (abs(vfrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(vfrc_ref(i,j,k))
        end if
        if (rel_error > max_error_v) max_error_v = rel_error
        if (rel_error > tolerance) error_count_v = error_count_v + 1
      end do
    end do
  end do

  ! Validate wfrc
  max_error_w = 0.0
  error_count_w = 0
  do k = 2, nk-1
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(wfrc(i,j,k) - wfrc_ref(i,j,k))
        if (abs(wfrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(wfrc_ref(i,j,k))
        end if
        if (rel_error > max_error_w) max_error_w = rel_error
        if (rel_error > tolerance) error_count_w = error_count_w + 1
      end do
    end do
  end do

  validation_passed = (error_count_u == 0) .and. (error_count_v == 0) .and. (error_count_w == 0)

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
  write(*,'(A,ES12.4)') ' Max rel error (ufrc): ', max_error_u
  write(*,'(A,ES12.4)') ' Max rel error (vfrc): ', max_error_v
  write(*,'(A,ES12.4)') ' Max rel error (wfrc): ', max_error_w
  write(*,'(A,ES12.4)') ' Tolerance:            ', tolerance
  write(*,'(A,I12)') ' Error count (ufrc):   ', error_count_u
  write(*,'(A,I12)') ' Error count (vfrc):   ', error_count_v
  write(*,'(A,I12)') ' Error count (wfrc):   ', error_count_w
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(jcb8u, jcb8v, jcb8w, ubr, vbr, rst8u, rst8v, rst8w)
  deallocate(u, v, w, ufrc, vfrc, wfrc, ufrc_init, vfrc_init, wfrc_init)
  deallocate(tmp1, tmp2, tmp3, tmp4, tmp5)
  deallocate(ufrc_ref, vfrc_ref, wfrc_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: smoo4uvw (OpenACC version)
  ! Perform 4th order numerical smoothing for velocity (u,v,w)
  !=====================================================================
  subroutine kernel_smoo4uvw(smtopt, iwest, ieast, jsouth, jnorth, &
       smhcoe, smvcoe, ni, nj, nk, &
       jcb8u, jcb8v, jcb8w, ubr, vbr, rst8u, rst8v, rst8w, &
       u, v, w, ufrc, vfrc, wfrc, tmp1, tmp2, tmp3, tmp4, tmp5)
    implicit none

    integer, intent(in) :: smtopt
    integer, intent(in) :: iwest, ieast, jsouth, jnorth
    real, intent(in) :: smhcoe, smvcoe
    integer, intent(in) :: ni, nj, nk

    real, intent(in) :: jcb8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ubr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: ufrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: wfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp3(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp4(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp5(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    integer :: nkm1, nkm2

    nkm1 = nk - 1
    nkm2 = nk - 2

    ! === U smoothing ===
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = jsouth, nj-jnorth
        !$acc loop independent
        do i = iwest, ni+1-ieast
          tmp4(i,j,k) = rst8u(i,j,k) * (u(i,j,k) - ubr(i,j,k)) / jcb8u(i,j,k)
          tmp5(i,j,k) = 2.0 * tmp4(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 1+iwest, ni-ieast
          tmp1(i,j,k) = (tmp4(i+1,j,k) + tmp4(i-1,j,k)) - tmp5(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 1+jsouth, nj-1-jnorth
        !$acc loop independent
        do i = 2, ni-1
          tmp2(i,j,k) = (tmp4(i,j+1,k) + tmp4(i,j-1,k)) - tmp5(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-1
          tmp3(i,j,k) = (tmp4(i,j,k+1) + tmp4(i,j,k-1)) - tmp5(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-1
          ufrc(i,j,k) = ufrc(i,j,k) &
               + smhcoe * (tmp1(i,j,k) + tmp2(i,j,k)) + smvcoe * tmp3(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    if (mod(smtopt, 10) == 2) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-2
        !$acc loop independent
        do j = 2+jsouth, nj-2-jnorth
          !$acc loop independent
          do i = 2+iwest, ni-1-ieast
            ufrc(i,j,k) = ufrc(i,j,k) &
                 + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
                 + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k))))
          end do
        end do
      end do
      !$acc end kernels
      !$acc wait
    else
      !$acc kernels
      !$acc loop independent
      do j = 2+jsouth, nj-2-jnorth
        !$acc loop independent
        do i = 2+iwest, ni-1-ieast
          tmp3(i,j,1) = tmp3(i,j,2)
          tmp3(i,j,nkm1) = tmp3(i,j,nkm2)
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
          do i = 2+iwest, ni-1-ieast
            ufrc(i,j,k) = ufrc(i,j,k) &
                 + (smvcoe * (tmp3(i,j,k) - (tmp3(i,j,k+1) + tmp3(i,j,k-1))) &
                 + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
                 + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k)))))
          end do
        end do
      end do
      !$acc end kernels
      !$acc wait
    end if

    ! === V smoothing ===
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = jsouth, nj+1-jnorth
        !$acc loop independent
        do i = iwest, ni-ieast
          tmp4(i,j,k) = rst8v(i,j,k) * (v(i,j,k) - vbr(i,j,k)) / jcb8v(i,j,k)
          tmp5(i,j,k) = 2.0 * tmp4(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-1
        !$acc loop independent
        do i = 1+iwest, ni-1-ieast
          tmp1(i,j,k) = (tmp4(i+1,j,k) + tmp4(i-1,j,k)) - tmp5(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 1+jsouth, nj-jnorth
        !$acc loop independent
        do i = 2, ni-2
          tmp2(i,j,k) = (tmp4(i,j+1,k) + tmp4(i,j-1,k)) - tmp5(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-1
        !$acc loop independent
        do i = 2, ni-2
          tmp3(i,j,k) = (tmp4(i,j,k+1) + tmp4(i,j,k-1)) - tmp5(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-2
      !$acc loop independent
      do j = 2, nj-1
        !$acc loop independent
        do i = 2, ni-2
          vfrc(i,j,k) = vfrc(i,j,k) &
               + smhcoe * (tmp1(i,j,k) + tmp2(i,j,k)) + smvcoe * tmp3(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    if (mod(smtopt, 10) == 2) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-2
        !$acc loop independent
        do j = 2+jsouth, nj-1-jnorth
          !$acc loop independent
          do i = 2+iwest, ni-2-ieast
            vfrc(i,j,k) = vfrc(i,j,k) &
                 + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
                 + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k))))
          end do
        end do
      end do
      !$acc end kernels
      !$acc wait
    else
      !$acc kernels
      !$acc loop independent
      do j = 2+jsouth, nj-1-jnorth
        !$acc loop independent
        do i = 2+iwest, ni-2-ieast
          tmp3(i,j,1) = tmp3(i,j,2)
          tmp3(i,j,nkm1) = tmp3(i,j,nkm2)
        end do
      end do
      !$acc end kernels
      !$acc wait

      !$acc kernels
      !$acc loop independent
      do k = 2, nk-2
        !$acc loop independent
        do j = 2+jsouth, nj-1-jnorth
          !$acc loop independent
          do i = 2+iwest, ni-2-ieast
            vfrc(i,j,k) = vfrc(i,j,k) &
                 + (smvcoe * (tmp3(i,j,k) - (tmp3(i,j,k+1) + tmp3(i,j,k-1))) &
                 + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
                 + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k)))))
          end do
        end do
      end do
      !$acc end kernels
      !$acc wait
    end if

    ! === W smoothing ===
    !$acc kernels
    !$acc loop independent
    do k = 1, nk
      !$acc loop independent
      do j = jsouth, nj-jnorth
        !$acc loop independent
        do i = iwest, ni-ieast
          tmp4(i,j,k) = rst8w(i,j,k) * w(i,j,k) / jcb8w(i,j,k)
          tmp5(i,j,k) = 2.0 * tmp4(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 1+iwest, ni-1-ieast
          tmp1(i,j,k) = (tmp4(i+1,j,k) + tmp4(i-1,j,k)) - tmp5(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      !$acc loop independent
      do j = 1+jsouth, nj-1-jnorth
        !$acc loop independent
        do i = 2, ni-2
          tmp2(i,j,k) = (tmp4(i,j+1,k) + tmp4(i,j-1,k)) - tmp5(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          tmp3(i,j,k) = (tmp4(i,j,k+1) + tmp4(i,j,k-1)) - tmp5(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    !$acc kernels
    !$acc loop independent
    do k = 2, nk-1
      !$acc loop independent
      do j = 2, nj-2
        !$acc loop independent
        do i = 2, ni-2
          wfrc(i,j,k) = wfrc(i,j,k) &
               + smhcoe * (tmp1(i,j,k) + tmp2(i,j,k)) + smvcoe * tmp3(i,j,k)
        end do
      end do
    end do
    !$acc end kernels
    !$acc wait

    if (mod(smtopt, 10) == 2) then
      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do j = 2+jsouth, nj-2-jnorth
          !$acc loop independent
          do i = 2+iwest, ni-2-ieast
            wfrc(i,j,k) = wfrc(i,j,k) &
                 + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
                 + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k))))
          end do
        end do
      end do
      !$acc end kernels
      !$acc wait
    else
      !$acc kernels
      !$acc loop independent
      do j = 2+jsouth, nj-2-jnorth
        !$acc loop independent
        do i = 2+iwest, ni-2-ieast
          tmp3(i,j,1) = tmp3(i,j,2)
          tmp3(i,j,nk) = tmp3(i,j,nkm1)
        end do
      end do
      !$acc end kernels
      !$acc wait

      !$acc kernels
      !$acc loop independent
      do k = 2, nk-1
        !$acc loop independent
        do j = 2+jsouth, nj-2-jnorth
          !$acc loop independent
          do i = 2+iwest, ni-2-ieast
            wfrc(i,j,k) = wfrc(i,j,k) &
                 + (smvcoe * (tmp3(i,j,k) - (tmp3(i,j,k+1) + tmp3(i,j,k-1))) &
                 + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
                 + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k)))))
          end do
        end do
      end do
      !$acc end kernels
      !$acc wait
    end if

  end subroutine kernel_smoo4uvw

  !=====================================================================
  ! Configuration reader
  !=====================================================================
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

  !=====================================================================
  ! Parameter reader
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk, smtopt, &
       iwest, ieast, jsouth, jnorth, smhcoe, smvcoe)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, smtopt
    integer, intent(out) :: iwest, ieast, jsouth, jnorth
    real, intent(out) :: smhcoe, smvcoe

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
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('smtopt')
            read(val, *) smtopt
          case ('iwest')
            read(val, *) iwest
          case ('ieast')
            read(val, *) ieast
          case ('jsouth')
            read(val, *) jsouth
          case ('jnorth')
            read(val, *) jnorth
          case ('smhcoe')
            read(val, *) smhcoe
          case ('smvcoe')
            read(val, *) smvcoe
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary array reader
  !=====================================================================
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

end program kernel_benchmark_gpu_smoo4uvw

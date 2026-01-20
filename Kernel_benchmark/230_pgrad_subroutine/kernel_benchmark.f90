!***********************************************************************
! Kernel Benchmark: pgrad (s_pgrad)
!***********************************************************************
!
! Source: Src/pgrad.f90
! Description: Calculates pressure gradient force for u, v, w equations
!              including divergence damping and terrain-following corrections.
!
! Note: This benchmark covers the OpenMP parallel section.
!       When divopt >= 1, tmp1 (divergence) is computed using diver3d
!       before running the pgrad kernel.
!
!***********************************************************************
program kernel_benchmark_pgrad
  use omp_lib
  implicit none

  ! Grid dimensions
  integer :: ni, nj, nk

  ! Physics options
  integer :: trnopt, mpopt, mfcopt, divopt

  ! Grid parameters
  real :: dx, dy, dz, dxiv, dyiv, dziv, dziv25
  real :: dts, divch, divcv
  real :: divndc  ! Divergence damping coefficient (read from params.txt)

  ! Arrays for pgrad
  real, allocatable :: j31(:,:,:), j32(:,:,:), jcb(:,:,:)
  real, allocatable :: mf8u(:,:), mf8v(:,:)
  real, allocatable :: pp(:,:,:)
  real, allocatable :: upg(:,:,:), vpg(:,:,:), wpg(:,:,:)
  real, allocatable :: upg_ref(:,:,:), vpg_ref(:,:,:), wpg_ref(:,:,:)
  real, allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:)

  ! Arrays for diver3d (to compute tmp1 when divopt >= 1)
  real, allocatable :: mf(:,:)
  real, allocatable :: rmf(:,:,:), rmf8u(:,:,:), rmf8v(:,:,:)
  real, allocatable :: rst8u(:,:,:), rst8v(:,:,:), rst8w(:,:,:)
  real, allocatable :: u(:,:,:), v(:,:,:), wc(:,:,:)
  real, allocatable :: div3d_tmp1(:,:,:), div3d_tmp2(:,:,:), div3d_tmp3(:,:,:)

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
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', &
       ni, nj, nk, trnopt, mpopt, mfcopt, divopt, &
       dx, dy, dz, dxiv, dyiv, dziv, dts, divndc)

  dziv25 = 0.25e0 * dziv
  divch = divndc / dts * dx * dy
  divcv = divndc / dts * dz * dz

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: pgrad'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' trnopt: ', trnopt
  write(*,'(A,I6)') ' mpopt:  ', mpopt
  write(*,'(A,I6)') ' mfcopt: ', mfcopt
  write(*,'(A,I6)') ' divopt: ', divopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  ! pgrad arrays
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf8u(0:ni+1, 0:nj+1))
  allocate(mf8v(0:ni+1, 0:nj+1))
  allocate(pp(0:ni+1, 0:nj+1, 1:nk))
  allocate(upg(0:ni+1, 0:nj+1, 1:nk))
  allocate(vpg(0:ni+1, 0:nj+1, 1:nk))
  allocate(wpg(0:ni+1, 0:nj+1, 1:nk))
  allocate(upg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vpg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(wpg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))

  ! diver3d arrays (for computing tmp1 when divopt >= 1)
  if (divopt >= 1) then
    allocate(mf(0:ni+1, 0:nj+1))
    allocate(rmf(0:ni+1, 0:nj+1, 1:4))
    allocate(rmf8u(0:ni+1, 0:nj+1, 1:3))
    allocate(rmf8v(0:ni+1, 0:nj+1, 1:3))
    allocate(rst8u(0:ni+1, 0:nj+1, 1:nk))
    allocate(rst8v(0:ni+1, 0:nj+1, 1:nk))
    allocate(rst8w(0:ni+1, 0:nj+1, 1:nk))
    allocate(u(0:ni+1, 0:nj+1, 1:nk))
    allocate(v(0:ni+1, 0:nj+1, 1:nk))
    allocate(wc(0:ni+1, 0:nj+1, 1:nk))
    allocate(div3d_tmp1(0:ni+1, 0:nj+1, 1:nk))
    allocate(div3d_tmp2(0:ni+1, 0:nj+1, 1:nk))
    allocate(div3d_tmp3(0:ni+1, 0:nj+1, 1:nk))
  end if

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/mf8u.bin', mf8u, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/mf8v.bin', mf8v, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/pp.bin', pp, 0, ni+1, 0, nj+1, 1, nk)

  ! Read diver3d inputs and compute tmp1 (divergence) when divopt >= 1
  tmp1 = 0.0
  if (divopt >= 1) then
    write(*,'(A)') ' Loading diver3d inputs for divergence computation...'
    call read_array_2d(trim(data_dir)//'/mf.bin', mf, 0, ni+1, 0, nj+1)
    call read_array_3d_4(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
    call read_array_3d_3(trim(data_dir)//'/rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
    call read_array_3d_3(trim(data_dir)//'/rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
    call read_array_3d(trim(data_dir)//'/rst8u.bin', rst8u, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/rst8v.bin', rst8v, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/rst8w.bin', rst8w, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/wc.bin', wc, 0, ni+1, 0, nj+1, 1, nk)

    ! Compute divergence using diver3d kernel
    write(*,'(A)') ' Computing divergence (tmp1) using diver3d...'
    div3d_tmp1 = 0.0; div3d_tmp2 = 0.0; div3d_tmp3 = 0.0
    call kernel_diver3d(mpopt, mfcopt, dxiv, dyiv, dziv, ni, nj, nk, &
         mf, rmf, rmf8u, rmf8v, rst8u, rst8v, rst8w, u, v, wc, &
         tmp1, div3d_tmp1, div3d_tmp2, div3d_tmp3)
  end if

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/upg_ref.bin', upg_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vpg_ref.bin', vpg_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wpg_ref.bin', wpg_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    upg = 0.0; vpg = 0.0; wpg = 0.0
    ! tmp1 already contains divergence from diver3d (for divopt >= 1)
    ! For trnopt = 0, tmp1 is only read, not modified
    tmp2 = 0.0; tmp3 = 0.0
    call kernel_pgrad(trnopt, mpopt, mfcopt, divopt, ni, nj, nk, &
         j31, j32, jcb, mf8u, mf8v, pp, upg, vpg, wpg, &
         tmp1, tmp2, tmp3, dxiv, dyiv, dziv, dziv25, divch, divcv)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    upg = 0.0; vpg = 0.0; wpg = 0.0
    ! tmp1 already contains divergence from diver3d (for divopt >= 1)
    ! For trnopt = 0, tmp1 is only read, not modified
    tmp2 = 0.0; tmp3 = 0.0

    t_start = omp_get_wtime()
    call kernel_pgrad(trnopt, mpopt, mfcopt, divopt, ni, nj, nk, &
         j31, j32, jcb, mf8u, mf8v, pp, upg, vpg, wpg, &
         tmp1, tmp2, tmp3, dxiv, dyiv, dziv, dziv25, divch, divcv)
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

  ! Validate upg
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-1
        rel_error = abs(upg(i,j,k) - upg_ref(i,j,k))
        if (abs(upg_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(upg_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate vpg
  do k = 2, nk-2
    do j = 2, nj-1
      do i = 2, ni-2
        rel_error = abs(vpg(i,j,k) - vpg_ref(i,j,k))
        if (abs(vpg_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(vpg_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate wpg
  do k = 2, nk-1
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(wpg(i,j,k) - wpg_ref(i,j,k))
        if (abs(wpg_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(wpg_ref(i,j,k))
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

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(j31, j32, jcb, mf8u, mf8v, pp)
  deallocate(upg, vpg, wpg, upg_ref, vpg_ref, wpg_ref)
  deallocate(tmp1, tmp2, tmp3, times)
  if (divopt >= 1) then
    deallocate(mf, rmf, rmf8u, rmf8v)
    deallocate(rst8u, rst8v, rst8w, u, v, wc)
    deallocate(div3d_tmp1, div3d_tmp2, div3d_tmp3)
  end if

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: pgrad (OpenMP parallel section only)
  !=====================================================================
  subroutine kernel_pgrad(trnopt, mpopt, mfcopt, divopt, ni, nj, nk, &
       j31, j32, jcb, mf8u, mf8v, pp, upg, vpg, wpg, &
       tmp1, tmp2, tmp3, dxiv, dyiv, dziv, dziv25, divch, divcv)
    implicit none

    integer, intent(in) :: trnopt, mpopt, mfcopt, divopt
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: j31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: mf8u(0:ni+1, 0:nj+1)
    real, intent(in) :: mf8v(0:ni+1, 0:nj+1)
    real, intent(in) :: pp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: upg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: vpg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: wpg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp3(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: dxiv, dyiv, dziv, dziv25, divch, divcv

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! Add divergence damping to pressure perturbation
    if(divopt.eq.0) then
      do k=1,nk-1
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=pp(i,j,k)
          end do
        end do
        !$omp end do
      end do
    else
      do k=1,nk-1
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=pp(i,j,k)+divcv*jcb(i,j,k)*tmp1(i,j,k)
          end do
        end do
        !$omp end do
      end do
    end if

    ! Vertical pressure gradient
    do k=2,nk-1
      !$omp do schedule(runtime) private(i,j)
      do j=2,nj-2
        do i=2,ni-2
          wpg(i,j,k)=(tmp3(i,j,k-1)-tmp3(i,j,k))*dziv
        end do
      end do
      !$omp end do
    end do

    ! Reset for anisotropic case
    if(divopt.eq.2) then
      do k=1,nk-1
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=pp(i,j,k)+divch*tmp1(i,j,k)/jcb(i,j,k)
          end do
        end do
        !$omp end do
      end do
    end if

    ! Horizontal pressure gradient
    if(trnopt.eq.0) then

      do k=2,nk-2
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=jcb(i,j,k)*tmp3(i,j,k)
          end do
        end do
        !$omp end do
      end do

      if(mfcopt.eq.0) then
        do k=2,nk-2
          !$omp do schedule(runtime) private(i,j)
          do j=2,nj-2
            do i=2,ni-1
              upg(i,j,k)=(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv
            end do
          end do
          !$omp end do

          !$omp do schedule(runtime) private(i,j)
          do j=2,nj-1
            do i=2,ni-2
              vpg(i,j,k)=(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv
            end do
          end do
          !$omp end do
        end do
      else
        if(mpopt.eq.0.or.mpopt.eq.10) then
          do k=2,nk-2
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=mf8u(i,j)*(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv
              end do
            end do
            !$omp end do
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv
              end do
            end do
            !$omp end do
          end do
        else if(mpopt.eq.5) then
          do k=2,nk-2
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv
              end do
            end do
            !$omp end do
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=mf8v(i,j)*(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv
              end do
            end do
            !$omp end do
          end do
        else
          do k=2,nk-2
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=mf8u(i,j)*(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv
              end do
            end do
            !$omp end do
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=mf8v(i,j)*(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv
              end do
            end do
            !$omp end do
          end do
        end if
      end if

    else if(trnopt.ge.1) then

      do k=2,nk-1
        !$omp do schedule(runtime) private(i,j)
        do j=2,nj-2
          do i=2,ni-1
            tmp1(i,j,k)=((tmp3(i-1,j,k-1)+tmp3(i,j,k-1)) &
                 +(tmp3(i-1,j,k)+tmp3(i,j,k)))*j31(i,j,k)
          end do
        end do
        !$omp end do

        !$omp do schedule(runtime) private(i,j)
        do j=2,nj-1
          do i=2,ni-2
            tmp2(i,j,k)=((tmp3(i,j-1,k-1)+tmp3(i,j,k-1)) &
                 +(tmp3(i,j-1,k)+tmp3(i,j,k)))*j32(i,j,k)
          end do
        end do
        !$omp end do
      end do

      do k=2,nk-2
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=jcb(i,j,k)*tmp3(i,j,k)
          end do
        end do
        !$omp end do
      end do

      if(mfcopt.eq.0) then
        do k=2,nk-2
          !$omp do schedule(runtime) private(i,j)
          do j=2,nj-2
            do i=2,ni-1
              upg(i,j,k)=(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv &
                   -(tmp1(i,j,k+1)-tmp1(i,j,k))*dziv25
            end do
          end do
          !$omp end do
          !$omp do schedule(runtime) private(i,j)
          do j=2,nj-1
            do i=2,ni-2
              vpg(i,j,k)=(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv &
                   -(tmp2(i,j,k+1)-tmp2(i,j,k))*dziv25
            end do
          end do
          !$omp end do
        end do
      else
        if(mpopt.eq.0.or.mpopt.eq.10) then
          do k=2,nk-2
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=mf8u(i,j)*((tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv &
                     -(tmp1(i,j,k+1)-tmp1(i,j,k))*dziv25)
              end do
            end do
            !$omp end do
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv &
                     -(tmp2(i,j,k+1)-tmp2(i,j,k))*dziv25
              end do
            end do
            !$omp end do
          end do
        else if(mpopt.eq.5) then
          do k=2,nk-2
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv &
                     -(tmp1(i,j,k+1)-tmp1(i,j,k))*dziv25
              end do
            end do
            !$omp end do
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=mf8v(i,j)*((tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv &
                     -(tmp2(i,j,k+1)-tmp2(i,j,k))*dziv25)
              end do
            end do
            !$omp end do
          end do
        else
          do k=2,nk-2
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=mf8u(i,j)*((tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv &
                     -(tmp1(i,j,k+1)-tmp1(i,j,k))*dziv25)
              end do
            end do
            !$omp end do
            !$omp do schedule(runtime) private(i,j)
            do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=mf8v(i,j)*((tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv &
                     -(tmp2(i,j,k+1)-tmp2(i,j,k))*dziv25)
              end do
            end do
            !$omp end do
          end do
        end if
      end if

    end if

    !$omp end parallel

  end subroutine kernel_pgrad

  !=====================================================================
  ! I/O subroutines
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

  subroutine read_parameters(filename, ni, nj, nk, trnopt, mpopt, mfcopt, divopt, &
       dx, dy, dz, dxiv, dyiv, dziv, dts, divndc)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, trnopt, mpopt, mfcopt, divopt
    real, intent(out) :: dx, dy, dz, dxiv, dyiv, dziv, dts, divndc
    character(len=256) :: line, name, val
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
        name = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))
        select case (trim(name))
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('trnopt')
            read(val, *) trnopt
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
          case ('divopt')
            read(val, *) divopt
          case ('dx')
            read(val, *) dx
          case ('dy')
            read(val, *) dy
          case ('dz')
            read(val, *) dz
          case ('dxiv')
            read(val, *) dxiv
          case ('dyiv')
            read(val, *) dyiv
          case ('dziv')
            read(val, *) dziv
          case ('dts')
            read(val, *) dts
          case ('divndc')
            read(val, *) divndc
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
    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
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
    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

  subroutine read_array_3d_optional(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(inout) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: ios
    logical :: exists
    inquire(file=filename, exist=exists)
    if (.not. exists) then
      write(*,'(A,A)') '  WARNING: File not found, using zero: ', trim(filename)
      return
    end if
    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') '  WARNING: Cannot open file, using zero: ', trim(filename)
      return
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d_optional

  subroutine read_array_3d_3(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: ios
    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d_3

  subroutine read_array_3d_4(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: ios
    open(unit=10, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d_4

  !=====================================================================
  ! Kernel: diver3d (3D negative divergence)
  ! Used to compute tmp1 (divergence) for pgrad when divopt >= 1
  !=====================================================================
  subroutine kernel_diver3d(mpopt, mfcopt, dxiv, dyiv, dziv, ni, nj, nk, &
       mf, rmf, rmf8u, rmf8v, var8u, var8v, var8w, u, v, wc, &
       div3d, tmp1, tmp2, tmp3)
    implicit none

    integer, intent(in) :: mpopt, mfcopt
    real, intent(in) :: dxiv, dyiv, dziv
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: mf(0:ni+1, 0:nj+1)
    real, intent(in) :: rmf(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: rmf8u(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: rmf8v(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: var8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: var8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: var8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: wc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: div3d(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp3(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! Optional variables at u, v and w points are multiplied by u, v and wc
    if (mfcopt == 0) then
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 1, nj-1
          do i = 1, ni
            tmp1(i,j,k) = var8u(i,j,k) * u(i,j,k)
          end do
        end do
        !$omp end do
        !$omp do schedule(runtime) private(i,j)
        do j = 1, nj
          do i = 1, ni-1
            tmp2(i,j,k) = var8v(i,j,k) * v(i,j,k)
          end do
        end do
        !$omp end do
      end do
    else
      if (mpopt == 0 .or. mpopt == 10) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
            do i = 1, ni
              tmp1(i,j,k) = var8u(i,j,k) * u(i,j,k)
            end do
          end do
          !$omp end do
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj
            do i = 1, ni-1
              tmp2(i,j,k) = rmf8v(i,j,2) * var8v(i,j,k) * v(i,j,k)
            end do
          end do
          !$omp end do
        end do
      else if (mpopt == 5) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
            do i = 1, ni
              tmp1(i,j,k) = rmf8u(i,j,2) * var8u(i,j,k) * u(i,j,k)
            end do
          end do
          !$omp end do
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj
            do i = 1, ni-1
              tmp2(i,j,k) = var8v(i,j,k) * v(i,j,k)
            end do
          end do
          !$omp end do
        end do
      else
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
            do i = 1, ni
              tmp1(i,j,k) = rmf8u(i,j,2) * var8u(i,j,k) * u(i,j,k)
            end do
          end do
          !$omp end do
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj
            do i = 1, ni-1
              tmp2(i,j,k) = rmf8v(i,j,2) * var8v(i,j,k) * v(i,j,k)
            end do
          end do
          !$omp end do
        end do
      end if
    end if

    do k = 1, nk
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          tmp3(i,j,k) = var8w(i,j,k) * wc(i,j,k)
        end do
      end do
      !$omp end do
    end do

    ! Calculate the 3 dimensional divergence
    if (mfcopt == 0) then
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 1, nj-1
          do i = 1, ni-1
            div3d(i,j,k) = (tmp3(i,j,k) - tmp3(i,j,k+1)) * dziv &
                 + ((tmp1(i,j,k) - tmp1(i+1,j,k)) * dxiv &
                 + (tmp2(i,j,k) - tmp2(i,j+1,k)) * dyiv)
          end do
        end do
        !$omp end do
      end do
    else
      if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
            do i = 1, ni-1
              div3d(i,j,k) = mf(i,j) * ((tmp1(i,j,k) - tmp1(i+1,j,k)) * dxiv &
                   + (tmp2(i,j,k) - tmp2(i,j+1,k)) * dyiv) &
                   + (tmp3(i,j,k) - tmp3(i,j,k+1)) * dziv
            end do
          end do
          !$omp end do
        end do
      else
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
            do i = 1, ni-1
              div3d(i,j,k) = rmf(i,j,1) * ((tmp1(i,j,k) - tmp1(i+1,j,k)) * dxiv &
                   + (tmp2(i,j,k) - tmp2(i,j+1,k)) * dyiv) &
                   + (tmp3(i,j,k) - tmp3(i,j,k+1)) * dziv
            end do
          end do
          !$omp end do
        end do
      end if
    end if

    !$omp end parallel

  end subroutine kernel_diver3d

end program kernel_benchmark_pgrad

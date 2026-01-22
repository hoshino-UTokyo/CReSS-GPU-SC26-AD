!***********************************************************************
! Kernel Benchmark: advuvw (s_advuvw) - GPU Version
!***********************************************************************
!
! Source: Src/advuvw.f90
! Description: Calculates velocity (u, v, w) advection using 2nd or 4th
!              order centered finite difference schemes.
!
!***********************************************************************
program kernel_benchmark_advuvw
  use omp_lib
  implicit none

  ! Math constants (from commath module)
  real, parameter :: oned24 = 1.0 / 24.0
  real, parameter :: fourd3 = 4.0 / 3.0

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: advopt
  integer :: iwest, ieast, jsouth, jnorth
  real :: dxv24, dxv25n, dxv48
  real :: dyv24, dyv25n, dyv48
  real :: dzv24, dzv25n, dzv48

  ! Input arrays
  real, allocatable :: rstxu(:,:,:)
  real, allocatable :: rstxv(:,:,:)
  real, allocatable :: rstxwc(:,:,:)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)
  real, allocatable :: w(:,:,:)

  ! Input/output arrays
  real, allocatable :: ufrc(:,:,:)
  real, allocatable :: vfrc(:,:,:)
  real, allocatable :: wfrc(:,:,:)
  real, allocatable :: hadv(:,:,:)
  real, allocatable :: vadv(:,:,:)
  real, allocatable :: tmp1(:,:,:)
  real, allocatable :: tmp2(:,:,:)
  real, allocatable :: tmp3(:,:,:)

  ! Arrays to hold initial input data
  real, allocatable :: ufrc_in(:,:,:)
  real, allocatable :: vfrc_in(:,:,:)
  real, allocatable :: wfrc_in(:,:,:)
  real, allocatable :: hadv_in(:,:,:)
  real, allocatable :: vadv_in(:,:,:)
  real, allocatable :: tmp1_in(:,:,:)
  real, allocatable :: tmp2_in(:,:,:)
  real, allocatable :: tmp3_in(:,:,:)

  ! Reference output for validation
  real, allocatable :: ufrc_ref(:,:,:)
  real, allocatable :: vfrc_ref(:,:,:)
  real, allocatable :: wfrc_ref(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', advopt, &
       iwest, ieast, jsouth, jnorth, &
       dxv24, dxv25n, dxv48, dyv24, dyv25n, dyv48, dzv24, dzv25n, dzv48, &
       ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: advuvw (GPU)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' advopt=', advopt
  write(*,'(A,I4,A,I4,A,I4,A,I4)') ' Boundaries: iwest=', iwest, ', ieast=', ieast, &
       ', jsouth=', jsouth, ', jnorth=', jnorth
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rstxu(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxv(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxwc(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(w(0:ni+1, 0:nj+1, 1:nk))

  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(hadv(0:ni+1, 0:nj+1, 1:nk))
  allocate(vadv(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))

  allocate(ufrc_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(hadv_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(vadv_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3_in(0:ni+1, 0:nj+1, 1:nk))

  allocate(ufrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc_ref(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rstxu.bin', rstxu, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rstxv.bin', rstxv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rstxwc.bin', rstxwc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ufrc_in.bin', ufrc_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_in.bin', vfrc_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wfrc_in.bin', wfrc_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/hadv_in.bin', hadv_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vadv_in.bin', vadv_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp1_in.bin', tmp1_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp2_in.bin', tmp2_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp3_in.bin', tmp3_in, 0, ni+1, 0, nj+1, 1, nk)

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
    ufrc = ufrc_in
    vfrc = vfrc_in
    wfrc = wfrc_in
    hadv = hadv_in
    vadv = vadv_in
    tmp1 = tmp1_in
    tmp2 = tmp2_in
    tmp3 = tmp3_in
    call kernel_advuvw(advopt, iwest, ieast, jsouth, jnorth, &
         dxv24, dxv25n, dxv48, dyv24, dyv25n, dyv48, dzv24, dzv25n, dzv48, &
         ni, nj, nk, rstxu, rstxv, rstxwc, &
         u, v, w, ufrc, vfrc, wfrc, hadv, vadv, tmp1, tmp2, tmp3)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ufrc = ufrc_in
    vfrc = vfrc_in
    wfrc = wfrc_in
    hadv = hadv_in
    vadv = vadv_in
    tmp1 = tmp1_in
    tmp2 = tmp2_in
    tmp3 = tmp3_in

    !$acc wait
    t_start = omp_get_wtime()

    call kernel_advuvw(advopt, iwest, ieast, jsouth, jnorth, &
         dxv24, dxv25n, dxv48, dyv24, dyv25n, dyv48, dzv24, dzv25n, dzv48, &
         ni, nj, nk, rstxu, rstxv, rstxwc, &
         u, v, w, ufrc, vfrc, wfrc, hadv, vadv, tmp1, tmp2, tmp3)

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

  ! Check ufrc
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

  ! Check vfrc
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

  ! Check wfrc
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

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(rstxu, rstxv, rstxwc, u, v, w)
  deallocate(ufrc, vfrc, wfrc, hadv, vadv, tmp1, tmp2, tmp3)
  deallocate(ufrc_in, vfrc_in, wfrc_in, hadv_in, vadv_in)
  deallocate(tmp1_in, tmp2_in, tmp3_in)
  deallocate(ufrc_ref, vfrc_ref, wfrc_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: advuvw (GPU version)
  ! Calculate velocity advection
  !=====================================================================
  subroutine kernel_advuvw(advopt, iwest, ieast, jsouth, jnorth, &
       dxv24, dxv25n, dxv48, dyv24, dyv25n, dyv48, dzv24, dzv25n, dzv48, &
       ni, nj, nk, rstxu, rstxv, rstxwc, &
       u, v, w, ufrc, vfrc, wfrc, hadv, vadv, tmp1, tmp2, tmp3)
    implicit none

    integer, intent(in) :: advopt
    integer, intent(in) :: iwest, ieast, jsouth, jnorth
    real, intent(in) :: dxv24, dxv25n, dxv48
    real, intent(in) :: dyv24, dyv25n, dyv48
    real, intent(in) :: dzv24, dzv25n, dzv48
    integer, intent(in) :: ni, nj, nk

    real, intent(in) :: rstxu(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rstxv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rstxwc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: w(0:ni+1, 0:nj+1, 1:nk)

    real, intent(inout) :: ufrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: wfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: hadv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vadv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp3(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !! Calculate the u advection (2nd order)
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 1, ni-1
          tmp1(i,j,k) = (rstxu(i,j,k) + rstxu(i+1,j,k)) * (u(i+1,j,k) - u(i,j,k)) * dxv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-1
        do i = 2, ni-1
          tmp2(i,j,k) = (rstxv(i-1,j,k) + rstxv(i,j,k)) * (u(i,j,k) - u(i,j-1,k)) * dyv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-2
        do i = 2, ni-1
          tmp3(i,j,k) = (rstxwc(i-1,j,k) + rstxwc(i,j,k)) * (u(i,j,k) - u(i,j,k-1)) * dzv25n
        end do
      end do
    end do
    !$acc end kernels

    if (advopt == 1) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2, nj-2
          do i = 2, ni-1
            ufrc(i,j,k) = ufrc(i,j,k) + ((tmp3(i,j,k) + tmp3(i,j,k+1)) &
                 + ((tmp1(i-1,j,k) + tmp1(i,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))))
          end do
        end do
      end do
      !$acc end kernels
    else
      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              hadv(i,j,k) = (tmp1(i-1,j,k) + tmp1(i,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              hadv(i,j,k) = (tmp1(i-1,j,k) + tmp1(i,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))
              vadv(i,j,k) = tmp3(i,j,k) + tmp3(i,j,k+1)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! 4th order u advection
    if (advopt == 2 .or. advopt == 3) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-2-jnorth
          do i = 1+iwest, ni-ieast
            tmp1(i,j,k) = (rstxu(i-1,j,k) + rstxu(i+1,j,k)) * (u(i+1,j,k) - u(i-1,j,k)) * dxv24
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 1+jsouth, nj-1-jnorth
          do i = 2+iwest, ni-1-ieast
            tmp2(i,j,k) = ((rstxv(i-1,j,k) + rstxv(i,j,k)) + (rstxv(i-1,j+1,k) + rstxv(i,j+1,k))) &
                 * (u(i,j+1,k) - u(i,j-1,k)) * dyv48
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-2-jnorth
          do i = 2+iwest, ni-1-ieast
            hadv(i,j,k) = fourd3 * hadv(i,j,k) &
                 + ((tmp1(i-1,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j+1,k)))
          end do
        end do
      end do
      !$acc end kernels

      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) + hadv(i,j,k) + (tmp3(i,j,k) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-1-ieast
              tmp3(i,j,k) = ((rstxwc(i-1,j,k) + rstxwc(i,j,k)) + (rstxwc(i-1,j,k+1) + rstxwc(i,j,k+1))) &
                   * (u(i,j,k+1) - u(i,j,k-1)) * dzv48
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 3, nk-3
          do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-1-ieast
              vadv(i,j,k) = fourd3 * vadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) + hadv(i,j,k) + vadv(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    !! Calculate the v advection (2nd order)
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-1
        do i = 2, ni-1
          tmp1(i,j,k) = (rstxu(i,j-1,k) + rstxu(i,j,k)) * (v(i,j,k) - v(i-1,j,k)) * dxv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 1, nj-1
        do i = 2, ni-2
          tmp2(i,j,k) = (rstxv(i,j,k) + rstxv(i,j+1,k)) * (v(i,j+1,k) - v(i,j,k)) * dyv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-1
        do i = 2, ni-2
          tmp3(i,j,k) = (rstxwc(i,j-1,k) + rstxwc(i,j,k)) * (v(i,j,k) - v(i,j,k-1)) * dzv25n
        end do
      end do
    end do
    !$acc end kernels

    if (advopt == 1) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2, nj-1
          do i = 2, ni-2
            vfrc(i,j,k) = vfrc(i,j,k) + ((tmp3(i,j,k) + tmp3(i,j,k+1)) &
                 + ((tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j,k))))
          end do
        end do
      end do
      !$acc end kernels
    else
      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j,k))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j,k))
              vadv(i,j,k) = tmp3(i,j,k) + tmp3(i,j,k+1)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! 4th order v advection
    if (advopt == 2 .or. advopt == 3) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-1-jnorth
          do i = 1+iwest, ni-1-ieast
            tmp1(i,j,k) = ((rstxu(i,j-1,k) + rstxu(i,j,k)) + (rstxu(i+1,j-1,k) + rstxu(i+1,j,k))) &
                 * (v(i+1,j,k) - v(i-1,j,k)) * dxv48
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 1+jsouth, nj-jnorth
          do i = 2+iwest, ni-2-ieast
            tmp2(i,j,k) = (rstxv(i,j-1,k) + rstxv(i,j+1,k)) * (v(i,j+1,k) - v(i,j-1,k)) * dyv24
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2+jsouth, nj-1-jnorth
          do i = 2+iwest, ni-2-ieast
            hadv(i,j,k) = fourd3 * hadv(i,j,k) &
                 + ((tmp1(i-1,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j+1,k)))
          end do
        end do
      end do
      !$acc end kernels

      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) + hadv(i,j,k) + (tmp3(i,j,k) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2+jsouth, nj-1-jnorth
            do i = 2+iwest, ni-2-ieast
              tmp3(i,j,k) = ((rstxwc(i,j-1,k) + rstxwc(i,j,k)) + (rstxwc(i,j-1,k+1) + rstxwc(i,j,k+1))) &
                   * (v(i,j,k+1) - v(i,j,k-1)) * dzv48
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 3, nk-3
          do j = 2+jsouth, nj-1-jnorth
            do i = 2+iwest, ni-2-ieast
              vadv(i,j,k) = fourd3 * vadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) + hadv(i,j,k) + vadv(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    !! Calculate the w advection (2nd order)
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-2
        do i = 2, ni-1
          tmp1(i,j,k) = (rstxu(i,j,k-1) + rstxu(i,j,k)) * (w(i,j,k) - w(i-1,j,k)) * dxv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-1
        do i = 2, ni-2
          tmp2(i,j,k) = (rstxv(i,j,k-1) + rstxv(i,j,k)) * (w(i,j,k) - w(i,j-1,k)) * dyv25n
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-1
      do j = 2, nj-2
        do i = 2, ni-2
          tmp3(i,j,k) = (rstxwc(i,j,k) + rstxwc(i,j,k+1)) * (w(i,j,k+1) - w(i,j,k)) * dzv25n
        end do
      end do
    end do
    !$acc end kernels

    if (advopt == 1) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-1
        do j = 2, nj-2
          do i = 2, ni-2
            wfrc(i,j,k) = wfrc(i,j,k) + ((tmp3(i,j,k-1) + tmp3(i,j,k)) &
                 + ((tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))))
          end do
        end do
      end do
      !$acc end kernels
    else
      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2, nj-2
            do i = 2, ni-2
              hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2, nj-2
            do i = 2, ni-2
              hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k))
              vadv(i,j,k) = tmp3(i,j,k-1) + tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! 4th order w advection
    if (advopt == 2 .or. advopt == 3) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-1
        do j = 2+jsouth, nj-2-jnorth
          do i = 1+iwest, ni-1-ieast
            tmp1(i,j,k) = ((rstxu(i,j,k-1) + rstxu(i,j,k)) + (rstxu(i+1,j,k-1) + rstxu(i+1,j,k))) &
                 * (w(i+1,j,k) - w(i-1,j,k)) * dxv48
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-1
        do j = 1+jsouth, nj-1-jnorth
          do i = 2+iwest, ni-2-ieast
            tmp2(i,j,k) = ((rstxv(i,j,k-1) + rstxv(i,j,k)) + (rstxv(i,j+1,k-1) + rstxv(i,j+1,k))) &
                 * (w(i,j+1,k) - w(i,j-1,k)) * dyv48
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-1
        do j = 2+jsouth, nj-2-jnorth
          do i = 2+iwest, ni-2-ieast
            hadv(i,j,k) = fourd3 * hadv(i,j,k) &
                 + ((tmp1(i-1,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j-1,k) + tmp2(i,j+1,k)))
          end do
        end do
      end do
      !$acc end kernels

      if (advopt == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2, nj-2
            do i = 2, ni-2
              wfrc(i,j,k) = wfrc(i,j,k) + hadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k))
            end do
          end do
        end do
        !$acc end kernels
      else if (advopt == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-2-ieast
              tmp3(i,j,k) = (rstxwc(i,j,k-1) + rstxwc(i,j,k+1)) * (w(i,j,k+1) - w(i,j,k-1)) * dzv24
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 3, nk-2
          do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-2-ieast
              vadv(i,j,k) = fourd3 * vadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-1
          do j = 2, nj-2
            do i = 2, ni-2
              wfrc(i,j,k) = wfrc(i,j,k) + hadv(i,j,k) + vadv(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_advuvw

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
  subroutine read_parameters(filename, advopt, iwest, ieast, jsouth, jnorth, &
       dxv24, dxv25n, dxv48, dyv24, dyv25n, dyv48, dzv24, dzv25n, dzv48, &
       ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: advopt, iwest, ieast, jsouth, jnorth
    real, intent(out) :: dxv24, dxv25n, dxv48
    real, intent(out) :: dyv24, dyv25n, dyv48
    real, intent(out) :: dzv24, dzv25n, dzv48
    integer, intent(out) :: ni, nj, nk

    integer :: ios
    character(len=256) :: line
    character(len=64) :: name
    character(len=64) :: val_str

    ! Defaults
    advopt = 1
    iwest = 0
    ieast = 0
    jsouth = 0
    jnorth = 0
    dxv24 = 0.0
    dxv25n = 0.0
    dxv48 = 0.0
    dyv24 = 0.0
    dyv25n = 0.0
    dyv48 = 0.0
    dzv24 = 0.0
    dzv25n = 0.0
    dzv48 = 0.0
    ni = 1
    nj = 1
    nk = 1

    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open parameter file: ', trim(filename)
      stop 1
    end if

    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit

      if (index(line, '=') > 0) then
        read(line, *) name
        val_str = line(index(line, '=')+1:)

        select case (trim(name))
        case ('advopt')
          read(val_str, *) advopt
        case ('iwest')
          read(val_str, *) iwest
        case ('ieast')
          read(val_str, *) ieast
        case ('jsouth')
          read(val_str, *) jsouth
        case ('jnorth')
          read(val_str, *) jnorth
        case ('dxv24')
          read(val_str, *) dxv24
        case ('dxv25n')
          read(val_str, *) dxv25n
        case ('dxv48')
          read(val_str, *) dxv48
        case ('dyv24')
          read(val_str, *) dyv24
        case ('dyv25n')
          read(val_str, *) dyv25n
        case ('dyv48')
          read(val_str, *) dyv48
        case ('dzv24')
          read(val_str, *) dzv24
        case ('dzv25n')
          read(val_str, *) dzv25n
        case ('dzv48')
          read(val_str, *) dzv48
        case ('ni')
          read(val_str, *) ni
        case ('nj')
          read(val_str, *) nj
        case ('nk')
          read(val_str, *) nk
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

end program kernel_benchmark_advuvw

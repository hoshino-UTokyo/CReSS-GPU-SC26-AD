!***********************************************************************
! Kernel Benchmark: turbtke (s_turbtke)
!***********************************************************************
!
! Source: Src/turbtke.f90
! Description: Calculate the turbulent kinetic energy mixing.
!              Computes divergence of turbulent fluxes with terrain
!              and map scale factor corrections.
!
!***********************************************************************
program kernel_benchmark_turbtke
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: trnopt, mpopt, mfcopt
  real :: dxiv2, dyiv2, dziv2

  ! Input arrays
  real, allocatable :: j31(:,:,:)
  real, allocatable :: j32(:,:,:)
  real, allocatable :: jcb8u(:,:,:)
  real, allocatable :: jcb8v(:,:,:)
  real, allocatable :: mf(:,:)
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: rmf8u(:,:,:)
  real, allocatable :: rmf8v(:,:,:)
  real, allocatable :: h1(:,:,:)
  real, allocatable :: h2(:,:,:)
  real, allocatable :: h3(:,:,:)

  ! Input/output arrays
  real, allocatable :: tkefrc(:,:,:)
  real, allocatable :: tmp1(:,:,:)
  real, allocatable :: tmp2(:,:,:)
  real, allocatable :: tmp3(:,:,:)

  ! Reference outputs for validation
  real, allocatable :: tkefrc_ref(:,:,:)
  real, allocatable :: tmp1_ref(:,:,:)
  real, allocatable :: tmp2_ref(:,:,:)
  real, allocatable :: tmp3_ref(:,:,:)

  ! Initial input (for resetting)
  real, allocatable :: tkefrc_in(:,:,:)
  real, allocatable :: tmp1_in(:,:,:)
  real, allocatable :: tmp2_in(:,:,:)
  real, allocatable :: tmp3_in(:,:,:)

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
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, &
                       trnopt, mpopt, mfcopt, dxiv2, dyiv2, dziv2)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: turbtke'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' trnopt: ', trnopt
  write(*,'(A,I6)') ' mpopt:  ', mpopt
  write(*,'(A,I6)') ' mfcopt: ', mfcopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rmf8u(0:ni+1, 0:nj+1, 1:3))
  allocate(rmf8v(0:ni+1, 0:nj+1, 1:3))
  allocate(h1(0:ni+1, 0:nj+1, 1:nk))
  allocate(h2(0:ni+1, 0:nj+1, 1:nk))
  allocate(h3(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8u.bin', jcb8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8v.bin', jcb8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/mf.bin', mf, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call read_array_3d(trim(data_dir)//'/rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call read_array_3d(trim(data_dir)//'/rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  call read_array_3d(trim(data_dir)//'/h1.bin', h1, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/h2.bin', h2, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/h3.bin', h3, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tkefrc_in.bin', tkefrc_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp1_in.bin', tmp1_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp2_in.bin', tmp2_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp3_in.bin', tmp3_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/tkefrc_ref.bin', tkefrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp1_ref.bin', tmp1_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp2_ref.bin', tmp2_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp3_ref.bin', tmp3_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    tkefrc = tkefrc_in
    tmp1 = tmp1_in
    tmp2 = tmp2_in
    tmp3 = tmp3_in
    call kernel_turbtke(trnopt, mpopt, mfcopt, dxiv2, dyiv2, dziv2, &
                        ni, nj, nk, j31, j32, jcb8u, jcb8v, &
                        mf, rmf, rmf8u, rmf8v, h1, h2, h3, &
                        tkefrc, tmp1, tmp2, tmp3)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    tkefrc = tkefrc_in
    tmp1 = tmp1_in
    tmp2 = tmp2_in
    tmp3 = tmp3_in

    t_start = omp_get_wtime()
    call kernel_turbtke(trnopt, mpopt, mfcopt, dxiv2, dyiv2, dziv2, &
                        ni, nj, nk, j31, j32, jcb8u, jcb8v, &
                        mf, rmf, rmf8u, rmf8v, h1, h2, h3, &
                        tkefrc, tmp1, tmp2, tmp3)
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output (interior points for tkefrc)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0

  ! Validate tkefrc - main output (interior points: 2:ni-2, 2:nj-2, 2:nk-2)
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(tkefrc(i,j,k) - tkefrc_ref(i,j,k))
        if (abs(tkefrc_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(tkefrc_ref(i,j,k))
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
  deallocate(j31, j32, jcb8u, jcb8v, mf, rmf, rmf8u, rmf8v)
  deallocate(h1, h2, h3, tkefrc, tmp1, tmp2, tmp3)
  deallocate(tkefrc_ref, tmp1_ref, tmp2_ref, tmp3_ref)
  deallocate(tkefrc_in, tmp1_in, tmp2_in, tmp3_in, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: turbtke
  ! Calculate the turbulent kinetic energy mixing
  !=====================================================================
  subroutine kernel_turbtke(trnopt, mpopt, mfcopt, dxiv2, dyiv2, dziv2, &
                            ni, nj, nk, j31, j32, jcb8u, jcb8v, &
                            mf, rmf, rmf8u, rmf8v, h1, h2, h3, &
                            tkefrc, tmp1, tmp2, tmp3)
    implicit none

    integer, intent(in) :: trnopt, mpopt, mfcopt
    real, intent(in) :: dxiv2, dyiv2, dziv2
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: j31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: mf(0:ni+1, 0:nj+1)
    real, intent(in) :: rmf(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: rmf8u(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: rmf8v(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: h1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: h2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: h3(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tkefrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp3(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    if (trnopt .eq. 0) then

      if (mfcopt .eq. 0) then

        do k = 2, nk-2
          !$omp do schedule(runtime) private(i,j)
          do j = 2, nj-2
            do i = 2, ni-1
              tmp1(i,j,k) = jcb8u(i,j,k) * h1(i,j,k)
            end do
          end do
          !$omp end do

          !$omp do schedule(runtime) private(i,j)
          do j = 2, nj-1
            do i = 2, ni-2
              tmp2(i,j,k) = jcb8v(i,j,k) * h2(i,j,k)
            end do
          end do
          !$omp end do
        end do

        do k = 2, nk-2
          !$omp do schedule(runtime) private(i,j)
          do j = 2, nj-2
            do i = 2, ni-2
              tkefrc(i,j,k) = tkefrc(i,j,k) &
                + ((h3(i,j,k+1) - h3(i,j,k)) * dziv2 &
                + ((tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv2 &
                + (tmp2(i,j+1,k) - tmp2(i,j,k)) * dyiv2))
            end do
          end do
          !$omp end do
        end do

      else

        if (mpopt .eq. 0 .or. mpopt .eq. 5 .or. mpopt .eq. 10) then

          if (mpopt .eq. 0 .or. mpopt .eq. 10) then

            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
                do i = 2, ni-1
                  tmp1(i,j,k) = jcb8u(i,j,k) * h1(i,j,k)
                end do
              end do
              !$omp end do

              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-1
                do i = 2, ni-2
                  tmp2(i,j,k) = rmf8v(i,j,2) * jcb8v(i,j,k) * h2(i,j,k)
                end do
              end do
              !$omp end do
            end do

          else

            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
                do i = 2, ni-1
                  tmp1(i,j,k) = rmf8u(i,j,2) * jcb8u(i,j,k) * h1(i,j,k)
                end do
              end do
              !$omp end do

              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-1
                do i = 2, ni-2
                  tmp2(i,j,k) = jcb8v(i,j,k) * h2(i,j,k)
                end do
              end do
              !$omp end do
            end do

          end if

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                tkefrc(i,j,k) = tkefrc(i,j,k) &
                  + ((h3(i,j,k+1) - h3(i,j,k)) * dziv2 &
                  + mf(i,j) * ((tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv2 &
                  + (tmp2(i,j+1,k) - tmp2(i,j,k)) * dyiv2))
              end do
            end do
            !$omp end do
          end do

        else

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-1
                tmp1(i,j,k) = rmf8u(i,j,2) * jcb8u(i,j,k) * h1(i,j,k)
              end do
            end do
            !$omp end do

            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-1
              do i = 2, ni-2
                tmp2(i,j,k) = rmf8v(i,j,2) * jcb8v(i,j,k) * h2(i,j,k)
              end do
            end do
            !$omp end do
          end do

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                tkefrc(i,j,k) = tkefrc(i,j,k) &
                  + ((h3(i,j,k+1) - h3(i,j,k)) * dziv2 &
                  + rmf(i,j,1) * ((tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv2 &
                  + (tmp2(i,j+1,k) - tmp2(i,j,k)) * dyiv2))
              end do
            end do
            !$omp end do
          end do

        end if

      end if

    else
      ! trnopt /= 0 case (terrain)

      do k = 2, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-1
            tmp1(i,j,k) = j31(i,j,k) * (h1(i,j,k-1) + h1(i,j,k))
          end do
        end do
        !$omp end do

        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-1
          do i = 2, ni-2
            tmp2(i,j,k) = j32(i,j,k) * (h2(i,j,k-1) + h2(i,j,k))
          end do
        end do
        !$omp end do
      end do

      do k = 2, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-2
            tmp3(i,j,k) = h3(i,j,k) + .25e0 &
              * ((tmp1(i,j,k) + tmp1(i+1,j,k)) + (tmp2(i,j,k) + tmp2(i,j+1,k)))
          end do
        end do
        !$omp end do
      end do

      if (mfcopt .eq. 0) then

        do k = 2, nk-2
          !$omp do schedule(runtime) private(i,j)
          do j = 2, nj-2
            do i = 2, ni-1
              tmp1(i,j,k) = jcb8u(i,j,k) * h1(i,j,k)
            end do
          end do
          !$omp end do

          !$omp do schedule(runtime) private(i,j)
          do j = 2, nj-1
            do i = 2, ni-2
              tmp2(i,j,k) = jcb8v(i,j,k) * h2(i,j,k)
            end do
          end do
          !$omp end do
        end do

        do k = 2, nk-2
          !$omp do schedule(runtime) private(i,j)
          do j = 2, nj-2
            do i = 2, ni-2
              tkefrc(i,j,k) = tkefrc(i,j,k) &
                + ((tmp3(i,j,k+1) - tmp3(i,j,k)) * dziv2 &
                + ((tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv2 &
                + (tmp2(i,j+1,k) - tmp2(i,j,k)) * dyiv2))
            end do
          end do
          !$omp end do
        end do

      else

        if (mpopt .eq. 0 .or. mpopt .eq. 5 .or. mpopt .eq. 10) then

          if (mpopt .eq. 0 .or. mpopt .eq. 10) then

            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
                do i = 2, ni-1
                  tmp1(i,j,k) = jcb8u(i,j,k) * h1(i,j,k)
                end do
              end do
              !$omp end do

              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-1
                do i = 2, ni-2
                  tmp2(i,j,k) = rmf8v(i,j,2) * jcb8v(i,j,k) * h2(i,j,k)
                end do
              end do
              !$omp end do
            end do

          else

            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
                do i = 2, ni-1
                  tmp1(i,j,k) = rmf8u(i,j,2) * jcb8u(i,j,k) * h1(i,j,k)
                end do
              end do
              !$omp end do

              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-1
                do i = 2, ni-2
                  tmp2(i,j,k) = jcb8v(i,j,k) * h2(i,j,k)
                end do
              end do
              !$omp end do
            end do

          end if

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                tkefrc(i,j,k) = tkefrc(i,j,k) &
                  + ((tmp3(i,j,k+1) - tmp3(i,j,k)) * dziv2 &
                  + mf(i,j) * ((tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv2 &
                  + (tmp2(i,j+1,k) - tmp2(i,j,k)) * dyiv2))
              end do
            end do
            !$omp end do
          end do

        else

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-1
                tmp1(i,j,k) = rmf8u(i,j,2) * jcb8u(i,j,k) * h1(i,j,k)
              end do
            end do
            !$omp end do

            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-1
              do i = 2, ni-2
                tmp2(i,j,k) = rmf8v(i,j,2) * jcb8v(i,j,k) * h2(i,j,k)
              end do
            end do
            !$omp end do
          end do

          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
              do i = 2, ni-2
                tkefrc(i,j,k) = tkefrc(i,j,k) &
                  + ((tmp3(i,j,k+1) - tmp3(i,j,k)) * dziv2 &
                  + rmf(i,j,1) * ((tmp1(i+1,j,k) - tmp1(i,j,k)) * dxiv2 &
                  + (tmp2(i,j+1,k) - tmp2(i,j,k)) * dyiv2))
              end do
            end do
            !$omp end do
          end do

        end if

      end if

    end if

    !$omp end parallel

  end subroutine kernel_turbtke

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
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

  !=====================================================================
  ! Parameter reader
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk, trnopt, mpopt, mfcopt, &
                             dxiv2, dyiv2, dziv2)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    integer, intent(out) :: trnopt, mpopt, mfcopt
    real, intent(out) :: dxiv2, dyiv2, dziv2

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    ni = 1
    nj = 1
    nk = 1
    trnopt = 0
    mpopt = 0
    mfcopt = 0
    dxiv2 = 1.0
    dyiv2 = 1.0
    dziv2 = 1.0

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
          case ('trnopt')
            read(val, *) trnopt
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
          case ('dxiv2')
            read(val, *) dxiv2
          case ('dyiv2')
            read(val, *) dyiv2
          case ('dziv2')
            read(val, *) dziv2
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary array readers
  !=====================================================================
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

end program kernel_benchmark_turbtke

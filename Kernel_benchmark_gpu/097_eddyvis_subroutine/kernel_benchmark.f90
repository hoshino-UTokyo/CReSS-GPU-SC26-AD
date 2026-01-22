!***********************************************************************
! Kernel Benchmark: eddyvis (s_eddyvis) - GPU Version
!***********************************************************************
!
! Source: Src/eddyvis.f90
! Description: Calculate the eddy viscosity and turbulent length scale
!              using Smagorinsky or Deardorff (TKE-based) formulations.
!
!***********************************************************************
program kernel_benchmark_eddyvis
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: mpopt, mfcopt, sfcopt, tubopt, isoopt
  real :: dtb, kappa, dx, dy, dz
  real :: ds3, cpriv, lnh, cslnh2, csnum2, ckmax, khmax, khmin

  ! Physical constants (from m_commath and m_comphy)
  real, parameter :: oned3 = 1.e0/3.e0
  real, parameter :: eps = 1.e-20
  real, parameter :: ckm = 0.1e0
  real, parameter :: ckmin = 1.e-6

  ! Input arrays
  real, allocatable :: zph(:,:,:)
  real, allocatable :: jcb(:,:,:)
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: rbr(:,:,:)
  real, allocatable :: tke(:,:,:)
  real, allocatable :: ssq(:,:,:)
  real, allocatable :: nsq8w(:,:,:)

  ! Output arrays
  real, allocatable :: rkh(:,:,:)
  real, allocatable :: rkv(:,:,:)
  real, allocatable :: priv(:,:,:)

  ! Reference output for validation
  real, allocatable :: rkh_ref(:,:,:)
  real, allocatable :: rkv_ref(:,:,:)
  real, allocatable :: priv_ref(:,:,:)

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
       ni, nj, nk, mpopt, mfcopt, sfcopt, tubopt, isoopt, &
       dtb, kappa, dx, dy, dz, ds3, cpriv, lnh, cslnh2, csnum2, &
       ckmax, khmax, khmin)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: eddyvis (GPU)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' mpopt=', mpopt
  write(*,'(A,I6)') ' mfcopt=', mfcopt
  write(*,'(A,I6)') ' sfcopt=', sfcopt
  write(*,'(A,I6)') ' tubopt=', tubopt
  write(*,'(A,I6)') ' isoopt=', isoopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(zph(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(ssq(0:ni+1, 0:nj+1, 1:nk))
  allocate(nsq8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkh(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkv(0:ni+1, 0:nj+1, 1:nk))
  allocate(priv(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkh_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkv_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(priv_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/zph.bin', zph, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tke.bin', tke, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ssq.bin', ssq, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/nsq8w.bin', nsq8w, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/rkh_ref.bin', rkh_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rkv_ref.bin', rkv_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/priv_ref.bin', priv_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_eddyvis(mpopt, mfcopt, sfcopt, tubopt, isoopt, &
         ni, nj, nk, ds3, cpriv, lnh, cslnh2, csnum2, ckmax, khmax, khmin, &
         kappa, zph, jcb, rmf, rbr, tke, ssq, nsq8w, rkh, rkv, priv)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    !$acc wait
    t_start = omp_get_wtime()

    call kernel_eddyvis(mpopt, mfcopt, sfcopt, tubopt, isoopt, &
         ni, nj, nk, ds3, cpriv, lnh, cslnh2, csnum2, ckmax, khmax, khmin, &
         kappa, zph, jcb, rmf, rbr, tke, ssq, nsq8w, rkh, rkv, priv)

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

  ! Validate rkh
  do k = 1, nk-1
    do j = 1, nj-1
      do i = 1, ni-1
        rel_error = abs(rkh(i,j,k) - rkh_ref(i,j,k))
        if (abs(rkh_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(rkh_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate rkv
  do k = 1, nk-1
    do j = 1, nj-1
      do i = 1, ni-1
        rel_error = abs(rkv(i,j,k) - rkv_ref(i,j,k))
        if (abs(rkv_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(rkv_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  ! Validate priv
  do k = 1, nk-1
    do j = 1, nj-1
      do i = 1, ni-1
        rel_error = abs(priv(i,j,k) - priv_ref(i,j,k))
        if (abs(priv_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(priv_ref(i,j,k))
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
  deallocate(zph, jcb, rmf, rbr, tke, ssq, nsq8w)
  deallocate(rkh, rkv, priv, rkh_ref, rkv_ref, priv_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: eddyvis (GPU version)
  ! Calculate eddy viscosity and turbulent length scale
  !=====================================================================
  subroutine kernel_eddyvis(mpopt, mfcopt, sfcopt, tubopt, isoopt, &
       ni, nj, nk, ds3, cpriv, lnh, cslnh2, csnum2, ckmax, khmax, khmin, &
       kappa, zph, jcb, rmf, rbr, tke, ssq, nsq8w, rkh, rkv, priv)
    implicit none

    integer, intent(in) :: mpopt, mfcopt, sfcopt, tubopt, isoopt
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: ds3, cpriv, lnh, cslnh2, csnum2, ckmax, khmax, khmin
    real, intent(in) :: kappa

    real, intent(in) :: zph(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rmf(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tke(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ssq(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: nsq8w(0:ni+1, 0:nj+1, 1:nk)

    real, intent(out) :: rkh(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: rkv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: priv(0:ni+1, 0:nj+1, 1:nk)

    ! Local variables
    integer :: i, j, k
    real :: ln, ln0, ln02, ln2, htskp, nsq, stabc, a

    ! Smagorinsky formulation (tubopt=1)
    if (tubopt == 1) then

      ! Isotropic case
      if (isoopt == 1) then

        if (mfcopt == 0) then

          !$acc kernels
          !$acc loop independent collapse(3) private(ln, ln2)
          do k = 1, nk-1
            do j = 1, nj-1
              do i = 1, ni-1
                ln = exp(oned3*log(ds3*jcb(i,j,k)))
                ln2 = ln*ln
                priv(i,j,k) = cpriv
                rkv(i,j,k) = csnum2*ln2*sqrt(max(ssq(i,j,k) &
                     -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv, 0.e0))
                rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln2)
                rkh(i,j,k) = rkv(i,j,k)
              end do
            end do
          end do
          !$acc end kernels

        else

          if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then

            !$acc kernels
            !$acc loop independent collapse(3) private(ln, ln2)
            do k = 1, nk-1
              do j = 1, nj-1
                do i = 1, ni-1
                  ln = exp(oned3*log(ds3*rmf(i,j,2)*jcb(i,j,k)))
                  ln2 = ln*ln
                  priv(i,j,k) = cpriv
                  rkv(i,j,k) = csnum2*ln2*sqrt(max(ssq(i,j,k) &
                       -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv, 0.e0))
                  rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln2)
                  rkh(i,j,k) = rkv(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

          else

            !$acc kernels
            !$acc loop independent collapse(3) private(ln, ln2)
            do k = 1, nk-1
              do j = 1, nj-1
                do i = 1, ni-1
                  ln = exp(oned3*log(ds3*rmf(i,j,3)*jcb(i,j,k)))
                  ln2 = ln*ln
                  priv(i,j,k) = cpriv
                  rkv(i,j,k) = csnum2*ln2*sqrt(max(ssq(i,j,k) &
                       -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv, 0.e0))
                  rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln2)
                  rkh(i,j,k) = rkv(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

          end if

        end if

      ! Anisotropic case (isoopt=2)
      else if (isoopt == 2) then

        if (mfcopt == 0) then

          !$acc kernels
          !$acc loop independent collapse(3) private(ln, ln2, stabc)
          do k = 1, nk-1
            do j = 1, nj-1
              do i = 1, ni-1
                ln = zph(i,j,k+1) - zph(i,j,k)
                ln2 = ln*ln
                priv(i,j,k) = cpriv
                stabc = sqrt(max(ssq(i,j,k) &
                     -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv, 0.e0))
                rkv(i,j,k) = rbr(i,j,k)*min(stabc*csnum2, ckmax)*ln2
                rkh(i,j,k) = rbr(i,j,k)*min(stabc*cslnh2, khmax)
              end do
            end do
          end do
          !$acc end kernels

        else

          if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then

            !$acc kernels
            !$acc loop independent collapse(3) private(ln, ln2, stabc)
            do k = 1, nk-1
              do j = 1, nj-1
                do i = 1, ni-1
                  ln = zph(i,j,k+1) - zph(i,j,k)
                  ln2 = ln*ln
                  priv(i,j,k) = cpriv
                  stabc = sqrt(max(ssq(i,j,k) &
                       -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv, 0.e0))
                  rkv(i,j,k) = rbr(i,j,k)*min(stabc*csnum2, ckmax)*ln2
                  rkh(i,j,k) = rbr(i,j,k)*min(stabc*cslnh2, khmax)*rmf(i,j,2)
                end do
              end do
            end do
            !$acc end kernels

          else

            !$acc kernels
            !$acc loop independent collapse(3) private(ln, ln2, stabc)
            do k = 1, nk-1
              do j = 1, nj-1
                do i = 1, ni-1
                  ln = zph(i,j,k+1) - zph(i,j,k)
                  ln2 = ln*ln
                  priv(i,j,k) = cpriv
                  stabc = sqrt(max(ssq(i,j,k) &
                       -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv, 0.e0))
                  rkv(i,j,k) = rbr(i,j,k)*min(stabc*csnum2, ckmax)*ln2
                  rkh(i,j,k) = rbr(i,j,k)*min(stabc*cslnh2, khmax)*rmf(i,j,3)
                end do
              end do
            end do
            !$acc end kernels

          end if

        end if

      end if

    ! Deardorff formulation (tubopt/=1)
    else

      ! No surface process (sfcopt=0)
      if (sfcopt == 0) then

        ! Isotropic case
        if (isoopt == 1) then

          if (mfcopt == 0) then

            !$acc kernels
            !$acc loop independent collapse(3) private(ln, ln0, ln02, nsq)
            do k = 1, nk-1
              do j = 1, nj-1
                do i = 1, ni-1
                  ln0 = exp(oned3*log(ds3*jcb(i,j,k)))
                  ln02 = ln0*ln0
                  nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                  if (nsq < 0.e0) then
                    ln = ln0
                  else
                    ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                  end if
                  priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                  rkv(i,j,k) = ckm*ln*sqrt(tke(i,j,k))
                  if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                    rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                  end if
                  rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                  rkh(i,j,k) = rkv(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

          else

            if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then

              !$acc kernels
              !$acc loop independent collapse(3) private(ln, ln0, ln02, nsq)
              do k = 1, nk-1
                do j = 1, nj-1
                  do i = 1, ni-1
                    ln0 = exp(oned3*log(ds3*rmf(i,j,2)*jcb(i,j,k)))
                    ln02 = ln0*ln0
                    nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                    if (nsq < 0.e0) then
                      ln = ln0
                    else
                      ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                    end if
                    priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                    rkv(i,j,k) = ckm*ln*sqrt(tke(i,j,k))
                    if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                      rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                    end if
                    rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                    rkh(i,j,k) = rkv(i,j,k)
                  end do
                end do
              end do
              !$acc end kernels

            else

              !$acc kernels
              !$acc loop independent collapse(3) private(ln, ln0, ln02, nsq)
              do k = 1, nk-1
                do j = 1, nj-1
                  do i = 1, ni-1
                    ln0 = exp(oned3*log(ds3*rmf(i,j,3)*jcb(i,j,k)))
                    ln02 = ln0*ln0
                    nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                    if (nsq < 0.e0) then
                      ln = ln0
                    else
                      ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                    end if
                    priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                    rkv(i,j,k) = ckm*ln*sqrt(tke(i,j,k))
                    if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                      rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                    end if
                    rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                    rkh(i,j,k) = rkv(i,j,k)
                  end do
                end do
              end do
              !$acc end kernels

            end if

          end if

        ! Anisotropic case (isoopt=2)
        else if (isoopt == 2) then

          if (mfcopt == 0) then

            !$acc kernels
            !$acc loop independent collapse(3) private(ln, ln0, ln02, nsq, a)
            do k = 1, nk-1
              do j = 1, nj-1
                do i = 1, ni-1
                  ln0 = zph(i,j,k+1) - zph(i,j,k)
                  ln02 = ln0*ln0
                  nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                  if (nsq < 0.e0) then
                    ln = ln0
                  else
                    ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                  end if
                  priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                  a = ckm*sqrt(tke(i,j,k))
                  rkv(i,j,k) = a*ln
                  rkh(i,j,k) = a*lnh
                  if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                    rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                    rkh(i,j,k) = max(rkh(i,j,k), khmin)
                  end if
                  rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                  rkh(i,j,k) = rbr(i,j,k)*min(rkh(i,j,k), khmax)
                end do
              end do
            end do
            !$acc end kernels

          else

            if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then

              !$acc kernels
              !$acc loop independent collapse(3) private(ln, ln0, ln02, nsq, a)
              do k = 1, nk-1
                do j = 1, nj-1
                  do i = 1, ni-1
                    ln0 = zph(i,j,k+1) - zph(i,j,k)
                    ln02 = ln0*ln0
                    nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                    if (nsq < 0.e0) then
                      ln = ln0
                    else
                      ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                    end if
                    priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                    a = ckm*sqrt(tke(i,j,k))
                    rkv(i,j,k) = a*ln
                    rkh(i,j,k) = a*lnh*rmf(i,j,4)
                    if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                      rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                      rkh(i,j,k) = max(rkh(i,j,k), khmin*rmf(i,j,2))
                    end if
                    rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                    rkh(i,j,k) = rbr(i,j,k)*min(rkh(i,j,k), khmax*rmf(i,j,2))
                  end do
                end do
              end do
              !$acc end kernels

            else

              !$acc kernels
              !$acc loop independent collapse(3) private(ln, ln0, ln02, nsq, a)
              do k = 1, nk-1
                do j = 1, nj-1
                  do i = 1, ni-1
                    ln0 = zph(i,j,k+1) - zph(i,j,k)
                    ln02 = ln0*ln0
                    nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                    if (nsq < 0.e0) then
                      ln = ln0
                    else
                      ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                    end if
                    priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                    a = ckm*sqrt(tke(i,j,k))
                    rkv(i,j,k) = a*ln
                    rkh(i,j,k) = a*lnh*rmf(i,j,2)
                    if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                      rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                      rkh(i,j,k) = max(rkh(i,j,k), khmin*rmf(i,j,3))
                    end if
                    rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                    rkh(i,j,k) = rbr(i,j,k)*min(rkh(i,j,k), khmax*rmf(i,j,3))
                  end do
                end do
              end do
              !$acc end kernels

            end if

          end if

        end if

      ! With surface process (sfcopt/=0)
      else

        ! Isotropic case
        if (isoopt == 1) then

          if (mfcopt == 0) then

            !$acc kernels
            !$acc loop independent collapse(3) private(ln, ln0, ln02, htskp, nsq)
            do k = 1, nk-1
              do j = 1, nj-1
                do i = 1, ni-1
                  ln0 = exp(oned3*log(ds3*jcb(i,j,k)))
                  ln02 = ln0*ln0
                  htskp = kappa*abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))
                  nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                  if (nsq < 0.e0) then
                    ln = ln0
                  else
                    ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                  end if
                  ln = ln*htskp/(htskp+ln)
                  priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                  rkv(i,j,k) = ckm*ln*sqrt(tke(i,j,k))
                  if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                    rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                  end if
                  rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                  rkh(i,j,k) = rkv(i,j,k)
                end do
              end do
            end do
            !$acc end kernels

          else

            if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then

              !$acc kernels
              !$acc loop independent collapse(3) private(ln, ln0, ln02, htskp, nsq)
              do k = 1, nk-1
                do j = 1, nj-1
                  do i = 1, ni-1
                    ln0 = exp(oned3*log(ds3*rmf(i,j,2)*jcb(i,j,k)))
                    ln02 = ln0*ln0
                    htskp = kappa*abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))
                    nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                    if (nsq < 0.e0) then
                      ln = ln0
                    else
                      ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                    end if
                    ln = ln*htskp/(htskp+ln)
                    priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                    rkv(i,j,k) = ckm*ln*sqrt(tke(i,j,k))
                    if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                      rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                    end if
                    rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                    rkh(i,j,k) = rkv(i,j,k)
                  end do
                end do
              end do
              !$acc end kernels

            else

              !$acc kernels
              !$acc loop independent collapse(3) private(ln, ln0, ln02, htskp, nsq)
              do k = 1, nk-1
                do j = 1, nj-1
                  do i = 1, ni-1
                    ln0 = exp(oned3*log(ds3*rmf(i,j,3)*jcb(i,j,k)))
                    ln02 = ln0*ln0
                    htskp = kappa*abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))
                    nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                    if (nsq < 0.e0) then
                      ln = ln0
                    else
                      ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                    end if
                    ln = ln*htskp/(htskp+ln)
                    priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                    rkv(i,j,k) = ckm*ln*sqrt(tke(i,j,k))
                    if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                      rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                    end if
                    rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                    rkh(i,j,k) = rkv(i,j,k)
                  end do
                end do
              end do
              !$acc end kernels

            end if

          end if

        ! Anisotropic case (isoopt=2)
        else if (isoopt == 2) then

          if (mfcopt == 0) then

            !$acc kernels
            !$acc loop independent collapse(3) private(ln, ln0, ln02, htskp, nsq, a)
            do k = 1, nk-1
              do j = 1, nj-1
                do i = 1, ni-1
                  ln0 = zph(i,j,k+1) - zph(i,j,k)
                  ln02 = ln0*ln0
                  htskp = kappa*abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))
                  nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                  if (nsq < 0.e0) then
                    ln = ln0
                  else
                    ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                  end if
                  ln = ln*htskp/(htskp+ln)
                  priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                  a = ckm*sqrt(tke(i,j,k))
                  rkv(i,j,k) = a*ln
                  rkh(i,j,k) = a*lnh
                  if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                    rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                    rkh(i,j,k) = max(rkh(i,j,k), khmin)
                  end if
                  rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                  rkh(i,j,k) = rbr(i,j,k)*min(rkh(i,j,k), khmax)
                end do
              end do
            end do
            !$acc end kernels

          else

            if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then

              !$acc kernels
              !$acc loop independent collapse(3) private(ln, ln0, ln02, htskp, nsq, a)
              do k = 1, nk-1
                do j = 1, nj-1
                  do i = 1, ni-1
                    ln0 = zph(i,j,k+1) - zph(i,j,k)
                    ln02 = ln0*ln0
                    htskp = kappa*abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))
                    nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                    if (nsq < 0.e0) then
                      ln = ln0
                    else
                      ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                    end if
                    ln = ln*htskp/(htskp+ln)
                    priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                    a = ckm*sqrt(tke(i,j,k))
                    rkv(i,j,k) = a*ln
                    rkh(i,j,k) = a*lnh*rmf(i,j,4)
                    if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                      rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                      rkh(i,j,k) = max(rkh(i,j,k), khmin*rmf(i,j,2))
                    end if
                    rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                    rkh(i,j,k) = rbr(i,j,k)*min(rkh(i,j,k), khmax*rmf(i,j,2))
                  end do
                end do
              end do
              !$acc end kernels

            else

              !$acc kernels
              !$acc loop independent collapse(3) private(ln, ln0, ln02, htskp, nsq, a)
              do k = 1, nk-1
                do j = 1, nj-1
                  do i = 1, ni-1
                    ln0 = zph(i,j,k+1) - zph(i,j,k)
                    ln02 = ln0*ln0
                    htskp = kappa*abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))
                    nsq = nsq8w(i,j,k) + nsq8w(i,j,k+1)
                    if (nsq < 0.e0) then
                      ln = ln0
                    else
                      ln = max(.1e0*ln0, min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)), ln0))
                    end if
                    ln = ln*htskp/(htskp+ln)
                    priv(i,j,k) = 1.e0 + 2.e0*ln/ln0
                    a = ckm*sqrt(tke(i,j,k))
                    rkv(i,j,k) = a*ln
                    rkh(i,j,k) = a*lnh*rmf(i,j,2)
                    if (priv(i,j,k)*nsq < ssq(i,j,k)) then
                      rkv(i,j,k) = max(rkv(i,j,k), ckmin*ln02)
                      rkh(i,j,k) = max(rkh(i,j,k), khmin*rmf(i,j,3))
                    end if
                    rkv(i,j,k) = rbr(i,j,k)*min(rkv(i,j,k), ckmax*ln02)
                    rkh(i,j,k) = rbr(i,j,k)*min(rkh(i,j,k), khmax*rmf(i,j,3))
                  end do
                end do
              end do
              !$acc end kernels

            end if

          end if

        end if

      end if

    end if

  end subroutine kernel_eddyvis

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

    ! Default values
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
  subroutine read_parameters(filename, ni, nj, nk, mpopt, mfcopt, sfcopt, &
       tubopt, isoopt, dtb, kappa, dx, dy, dz, ds3, cpriv, lnh, cslnh2, &
       csnum2, ckmax, khmax, khmin)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk, mpopt, mfcopt, sfcopt, tubopt, isoopt
    real, intent(out) :: dtb, kappa, dx, dy, dz, ds3, cpriv, lnh
    real, intent(out) :: cslnh2, csnum2, ckmax, khmax, khmin

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
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
          case ('sfcopt')
            read(val, *) sfcopt
          case ('tubopt')
            read(val, *) tubopt
          case ('isoopt')
            read(val, *) isoopt
          case ('dtb')
            read(val, *) dtb
          case ('kappa')
            read(val, *) kappa
          case ('dx')
            read(val, *) dx
          case ('dy')
            read(val, *) dy
          case ('dz')
            read(val, *) dz
          case ('ds3')
            read(val, *) ds3
          case ('cpriv')
            read(val, *) cpriv
          case ('lnh')
            read(val, *) lnh
          case ('cslnh2')
            read(val, *) cslnh2
          case ('csnum2')
            read(val, *) csnum2
          case ('ckmax')
            read(val, *) ckmax
          case ('khmax')
            read(val, *) khmax
          case ('khmin')
            read(val, *) khmin
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

end program kernel_benchmark_eddyvis

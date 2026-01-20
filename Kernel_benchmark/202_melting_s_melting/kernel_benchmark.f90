!***********************************************************************
! Kernel Benchmark: melting (s_melting)
!***********************************************************************
!
! Source: Src/melting.f90
! Description: Calculate the melting rate from the cloud ice to the cloud water
!              and from the snow and graupel to the rain water.
!
! GPU Difficulty: Easy
!   - Uses intrinsic max() function - GPU compatible
!   - Conditional logic based on thresq threshold and t0cel temperature
!   - No synchronization constructs
!
!***********************************************************************
program kernel_benchmark_melting
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  real :: dtb, thresq, cw, t0cel, cc

  ! Derived parameter
  real :: cc2dt

  ! Input arrays
  real, allocatable :: rbr(:,:,:)
  real, allocatable :: rbv(:,:,:)
  real, allocatable :: qi(:,:,:)
  real, allocatable :: qs(:,:,:)
  real, allocatable :: qg(:,:,:)
  real, allocatable :: tcel(:,:,:)
  real, allocatable :: qvsst0(:,:,:)
  real, allocatable :: lv(:,:,:)
  real, allocatable :: lf(:,:,:)
  real, allocatable :: kp(:,:,:)
  real, allocatable :: dv(:,:,:)
  real, allocatable :: vnts(:,:,:)
  real, allocatable :: vntg(:,:,:)
  real, allocatable :: clcs(:,:,:)
  real, allocatable :: clcg(:,:,:)
  real, allocatable :: clrs(:,:,:)
  real, allocatable :: clrg(:,:,:)

  ! Output arrays
  real, allocatable :: mlic(:,:,:)
  real, allocatable :: mlsr(:,:,:)
  real, allocatable :: mlgr(:,:,:)

  ! Reference output for validation
  real, allocatable :: mlic_ref(:,:,:)
  real, allocatable :: mlsr_ref(:,:,:)
  real, allocatable :: mlgr_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
  integer :: error_count, total_errors
  logical :: validation_passed

  ! Loop variables
  integer :: iter

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, dtb, thresq, cw, t0cel, cc)

  ! Compute derived parameter
  cc2dt = 2.0e0 * cc * dtb

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: melting'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,ES12.4)') ' dtb=', dtb
  write(*,'(A,ES12.4)') ' thresq=', thresq
  write(*,'(A,ES12.4)') ' cw=', cw
  write(*,'(A,ES12.4)') ' t0cel=', t0cel
  write(*,'(A,ES12.4)') ' cc=', cc
  write(*,'(A,ES12.4)') ' cc2dt=', cc2dt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(qg(0:ni+1, 0:nj+1, 1:nk))
  allocate(tcel(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsst0(0:ni+1, 0:nj+1, 1:nk))
  allocate(lv(0:ni+1, 0:nj+1, 1:nk))
  allocate(lf(0:ni+1, 0:nj+1, 1:nk))
  allocate(kp(0:ni+1, 0:nj+1, 1:nk))
  allocate(dv(0:ni+1, 0:nj+1, 1:nk))
  allocate(vnts(0:ni+1, 0:nj+1, 1:nk))
  allocate(vntg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrs(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrg(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlic(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlsr(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlgr(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlic_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlsr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(mlgr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rbv.bin', rbv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qi.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qg.bin', qg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tcel.bin', tcel, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvsst0.bin', qvsst0, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/lv.bin', lv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/lf.bin', lf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/kp.bin', kp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/dv.bin', dv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vnts.bin', vnts, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vntg.bin', vntg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clcs.bin', clcs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clcg.bin', clcg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clrs.bin', clrs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clrg.bin', clrg, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/mlic_ref.bin', mlic_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/mlsr_ref.bin', mlsr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/mlgr_ref.bin', mlgr_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_melting(cc2dt, thresq, cw, t0cel, ni, nj, nk, &
         rbr, rbv, qi, qs, qg, tcel, qvsst0, lv, lf, kp, dv, &
         vnts, vntg, clcs, clcg, clrs, clrg, mlic, mlsr, mlgr)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    t_start = omp_get_wtime()

    call kernel_melting(cc2dt, thresq, cw, t0cel, ni, nj, nk, &
         rbr, rbv, qi, qs, qg, tcel, qvsst0, lv, lf, kp, dv, &
         vnts, vntg, clcs, clcg, clrs, clrg, mlic, mlsr, mlgr)

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
  max_error = 0.0
  total_errors = 0

  call validate_output(mlic, mlic_ref, ni, nj, nk, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_output(mlsr, mlsr_ref, ni, nj, nk, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_output(mlgr, mlgr_ref, ni, nj, nk, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  validation_passed = (total_errors == 0)

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
  write(*,'(A,I12)')    ' Error count:        ', total_errors
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(rbr, rbv, qi, qs, qg, tcel, qvsst0, lv, lf, kp, dv)
  deallocate(vnts, vntg, clcs, clcg, clrs, clrg)
  deallocate(mlic, mlsr, mlgr)
  deallocate(mlic_ref, mlsr_ref, mlgr_ref)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: melting - calculate melting rate
  !-------------------------------------------------------------------
  subroutine kernel_melting(cc2dt, thresq, cw, t0cel, ni, nj, nk, &
       rbr, rbv, qi, qs, qg, tcel, qvsst0, lv, lf, kp, dv, &
       vnts, vntg, clcs, clcg, clrs, clrg, mlic, mlsr, mlgr)
    implicit none

    real, intent(in) :: cc2dt, thresq, cw, t0cel
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: rbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rbv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qi(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tcel(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: qvsst0(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: lv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: lf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: kp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: dv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vnts(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vntg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clcs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clcg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clrs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: clrg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: mlic(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: mlsr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(out) :: mlgr(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    real :: cmlxr1, cmlxr2

    !$omp parallel default(shared) private(k)

    ! In the case nk = 1.
    if (nk == 1) then

      !$omp do schedule(runtime) private(i,j,cmlxr1,cmlxr2)
      do j = 1, nj-1
        do i = 1, ni-1

          ! Set the common used variables.
          cmlxr1 = cc2dt * rbv(i,j,1) * (kp(i,j,1) * tcel(i,j,1) &
               + lv(i,j,1) * dv(i,j,1) * rbr(i,j,1) * qvsst0(i,j,1))
          cmlxr2 = cw * tcel(i,j,1)

          ! Calculate the melting rate from the cloud ice to the cloud water.
          if (qi(i,j,1) > thresq) then
            if (tcel(i,j,1) >= t0cel) then
              mlic(i,j,1) = qi(i,j,1)
            else
              mlic(i,j,1) = 0.0e0
            end if
          else
            mlic(i,j,1) = 0.0e0
          end if

          ! Calculate the melting rate from the snow to the rain water.
          if (qs(i,j,1) > thresq) then
            if (tcel(i,j,1) >= t0cel) then
              mlsr(i,j,1) = max((cmlxr1 * vnts(i,j,1) &
                   + cmlxr2 * (clcs(i,j,1) + clrs(i,j,1))) / lf(i,j,1), 0.0e0)
            else
              mlsr(i,j,1) = 0.0e0
            end if
          else
            mlsr(i,j,1) = 0.0e0
          end if

          ! Calculate the melting rate from the graupel to the rain water.
          if (qg(i,j,1) > thresq) then
            if (tcel(i,j,1) >= t0cel) then
              mlgr(i,j,1) = max((cmlxr1 * vntg(i,j,1) &
                   + cmlxr2 * (clcg(i,j,1) + clrg(i,j,1))) / lf(i,j,1), 0.0e0)
            else
              mlgr(i,j,1) = 0.0e0
            end if
          else
            mlgr(i,j,1) = 0.0e0
          end if

        end do
      end do
      !$omp end do

    ! In the case nk > 1.
    else

      do k = 1, nk-1

        !$omp do schedule(runtime) private(i,j,cmlxr1,cmlxr2)
        do j = 1, nj-1
          do i = 1, ni-1

            ! Set the common used variables.
            cmlxr1 = cc2dt * rbv(i,j,k) * (kp(i,j,k) * tcel(i,j,k) &
                 + lv(i,j,k) * dv(i,j,k) * rbr(i,j,k) * qvsst0(i,j,k))
            cmlxr2 = cw * tcel(i,j,k)

            ! Calculate the melting rate from the cloud ice to the cloud water.
            if (qi(i,j,k) > thresq) then
              if (tcel(i,j,k) >= t0cel) then
                mlic(i,j,k) = qi(i,j,k)
              else
                mlic(i,j,k) = 0.0e0
              end if
            else
              mlic(i,j,k) = 0.0e0
            end if

            ! Calculate the melting rate from the snow to the rain water.
            if (qs(i,j,k) > thresq) then
              if (tcel(i,j,k) >= t0cel) then
                mlsr(i,j,k) = max((cmlxr1 * vnts(i,j,k) &
                     + cmlxr2 * (clcs(i,j,k) + clrs(i,j,k))) / lf(i,j,k), 0.0e0)
              else
                mlsr(i,j,k) = 0.0e0
              end if
            else
              mlsr(i,j,k) = 0.0e0
            end if

            ! Calculate the melting rate from the graupel to the rain water.
            if (qg(i,j,k) > thresq) then
              if (tcel(i,j,k) >= t0cel) then
                mlgr(i,j,k) = max((cmlxr1 * vntg(i,j,k) &
                     + cmlxr2 * (clcg(i,j,k) + clrg(i,j,k))) / lf(i,j,k), 0.0e0)
              else
                mlgr(i,j,k) = 0.0e0
              end if
            else
              mlgr(i,j,k) = 0.0e0
            end if

          end do
        end do
        !$omp end do

      end do

    end if

    !$omp end parallel

  end subroutine kernel_melting

  !-------------------------------------------------------------------
  ! Read benchmark configuration
  !-------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios

    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      data_dir = './data'
      num_iter = 10
      warmup_iter = 2
      tol = 1.0e-5
      return
    end if
    read(10, '(A)') data_dir
    read(10, *) num_iter
    read(10, *) warmup_iter
    read(10, *) tol
    close(10)
  end subroutine read_config

  !-------------------------------------------------------------------
  ! Read parameters (key = value format)
  !-------------------------------------------------------------------
  subroutine read_parameters(filename, ni, nj, nk, dtb, thresq, cw, t0cel, cc)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    real, intent(out) :: dtb, thresq, cw, t0cel, cc

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
          case ('dtb')
            read(val, *) dtb
          case ('thresq')
            read(val, *) thresq
          case ('cw')
            read(val, *) cw
          case ('t0cel')
            read(val, *) t0cel
          case ('cc')
            read(val, *) cc
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

    open(unit=10, file=filename, access='stream', form='unformatted', &
         status='old', iostat=ios, convert='big_endian')
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
  subroutine validate_output(output, reference, ni, nj, nk, tol, &
       max_err, err_count)
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

    do k = 1, nk
      do j = 1, nj-1
        do i = 1, ni-1
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

end program kernel_benchmark_melting

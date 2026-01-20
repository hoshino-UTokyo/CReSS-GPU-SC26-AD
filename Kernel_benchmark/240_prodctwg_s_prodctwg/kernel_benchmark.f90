!***********************************************************************
! Kernel Benchmark: prodctwg
!***********************************************************************
! Source: Src/prodctwg.f90
! Subroutine: s_prodctwg
! Profile ID: 240
! Description: Calculate graupel production rate for wet/dry growth
!***********************************************************************
program kernel_benchmark_prodctwg
  use omp_lib
  implicit none

  ! === Grid dimensions ===
  integer :: cphopt, ni, nj, nk

  ! === Scalar parameters ===
  real :: dtb, thresq, ci, cw, cc2dtn, eigiv, esgiv

  ! === Physical constants ===
  real, parameter :: t0cel = 0.e0   ! Melting point in Celsius

  ! === Arrays (allocatable) ===
  real, allocatable :: rbr(:,:,:)
  real, allocatable :: qi(:,:,:), qs(:,:,:), qg(:,:,:), ncs(:,:,:)
  real, allocatable :: tcel(:,:,:), qvsst0(:,:,:)
  real, allocatable :: lv(:,:,:), lf(:,:,:), kp(:,:,:), dv(:,:,:)
  real, allocatable :: vntg(:,:,:), clcg(:,:,:), clrg(:,:,:)
  ! Input/output arrays
  real, allocatable :: clir(:,:,:), clis(:,:,:), clig(:,:,:)
  real, allocatable :: clsr(:,:,:), clsg(:,:,:)
  real, allocatable :: clsrn(:,:,:), clsgn(:,:,:)
  real, allocatable :: pgwet(:,:,:)
  ! Reference arrays for validation
  real, allocatable :: clir_ref(:,:,:), clis_ref(:,:,:), clig_ref(:,:,:)
  real, allocatable :: clsr_ref(:,:,:), clsg_ref(:,:,:)
  real, allocatable :: clsrn_ref(:,:,:), clsgn_ref(:,:,:)
  real, allocatable :: pgwet_ref(:,:,:)
  ! Backup arrays for iterations
  real, allocatable :: clir_in(:,:,:), clis_in(:,:,:), clig_in(:,:,:)
  real, allocatable :: clsr_in(:,:,:), clsg_in(:,:,:)
  real, allocatable :: clsrn_in(:,:,:), clsgn_in(:,:,:)

  ! === Benchmark control ===
  integer :: num_iterations, warmup_iterations, iter
  real :: tolerance
  character(len=512) :: data_dir

  ! === Timing ===
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max

  ! === Validation ===
  logical :: passed
  real :: max_err, max_err_total
  integer :: err_count, err_count_total

  ! === Main program flow ===

  ! 1. Read configuration
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! 2. Read parameters
  call read_parameters(trim(data_dir)//'/params.txt', cphopt, ni, nj, nk, &
       dtb, thresq, ci, cw, cc2dtn, eigiv, esgiv)

  ! 3. Allocate arrays
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(qg(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncs(0:ni+1, 0:nj+1, 1:nk))
  allocate(tcel(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsst0(0:ni+1, 0:nj+1, 1:nk))
  allocate(lv(0:ni+1, 0:nj+1, 1:nk))
  allocate(lf(0:ni+1, 0:nj+1, 1:nk))
  allocate(kp(0:ni+1, 0:nj+1, 1:nk))
  allocate(dv(0:ni+1, 0:nj+1, 1:nk))
  allocate(vntg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clcg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clrg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clir(0:ni+1, 0:nj+1, 1:nk))
  allocate(clis(0:ni+1, 0:nj+1, 1:nk))
  allocate(clig(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsr(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsg(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsrn(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsgn(0:ni+1, 0:nj+1, 1:nk))
  allocate(pgwet(0:ni+1, 0:nj+1, 1:nk))
  ! Reference and backup
  allocate(clir_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clis_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clig_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsrn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsgn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(pgwet_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(clir_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(clis_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(clig_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsr_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsg_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsrn_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(clsgn_in(0:ni+1, 0:nj+1, 1:nk))

  ! 4. Read input data
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qi.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qg.bin', qg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncs.bin', ncs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tcel.bin', tcel, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvsst0.bin', qvsst0, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/lv.bin', lv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/lf.bin', lf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/kp.bin', kp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/dv.bin', dv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vntg.bin', vntg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clcg.bin', clcg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clrg.bin', clrg, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clir_in.bin', clir_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clis_in.bin', clis_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clig_in.bin', clig_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsr_in.bin', clsr_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsg_in.bin', clsg_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsrn_in.bin', clsrn_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsgn_in.bin', clsgn_in, 0, ni+1, 0, nj+1, 1, nk)

  ! 5. Read reference output
  call read_array_3d(trim(data_dir)//'/pgwet_ref.bin', pgwet_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clir_ref.bin', clir_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clis_ref.bin', clis_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clig_ref.bin', clig_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsr_ref.bin', clsr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsg_ref.bin', clsg_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsrn_ref.bin', clsrn_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/clsgn_ref.bin', clsgn_ref, 0, ni+1, 0, nj+1, 1, nk)

  ! Print header
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: prodctwg'
  write(*,'(A)') '=================================================='
  write(*,'(A,I8,A,I8,A,I8)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I8)') ' cphopt: ', cphopt
  write(*,'(A,I8)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  ! 6. Warmup iterations
  do iter = 1, warmup_iterations
    ! Reset inout arrays
    clir = clir_in; clis = clis_in; clig = clig_in
    clsr = clsr_in; clsg = clsg_in
    clsrn = clsrn_in; clsgn = clsgn_in
    call kernel_prodctwg(cphopt, ni, nj, nk, dtb, thresq, ci, cw, cc2dtn, eigiv, esgiv, &
         rbr, qi, qs, qg, ncs, tcel, qvsst0, lv, lf, kp, dv, vntg, clcg, clrg, &
         clir, clis, clig, clsr, clsg, clsrn, clsgn, pgwet)
  end do

  ! 7. Benchmark iterations (timed)
  t_total = 0.0d0
  t_min = huge(1.0d0)
  t_max = 0.0d0

  do iter = 1, num_iterations
    ! Reset inout arrays
    clir = clir_in; clis = clis_in; clig = clig_in
    clsr = clsr_in; clsg = clsg_in
    clsrn = clsrn_in; clsgn = clsgn_in

    t_start = omp_get_wtime()
    call kernel_prodctwg(cphopt, ni, nj, nk, dtb, thresq, ci, cw, cc2dtn, eigiv, esgiv, &
         rbr, qi, qs, qg, ncs, tcel, qvsst0, lv, lf, kp, dv, vntg, clcg, clrg, &
         clir, clis, clig, clsr, clsg, clsrn, clsgn, pgwet)
    t_end = omp_get_wtime()

    t_total = t_total + (t_end - t_start)
    if (t_end - t_start < t_min) t_min = t_end - t_start
    if (t_end - t_start > t_max) t_max = t_end - t_start
  end do
  t_avg = t_total / dble(num_iterations)

  ! 8. Validate output (only interior points processed by kernel: 1:ni-1, 1:nj-1, 1:nk-1 or 1 for nk=1)
  max_err_total = 0.0
  err_count_total = 0

  call validate_output_3d(pgwet, pgwet_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(clir, clir_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(clis, clis_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(clig, clig_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(clsr, clsr_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(clsg, clsg_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(clsrn, clsrn_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  call validate_output_3d(clsgn, clsgn_ref, tolerance, max_err, err_count, 1, ni-1, 1, nj-1, 1, max(nk-1,1))
  if (max_err > max_err_total) max_err_total = max_err
  err_count_total = err_count_total + err_count

  passed = (err_count_total == 0)

  ! 9. Report results
  write(*,'(A,F12.3,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:     ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:     ', t_max * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_err_total
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Total error count:  ', err_count_total
  if (passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  ! 10. Cleanup
  deallocate(rbr, qi, qs, qg, ncs, tcel, qvsst0, lv, lf, kp, dv, vntg, clcg, clrg)
  deallocate(clir, clis, clig, clsr, clsg, clsrn, clsgn, pgwet)
  deallocate(clir_ref, clis_ref, clig_ref, clsr_ref, clsg_ref, clsrn_ref, clsgn_ref, pgwet_ref)
  deallocate(clir_in, clis_in, clig_in, clsr_in, clsg_in, clsrn_in, clsgn_in)

contains

  !-----------------------------------------------------------------
  ! The kernel (extracted from original source)
  !-----------------------------------------------------------------
  subroutine kernel_prodctwg(cphopt, ni, nj, nk, dtb, thresq, ci, cw, cc2dtn, eigiv, esgiv, &
       rbr, qi, qs, qg, ncs, tcel, qvsst0, lv, lf, kp, dv, vntg, clcg, clrg, &
       clir, clis, clig, clsr, clsg, clsrn, clsgn, pgwet)
    implicit none

    integer, intent(in) :: cphopt, ni, nj, nk
    real, intent(in) :: dtb, thresq, ci, cw, cc2dtn, eigiv, esgiv
    real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: tcel(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qvsst0(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: lv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: lf(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: kp(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: dv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: vntg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: clcg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: clrg(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: clir(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: clis(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: clig(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: clsr(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: clsg(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: clsrn(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: clsgn(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: pgwet(0:ni+1,0:nj+1,1:nk)

    integer :: i, j, k
    real :: pgdry, lfice, cligw, clsgw, sink, a

    !$omp parallel default(shared) private(k)

    if (nk == 1) then

      if (abs(cphopt) == 2) then
        !$omp do schedule(runtime) private(i,j,pgdry,lfice,cligw,clsgw,sink,a)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qg(i,j,1) > thresq) then
              if (tcel(i,j,1) < t0cel) then
                pgdry = clcg(i,j,1) + clrg(i,j,1) + clig(i,j,1) + clsg(i,j,1)
                lfice = lf(i,j,1) + cw*tcel(i,j,1)
                cligw = eigiv*clig(i,j,1)
                clsgw = esgiv*clsg(i,j,1)
                pgwet(i,j,1) = cc2dtn*vntg(i,j,1) &
                     *(lv(i,j,1)*dv(i,j,1)*rbr(i,j,1)*qvsst0(i,j,1) &
                     +kp(i,j,1)*tcel(i,j,1))/(lfice*rbr(i,j,1)) &
                     +(cligw+clsgw)*(1.e0-ci*tcel(i,j,1)/lfice)
                if (pgwet(i,j,1) > 0.e0 .and. pgwet(i,j,1) < pgdry) then
                  clig(i,j,1) = cligw
                  clsg(i,j,1) = clsgw
                  sink = clir(i,j,1) + clis(i,j,1) + clig(i,j,1)
                  if (qi(i,j,1) < sink) then
                    a = qi(i,j,1)/sink
                    clir(i,j,1) = clir(i,j,1)*a
                    clis(i,j,1) = clis(i,j,1)*a
                    clig(i,j,1) = clig(i,j,1)*a
                  end if
                  sink = clsr(i,j,1) + clsg(i,j,1)
                  if (qs(i,j,1) < sink) then
                    a = qs(i,j,1)/sink
                    clsr(i,j,1) = clsr(i,j,1)*a
                    clsg(i,j,1) = clsg(i,j,1)*a
                  end if
                else
                  pgwet(i,j,1) = -1.e0
                end if
              else
                pgwet(i,j,1) = -1.e0
              end if
            else
              pgwet(i,j,1) = -1.e0
            end if
          end do
        end do
        !$omp end do

      else if (abs(cphopt) >= 3) then
        !$omp do schedule(runtime) private(i,j,pgdry,lfice,cligw,clsgw,sink,a)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qg(i,j,1) > thresq) then
              if (tcel(i,j,1) < t0cel) then
                pgdry = clcg(i,j,1) + clrg(i,j,1) + clig(i,j,1) + clsg(i,j,1)
                lfice = lf(i,j,1) + cw*tcel(i,j,1)
                cligw = eigiv*clig(i,j,1)
                clsgw = esgiv*clsg(i,j,1)
                pgwet(i,j,1) = cc2dtn*vntg(i,j,1) &
                     *(lv(i,j,1)*dv(i,j,1)*rbr(i,j,1)*qvsst0(i,j,1) &
                     +kp(i,j,1)*tcel(i,j,1))/(lfice*rbr(i,j,1)) &
                     +(cligw+clsgw)*(1.e0-ci*tcel(i,j,1)/lfice)
                if (pgwet(i,j,1) > 0.e0 .and. pgwet(i,j,1) < pgdry) then
                  clig(i,j,1) = cligw
                  clsg(i,j,1) = clsgw
                  clsgn(i,j,1) = esgiv*clsgn(i,j,1)
                  sink = clir(i,j,1) + clis(i,j,1) + clig(i,j,1)
                  if (qi(i,j,1) < sink) then
                    a = qi(i,j,1)/sink
                    clir(i,j,1) = clir(i,j,1)*a
                    clis(i,j,1) = clis(i,j,1)*a
                    clig(i,j,1) = clig(i,j,1)*a
                  end if
                  sink = clsr(i,j,1) + clsg(i,j,1)
                  if (qs(i,j,1) < sink) then
                    a = qs(i,j,1)/sink
                    clsr(i,j,1) = clsr(i,j,1)*a
                    clsg(i,j,1) = clsg(i,j,1)*a
                  end if
                  sink = clsrn(i,j,1) + clsgn(i,j,1)
                  if (ncs(i,j,1) < sink) then
                    a = ncs(i,j,1)/sink
                    clsrn(i,j,1) = clsrn(i,j,1)*a
                    clsgn(i,j,1) = clsgn(i,j,1)*a
                  end if
                else
                  pgwet(i,j,1) = -1.e0
                end if
              else
                pgwet(i,j,1) = -1.e0
              end if
            else
              pgwet(i,j,1) = -1.e0
            end if
          end do
        end do
        !$omp end do
      end if

    else
      ! nk > 1

      if (abs(cphopt) == 2) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j,pgdry,lfice,cligw,clsgw,sink,a)
          do j = 1, nj-1
            do i = 1, ni-1
              if (qg(i,j,k) > thresq) then
                if (tcel(i,j,k) < t0cel) then
                  pgdry = clcg(i,j,k) + clrg(i,j,k) + clig(i,j,k) + clsg(i,j,k)
                  lfice = lf(i,j,k) + cw*tcel(i,j,k)
                  cligw = eigiv*clig(i,j,k)
                  clsgw = esgiv*clsg(i,j,k)
                  pgwet(i,j,k) = cc2dtn*vntg(i,j,k) &
                       *(lv(i,j,k)*dv(i,j,k)*rbr(i,j,k)*qvsst0(i,j,k) &
                       +kp(i,j,k)*tcel(i,j,k))/(lfice*rbr(i,j,k)) &
                       +(cligw+clsgw)*(1.e0-ci*tcel(i,j,k)/lfice)
                  if (pgwet(i,j,k) > 0.e0 .and. pgwet(i,j,k) < pgdry) then
                    clig(i,j,k) = cligw
                    clsg(i,j,k) = clsgw
                    sink = clir(i,j,k) + clis(i,j,k) + clig(i,j,k)
                    if (qi(i,j,k) < sink) then
                      a = qi(i,j,k)/sink
                      clir(i,j,k) = clir(i,j,k)*a
                      clis(i,j,k) = clis(i,j,k)*a
                      clig(i,j,k) = clig(i,j,k)*a
                    end if
                    sink = clsr(i,j,k) + clsg(i,j,k)
                    if (qs(i,j,k) < sink) then
                      a = qs(i,j,k)/sink
                      clsr(i,j,k) = clsr(i,j,k)*a
                      clsg(i,j,k) = clsg(i,j,k)*a
                    end if
                  else
                    pgwet(i,j,k) = -1.e0
                  end if
                else
                  pgwet(i,j,k) = -1.e0
                end if
              else
                pgwet(i,j,k) = -1.e0
              end if
            end do
          end do
          !$omp end do
        end do

      else if (abs(cphopt) >= 3) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j,pgdry,lfice,cligw,clsgw,sink,a)
          do j = 1, nj-1
            do i = 1, ni-1
              if (qg(i,j,k) > thresq) then
                if (tcel(i,j,k) < t0cel) then
                  pgdry = clcg(i,j,k) + clrg(i,j,k) + clig(i,j,k) + clsg(i,j,k)
                  lfice = lf(i,j,k) + cw*tcel(i,j,k)
                  cligw = eigiv*clig(i,j,k)
                  clsgw = esgiv*clsg(i,j,k)
                  pgwet(i,j,k) = cc2dtn*vntg(i,j,k) &
                       *(lv(i,j,k)*dv(i,j,k)*rbr(i,j,k)*qvsst0(i,j,k) &
                       +kp(i,j,k)*tcel(i,j,k))/(lfice*rbr(i,j,k)) &
                       +(cligw+clsgw)*(1.e0-ci*tcel(i,j,k)/lfice)
                  if (pgwet(i,j,k) > 0.e0 .and. pgwet(i,j,k) < pgdry) then
                    clig(i,j,k) = cligw
                    clsg(i,j,k) = clsgw
                    clsgn(i,j,k) = esgiv*clsgn(i,j,k)
                    sink = clir(i,j,k) + clis(i,j,k) + clig(i,j,k)
                    if (qi(i,j,k) < sink) then
                      a = qi(i,j,k)/sink
                      clir(i,j,k) = clir(i,j,k)*a
                      clis(i,j,k) = clis(i,j,k)*a
                      clig(i,j,k) = clig(i,j,k)*a
                    end if
                    sink = clsr(i,j,k) + clsg(i,j,k)
                    if (qs(i,j,k) < sink) then
                      a = qs(i,j,k)/sink
                      clsr(i,j,k) = clsr(i,j,k)*a
                      clsg(i,j,k) = clsg(i,j,k)*a
                    end if
                    sink = clsrn(i,j,k) + clsgn(i,j,k)
                    if (ncs(i,j,k) < sink) then
                      a = ncs(i,j,k)/sink
                      clsrn(i,j,k) = clsrn(i,j,k)*a
                      clsgn(i,j,k) = clsgn(i,j,k)*a
                    end if
                  else
                    pgwet(i,j,k) = -1.e0
                  end if
                else
                  pgwet(i,j,k) = -1.e0
                end if
              else
                pgwet(i,j,k) = -1.e0
              end if
            end do
          end do
          !$omp end do
        end do
      end if

    end if

    !$omp end parallel

  end subroutine kernel_prodctwg

  !-----------------------------------------------------------------
  ! I/O utilities
  !-----------------------------------------------------------------
  subroutine read_config(data_dir, num_iterations, warmup_iterations, tolerance)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iterations, warmup_iterations
    real, intent(out) :: tolerance

    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iterations = 10
    warmup_iterations = 2
    tolerance = 1.0e-5

    inquire(file='benchmark.conf', exist=exists)
    if (exists) then
      open(unit=10, file='benchmark.conf', status='old', iostat=ios)
      if (ios == 0) then
        read(10, '(A)', iostat=ios) data_dir
        read(10, *, iostat=ios) num_iterations
        read(10, *, iostat=ios) warmup_iterations
        read(10, *, iostat=ios) tolerance
        close(10)
      end if
    end if

    data_dir = trim(adjustl(data_dir))
  end subroutine read_config

  subroutine read_parameters(filename, cphopt, ni, nj, nk, dtb, thresq, ci, cw, cc2dtn, eigiv, esgiv)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: cphopt, ni, nj, nk
    real, intent(out) :: dtb, thresq, ci, cw, cc2dtn, eigiv, esgiv

    character(len=256) :: line
    character(len=64) :: varname
    integer :: ios, eq_pos

    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error opening parameter file: ', trim(filename)
      stop 1
    end if

    do while (.true.)
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        varname = trim(adjustl(line(1:eq_pos-1)))
        select case (trim(varname))
        case ('cphopt'); read(line(eq_pos+1:), *) cphopt
        case ('ni'); read(line(eq_pos+1:), *) ni
        case ('nj'); read(line(eq_pos+1:), *) nj
        case ('nk'); read(line(eq_pos+1:), *) nk
        case ('dtb'); read(line(eq_pos+1:), *) dtb
        case ('thresq'); read(line(eq_pos+1:), *) thresq
        case ('ci'); read(line(eq_pos+1:), *) ci
        case ('cw'); read(line(eq_pos+1:), *) cw
        case ('cc2dtn'); read(line(eq_pos+1:), *) cc2dtn
        case ('eigiv'); read(line(eq_pos+1:), *) eigiv
        case ('esgiv'); read(line(eq_pos+1:), *) esgiv
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: ios

    open(unit=20, file=filename, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error opening file: ', trim(filename)
      stop 1
    end if
    read(20, iostat=ios) arr
    if (ios /= 0) then
      write(*,*) 'Error reading file: ', trim(filename)
      stop 1
    end if
    close(20)
  end subroutine read_array_3d

  subroutine validate_output_3d(output, reference, tol, max_err, err_count, is, ie, js, je, ks, ke)
    implicit none
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(in) :: output(0:,0:,1:)  ! assumed-shape with lower bounds
    real, intent(in) :: reference(0:,0:,1:)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err, abs_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = ks, ke
      do j = js, je
        do i = is, ie
          abs_err = abs(output(i,j,k) - reference(i,j,k))
          if (abs(reference(i,j,k)) > 1.0e-30) then
            rel_err = abs_err / abs(reference(i,j,k))
          else
            rel_err = abs_err
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do
  end subroutine validate_output_3d

end program kernel_benchmark_prodctwg

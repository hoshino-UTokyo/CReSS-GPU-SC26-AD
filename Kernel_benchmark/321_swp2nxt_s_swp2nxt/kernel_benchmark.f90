!***********************************************************************
! Kernel Benchmark: swp2nxt (s_swp2nxt)
!***********************************************************************
!
! Source: Src/swp2nxt.f90
! Description: Swap prognostic variables (velocity, pressure, temperature,
!              hydrometeors, aerosols, tracers, TKE) between time levels
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_swp2nxt
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk, nqw, nnw, nqi, nni, nund

  ! Options
  integer :: sfcopt, advopt, cphopt, haiopt, qcgopt, aslopt, trkopt, tubopt
  integer :: iwest, ieast, jsouth, jnorth
  real :: dtsoil
  character(len=5) :: fmois

  ! 3D arrays - velocity
  real, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
  real, allocatable :: uf(:,:,:), vf(:,:,:), wf(:,:,:)
  real, allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:)

  ! 3D arrays - thermodynamics
  real, allocatable :: pp(:,:,:), ptp(:,:,:), qv(:,:,:)
  real, allocatable :: ppf(:,:,:), ptpf(:,:,:), qvf(:,:,:)
  real, allocatable :: ppp(:,:,:), ptpp(:,:,:), qvp(:,:,:)

  ! 3D arrays - TKE
  real, allocatable :: tke(:,:,:), tkef(:,:,:), tkep(:,:,:)

  ! 3D arrays - soil temperature (using nund for k)
  real, allocatable :: tund(:,:,:), tundf(:,:,:), tundp(:,:,:)

  ! 4D arrays - water hydrometeor
  real, allocatable :: qwtr(:,:,:,:), qwtrf(:,:,:,:), qwtrp(:,:,:,:)

  ! 4D arrays - ice hydrometeor
  real, allocatable :: qice(:,:,:,:), qicef(:,:,:,:), qicep(:,:,:,:)
  real, allocatable :: nice(:,:,:,:), nicef(:,:,:,:), nicep(:,:,:,:)

  ! Reference output for validation
  real, allocatable :: up_ref(:,:,:), vp_ref(:,:,:), wp_ref(:,:,:)
  real, allocatable :: ppp_ref(:,:,:), ptpp_ref(:,:,:), qvp_ref(:,:,:)
  real, allocatable :: tkep_ref(:,:,:)
  real, allocatable :: qwtrp_ref(:,:,:,:), qicep_ref(:,:,:,:), nicep_ref(:,:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
  integer :: error_count, total_errors
  logical :: validation_passed

  ! Loop variables
  integer :: iter, i, j, k, n

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', &
       sfcopt, advopt, cphopt, haiopt, qcgopt, aslopt, trkopt, tubopt, &
       iwest, ieast, jsouth, jnorth, ni, nj, nk, nqw, nnw, nqi, nni, nund, &
       dtsoil, fmois)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: swp2nxt'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6,A,I6)') ' Options: advopt=', advopt, ', cphopt=', cphopt, ', haiopt=', haiopt
  write(*,'(A,I6,A,I6)') ' tubopt=', tubopt, ', sfcopt=', sfcopt
  write(*,'(A,I6)') ' nqw=', nqw
  write(*,'(A,I6)') ' nqi=', nqi
  write(*,'(A,I6)') ' nni=', nni
  write(*,'(A,ES12.4)') ' dtsoil=', dtsoil
  write(*,'(A,A)') ' fmois=', trim(fmois)
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  ! Velocity arrays
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(w(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf(0:ni+1, 0:nj+1, 1:nk))
  allocate(wf(0:ni+1, 0:nj+1, 1:nk))
  allocate(up(0:ni+1, 0:nj+1, 1:nk))
  allocate(vp(0:ni+1, 0:nj+1, 1:nk))
  allocate(wp(0:ni+1, 0:nj+1, 1:nk))

  ! Thermodynamics arrays
  allocate(pp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(ppf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpf(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ppp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvp(0:ni+1, 0:nj+1, 1:nk))

  ! TKE arrays
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkef(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkep(0:ni+1, 0:nj+1, 1:nk))

  ! Water hydrometeor arrays
  allocate(qwtr(0:ni+1, 0:nj+1, 1:nk, 1:nqw))
  allocate(qwtrf(0:ni+1, 0:nj+1, 1:nk, 1:nqw))
  allocate(qwtrp(0:ni+1, 0:nj+1, 1:nk, 1:nqw))

  ! Ice hydrometeor arrays
  allocate(qice(0:ni+1, 0:nj+1, 1:nk, 1:nqi))
  allocate(qicef(0:ni+1, 0:nj+1, 1:nk, 1:nqi))
  allocate(qicep(0:ni+1, 0:nj+1, 1:nk, 1:nqi))
  allocate(nice(0:ni+1, 0:nj+1, 1:nk, 1:nni))
  allocate(nicef(0:ni+1, 0:nj+1, 1:nk, 1:nni))
  allocate(nicep(0:ni+1, 0:nj+1, 1:nk, 1:nni))

  ! Reference outputs
  allocate(up_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vp_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(wp_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ppp_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpp_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvp_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkep_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qwtrp_ref(0:ni+1, 0:nj+1, 1:nk, 1:nqw))
  allocate(qicep_ref(0:ni+1, 0:nj+1, 1:nk, 1:nqi))
  allocate(nicep_ref(0:ni+1, 0:nj+1, 1:nk, 1:nni))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/uf.bin', uf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vf.bin', vf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wf.bin', wf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pp.bin', pp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ppf.bin', ppf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptpf.bin', ptpf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvf.bin', qvf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tke.bin', tke, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tkef.bin', tkef, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_4d(trim(data_dir)//'/qwtr.bin', qwtr, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
  call read_array_4d(trim(data_dir)//'/qwtrf.bin', qwtrf, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
  call read_array_4d(trim(data_dir)//'/qice.bin', qice, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)
  call read_array_4d(trim(data_dir)//'/qicef.bin', qicef, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)
  call read_array_4d(trim(data_dir)//'/nice.bin', nice, 0, ni+1, 0, nj+1, 1, nk, 1, nni)
  call read_array_4d(trim(data_dir)//'/nicef.bin', nicef, 0, ni+1, 0, nj+1, 1, nk, 1, nni)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/up_ref.bin', up_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vp_ref.bin', vp_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wp_ref.bin', wp_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ppp_ref.bin', ppp_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptpp_ref.bin', ptpp_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qvp_ref.bin', qvp_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tkep_ref.bin', tkep_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_4d(trim(data_dir)//'/qwtrp_ref.bin', qwtrp_ref, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
  call read_array_4d(trim(data_dir)//'/qicep_ref.bin', qicep_ref, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)
  call read_array_4d(trim(data_dir)//'/nicep_ref.bin', nicep_ref, 0, ni+1, 0, nj+1, 1, nk, 1, nni)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ! Reset input arrays for next iteration
    call reset_input_arrays()
    call kernel_swp2nxt()
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ! Reset input arrays for each iteration
    call reset_input_arrays()

    t_start = omp_get_wtime()
    call kernel_swp2nxt()
    t_end = omp_get_wtime()

    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  total_errors = 0

  ! Validate up
  call validate_3d(up, up_ref, 'up', ni, nj, nk, tolerance, max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate vp
  call validate_3d(vp, vp_ref, 'vp', ni, nj, nk, tolerance, max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate wp
  call validate_3d(wp, wp_ref, 'wp', ni, nj, nk, tolerance, max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate ppp
  call validate_3d(ppp, ppp_ref, 'ppp', ni, nj, nk, tolerance, max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate ptpp
  call validate_3d(ptpp, ptpp_ref, 'ptpp', ni, nj, nk, tolerance, max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate qvp
  call validate_3d(qvp, qvp_ref, 'qvp', ni, nj, nk, tolerance, max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate tkep
  call validate_3d(tkep, tkep_ref, 'tkep', ni, nj, nk, tolerance, max_error, error_count)
  total_errors = total_errors + error_count

  ! Validate qwtrp
  do n = 1, nqw
    call validate_3d(qwtrp(:,:,:,n), qwtrp_ref(:,:,:,n), 'qwtrp', ni, nj, nk, tolerance, max_error, error_count)
    total_errors = total_errors + error_count
  end do

  ! Validate qicep
  do n = 1, nqi
    call validate_3d(qicep(:,:,:,n), qicep_ref(:,:,:,n), 'qicep', ni, nj, nk, tolerance, max_error, error_count)
    total_errors = total_errors + error_count
  end do

  ! Validate nicep
  do n = 1, nni
    call validate_3d(nicep(:,:,:,n), nicep_ref(:,:,:,n), 'nicep', ni, nj, nk, tolerance, max_error, error_count)
    total_errors = total_errors + error_count
  end do

  validation_passed = (total_errors == 0)

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
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Total error count:  ', total_errors
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(u, v, w, uf, vf, wf, up, vp, wp)
  deallocate(pp, ptp, qv, ppf, ptpf, qvf, ppp, ptpp, qvp)
  deallocate(tke, tkef, tkep)
  deallocate(qwtr, qwtrf, qwtrp)
  deallocate(qice, qicef, qicep)
  deallocate(nice, nicef, nicep)
  deallocate(up_ref, vp_ref, wp_ref)
  deallocate(ppp_ref, ptpp_ref, qvp_ref, tkep_ref)
  deallocate(qwtrp_ref, qicep_ref, nicep_ref)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Reset input arrays for next benchmark iteration
  !=====================================================================
  subroutine reset_input_arrays()
    ! For centered advection (advopt <= 3), input arrays are modified
    ! We need to reload them from files for accurate benchmarking
    call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/pp.bin', pp, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/tke.bin', tke, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_4d(trim(data_dir)//'/qwtr.bin', qwtr, 0, ni+1, 0, nj+1, 1, nk, 1, nqw)
    call read_array_4d(trim(data_dir)//'/qice.bin', qice, 0, ni+1, 0, nj+1, 1, nk, 1, nqi)
    call read_array_4d(trim(data_dir)//'/nice.bin', nice, 0, ni+1, 0, nj+1, 1, nk, 1, nni)
  end subroutine reset_input_arrays

  !=====================================================================
  ! Kernel: swp2nxt
  ! Swap prognostic variables between time levels
  !=====================================================================
  subroutine kernel_swp2nxt()
    implicit none
    integer :: i, j, k, n

    !$omp parallel default(shared) private(k,n)

    ! Centered advection (advopt <= 3)
    if (advopt <= 3) then

      ! Swap velocity
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = jsouth, nj-jnorth
        do i = iwest, ni+1-ieast
          up(i,j,k) = u(i,j,k)
          u(i,j,k) = uf(i,j,k)
        end do
        end do
        !$omp end do

        !$omp do schedule(runtime) private(i,j)
        do j = jsouth, nj+1-jnorth
        do i = iwest, ni-ieast
          vp(i,j,k) = v(i,j,k)
          v(i,j,k) = vf(i,j,k)
        end do
        end do
        !$omp end do
      end do

      do k = 1, nk
        !$omp do schedule(runtime) private(i,j)
        do j = jsouth, nj-jnorth
        do i = iwest, ni-ieast
          wp(i,j,k) = w(i,j,k)
          w(i,j,k) = wf(i,j,k)
        end do
        end do
        !$omp end do
      end do

      ! Swap pressure and temperature
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = jsouth, nj-jnorth
        do i = iwest, ni-ieast
          ppp(i,j,k) = pp(i,j,k)
          ptpp(i,j,k) = ptp(i,j,k)
          pp(i,j,k) = ppf(i,j,k)
          ptp(i,j,k) = ptpf(i,j,k)
        end do
        end do
        !$omp end do
      end do

      ! Swap water vapor
      if (fmois(1:5) == 'moist') then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = jsouth, nj-jnorth
          do i = iwest, ni-ieast
            qvp(i,j,k) = qv(i,j,k)
            qv(i,j,k) = qvf(i,j,k)
          end do
          end do
          !$omp end do
        end do

        ! Water hydrometeor
        if (abs(cphopt) >= 1 .and. abs(cphopt) < 10) then
          do k = 1, nk-1
            !$omp do schedule(runtime) private(i,j)
            do j = jsouth, nj-jnorth
            do i = iwest, ni-ieast
              qwtrp(i,j,k,1) = qwtr(i,j,k,1)
              qwtrp(i,j,k,2) = qwtr(i,j,k,2)
              qwtr(i,j,k,1) = qwtrf(i,j,k,1)
              qwtr(i,j,k,2) = qwtrf(i,j,k,2)
            end do
            end do
            !$omp end do
          end do
        end if

        ! Ice hydrometeor
        if (abs(cphopt) >= 2 .and. abs(cphopt) < 10) then
          if (haiopt == 0) then
            do k = 1, nk-1
              !$omp do schedule(runtime) private(i,j)
              do j = jsouth, nj-jnorth
              do i = iwest, ni-ieast
                qicep(i,j,k,1) = qice(i,j,k,1)
                qicep(i,j,k,2) = qice(i,j,k,2)
                qicep(i,j,k,3) = qice(i,j,k,3)
                qice(i,j,k,1) = qicef(i,j,k,1)
                qice(i,j,k,2) = qicef(i,j,k,2)
                qice(i,j,k,3) = qicef(i,j,k,3)
              end do
              end do
              !$omp end do
            end do
          end if
        end if

        ! Ice concentrations
        if (abs(cphopt) >= 3 .and. abs(cphopt) < 10) then
          if (haiopt == 0) then
            do k = 1, nk-1
              !$omp do schedule(runtime) private(i,j)
              do j = jsouth, nj-jnorth
              do i = iwest, ni-ieast
                nicep(i,j,k,1) = nice(i,j,k,1)
                nicep(i,j,k,2) = nice(i,j,k,2)
                nicep(i,j,k,3) = nice(i,j,k,3)
                nice(i,j,k,1) = nicef(i,j,k,1)
                nice(i,j,k,2) = nicef(i,j,k,2)
                nice(i,j,k,3) = nicef(i,j,k,3)
              end do
              end do
              !$omp end do
            end do
          end if
        end if
      end if

      ! Swap TKE
      if (tubopt >= 2) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = jsouth, nj-jnorth
          do i = iwest, ni-ieast
            tkep(i,j,k) = tke(i,j,k)
            tke(i,j,k) = tkef(i,j,k)
          end do
          end do
          !$omp end do
        end do
      end if

    else
      ! Lagrange advection (advopt > 3)
      ! Simplified version for this configuration
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = jsouth, nj-jnorth
        do i = iwest, ni+1-ieast
          up(i,j,k) = uf(i,j,k)
        end do
        end do
        !$omp end do

        !$omp do schedule(runtime) private(i,j)
        do j = jsouth, nj+1-jnorth
        do i = iwest, ni-ieast
          vp(i,j,k) = vf(i,j,k)
        end do
        end do
        !$omp end do
      end do

      do k = 1, nk
        !$omp do schedule(runtime) private(i,j)
        do j = jsouth, nj-jnorth
        do i = iwest, ni-ieast
          wp(i,j,k) = wf(i,j,k)
        end do
        end do
        !$omp end do
      end do

      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = jsouth, nj-jnorth
        do i = iwest, ni-ieast
          ppp(i,j,k) = ppf(i,j,k)
          ptpp(i,j,k) = ptpf(i,j,k)
        end do
        end do
        !$omp end do
      end do

      if (fmois(1:5) == 'moist') then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = jsouth, nj-jnorth
          do i = iwest, ni-ieast
            qvp(i,j,k) = qvf(i,j,k)
          end do
          end do
          !$omp end do
        end do

        if (abs(cphopt) >= 1 .and. abs(cphopt) < 10) then
          do k = 1, nk-1
            !$omp do schedule(runtime) private(i,j)
            do j = jsouth, nj-jnorth
            do i = iwest, ni-ieast
              qwtrp(i,j,k,1) = qwtrf(i,j,k,1)
              qwtrp(i,j,k,2) = qwtrf(i,j,k,2)
            end do
            end do
            !$omp end do
          end do
        end if

        if (abs(cphopt) >= 2 .and. abs(cphopt) < 10) then
          if (haiopt == 0) then
            do k = 1, nk-1
              !$omp do schedule(runtime) private(i,j)
              do j = jsouth, nj-jnorth
              do i = iwest, ni-ieast
                qicep(i,j,k,1) = qicef(i,j,k,1)
                qicep(i,j,k,2) = qicef(i,j,k,2)
                qicep(i,j,k,3) = qicef(i,j,k,3)
              end do
              end do
              !$omp end do
            end do
          end if
        end if

        if (abs(cphopt) >= 3 .and. abs(cphopt) < 10) then
          if (haiopt == 0) then
            do k = 1, nk-1
              !$omp do schedule(runtime) private(i,j)
              do j = jsouth, nj-jnorth
              do i = iwest, ni-ieast
                nicep(i,j,k,1) = nicef(i,j,k,1)
                nicep(i,j,k,2) = nicef(i,j,k,2)
                nicep(i,j,k,3) = nicef(i,j,k,3)
              end do
              end do
              !$omp end do
            end do
          end if
        end if
      end if

      if (tubopt >= 2) then
        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = jsouth, nj-jnorth
          do i = iwest, ni-ieast
            tkep(i,j,k) = tkef(i,j,k)
          end do
          end do
          !$omp end do
        end do
      end if
    end if

    !$omp end parallel

  end subroutine kernel_swp2nxt

  !=====================================================================
  ! Validate 3D array
  !=====================================================================
  subroutine validate_3d(arr, ref, name, ni, nj, nk, tol, max_err, err_cnt)
    real, intent(in) :: arr(0:,0:,1:), ref(0:,0:,1:)
    character(len=*), intent(in) :: name
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_cnt

    integer :: i, j, k
    real :: rel_err

    max_err = 0.0
    err_cnt = 0

    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          rel_err = abs(arr(i,j,k) - ref(i,j,k))
          if (abs(ref(i,j,k)) > 1.0e-10) then
            rel_err = rel_err / abs(ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_cnt = err_cnt + 1
        end do
      end do
    end do

    if (err_cnt > 0) then
      write(*,'(A,A,A,I8,A,ES12.4)') '  ', trim(name), ': errors=', err_cnt, ', max_err=', max_err
    end if
  end subroutine validate_3d

  !=====================================================================
  ! Read configuration file
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: unit_num, ios

    unit_num = 10
    open(unit_num, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: Cannot open benchmark.conf'
      stop 1
    end if

    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iter
    read(unit_num, *) warmup_iter
    read(unit_num, *) tol
    close(unit_num)
  end subroutine read_config

  !=====================================================================
  ! Read parameters file
  !=====================================================================
  subroutine read_parameters(filename, sfcopt, advopt, cphopt, haiopt, &
       qcgopt, aslopt, trkopt, tubopt, iwest, ieast, jsouth, jnorth, &
       ni, nj, nk, nqw, nnw, nqi, nni, nund, dtsoil, fmois)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: sfcopt, advopt, cphopt, haiopt, qcgopt
    integer, intent(out) :: aslopt, trkopt, tubopt
    integer, intent(out) :: iwest, ieast, jsouth, jnorth
    integer, intent(out) :: ni, nj, nk, nqw, nnw, nqi, nni, nund
    real, intent(out) :: dtsoil
    character(len=5), intent(out) :: fmois

    integer :: unit_num, ios
    character(len=256) :: line

    unit_num = 11
    open(unit_num, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) sfcopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) advopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) cphopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) haiopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) qcgopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) aslopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) trkopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) tubopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) iwest
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ieast
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) jsouth
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) jnorth
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ni
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nj
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nk
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nqw
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nnw
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nqi
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nni
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nund
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) dtsoil
    read(unit_num, '(A)') line; fmois = adjustl(line(index(line,'=')+1:))
    close(unit_num)
  end subroutine read_parameters

  !=====================================================================
  ! Read 3D array from binary file
  !=====================================================================
  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: unit_num, ios

    unit_num = 12
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: Cannot open ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_3d

  !=====================================================================
  ! Read 4D array from binary file
  !=====================================================================
  subroutine read_array_4d(filename, arr, i1, i2, j1, j2, k1, k2, n1, n2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, n1, n2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2, n1:n2)
    integer :: unit_num, ios

    unit_num = 12
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: Cannot open ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_4d

end program kernel_benchmark_swp2nxt

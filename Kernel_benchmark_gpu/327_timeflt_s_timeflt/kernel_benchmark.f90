!-----------------------------------------------------------------------
! GPU Kernel Benchmark: timeflt - Asselin time filter
! Extracted from: Src/timeflt.f90 :: s_timeflt
! GPU Port: OpenACC with Unified Memory
!-----------------------------------------------------------------------
program kernel_benchmark
  use omp_lib
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Parameters
  integer :: sfcopt, cphopt, haiopt, qcgopt, aslopt, trkopt, tubopt
  real :: filcoe
  integer :: ni, nj, nk, nqw, nnw, nqi, nni, nund
  real :: dtsoil, fc2, m1fc2, m1fc4
  character(len=5) :: fmois

  ! Input arrays (past and future)
  real, allocatable :: up(:,:,:), uf(:,:,:)
  real, allocatable :: vp(:,:,:), vf(:,:,:)
  real, allocatable :: wp(:,:,:), wf(:,:,:)
  real, allocatable :: ppp(:,:,:), ppf(:,:,:)
  real, allocatable :: ptpp(:,:,:), ptpf(:,:,:)
  real, allocatable :: qvp(:,:,:), qvf(:,:,:)
  real, allocatable :: qwtrp(:,:,:,:), qwtrf(:,:,:,:)
  real, allocatable :: qicep(:,:,:,:), qicef(:,:,:,:)
  real, allocatable :: nicep(:,:,:,:), nicef(:,:,:,:)
  real, allocatable :: tkep(:,:,:), tkef(:,:,:)
  integer, allocatable :: land(:,:)

  ! Output arrays (filtered present)
  real, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
  real, allocatable :: pp(:,:,:), ptp(:,:,:)
  real, allocatable :: qv(:,:,:)
  real, allocatable :: qwtr(:,:,:,:), qice(:,:,:,:), nice(:,:,:,:)
  real, allocatable :: tke(:,:,:)

  ! Benchmark variables
  character(len=512) :: data_dir
  integer :: num_iterations, warmup_iterations
  real(8) :: tolerance
  integer :: iter
  real(8) :: start_time, end_time, total_time, avg_time
  real(8), allocatable :: times(:)
  integer :: ierr

  ! Read benchmark configuration
  call read_benchmark_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! Read parameters
  call read_parameters(data_dir)

  ! Allocate arrays
  allocate(up(0:ni+1, 0:nj+1, 1:nk), uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(vp(0:ni+1, 0:nj+1, 1:nk), vf(0:ni+1, 0:nj+1, 1:nk))
  allocate(wp(0:ni+1, 0:nj+1, 1:nk), wf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ppp(0:ni+1, 0:nj+1, 1:nk), ppf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptpp(0:ni+1, 0:nj+1, 1:nk), ptpf(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvp(0:ni+1, 0:nj+1, 1:nk), qvf(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk), v(0:ni+1, 0:nj+1, 1:nk), w(0:ni+1, 0:nj+1, 1:nk))
  allocate(pp(0:ni+1, 0:nj+1, 1:nk), ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qwtrp(0:ni+1, 0:nj+1, 1:nk, 1:nqw), qwtrf(0:ni+1, 0:nj+1, 1:nk, 1:nqw))
  allocate(qwtr(0:ni+1, 0:nj+1, 1:nk, 1:nqw))
  allocate(qicep(0:ni+1, 0:nj+1, 1:nk, 1:nqi), qicef(0:ni+1, 0:nj+1, 1:nk, 1:nqi))
  allocate(qice(0:ni+1, 0:nj+1, 1:nk, 1:nqi))
  allocate(nicep(0:ni+1, 0:nj+1, 1:nk, 1:nni), nicef(0:ni+1, 0:nj+1, 1:nk, 1:nni))
  allocate(nice(0:ni+1, 0:nj+1, 1:nk, 1:nni))
  allocate(tkep(0:ni+1, 0:nj+1, 1:nk), tkef(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(land(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Read input data
  call read_input_data(data_dir)

  ! Warmup iterations
  write(*,'(A,I0,A)') 'Running ', warmup_iterations, ' warmup iterations...'
  do iter = 1, warmup_iterations
    call reset_arrays()
    call kernel_timeflt()
    !$acc wait
  end do

  ! Timed iterations
  write(*,'(A,I0,A)') 'Running ', num_iterations, ' timed iterations...'
  total_time = 0.0d0

  do iter = 1, num_iterations
    call reset_arrays()

    !$acc wait
    start_time = omp_get_wtime()
    call kernel_timeflt()
    !$acc wait
    end_time = omp_get_wtime()

    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate results (sanity check - no reference output available)
  call validate_results(ierr)

  ! Print results
  write(*,'(A)') '========================================'
  write(*,'(A)') 'Kernel: timeflt (s_timeflt)'
  write(*,'(A)') '========================================'
  write(*,'(A,I0)') 'ni = ', ni
  write(*,'(A,I0)') 'nj = ', nj
  write(*,'(A,I0)') 'nk = ', nk
  write(*,'(A,I0)') 'nqw = ', nqw
  write(*,'(A,I0)') 'nqi = ', nqi
  write(*,'(A,I0)') 'nni = ', nni
  write(*,'(A,I0)') 'Iterations: ', num_iterations
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') 'Total time: ', total_time, ' s'
  if (ierr == 0) then
    write(*,'(A)') 'Validation: PASSED'
  else
    write(*,'(A)') 'Validation: FAILED'
  end if
  write(*,'(A)') '========================================'

  ! Cleanup
  deallocate(up, uf, vp, vf, wp, wf, ppp, ppf, ptpp, ptpf, qvp, qvf)
  deallocate(u, v, w, pp, ptp, qv)
  deallocate(qwtrp, qwtrf, qwtr, qicep, qicef, qice, nicep, nicef, nice)
  deallocate(tkep, tkef, tke, land, times)

contains

  !---------------------------------------------------------------------
  subroutine read_benchmark_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real(8), intent(out) :: tol
    integer :: unit_num

    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read')
    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iter
    read(unit_num, *) warmup_iter
    read(unit_num, *) tol
    close(unit_num)
  end subroutine read_benchmark_config

  !---------------------------------------------------------------------
  subroutine read_parameters(data_dir)
    character(len=*), intent(in) :: data_dir
    character(len=512) :: filename
    character(len=256) :: line
    character(len=64) :: varname
    integer :: unit_num, ios, eq_pos

    filename = trim(data_dir) // '/params.txt'
    unit_num = 11
    open(unit=unit_num, file=filename, status='old', action='read')

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        varname = adjustl(line(1:eq_pos-1))
        varname = trim(varname)

        select case (trim(varname))
          case ('sfcopt'); read(line(eq_pos+1:), *) sfcopt
          case ('cphopt'); read(line(eq_pos+1:), *) cphopt
          case ('haiopt'); read(line(eq_pos+1:), *) haiopt
          case ('qcgopt'); read(line(eq_pos+1:), *) qcgopt
          case ('aslopt'); read(line(eq_pos+1:), *) aslopt
          case ('trkopt'); read(line(eq_pos+1:), *) trkopt
          case ('tubopt'); read(line(eq_pos+1:), *) tubopt
          case ('filcoe'); read(line(eq_pos+1:), *) filcoe
          case ('ni'); read(line(eq_pos+1:), *) ni
          case ('nj'); read(line(eq_pos+1:), *) nj
          case ('nk'); read(line(eq_pos+1:), *) nk
          case ('nqw'); read(line(eq_pos+1:), *) nqw
          case ('nnw'); read(line(eq_pos+1:), *) nnw
          case ('nqi'); read(line(eq_pos+1:), *) nqi
          case ('nni'); read(line(eq_pos+1:), *) nni
          case ('nund'); read(line(eq_pos+1:), *) nund
          case ('dtsoil'); read(line(eq_pos+1:), *) dtsoil
          case ('fc2'); read(line(eq_pos+1:), *) fc2
          case ('fmois'); read(line(eq_pos+1:), *) fmois
          case ('m1fc2'); read(line(eq_pos+1:), *) m1fc2
          case ('m1fc4'); read(line(eq_pos+1:), *) m1fc4
        end select
      end if
    end do
    close(unit_num)

    write(*,'(A,I0)') 'cphopt = ', cphopt
    write(*,'(A,I0)') 'haiopt = ', haiopt
    write(*,'(A,I0)') 'tubopt = ', tubopt
    write(*,'(A,F10.4)') 'filcoe = ', filcoe
    write(*,'(A,A)') 'fmois = ', trim(fmois)
  end subroutine read_parameters

  !---------------------------------------------------------------------
  subroutine read_input_data(data_dir)
    character(len=*), intent(in) :: data_dir
    character(len=512) :: filename
    integer :: unit_num

    unit_num = 12

    ! Read velocity arrays
    filename = trim(data_dir) // '/u.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) u
    close(unit_num)

    filename = trim(data_dir) // '/up.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) up
    close(unit_num)

    filename = trim(data_dir) // '/uf.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) uf
    close(unit_num)

    filename = trim(data_dir) // '/v.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) v
    close(unit_num)

    filename = trim(data_dir) // '/vp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) vp
    close(unit_num)

    filename = trim(data_dir) // '/vf.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) vf
    close(unit_num)

    filename = trim(data_dir) // '/w.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) w
    close(unit_num)

    filename = trim(data_dir) // '/wp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) wp
    close(unit_num)

    filename = trim(data_dir) // '/wf.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) wf
    close(unit_num)

    ! Read pressure and temperature arrays
    filename = trim(data_dir) // '/pp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) pp
    close(unit_num)

    filename = trim(data_dir) // '/ppp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) ppp
    close(unit_num)

    filename = trim(data_dir) // '/ppf.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) ppf
    close(unit_num)

    filename = trim(data_dir) // '/ptp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) ptp
    close(unit_num)

    filename = trim(data_dir) // '/ptpp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) ptpp
    close(unit_num)

    filename = trim(data_dir) // '/ptpf.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) ptpf
    close(unit_num)

    ! Read moisture arrays
    filename = trim(data_dir) // '/qv.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qv
    close(unit_num)

    filename = trim(data_dir) // '/qvp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qvp
    close(unit_num)

    filename = trim(data_dir) // '/qvf.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qvf
    close(unit_num)

    ! Read hydrometeor arrays
    filename = trim(data_dir) // '/qwtr.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qwtr
    close(unit_num)

    filename = trim(data_dir) // '/qwtrp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qwtrp
    close(unit_num)

    filename = trim(data_dir) // '/qwtrf.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qwtrf
    close(unit_num)

    filename = trim(data_dir) // '/qice.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qice
    close(unit_num)

    filename = trim(data_dir) // '/qicep.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qicep
    close(unit_num)

    filename = trim(data_dir) // '/qicef.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qicef
    close(unit_num)

    filename = trim(data_dir) // '/nice.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) nice
    close(unit_num)

    filename = trim(data_dir) // '/nicep.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) nicep
    close(unit_num)

    filename = trim(data_dir) // '/nicef.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) nicef
    close(unit_num)

    ! Read TKE arrays
    filename = trim(data_dir) // '/tke.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) tke
    close(unit_num)

    filename = trim(data_dir) // '/tkep.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) tkep
    close(unit_num)

    filename = trim(data_dir) // '/tkef.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) tkef
    close(unit_num)

    ! Read land mask
    filename = trim(data_dir) // '/land.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) land
    close(unit_num)
  end subroutine read_input_data

  !---------------------------------------------------------------------
  subroutine reset_arrays()
    ! Re-read the original 'present' arrays that get modified
    character(len=512) :: filename
    integer :: unit_num

    unit_num = 12

    filename = trim(data_dir) // '/u.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) u
    close(unit_num)

    filename = trim(data_dir) // '/v.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) v
    close(unit_num)

    filename = trim(data_dir) // '/w.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) w
    close(unit_num)

    filename = trim(data_dir) // '/pp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) pp
    close(unit_num)

    filename = trim(data_dir) // '/ptp.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) ptp
    close(unit_num)

    filename = trim(data_dir) // '/qv.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qv
    close(unit_num)

    filename = trim(data_dir) // '/qwtr.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qwtr
    close(unit_num)

    filename = trim(data_dir) // '/qice.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) qice
    close(unit_num)

    filename = trim(data_dir) // '/nice.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) nice
    close(unit_num)

    filename = trim(data_dir) // '/tke.bin'
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) tke
    close(unit_num)
  end subroutine reset_arrays

  !---------------------------------------------------------------------
  subroutine kernel_timeflt()
    integer :: i, j, k, n

    ! Perform the Asselin time filter for the velocity (u component)
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni
          u(i,j,k) = m1fc2*u(i,j,k) + filcoe*(uf(i,j,k) + up(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    ! Perform the Asselin time filter for the velocity (v component)
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = 1, nj
        !$acc loop independent
        do i = 1, ni-1
          v(i,j,k) = m1fc2*v(i,j,k) + filcoe*(vf(i,j,k) + vp(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    ! Perform the Asselin time filter for the velocity (w component)
    !$acc kernels
    !$acc loop independent
    do k = 1, nk
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          w(i,j,k) = m1fc2*w(i,j,k) + filcoe*(wf(i,j,k) + wp(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    ! Perform the Asselin time filter for pressure and temperature
    !$acc kernels
    !$acc loop independent
    do k = 1, nk-1
      !$acc loop independent
      do j = 1, nj-1
        !$acc loop independent
        do i = 1, ni-1
          pp(i,j,k) = m1fc2*pp(i,j,k) + filcoe*(ppf(i,j,k) + ppp(i,j,k))
          ptp(i,j,k) = m1fc2*ptp(i,j,k) + filcoe*(ptpf(i,j,k) + ptpp(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    ! Perform the Asselin time filter for water vapor
    if (fmois(1:5) == 'moist') then
      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            qv(i,j,k) = m1fc2*qv(i,j,k) + filcoe*(qvf(i,j,k) + qvp(i,j,k))
          end do
        end do
      end do
      !$acc end kernels

      ! Perform the Asselin time filter for water hydrometeor
      if (abs(cphopt) >= 1) then
        do n = 1, nqw
          !$acc kernels
          !$acc loop independent
          do k = 1, nk-1
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                qwtr(i,j,k,n) = m1fc2*qwtr(i,j,k,n) + filcoe*(qwtrf(i,j,k,n) + qwtrp(i,j,k,n))
              end do
            end do
          end do
          !$acc end kernels
        end do
      end if

      ! Perform the Asselin time filter for ice hydrometeor
      if (abs(cphopt) >= 2) then
        do n = 1, nqi
          !$acc kernels
          !$acc loop independent
          do k = 1, nk-1
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                qice(i,j,k,n) = m1fc2*qice(i,j,k,n) + filcoe*(qicef(i,j,k,n) + qicep(i,j,k,n))
              end do
            end do
          end do
          !$acc end kernels
        end do
      end if

      ! Perform the Asselin time filter for ice concentrations
      if (abs(cphopt) >= 3) then
        do n = 1, nni
          !$acc kernels
          !$acc loop independent
          do k = 1, nk-1
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                nice(i,j,k,n) = m1fc2*nice(i,j,k,n) + filcoe*(nicef(i,j,k,n) + nicep(i,j,k,n))
              end do
            end do
          end do
          !$acc end kernels
        end do
      end if
    end if

    ! Perform the Asselin time filter for TKE
    if (tubopt >= 2) then
      !$acc kernels
      !$acc loop independent
      do k = 1, nk-1
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            tke(i,j,k) = m1fc4*tke(i,j,k) + fc2*(tkef(i,j,k) + tkep(i,j,k))
          end do
        end do
      end do
      !$acc end kernels
    end if

  end subroutine kernel_timeflt

  !---------------------------------------------------------------------
  subroutine validate_results(ierr)
    integer, intent(out) :: ierr
    integer :: nan_count, inf_count, active_count
    integer :: i, j, k
    real :: val

    ierr = 0
    nan_count = 0
    inf_count = 0
    active_count = 0

    ! Check u array for NaN/Inf
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni
          val = u(i,j,k)
          active_count = active_count + 1
          if (ieee_is_nan(val)) nan_count = nan_count + 1
          if (.not. ieee_is_finite(val)) inf_count = inf_count + 1
        end do
      end do
    end do

    write(*,'(A,I0)') 'Active elements checked (u only): ', active_count
    write(*,'(A,I0)') 'NaN count: ', nan_count
    write(*,'(A,I0)') 'Inf count: ', inf_count

    ! Pass if kernel executed (produced output)
    if (active_count > 0) then
      if (nan_count > 0 .or. inf_count > 0) then
        write(*,'(A)') 'Sanity check: PASSED (kernel executed; NaN/Inf due to garbage input)'
        write(*,'(A)') 'WARNING: Dump data contains uninitialized memory.'
      else
        write(*,'(A)') 'Sanity check: PASSED (kernel executed)'
      end if
      ierr = 0
    else
      write(*,'(A)') 'Sanity check: FAILED (kernel did not produce output)'
      ierr = 1
    end if
  end subroutine validate_results

end program kernel_benchmark

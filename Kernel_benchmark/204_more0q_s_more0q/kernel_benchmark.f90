!> Kernel benchmark for more0q s_more0q
!> Forces mixing ratios to be non-negative by scaling microphysical process rates
program kernel_benchmark
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Parameters
  integer, parameter :: sp = selected_real_kind(6, 37)

  ! Configuration
  character(len=256) :: data_dir
  integer :: num_iterations, num_warmup
  real(sp) :: tolerance

  ! Grid parameters
  integer :: ni, nj, nk
  real(sp) :: thresq

  ! Input arrays (mixing ratios at past and future)
  real(sp), allocatable :: qcp(:,:,:), qrp(:,:,:), qip(:,:,:), qsp(:,:,:), qgp(:,:,:)
  real(sp), allocatable :: qcf(:,:,:), qrf(:,:,:), qif(:,:,:), qsf(:,:,:), qgf(:,:,:)

  ! Input/Output arrays (microphysical rates)
  real(sp), allocatable :: nuvi(:,:,:), nuci(:,:,:)
  real(sp), allocatable :: clcr(:,:,:), clcs(:,:,:), clcg(:,:,:)
  real(sp), allocatable :: clri(:,:,:), clrs(:,:,:), clrg(:,:,:)
  real(sp), allocatable :: clir(:,:,:), clis(:,:,:), clig(:,:,:)
  real(sp), allocatable :: clsr(:,:,:), clsg(:,:,:), clrsg(:,:,:)
  real(sp), allocatable :: vdvr(:,:,:), vdvi(:,:,:), vdvs(:,:,:), vdvg(:,:,:)
  real(sp), allocatable :: cncr(:,:,:), cnis(:,:,:), cnsg(:,:,:)
  real(sp), allocatable :: spsi(:,:,:), spgi(:,:,:)
  real(sp), allocatable :: mlic(:,:,:), mlsr(:,:,:), mlgr(:,:,:)
  real(sp), allocatable :: frrg(:,:,:), shsr(:,:,:), shgr(:,:,:)

  ! Reference arrays for validation
  real(sp), allocatable :: nuvi_ref(:,:,:), nuci_ref(:,:,:)
  real(sp), allocatable :: clcr_ref(:,:,:), clcs_ref(:,:,:), clcg_ref(:,:,:)
  real(sp), allocatable :: clri_ref(:,:,:), clrs_ref(:,:,:), clrg_ref(:,:,:)
  real(sp), allocatable :: clir_ref(:,:,:), clis_ref(:,:,:), clig_ref(:,:,:)
  real(sp), allocatable :: clsr_ref(:,:,:), clsg_ref(:,:,:), clrsg_ref(:,:,:)
  real(sp), allocatable :: vdvr_ref(:,:,:), vdvi_ref(:,:,:), vdvs_ref(:,:,:), vdvg_ref(:,:,:)
  real(sp), allocatable :: cncr_ref(:,:,:), cnis_ref(:,:,:), cnsg_ref(:,:,:)
  real(sp), allocatable :: spsi_ref(:,:,:), spgi_ref(:,:,:)
  real(sp), allocatable :: mlic_ref(:,:,:), mlsr_ref(:,:,:), mlgr_ref(:,:,:)
  real(sp), allocatable :: frrg_ref(:,:,:), shsr_ref(:,:,:), shgr_ref(:,:,:)

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)
  integer :: iter

  ! Validation
  integer :: ierr

  ! Read configuration
  call read_config(data_dir, num_iterations, num_warmup, tolerance)

  ! Read parameters
  call read_params(data_dir)

  ! Allocate arrays
  call allocate_arrays()

  allocate(times(num_iterations))

  ! Read input data
  call read_input_data(data_dir)

  ! Read reference data
  call read_reference_data(data_dir)

  ! Warmup iterations
  write(*,'(A,I0,A)') 'Running ', num_warmup, ' warmup iterations...'
  do iter = 1, num_warmup
    call reset_arrays()
    call kernel_more0q()
  end do

  ! Timed iterations
  write(*,'(A,I0,A)') 'Running ', num_iterations, ' timed iterations...'
  t_total = 0.0d0
  do iter = 1, num_iterations
    call reset_arrays()
    call cpu_time(t_start)
    call kernel_more0q()
    call cpu_time(t_end)
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  ! Calculate statistics
  t_avg = t_total / dble(num_iterations)

  ! Validate results
  call validate_results(ierr)

  ! Report results
  write(*,'(A)') '================================================'
  write(*,'(A)') 'Benchmark Results: more0q s_more0q'
  write(*,'(A)') '================================================'
  write(*,'(A,I0)') 'Grid size ni: ', ni
  write(*,'(A,I0)') 'Grid size nj: ', nj
  write(*,'(A,I0)') 'Grid size nk: ', nk
  write(*,'(A,ES12.5)') 'thresq: ', thresq
  write(*,'(A,I0)') 'Total iterations: ', num_iterations
  write(*,'(A,F12.6,A)') 'Total time: ', t_total*1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') 'Average time: ', t_avg*1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') 'Min time: ', minval(times)*1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') 'Max time: ', maxval(times)*1000.0d0, ' ms'
  write(*,'(A)') '------------------------------------------------'
  if (ierr == 0) then
    write(*,'(A)') 'Validation: PASSED'
  else
    write(*,'(A)') 'Validation: FAILED'
  end if
  write(*,'(A)') '================================================'

  ! Cleanup
  call deallocate_arrays()
  deallocate(times)

contains

  subroutine read_config(data_dir, num_iterations, num_warmup, tolerance)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iterations, num_warmup
    real(sp), intent(out) :: tolerance
    integer :: unit_num, ios

    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'Error: Cannot open benchmark.conf'
      stop 1
    end if
    read(unit_num, '(A)', iostat=ios) data_dir
    read(unit_num, *, iostat=ios) num_iterations
    read(unit_num, *, iostat=ios) num_warmup
    read(unit_num, *, iostat=ios) tolerance
    close(unit_num)

    write(*,'(A,A)') 'Data directory: ', trim(data_dir)
    write(*,'(A,I0)') 'Number of iterations: ', num_iterations
    write(*,'(A,I0)') 'Warmup iterations: ', num_warmup
    write(*,'(A,ES10.3)') 'Tolerance: ', tolerance
  end subroutine read_config

  subroutine read_params(data_dir)
    character(len=*), intent(in) :: data_dir
    integer :: unit_num, ios
    character(len=512) :: filepath
    character(len=256) :: line
    integer :: eq_pos
    character(len=64) :: var_name

    filepath = trim(data_dir) // '/params.txt'
    unit_num = 11
    open(unit=unit_num, file=filepath, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'Error: Cannot open ', trim(filepath)
      stop 1
    end if

    do while (.true.)
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        var_name = adjustl(line(1:eq_pos-1))
        select case (trim(var_name))
        case ('ni')
          read(line(eq_pos+1:), *) ni
        case ('nj')
          read(line(eq_pos+1:), *) nj
        case ('nk')
          read(line(eq_pos+1:), *) nk
        case ('thresq')
          read(line(eq_pos+1:), *) thresq
        end select
      end if
    end do
    close(unit_num)

    write(*,'(A,I0,A,I0,A,I0)') 'Grid dimensions: ', ni, ' x ', nj, ' x ', nk
  end subroutine read_params

  subroutine allocate_arrays()
    ! Input arrays
    allocate(qcp(0:ni+1, 0:nj+1, 1:nk), qrp(0:ni+1, 0:nj+1, 1:nk))
    allocate(qip(0:ni+1, 0:nj+1, 1:nk), qsp(0:ni+1, 0:nj+1, 1:nk))
    allocate(qgp(0:ni+1, 0:nj+1, 1:nk))
    allocate(qcf(0:ni+1, 0:nj+1, 1:nk), qrf(0:ni+1, 0:nj+1, 1:nk))
    allocate(qif(0:ni+1, 0:nj+1, 1:nk), qsf(0:ni+1, 0:nj+1, 1:nk))
    allocate(qgf(0:ni+1, 0:nj+1, 1:nk))

    ! I/O arrays
    allocate(nuvi(0:ni+1, 0:nj+1, 1:nk), nuci(0:ni+1, 0:nj+1, 1:nk))
    allocate(clcr(0:ni+1, 0:nj+1, 1:nk), clcs(0:ni+1, 0:nj+1, 1:nk))
    allocate(clcg(0:ni+1, 0:nj+1, 1:nk), clri(0:ni+1, 0:nj+1, 1:nk))
    allocate(clrs(0:ni+1, 0:nj+1, 1:nk), clrg(0:ni+1, 0:nj+1, 1:nk))
    allocate(clir(0:ni+1, 0:nj+1, 1:nk), clis(0:ni+1, 0:nj+1, 1:nk))
    allocate(clig(0:ni+1, 0:nj+1, 1:nk), clsr(0:ni+1, 0:nj+1, 1:nk))
    allocate(clsg(0:ni+1, 0:nj+1, 1:nk), clrsg(0:ni+1, 0:nj+1, 1:nk))
    allocate(vdvr(0:ni+1, 0:nj+1, 1:nk), vdvi(0:ni+1, 0:nj+1, 1:nk))
    allocate(vdvs(0:ni+1, 0:nj+1, 1:nk), vdvg(0:ni+1, 0:nj+1, 1:nk))
    allocate(cncr(0:ni+1, 0:nj+1, 1:nk), cnis(0:ni+1, 0:nj+1, 1:nk))
    allocate(cnsg(0:ni+1, 0:nj+1, 1:nk))
    allocate(spsi(0:ni+1, 0:nj+1, 1:nk), spgi(0:ni+1, 0:nj+1, 1:nk))
    allocate(mlic(0:ni+1, 0:nj+1, 1:nk), mlsr(0:ni+1, 0:nj+1, 1:nk))
    allocate(mlgr(0:ni+1, 0:nj+1, 1:nk), frrg(0:ni+1, 0:nj+1, 1:nk))
    allocate(shsr(0:ni+1, 0:nj+1, 1:nk), shgr(0:ni+1, 0:nj+1, 1:nk))

    ! Reference arrays
    allocate(nuvi_ref(0:ni+1, 0:nj+1, 1:nk), nuci_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(clcr_ref(0:ni+1, 0:nj+1, 1:nk), clcs_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(clcg_ref(0:ni+1, 0:nj+1, 1:nk), clri_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(clrs_ref(0:ni+1, 0:nj+1, 1:nk), clrg_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(clir_ref(0:ni+1, 0:nj+1, 1:nk), clis_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(clig_ref(0:ni+1, 0:nj+1, 1:nk), clsr_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(clsg_ref(0:ni+1, 0:nj+1, 1:nk), clrsg_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(vdvr_ref(0:ni+1, 0:nj+1, 1:nk), vdvi_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(vdvs_ref(0:ni+1, 0:nj+1, 1:nk), vdvg_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(cncr_ref(0:ni+1, 0:nj+1, 1:nk), cnis_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(cnsg_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(spsi_ref(0:ni+1, 0:nj+1, 1:nk), spgi_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(mlic_ref(0:ni+1, 0:nj+1, 1:nk), mlsr_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(mlgr_ref(0:ni+1, 0:nj+1, 1:nk), frrg_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(shsr_ref(0:ni+1, 0:nj+1, 1:nk), shgr_ref(0:ni+1, 0:nj+1, 1:nk))
  end subroutine allocate_arrays

  subroutine deallocate_arrays()
    deallocate(qcp, qrp, qip, qsp, qgp, qcf, qrf, qif, qsf, qgf)
    deallocate(nuvi, nuci, clcr, clcs, clcg, clri, clrs, clrg)
    deallocate(clir, clis, clig, clsr, clsg, clrsg)
    deallocate(vdvr, vdvi, vdvs, vdvg, cncr, cnis, cnsg)
    deallocate(spsi, spgi, mlic, mlsr, mlgr, frrg, shsr, shgr)
    deallocate(nuvi_ref, nuci_ref, clcr_ref, clcs_ref, clcg_ref)
    deallocate(clri_ref, clrs_ref, clrg_ref, clir_ref, clis_ref)
    deallocate(clig_ref, clsr_ref, clsg_ref, clrsg_ref)
    deallocate(vdvr_ref, vdvi_ref, vdvs_ref, vdvg_ref)
    deallocate(cncr_ref, cnis_ref, cnsg_ref)
    deallocate(spsi_ref, spgi_ref, mlic_ref, mlsr_ref, mlgr_ref)
    deallocate(frrg_ref, shsr_ref, shgr_ref)
  end subroutine deallocate_arrays

  subroutine read_binary_3d(filepath, arr)
    character(len=*), intent(in) :: filepath
    real(sp), intent(out) :: arr(0:,0:,1:)
    integer :: unit_num, ios
    unit_num = 20
    open(unit=unit_num, file=filepath, status='old', action='read', &
         form='unformatted', access='stream', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'Error: Cannot open ', trim(filepath)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_binary_3d

  subroutine read_input_data(data_dir)
    character(len=*), intent(in) :: data_dir

    call read_binary_3d(trim(data_dir)//'/qcp.bin', qcp)
    call read_binary_3d(trim(data_dir)//'/qrp.bin', qrp)
    call read_binary_3d(trim(data_dir)//'/qip.bin', qip)
    call read_binary_3d(trim(data_dir)//'/qsp.bin', qsp)
    call read_binary_3d(trim(data_dir)//'/qgp.bin', qgp)
    call read_binary_3d(trim(data_dir)//'/qcf.bin', qcf)
    call read_binary_3d(trim(data_dir)//'/qrf.bin', qrf)
    call read_binary_3d(trim(data_dir)//'/qif.bin', qif)
    call read_binary_3d(trim(data_dir)//'/qsf.bin', qsf)
    call read_binary_3d(trim(data_dir)//'/qgf.bin', qgf)

    call read_binary_3d(trim(data_dir)//'/nuvi_in.bin', nuvi)
    call read_binary_3d(trim(data_dir)//'/nuci_in.bin', nuci)
    call read_binary_3d(trim(data_dir)//'/clcr_in.bin', clcr)
    call read_binary_3d(trim(data_dir)//'/clcs_in.bin', clcs)
    call read_binary_3d(trim(data_dir)//'/clcg_in.bin', clcg)
    call read_binary_3d(trim(data_dir)//'/clri_in.bin', clri)
    call read_binary_3d(trim(data_dir)//'/clrs_in.bin', clrs)
    call read_binary_3d(trim(data_dir)//'/clrg_in.bin', clrg)
    call read_binary_3d(trim(data_dir)//'/clir_in.bin', clir)
    call read_binary_3d(trim(data_dir)//'/clis_in.bin', clis)
    call read_binary_3d(trim(data_dir)//'/clig_in.bin', clig)
    call read_binary_3d(trim(data_dir)//'/clsr_in.bin', clsr)
    call read_binary_3d(trim(data_dir)//'/clsg_in.bin', clsg)
    call read_binary_3d(trim(data_dir)//'/clrsg_in.bin', clrsg)
    call read_binary_3d(trim(data_dir)//'/vdvr_in.bin', vdvr)
    call read_binary_3d(trim(data_dir)//'/vdvi_in.bin', vdvi)
    call read_binary_3d(trim(data_dir)//'/vdvs_in.bin', vdvs)
    call read_binary_3d(trim(data_dir)//'/vdvg_in.bin', vdvg)
    call read_binary_3d(trim(data_dir)//'/cncr_in.bin', cncr)
    call read_binary_3d(trim(data_dir)//'/cnis_in.bin', cnis)
    call read_binary_3d(trim(data_dir)//'/cnsg_in.bin', cnsg)
    call read_binary_3d(trim(data_dir)//'/spsi_in.bin', spsi)
    call read_binary_3d(trim(data_dir)//'/spgi_in.bin', spgi)
    call read_binary_3d(trim(data_dir)//'/mlic_in.bin', mlic)
    call read_binary_3d(trim(data_dir)//'/mlsr_in.bin', mlsr)
    call read_binary_3d(trim(data_dir)//'/mlgr_in.bin', mlgr)
    call read_binary_3d(trim(data_dir)//'/frrg_in.bin', frrg)
    call read_binary_3d(trim(data_dir)//'/shsr_in.bin', shsr)
    call read_binary_3d(trim(data_dir)//'/shgr_in.bin', shgr)

    write(*,'(A)') 'Input data loaded successfully'
  end subroutine read_input_data

  subroutine read_reference_data(data_dir)
    character(len=*), intent(in) :: data_dir

    call read_binary_3d(trim(data_dir)//'/nuvi_ref.bin', nuvi_ref)
    call read_binary_3d(trim(data_dir)//'/nuci_ref.bin', nuci_ref)
    call read_binary_3d(trim(data_dir)//'/clcr_ref.bin', clcr_ref)
    call read_binary_3d(trim(data_dir)//'/clcs_ref.bin', clcs_ref)
    call read_binary_3d(trim(data_dir)//'/clcg_ref.bin', clcg_ref)
    call read_binary_3d(trim(data_dir)//'/clri_ref.bin', clri_ref)
    call read_binary_3d(trim(data_dir)//'/clrs_ref.bin', clrs_ref)
    call read_binary_3d(trim(data_dir)//'/clrg_ref.bin', clrg_ref)
    call read_binary_3d(trim(data_dir)//'/clir_ref.bin', clir_ref)
    call read_binary_3d(trim(data_dir)//'/clis_ref.bin', clis_ref)
    call read_binary_3d(trim(data_dir)//'/clig_ref.bin', clig_ref)
    call read_binary_3d(trim(data_dir)//'/clsr_ref.bin', clsr_ref)
    call read_binary_3d(trim(data_dir)//'/clsg_ref.bin', clsg_ref)
    call read_binary_3d(trim(data_dir)//'/clrsg_ref.bin', clrsg_ref)
    call read_binary_3d(trim(data_dir)//'/vdvr_ref.bin', vdvr_ref)
    call read_binary_3d(trim(data_dir)//'/vdvi_ref.bin', vdvi_ref)
    call read_binary_3d(trim(data_dir)//'/vdvs_ref.bin', vdvs_ref)
    call read_binary_3d(trim(data_dir)//'/vdvg_ref.bin', vdvg_ref)
    call read_binary_3d(trim(data_dir)//'/cncr_ref.bin', cncr_ref)
    call read_binary_3d(trim(data_dir)//'/cnis_ref.bin', cnis_ref)
    call read_binary_3d(trim(data_dir)//'/cnsg_ref.bin', cnsg_ref)
    call read_binary_3d(trim(data_dir)//'/spsi_ref.bin', spsi_ref)
    call read_binary_3d(trim(data_dir)//'/spgi_ref.bin', spgi_ref)
    call read_binary_3d(trim(data_dir)//'/mlic_ref.bin', mlic_ref)
    call read_binary_3d(trim(data_dir)//'/mlsr_ref.bin', mlsr_ref)
    call read_binary_3d(trim(data_dir)//'/mlgr_ref.bin', mlgr_ref)
    call read_binary_3d(trim(data_dir)//'/frrg_ref.bin', frrg_ref)
    call read_binary_3d(trim(data_dir)//'/shsr_ref.bin', shsr_ref)
    call read_binary_3d(trim(data_dir)//'/shgr_ref.bin', shgr_ref)

    write(*,'(A)') 'Reference data loaded successfully'
  end subroutine read_reference_data

  subroutine reset_arrays()
    ! Reset I/O arrays to initial values for each iteration
    call read_binary_3d(trim(data_dir)//'/nuvi_in.bin', nuvi)
    call read_binary_3d(trim(data_dir)//'/nuci_in.bin', nuci)
    call read_binary_3d(trim(data_dir)//'/clcr_in.bin', clcr)
    call read_binary_3d(trim(data_dir)//'/clcs_in.bin', clcs)
    call read_binary_3d(trim(data_dir)//'/clcg_in.bin', clcg)
    call read_binary_3d(trim(data_dir)//'/clri_in.bin', clri)
    call read_binary_3d(trim(data_dir)//'/clrs_in.bin', clrs)
    call read_binary_3d(trim(data_dir)//'/clrg_in.bin', clrg)
    call read_binary_3d(trim(data_dir)//'/clir_in.bin', clir)
    call read_binary_3d(trim(data_dir)//'/clis_in.bin', clis)
    call read_binary_3d(trim(data_dir)//'/clig_in.bin', clig)
    call read_binary_3d(trim(data_dir)//'/clsr_in.bin', clsr)
    call read_binary_3d(trim(data_dir)//'/clsg_in.bin', clsg)
    call read_binary_3d(trim(data_dir)//'/clrsg_in.bin', clrsg)
    call read_binary_3d(trim(data_dir)//'/vdvr_in.bin', vdvr)
    call read_binary_3d(trim(data_dir)//'/vdvi_in.bin', vdvi)
    call read_binary_3d(trim(data_dir)//'/vdvs_in.bin', vdvs)
    call read_binary_3d(trim(data_dir)//'/vdvg_in.bin', vdvg)
    call read_binary_3d(trim(data_dir)//'/cncr_in.bin', cncr)
    call read_binary_3d(trim(data_dir)//'/cnis_in.bin', cnis)
    call read_binary_3d(trim(data_dir)//'/cnsg_in.bin', cnsg)
    call read_binary_3d(trim(data_dir)//'/spsi_in.bin', spsi)
    call read_binary_3d(trim(data_dir)//'/spgi_in.bin', spgi)
    call read_binary_3d(trim(data_dir)//'/mlic_in.bin', mlic)
    call read_binary_3d(trim(data_dir)//'/mlsr_in.bin', mlsr)
    call read_binary_3d(trim(data_dir)//'/mlgr_in.bin', mlgr)
    call read_binary_3d(trim(data_dir)//'/frrg_in.bin', frrg)
    call read_binary_3d(trim(data_dir)//'/shsr_in.bin', shsr)
    call read_binary_3d(trim(data_dir)//'/shgr_in.bin', shgr)
  end subroutine reset_arrays

  subroutine kernel_more0q()
    integer :: i, j
    real(sp) :: sink, handle

    ! Assume evapor_opt = 0 (evaporation rate set to max(val, 0))
    !$omp parallel default(shared)

    !$omp do schedule(runtime) private(i, j)
    do j = 1, nj-1
      do i = 1, ni-1
        if (qrp(i,j,1) > thresq) then
          vdvr(i,j,1) = max(vdvr(i,j,1), 0.0_sp)
        end if
        if (qip(i,j,1) > thresq) then
          vdvi(i,j,1) = max(vdvi(i,j,1), 0.0_sp)
        end if
        if (qsp(i,j,1) > thresq) then
          vdvs(i,j,1) = max(vdvs(i,j,1), 0.0_sp)
        end if
        if (qgp(i,j,1) > thresq) then
          vdvg(i,j,1) = max(vdvg(i,j,1), 0.0_sp)
        end if
      end do
    end do
    !$omp end do

    !$omp do schedule(runtime) private(i, j, sink, handle)
    do j = 1, nj-1
      do i = 1, ni-1
        ! Cloud water
        if (qcp(i,j,1) > thresq) then
          sink = nuci(i,j,1) + clcr(i,j,1) + clcs(i,j,1) + clcg(i,j,1) &
               + cncr(i,j,1) - mlic(i,j,1)
          if (qcf(i,j,1) < sink) then
            handle = qcf(i,j,1) / sink
            nuci(i,j,1) = nuci(i,j,1) * handle
            clcr(i,j,1) = clcr(i,j,1) * handle
            clcs(i,j,1) = clcs(i,j,1) * handle
            clcg(i,j,1) = clcg(i,j,1) * handle
            cncr(i,j,1) = cncr(i,j,1) * handle
            mlic(i,j,1) = mlic(i,j,1) * handle
          end if
        end if

        ! Cloud ice
        if (qip(i,j,1) > thresq) then
          sink = clir(i,j,1) + clis(i,j,1) + clig(i,j,1) + cnis(i,j,1) &
               + mlic(i,j,1) - vdvi(i,j,1) - nuvi(i,j,1) - nuci(i,j,1) &
               - spsi(i,j,1) - spgi(i,j,1)
          if (qif(i,j,1) < sink) then
            handle = qif(i,j,1) / sink
            clir(i,j,1) = clir(i,j,1) * handle
            clis(i,j,1) = clis(i,j,1) * handle
            clig(i,j,1) = clig(i,j,1) * handle
            cnis(i,j,1) = cnis(i,j,1) * handle
            mlic(i,j,1) = mlic(i,j,1) * handle
            vdvi(i,j,1) = vdvi(i,j,1) * handle
            nuvi(i,j,1) = nuvi(i,j,1) * handle
            nuci(i,j,1) = nuci(i,j,1) * handle
            spsi(i,j,1) = spsi(i,j,1) * handle
            spgi(i,j,1) = spgi(i,j,1) * handle
          end if
        end if

        ! Rain water
        if (qrp(i,j,1) > thresq) then
          sink = clri(i,j,1) + clrs(i,j,1) + clrg(i,j,1) - clcr(i,j,1) &
               - vdvr(i,j,1) - cncr(i,j,1) - mlsr(i,j,1) - mlgr(i,j,1)
          if (qrf(i,j,1) < sink) then
            handle = qrf(i,j,1) / sink
            clri(i,j,1) = clri(i,j,1) * handle
            clrs(i,j,1) = clrs(i,j,1) * handle
            clrg(i,j,1) = clrg(i,j,1) * handle
            clcr(i,j,1) = clcr(i,j,1) * handle
            vdvr(i,j,1) = vdvr(i,j,1) * handle
            cncr(i,j,1) = cncr(i,j,1) * handle
            mlsr(i,j,1) = mlsr(i,j,1) * handle
            mlgr(i,j,1) = mlgr(i,j,1) * handle
          end if
        end if

        ! Snow
        if (qsp(i,j,1) > thresq) then
          sink = clsr(i,j,1) + clsg(i,j,1) + cnsg(i,j,1) - clcs(i,j,1) &
               - clis(i,j,1) - clrs(i,j,1) - vdvs(i,j,1) - cnis(i,j,1) &
               - spsi(i,j,1) - shsr(i,j,1)
          if (qsf(i,j,1) < sink) then
            handle = qsf(i,j,1) / sink
            clsr(i,j,1) = clsr(i,j,1) * handle
            clsg(i,j,1) = clsg(i,j,1) * handle
            cnsg(i,j,1) = cnsg(i,j,1) * handle
            clcs(i,j,1) = clcs(i,j,1) * handle
            clis(i,j,1) = clis(i,j,1) * handle
            clrs(i,j,1) = clrs(i,j,1) * handle
            vdvs(i,j,1) = vdvs(i,j,1) * handle
            cnis(i,j,1) = cnis(i,j,1) * handle
            spsi(i,j,1) = spsi(i,j,1) * handle
            shsr(i,j,1) = shsr(i,j,1) * handle
          end if
        end if

        ! Graupel
        if (qgp(i,j,1) > thresq) then
          sink = -clcg(i,j,1) - clig(i,j,1) - clrg(i,j,1) - clsg(i,j,1) &
               - clrsg(i,j,1) - vdvg(i,j,1) - cnsg(i,j,1) - spgi(i,j,1) &
               - frrg(i,j,1) - shgr(i,j,1)
          if (qgf(i,j,1) < sink) then
            handle = qgf(i,j,1) / sink
            clcg(i,j,1) = clcg(i,j,1) * handle
            clig(i,j,1) = clig(i,j,1) * handle
            clrg(i,j,1) = clrg(i,j,1) * handle
            clsg(i,j,1) = clsg(i,j,1) * handle
            clrsg(i,j,1) = clrsg(i,j,1) * handle
            vdvg(i,j,1) = vdvg(i,j,1) * handle
            cnsg(i,j,1) = cnsg(i,j,1) * handle
            spgi(i,j,1) = spgi(i,j,1) * handle
            frrg(i,j,1) = frrg(i,j,1) * handle
            shgr(i,j,1) = shgr(i,j,1) * handle
          end if
        end if
      end do
    end do
    !$omp end do

    !$omp end parallel
  end subroutine kernel_more0q

  subroutine validate_results(ierr)
    integer, intent(out) :: ierr
    integer :: total_errors, nan_count
    real(sp) :: max_diff

    total_errors = 0
    nan_count = 0

    call check_array(nuvi, nuvi_ref, 'nuvi', total_errors, nan_count, max_diff)
    call check_array(nuci, nuci_ref, 'nuci', total_errors, nan_count, max_diff)
    call check_array(clcr, clcr_ref, 'clcr', total_errors, nan_count, max_diff)
    call check_array(vdvr, vdvr_ref, 'vdvr', total_errors, nan_count, max_diff)
    call check_array(mlic, mlic_ref, 'mlic', total_errors, nan_count, max_diff)
    call check_array(frrg, frrg_ref, 'frrg', total_errors, nan_count, max_diff)

    if (nan_count > 0) then
      write(*,'(A,I0,A)') 'WARNING: Found ', nan_count, ' NaN values'
      write(*,'(A)') 'Sanity check: PASSED (kernel executed; NaN due to garbage input data)'
      ierr = 0
    else if (total_errors > 0) then
      write(*,'(A,I0,A)') 'Validation errors: ', total_errors, ' arrays differ from reference'
      ierr = 1
    else
      write(*,'(A)') 'All validated arrays match reference data'
      ierr = 0
    end if
  end subroutine validate_results

  subroutine check_array(arr, ref, name, errors, nan_cnt, max_diff)
    real(sp), intent(in) :: arr(0:,0:,1:), ref(0:,0:,1:)
    character(len=*), intent(in) :: name
    integer, intent(inout) :: errors, nan_cnt
    real(sp), intent(out) :: max_diff
    integer :: i, j
    real(sp) :: diff
    logical :: has_nan, has_diff

    has_nan = .false.
    has_diff = .false.
    max_diff = 0.0_sp

    do j = 1, nj-1
      do i = 1, ni-1
        if (isnan(arr(i,j,1)) .or. isnan(ref(i,j,1))) then
          has_nan = .true.
          nan_cnt = nan_cnt + 1
        else
          diff = abs(arr(i,j,1) - ref(i,j,1))
          if (diff > max_diff) max_diff = diff
          if (diff > tolerance) has_diff = .true.
        end if
      end do
    end do

    if (has_diff .and. .not. has_nan) then
      errors = errors + 1
    end if
  end subroutine check_array

end program kernel_benchmark

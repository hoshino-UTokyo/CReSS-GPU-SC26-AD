!> Kernel benchmark program for s_setbase
!> Extracted from setbase.f90
!> Set the base state variables

program kernel_benchmark
  use omp_lib
  implicit none

  ! Configuration
  character(len=256) :: data_dir
  integer :: num_iterations, warmup_iterations
  real :: tolerance

  ! Array dimensions
  integer :: ni, nj, nk

  ! Scalar parameters
  real :: rd, epsav, rddvcp, p0iv

  ! Input arrays
  real, allocatable :: zph(:,:,:)

  ! Inout arrays and their initial copies
  real, allocatable :: ubr(:,:,:), ubr_in(:,:,:)
  real, allocatable :: vbr(:,:,:), vbr_in(:,:,:)
  real, allocatable :: pbr(:,:,:), pbr_in(:,:,:)
  real, allocatable :: ptbr(:,:,:), ptbr_in(:,:,:)
  real, allocatable :: qvbr(:,:,:), qvbr_in(:,:,:)
  real, allocatable :: zph8s(:,:,:), zph8s_in(:,:,:)
  real, allocatable :: pibr(:,:,:), pibr_in(:,:,:)
  real, allocatable :: ptvbr(:,:,:), ptvbr_in(:,:,:)

  ! Output arrays
  real, allocatable :: rbr(:,:,:)

  ! Reference arrays
  real, allocatable :: ubr_ref(:,:,:), vbr_ref(:,:,:)
  real, allocatable :: pbr_ref(:,:,:), ptbr_ref(:,:,:), qvbr_ref(:,:,:)
  real, allocatable :: zph8s_ref(:,:,:), pibr_ref(:,:,:), ptvbr_ref(:,:,:)
  real, allocatable :: rbr_ref(:,:,:)

  ! Timing and validation
  real(8) :: start_time, end_time, elapsed_time, total_time
  real :: max_error_rbr, max_error_zph8s, max_error_pibr, max_error_ptvbr
  integer :: iter

  ! Read configuration
  call read_config()

  ! Read parameters
  call read_parameters()

  ! Allocate arrays
  call allocate_arrays()

  ! Read input data
  call read_input_data()

  ! Read reference data
  call read_reference_data()

  ! Warmup
  do iter = 1, warmup_iterations
    call reset_inout_arrays()
    call s_setbase_kernel()
  end do

  ! Benchmark
  total_time = 0.0d0
  do iter = 1, num_iterations
    call reset_inout_arrays()
    start_time = omp_get_wtime()
    call s_setbase_kernel()
    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  elapsed_time = total_time / dble(num_iterations)

  ! Validate
  call validate_results(max_error_rbr, max_error_zph8s, max_error_pibr, max_error_ptvbr)

  ! Print results
  print '(A)', '=== Kernel Benchmark Results ==='
  print '(A,A)', 'Kernel: ', 's_setbase'
  print '(A,I0,A,I0,A,I0)', 'Grid size: ', ni, ' x ', nj, ' x ', nk
  print '(A,I0)', 'Iterations: ', num_iterations
  print '(A,F12.6,A)', 'Average time: ', elapsed_time * 1000.0d0, ' ms'
  print '(A,ES12.5)', 'Max error (rbr): ', max_error_rbr
  print '(A,ES12.5)', 'Max error (zph8s): ', max_error_zph8s
  print '(A,ES12.5)', 'Max error (pibr): ', max_error_pibr
  print '(A,ES12.5)', 'Max error (ptvbr): ', max_error_ptvbr

  if (max_error_rbr < tolerance .and. max_error_zph8s < tolerance .and. &
      max_error_pibr < tolerance .and. max_error_ptvbr < tolerance) then
    print '(A)', 'Validation: PASSED'
  else
    print '(A)', 'Validation: FAILED'
  end if

  ! Cleanup
  call deallocate_arrays()

contains

  subroutine read_config()
    integer :: unit_num, ios

    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open benchmark.conf'
      stop 1
    end if

    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iterations
    read(unit_num, *) warmup_iterations
    read(unit_num, *) tolerance

    close(unit_num)
  end subroutine read_config

  subroutine read_parameters()
    integer :: unit_num, ios
    character(len=256) :: line, key, value_str
    integer :: eq_pos

    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open params.txt'
      stop 1
    end if

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit

      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      eq_pos = index(line, '=')
      if (eq_pos == 0) cycle

      key = adjustl(trim(line(1:eq_pos-1)))
      value_str = adjustl(trim(line(eq_pos+1:)))

      select case (trim(key))
        case ('ni')
          read(value_str, *) ni
        case ('nj')
          read(value_str, *) nj
        case ('nk')
          read(value_str, *) nk
        case ('rd')
          read(value_str, *) rd
        case ('epsav')
          read(value_str, *) epsav
        case ('rddvcp')
          read(value_str, *) rddvcp
        case ('p0iv')
          read(value_str, *) p0iv
      end select
    end do

    close(unit_num)
  end subroutine read_parameters

  subroutine allocate_arrays()
    allocate(zph(0:ni+1, 0:nj+1, 1:nk))

    allocate(ubr(0:ni+1, 0:nj+1, 1:nk))
    allocate(ubr_in(0:ni+1, 0:nj+1, 1:nk))
    allocate(vbr(0:ni+1, 0:nj+1, 1:nk))
    allocate(vbr_in(0:ni+1, 0:nj+1, 1:nk))
    allocate(pbr(0:ni+1, 0:nj+1, 1:nk))
    allocate(pbr_in(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptbr_in(0:ni+1, 0:nj+1, 1:nk))
    allocate(qvbr(0:ni+1, 0:nj+1, 1:nk))
    allocate(qvbr_in(0:ni+1, 0:nj+1, 1:nk))
    allocate(zph8s(0:ni+1, 0:nj+1, 1:nk))
    allocate(zph8s_in(0:ni+1, 0:nj+1, 1:nk))
    allocate(pibr(0:ni+1, 0:nj+1, 1:nk))
    allocate(pibr_in(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptvbr(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptvbr_in(0:ni+1, 0:nj+1, 1:nk))

    allocate(rbr(0:ni+1, 0:nj+1, 1:nk))

    allocate(ubr_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(vbr_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(pbr_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptbr_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(qvbr_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(zph8s_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(pibr_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptvbr_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(rbr_ref(0:ni+1, 0:nj+1, 1:nk))
  end subroutine allocate_arrays

  subroutine read_input_data()
    call read_3d_real(trim(data_dir)//'/zph.bin', zph, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/ubr_in.bin', ubr_in, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/vbr_in.bin', vbr_in, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/pbr_in.bin', pbr_in, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/ptbr_in.bin', ptbr_in, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/qvbr_in.bin', qvbr_in, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/zph8s_in.bin', zph8s_in, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/pibr_in.bin', pibr_in, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/ptvbr_in.bin', ptvbr_in, 0, ni+1, 0, nj+1, 1, nk)
  end subroutine read_input_data

  subroutine read_reference_data()
    call read_3d_real(trim(data_dir)//'/rbr_ref.bin', rbr_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/ubr_ref.bin', ubr_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/vbr_ref.bin', vbr_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/pbr_ref.bin', pbr_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/ptbr_ref.bin', ptbr_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/qvbr_ref.bin', qvbr_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/zph8s_ref.bin', zph8s_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/pibr_ref.bin', pibr_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_3d_real(trim(data_dir)//'/ptvbr_ref.bin', ptvbr_ref, 0, ni+1, 0, nj+1, 1, nk)
  end subroutine read_reference_data

  subroutine reset_inout_arrays()
    ubr(:,:,:) = ubr_in(:,:,:)
    vbr(:,:,:) = vbr_in(:,:,:)
    pbr(:,:,:) = pbr_in(:,:,:)
    ptbr(:,:,:) = ptbr_in(:,:,:)
    qvbr(:,:,:) = qvbr_in(:,:,:)
    zph8s(:,:,:) = zph8s_in(:,:,:)
    pibr(:,:,:) = pibr_in(:,:,:)
    ptvbr(:,:,:) = ptvbr_in(:,:,:)
  end subroutine reset_inout_arrays

  subroutine read_3d_real(filename, array, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: array(is:ie, js:je, ks:ke)
    integer :: unit_num, ios

    unit_num = 20
    open(unit=unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if
    read(unit_num) array
    close(unit_num)
  end subroutine read_3d_real

  subroutine s_setbase_kernel()
    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! Calculate z physical coordinates at scalar points
    do k = 1, nk-1
      !$omp do schedule(runtime) private(i, j)
      do j = 0, nj
        do i = 0, ni
          zph8s(i,j,k) = 0.5e0 * (zph(i,j,k+1) + zph(i,j,k))
        end do
      end do
      !$omp end do
    end do

    ! Get base state Exner function and density
    do k = 2, nk-2
      !$omp do schedule(runtime) private(i, j)
      do j = 0, nj
        do i = 0, ni
          ptvbr(i,j,k) = ptbr(i,j,k) * (1.e0 + epsav * qvbr(i,j,k)) / (1.e0 + qvbr(i,j,k))
          pibr(i,j,k) = exp(rddvcp * log(p0iv * pbr(i,j,k)))
          rbr(i,j,k) = pbr(i,j,k) / (rd * ptvbr(i,j,k) * pibr(i,j,k))
        end do
      end do
      !$omp end do
    end do

    !$omp end parallel
  end subroutine s_setbase_kernel

  subroutine validate_results(max_err_rbr, max_err_zph8s, max_err_pibr, max_err_ptvbr)
    real, intent(out) :: max_err_rbr, max_err_zph8s, max_err_pibr, max_err_ptvbr
    integer :: i, j, k
    real :: err, ref_val

    max_err_rbr = 0.0
    max_err_zph8s = 0.0
    max_err_pibr = 0.0
    max_err_ptvbr = 0.0

    ! Validate zph8s for k=1 to nk-1
    do k = 1, nk-1
      do j = 0, nj
        do i = 0, ni
          ref_val = abs(zph8s_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(zph8s(i,j,k) - zph8s_ref(i,j,k)) / ref_val
            max_err_zph8s = max(max_err_zph8s, err)
          end if
        end do
      end do
    end do

    ! Validate ptvbr, pibr, rbr for k=2 to nk-2
    do k = 2, nk-2
      do j = 0, nj
        do i = 0, ni
          ref_val = abs(rbr_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(rbr(i,j,k) - rbr_ref(i,j,k)) / ref_val
            max_err_rbr = max(max_err_rbr, err)
          end if

          ref_val = abs(pibr_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(pibr(i,j,k) - pibr_ref(i,j,k)) / ref_val
            max_err_pibr = max(max_err_pibr, err)
          end if

          ref_val = abs(ptvbr_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(ptvbr(i,j,k) - ptvbr_ref(i,j,k)) / ref_val
            max_err_ptvbr = max(max_err_ptvbr, err)
          end if
        end do
      end do
    end do
  end subroutine validate_results

  subroutine deallocate_arrays()
    deallocate(zph)
    deallocate(ubr, ubr_in, vbr, vbr_in, pbr, pbr_in)
    deallocate(ptbr, ptbr_in, qvbr, qvbr_in)
    deallocate(zph8s, zph8s_in, pibr, pibr_in, ptvbr, ptvbr_in)
    deallocate(rbr)
    deallocate(ubr_ref, vbr_ref, pbr_ref, ptbr_ref, qvbr_ref)
    deallocate(zph8s_ref, pibr_ref, ptvbr_ref, rbr_ref)
  end subroutine deallocate_arrays

end program kernel_benchmark

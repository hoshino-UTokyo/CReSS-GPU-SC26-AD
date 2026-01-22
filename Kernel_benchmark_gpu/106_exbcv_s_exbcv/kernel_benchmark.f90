!> GPU Kernel benchmark program for exbcv
!> Force lateral boundary values to external GPV values for v velocity
!> OpenACC version with Unified Memory
program kernel_benchmark
  use omp_lib
  implicit none

  ! Grid dimensions and parameters
  integer :: ni, nj, nk, nim1, nim2, njm1
  integer :: wbc, ebc
  integer :: ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub
  real :: dts, tdmpdt, ndmpdt, tpdt
  character(len=108) :: exbvar

  ! Input arrays
  real, allocatable :: vcpx(:,:,:)    ! Phase speed on west/east boundary
  real, allocatable :: vcpy(:,:,:)    ! Phase speed on south/north boundary
  real, allocatable :: vgpv(:,:,:)    ! GPV v velocity
  real, allocatable :: vtd(:,:,:)     ! Time tendency of GPV data

  ! Input/output array
  real, allocatable :: v(:,:,:)       ! v velocity
  real, allocatable :: v_init(:,:,:)
  real, allocatable :: v_ref(:,:,:)

  ! Benchmark parameters
  character(len=256) :: data_dir
  integer :: warmup_iterations, benchmark_iterations
  real :: tolerance

  ! Timing variables
  real(8), allocatable :: times(:)
  real(8) :: start_time, end_time
  real(8) :: avg_time, min_time, max_time, total_time

  ! Validation variables
  real :: max_rel_error
  integer :: error_count, total_elements

  integer :: iter

  ! Read benchmark configuration
  call read_config(data_dir, warmup_iterations, benchmark_iterations, tolerance)

  ! Read parameters
  call read_params(data_dir, ni, nj, nk, nim1, nim2, njm1, &
                   wbc, ebc, ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
                   dts, tdmpdt, ndmpdt, tpdt, exbvar)

  ! Print benchmark info
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: exbcv'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I3,A,I3)') ' wbc=', wbc, ', ebc=', ebc
  write(*,'(A,A)') ' exbvar=', trim(exbvar)
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', benchmark_iterations
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(vcpx(1:nj,1:nk,1:2))
  allocate(vcpy(1:ni,1:nk,1:2))
  allocate(vgpv(0:ni+1,0:nj+1,1:nk))
  allocate(vtd(0:ni+1,0:nj+1,1:nk))
  allocate(v(0:ni+1,0:nj+1,1:nk))
  allocate(v_init(0:ni+1,0:nj+1,1:nk))
  allocate(v_ref(0:ni+1,0:nj+1,1:nk))
  allocate(times(benchmark_iterations))

  ! Load input data
  write(*,'(A)') ' Loading input data...'
  call read_3d_array(trim(data_dir)//'/vcpx.bin', vcpx, 1, nj, 1, nk, 1, 2)
  call read_3d_array(trim(data_dir)//'/vcpy.bin', vcpy, 1, ni, 1, nk, 1, 2)
  call read_3d_array(trim(data_dir)//'/vgpv.bin', vgpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/vtd.bin', vtd, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/v_in.bin', v_init, 0, ni+1, 0, nj+1, 1, nk)

  ! Load reference output
  write(*,'(A)') ' Loading reference output...'
  call read_3d_array(trim(data_dir)//'/v_ref.bin', v_ref, 0, ni+1, 0, nj+1, 1, nk)

  ! Warmup iterations
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    v = v_init
    call kernel_exbcv()
    !$acc wait
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  do iter = 1, benchmark_iterations
    v = v_init
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_exbcv()
    !$acc wait
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
  end do

  ! Validate results
  write(*,'(A)') ' Validating output...'
  max_rel_error = 0.0
  error_count = 0
  total_elements = 0
  call validate_3d_array(v, v_ref, 0, ni+1, 0, nj+1, 1, nk, &
                         tolerance, max_rel_error, error_count, total_elements)

  ! Calculate timing statistics
  total_time = sum(times) * 1000.0d0
  avg_time = total_time / benchmark_iterations
  min_time = minval(times) * 1000.0d0
  max_time = maxval(times) * 1000.0d0

  ! Print results
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:    ', avg_time, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:        ', min_time, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:        ', max_time, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:      ', total_time, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,E12.4)') ' Max relative error: ', max_rel_error
  write(*,'(A,E12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', error_count
  if (error_count == 0) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  deallocate(vcpx, vcpy, vgpv, vtd, v, v_init, v_ref, times)

contains

  subroutine kernel_exbcv()
    integer :: i, j, k
    real :: vb1, vb2

    ! Force south boundary
    if (ebs == 1 .and. jsub == 0) then
      if (exbvar(2:2) == '-') then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent private(vb1,vb2)
          do i = 1, ni-1
            vb1 = vgpv(i,1,k) + vtd(i,1,k)*tpdt
            vb2 = vgpv(i,2,k) + vtd(i,2,k)*tpdt
            v(i,1,k) = v(i,1,k) + vtd(i,1,k)*dts &
              - vcpy(i,k,1)*((v(i,2,k)-v(i,1,k))-(vb2-vb1)) &
              - ndmpdt*(v(i,1,k)-vb1)
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 2, ni-2
            v(i,1,k) = v(i,1,k) + vtd(i,1,k)*dts
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! Force north boundary
    if (ebn == 1 .and. jsub == njsub-1) then
      if (exbvar(2:2) == '-') then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent private(vb1,vb2)
          do i = 1, ni-1
            vb1 = vgpv(i,nj,k) + vtd(i,nj,k)*tpdt
            vb2 = vgpv(i,njm1,k) + vtd(i,njm1,k)*tpdt
            v(i,nj,k) = v(i,nj,k) + vtd(i,nj,k)*dts &
              + vcpy(i,k,2)*((v(i,njm1,k)-v(i,nj,k))-(vb2-vb1)) &
              - ndmpdt*(v(i,nj,k)-vb1)
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 2, ni-2
            v(i,nj,k) = v(i,nj,k) + vtd(i,nj,k)*dts
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! Force west boundary
    if (ebw == 1 .and. isub == 0) then
      if (abs(wbc) /= 1) then
        if (exbvar(2:2) == '-' .or. exbvar(2:2) == '+') then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent private(vb1,vb2)
            do j = 2, nj-1
              vb1 = vgpv(1,j,k) + vtd(1,j,k)*tpdt
              vb2 = vgpv(2,j,k) + vtd(2,j,k)*tpdt
              v(1,j,k) = v(1,j,k) + vtd(1,j,k)*dts &
                - vcpx(j,k,1)*((v(2,j,k)-v(1,j,k))-(vb2-vb1)) &
                - tdmpdt*(v(1,j,k)-vb1)
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 1, nj
              v(1,j,k) = v(1,j,k) + vtd(1,j,k)*dts
            end do
          end do
          !$acc end kernels
        end if
      end if
    end if

    ! Force east boundary
    if (ebe == 1 .and. isub == nisub-1) then
      if (abs(ebc) /= 1) then
        if (exbvar(2:2) == '-' .or. exbvar(2:2) == '+') then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent private(vb1,vb2)
            do j = 2, nj-1
              vb1 = vgpv(nim1,j,k) + vtd(nim1,j,k)*tpdt
              vb2 = vgpv(nim2,j,k) + vtd(nim2,j,k)*tpdt
              v(nim1,j,k) = v(nim1,j,k) + vtd(nim1,j,k)*dts &
                + vcpx(j,k,2)*((v(nim2,j,k)-v(nim1,j,k))-(vb2-vb1)) &
                - tdmpdt*(v(nim1,j,k)-vb1)
            end do
          end do
          !$acc end kernels
        else
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            !$acc loop independent
            do j = 1, nj
              v(nim1,j,k) = v(nim1,j,k) + vtd(nim1,j,k)*dts
            end do
          end do
          !$acc end kernels
        end if
      end if
    end if

  end subroutine kernel_exbcv

  subroutine read_config(data_dir, warmup_iters, bench_iters, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: warmup_iters, bench_iters
    real, intent(out) :: tol
    integer :: unit_num

    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read')
    read(unit_num, '(A)') data_dir
    read(unit_num, *) bench_iters
    read(unit_num, *) warmup_iters
    read(unit_num, *) tol
    close(unit_num)
  end subroutine read_config

  subroutine read_params(data_dir, ni, nj, nk, nim1, nim2, njm1, &
                         wbc, ebc, ebw, ebe, ebs, ebn, isub, jsub, &
                         nisub, njsub, dts, tdmpdt, ndmpdt, tpdt, exbvar)
    character(len=*), intent(in) :: data_dir
    integer, intent(out) :: ni, nj, nk, nim1, nim2, njm1
    integer, intent(out) :: wbc, ebc
    integer, intent(out) :: ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub
    real, intent(out) :: dts, tdmpdt, ndmpdt, tpdt
    character(len=*), intent(out) :: exbvar
    character(len=256) :: line, key, sval
    integer :: unit_num, ios, eq_pos
    real :: val

    exbvar = ''
    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', action='read')
    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        sval = adjustl(line(eq_pos+1:))
        select case(trim(key))
          case('exbvar')
            exbvar = trim(sval)
          case default
            read(sval, *) val
            select case(trim(key))
              case('ni'); ni = int(val)
              case('nj'); nj = int(val)
              case('nk'); nk = int(val)
              case('nim1'); nim1 = int(val)
              case('nim2'); nim2 = int(val)
              case('njm1'); njm1 = int(val)
              case('wbc'); wbc = int(val)
              case('ebc'); ebc = int(val)
              case('ebw'); ebw = int(val)
              case('ebe'); ebe = int(val)
              case('ebs'); ebs = int(val)
              case('ebn'); ebn = int(val)
              case('isub'); isub = int(val)
              case('jsub'); jsub = int(val)
              case('nisub'); nisub = int(val)
              case('njsub'); njsub = int(val)
              case('dts'); dts = val
              case('tdmpdt'); tdmpdt = val
              case('ndmpdt'); ndmpdt = val
              case('tpdt'); tpdt = val
            end select
        end select
      end if
    end do
    close(unit_num)
  end subroutine read_params

  subroutine read_3d_array(filename, arr, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: unit_num

    unit_num = 12
    open(unit=unit_num, file=filename, form='unformatted', access='stream', status='old')
    read(unit_num) arr
    close(unit_num)
  end subroutine read_3d_array

  subroutine validate_3d_array(computed, reference, is, ie, js, je, ks, ke, &
                               tol, max_err, err_count, total_count)
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(in) :: computed(is:ie, js:je, ks:ke)
    real, intent(in) :: reference(is:ie, js:je, ks:ke)
    real, intent(in) :: tol
    real, intent(inout) :: max_err
    integer, intent(inout) :: err_count, total_count
    real :: rel_err, abs_ref
    integer :: i, j, k

    do k = ks, ke
      do j = js, je
        do i = is, ie
          total_count = total_count + 1
          abs_ref = abs(reference(i,j,k))
          if (abs_ref > 1.0e-30) then
            rel_err = abs(computed(i,j,k) - reference(i,j,k)) / abs_ref
          else
            rel_err = abs(computed(i,j,k) - reference(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do
  end subroutine validate_3d_array

  function get_time() result(t)
    real(8) :: t
    integer(8) :: count, count_rate
    call system_clock(count, count_rate)
    t = dble(count) / dble(count_rate)
  end function get_time

end program kernel_benchmark

!> Kernel benchmark program for exbcu
!> Force lateral boundary values to external GPV values for u velocity
program kernel_benchmark
  implicit none

  ! Grid dimensions and parameters
  integer :: ni, nj, nk, nim1, njm1, njm2
  integer :: wbc, ebc
  integer :: ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub
  real :: dts, tdmpdt, ndmpdt, tpdt
  character(len=108) :: exbvar

  ! Input arrays
  real, allocatable :: ucpx(:,:,:)    ! Phase speed on west/east boundary
  real, allocatable :: ucpy(:,:,:)    ! Phase speed on south/north boundary
  real, allocatable :: ugpv(:,:,:)    ! GPV u velocity
  real, allocatable :: utd(:,:,:)     ! Time tendency of GPV data

  ! Input/output array
  real, allocatable :: u(:,:,:)       ! u velocity
  real, allocatable :: u_init(:,:,:)
  real, allocatable :: u_ref(:,:,:)

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
  call read_params(data_dir, ni, nj, nk, nim1, njm1, njm2, &
                   wbc, ebc, ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
                   dts, tdmpdt, ndmpdt, tpdt, exbvar)

  ! Print benchmark info
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: exbcu'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I3,A,I3)') ' wbc=', wbc, ', ebc=', ebc
  write(*,'(A,A)') ' exbvar=', trim(exbvar)
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', benchmark_iterations
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(ucpx(1:nj,1:nk,1:2))
  allocate(ucpy(1:ni,1:nk,1:2))
  allocate(ugpv(0:ni+1,0:nj+1,1:nk))
  allocate(utd(0:ni+1,0:nj+1,1:nk))
  allocate(u(0:ni+1,0:nj+1,1:nk))
  allocate(u_init(0:ni+1,0:nj+1,1:nk))
  allocate(u_ref(0:ni+1,0:nj+1,1:nk))
  allocate(times(benchmark_iterations))

  ! Load input data
  write(*,'(A)') ' Loading input data...'
  call read_3d_array(trim(data_dir)//'/ucpx.bin', ucpx, 1, nj, 1, nk, 1, 2)
  call read_3d_array(trim(data_dir)//'/ucpy.bin', ucpy, 1, ni, 1, nk, 1, 2)
  call read_3d_array(trim(data_dir)//'/ugpv.bin', ugpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/utd.bin', utd, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/u_in.bin', u_init, 0, ni+1, 0, nj+1, 1, nk)

  ! Load reference output
  write(*,'(A)') ' Loading reference output...'
  call read_3d_array(trim(data_dir)//'/u_ref.bin', u_ref, 0, ni+1, 0, nj+1, 1, nk)

  ! Warmup iterations
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    u = u_init
    call kernel_exbcu()
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  do iter = 1, benchmark_iterations
    u = u_init
    start_time = get_time()
    call kernel_exbcu()
    end_time = get_time()
    times(iter) = end_time - start_time
  end do

  ! Validate results
  write(*,'(A)') ' Validating output...'
  max_rel_error = 0.0
  error_count = 0
  total_elements = 0
  call validate_3d_array(u, u_ref, 0, ni+1, 0, nj+1, 1, nk, &
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

  deallocate(ucpx, ucpy, ugpv, utd, u, u_init, u_ref, times)

contains

  subroutine kernel_exbcu()
    integer :: i, j, k
    real :: ub1, ub2

    !$omp parallel default(shared)

    ! Force west boundary
    if (ebw == 1 .and. isub == 0) then
      if (abs(wbc) /= 1) then
        if (exbvar(1:1) == '-') then
          !$omp do schedule(runtime) private(j,k,ub1,ub2)
          do k = 2, nk-2
            do j = 1, nj-1
              ub1 = ugpv(1,j,k) + utd(1,j,k)*tpdt
              ub2 = ugpv(2,j,k) + utd(2,j,k)*tpdt
              u(1,j,k) = u(1,j,k) + utd(1,j,k)*dts &
                - ucpx(j,k,1)*((u(2,j,k)-u(1,j,k))-(ub2-ub1)) &
                - ndmpdt*(u(1,j,k)-ub1)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(j,k)
          do k = 2, nk-2
            do j = 2, nj-2
              u(1,j,k) = u(1,j,k) + utd(1,j,k)*dts
            end do
          end do
          !$omp end do
        end if
      end if
    end if

    ! Force east boundary
    if (ebe == 1 .and. isub == nisub-1) then
      if (abs(ebc) /= 1) then
        if (exbvar(1:1) == '-') then
          !$omp do schedule(runtime) private(j,k,ub1,ub2)
          do k = 2, nk-2
            do j = 1, nj-1
              ub1 = ugpv(ni,j,k) + utd(ni,j,k)*tpdt
              ub2 = ugpv(nim1,j,k) + utd(nim1,j,k)*tpdt
              u(ni,j,k) = u(ni,j,k) + utd(ni,j,k)*dts &
                + ucpx(j,k,2)*((u(nim1,j,k)-u(ni,j,k))-(ub2-ub1)) &
                - ndmpdt*(u(ni,j,k)-ub1)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(j,k)
          do k = 2, nk-2
            do j = 2, nj-2
              u(ni,j,k) = u(ni,j,k) + utd(ni,j,k)*dts
            end do
          end do
          !$omp end do
        end if
      end if
    end if

    ! Force south boundary
    if (ebs == 1 .and. jsub == 0) then
      if (exbvar(1:1) == '-' .or. exbvar(1:1) == '+') then
        !$omp do schedule(runtime) private(i,k,ub1,ub2)
        do k = 2, nk-2
          do i = 2, ni-1
            ub1 = ugpv(i,1,k) + utd(i,1,k)*tpdt
            ub2 = ugpv(i,2,k) + utd(i,2,k)*tpdt
            u(i,1,k) = u(i,1,k) + utd(i,1,k)*dts &
              - ucpy(i,k,1)*((u(i,2,k)-u(i,1,k))-(ub2-ub1)) &
              - tdmpdt*(u(i,1,k)-ub1)
          end do
        end do
        !$omp end do
      else
        !$omp do schedule(runtime) private(i,k)
        do k = 2, nk-2
          do i = 1, ni
            u(i,1,k) = u(i,1,k) + utd(i,1,k)*dts
          end do
        end do
        !$omp end do
      end if
    end if

    ! Force north boundary
    if (ebn == 1 .and. jsub == njsub-1) then
      if (exbvar(1:1) == '-' .or. exbvar(1:1) == '+') then
        !$omp do schedule(runtime) private(i,k,ub1,ub2)
        do k = 2, nk-2
          do i = 2, ni-1
            ub1 = ugpv(i,njm1,k) + utd(i,njm1,k)*tpdt
            ub2 = ugpv(i,njm2,k) + utd(i,njm2,k)*tpdt
            u(i,njm1,k) = u(i,njm1,k) + utd(i,njm1,k)*dts &
              + ucpy(i,k,2)*((u(i,njm2,k)-u(i,njm1,k))-(ub2-ub1)) &
              - tdmpdt*(u(i,njm1,k)-ub1)
          end do
        end do
        !$omp end do
      else
        !$omp do schedule(runtime) private(i,k)
        do k = 2, nk-2
          do i = 1, ni
            u(i,njm1,k) = u(i,njm1,k) + utd(i,njm1,k)*dts
          end do
        end do
        !$omp end do
      end if
    end if

    !$omp end parallel
  end subroutine kernel_exbcu

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

  subroutine read_params(data_dir, ni, nj, nk, nim1, njm1, njm2, &
                         wbc, ebc, ebw, ebe, ebs, ebn, isub, jsub, &
                         nisub, njsub, dts, tdmpdt, ndmpdt, tpdt, exbvar)
    character(len=*), intent(in) :: data_dir
    integer, intent(out) :: ni, nj, nk, nim1, njm1, njm2
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
              case('njm1'); njm1 = int(val)
              case('njm2'); njm2 = int(val)
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

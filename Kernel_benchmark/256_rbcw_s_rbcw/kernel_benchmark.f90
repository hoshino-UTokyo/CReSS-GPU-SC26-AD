!> Kernel benchmark program for rbcw
!> Radiative lateral boundary conditions for w velocity
program kernel_benchmark
  implicit none

  ! Grid dimensions and parameters
  integer :: ni, nj, nk, nim1, nim2, njm1, njm2
  integer :: wbc, ebc, sbc, nbc
  integer :: nggopt, lspopt, vspopt
  integer :: ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub
  real :: gtinc, dmpdt, tpdt
  character(len=108) :: gpvvar

  ! Input arrays
  real, allocatable :: wcpx(:,:,:)    ! Phase speed on west/east boundary
  real, allocatable :: wcpy(:,:,:)    ! Phase speed on south/north boundary
  real, allocatable :: wgpv(:,:,:)    ! GPV w velocity
  real, allocatable :: wtd(:,:,:)     ! Time tendency of GPV data

  ! Input/output array
  real, allocatable :: w(:,:,:)       ! w velocity
  real, allocatable :: w_init(:,:,:)
  real, allocatable :: w_ref(:,:,:)

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
  call read_params(data_dir, ni, nj, nk, nim1, nim2, njm1, njm2, &
                   wbc, ebc, sbc, nbc, nggopt, lspopt, vspopt, &
                   ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
                   gtinc, dmpdt, tpdt, gpvvar)

  ! Print benchmark info
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: rbcw'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I3,A,I3,A,I3,A,I3)') ' wbc=', wbc, ', ebc=', ebc, ', sbc=', sbc, ', nbc=', nbc
  write(*,'(A,A)') ' gpvvar=', trim(gpvvar)
  write(*,'(A,I3)') ' vspopt=', vspopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', benchmark_iterations
  write(*,'(A)') '=================================================='

  ! Allocate arrays
  allocate(wcpx(1:nj,1:nk,1:2))
  allocate(wcpy(1:ni,1:nk,1:2))
  allocate(wgpv(0:ni+1,0:nj+1,1:nk))
  allocate(wtd(0:ni+1,0:nj+1,1:nk))
  allocate(w(0:ni+1,0:nj+1,1:nk))
  allocate(w_init(0:ni+1,0:nj+1,1:nk))
  allocate(w_ref(0:ni+1,0:nj+1,1:nk))
  allocate(times(benchmark_iterations))

  ! Load input data
  write(*,'(A)') ' Loading input data...'
  call read_3d_array(trim(data_dir)//'/wcpx.bin', wcpx, 1, nj, 1, nk, 1, 2)
  call read_3d_array(trim(data_dir)//'/wcpy.bin', wcpy, 1, ni, 1, nk, 1, 2)
  call read_3d_array(trim(data_dir)//'/wgpv.bin', wgpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/wtd.bin', wtd, 0, ni+1, 0, nj+1, 1, nk)
  call read_3d_array(trim(data_dir)//'/w_in.bin', w_init, 0, ni+1, 0, nj+1, 1, nk)

  ! Load reference output
  write(*,'(A)') ' Loading reference output...'
  call read_3d_array(trim(data_dir)//'/w_ref.bin', w_ref, 0, ni+1, 0, nj+1, 1, nk)

  ! Warmup iterations
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    w = w_init
    call kernel_rbcw()
  end do

  ! Benchmark iterations
  write(*,'(A)') ' Running benchmark iterations...'
  do iter = 1, benchmark_iterations
    w = w_init
    start_time = get_time()
    call kernel_rbcw()
    end_time = get_time()
    times(iter) = end_time - start_time
  end do

  ! Validate results
  write(*,'(A)') ' Validating output...'
  max_rel_error = 0.0
  error_count = 0
  total_elements = 0
  call validate_3d_array(w, w_ref, 0, ni+1, 0, nj+1, 1, nk, &
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

  deallocate(wcpx, wcpy, wgpv, wtd, w, w_init, w_ref, times)

contains

  subroutine kernel_rbcw()
    integer :: i, j, k
    real :: radwe, radsn
    logical :: use_gpv

    use_gpv = (gpvvar(1:1) == 'o') .and. &
              (nggopt == 1 .or. mod(lspopt,10) == 1 .or. vspopt == 1)

    !$omp parallel default(shared)

    ! Four corners
    if (ebs == 1 .and. jsub == 0) then
      if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. sbc >= 4) then
        if (use_gpv) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-1
            radwe = (w(2,1,k) - w(1,1,k)) * wcpx(1,k,1)
            radsn = (w(1,2,k) - w(1,1,k)) * wcpy(1,k,1)
            w(1,1,k) = w(1,1,k) - (radwe+radsn) &
              - dmpdt * (w(1,1,k) - (wgpv(1,1,k) + wtd(1,1,k)*tpdt))
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-1
            radwe = (w(2,1,k) - w(1,1,k)) * wcpx(1,k,1)
            radsn = (w(1,2,k) - w(1,1,k)) * wcpy(1,k,1)
            w(1,1,k) = w(1,1,k) - (radwe+radsn) - dmpdt*w(1,1,k)
          end do
          !$omp end do
        end if
      end if
      if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. sbc >= 4) then
        if (use_gpv) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-1
            radwe = (w(nim2,1,k) - w(nim1,1,k)) * wcpx(1,k,2)
            radsn = (w(nim1,2,k) - w(nim1,1,k)) * wcpy(nim1,k,1)
            w(nim1,1,k) = w(nim1,1,k) + (radwe-radsn) &
              - dmpdt * (w(nim1,1,k) - (wgpv(nim1,1,k) + wtd(nim1,1,k)*tpdt))
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-1
            radwe = (w(nim2,1,k) - w(nim1,1,k)) * wcpx(1,k,2)
            radsn = (w(nim1,2,k) - w(nim1,1,k)) * wcpy(nim1,k,1)
            w(nim1,1,k) = w(nim1,1,k) + (radwe-radsn) - dmpdt*w(nim1,1,k)
          end do
          !$omp end do
        end if
      end if
    end if

    if (ebn == 1 .and. jsub == njsub-1) then
      if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. nbc >= 4) then
        if (use_gpv) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-1
            radwe = (w(2,njm1,k) - w(1,njm1,k)) * wcpx(njm1,k,1)
            radsn = (w(1,njm2,k) - w(1,njm1,k)) * wcpy(1,k,2)
            w(1,njm1,k) = w(1,njm1,k) - (radwe-radsn) &
              - dmpdt * (w(1,njm1,k) - (wgpv(1,njm1,k) + wtd(1,njm1,k)*tpdt))
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-1
            radwe = (w(2,njm1,k) - w(1,njm1,k)) * wcpx(njm1,k,1)
            radsn = (w(1,njm2,k) - w(1,njm1,k)) * wcpy(1,k,2)
            w(1,njm1,k) = w(1,njm1,k) - (radwe-radsn) - dmpdt*w(1,njm1,k)
          end do
          !$omp end do
        end if
      end if
      if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. nbc >= 4) then
        if (use_gpv) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-1
            radwe = (w(nim2,njm1,k) - w(nim1,njm1,k)) * wcpx(njm1,k,2)
            radsn = (w(nim1,njm2,k) - w(nim1,njm1,k)) * wcpy(nim1,k,2)
            w(nim1,njm1,k) = w(nim1,njm1,k) + (radwe+radsn) &
              - dmpdt * (w(nim1,njm1,k) - (wgpv(nim1,njm1,k) + wtd(nim1,njm1,k)*tpdt))
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-1
            radwe = (w(nim2,njm1,k) - w(nim1,njm1,k)) * wcpx(njm1,k,2)
            radsn = (w(nim1,njm2,k) - w(nim1,njm1,k)) * wcpy(nim1,k,2)
            w(nim1,njm1,k) = w(nim1,njm1,k) + (radwe+radsn) - dmpdt*w(nim1,njm1,k)
          end do
          !$omp end do
        end if
      end if
    end if

    ! West boundary
    if (ebw == 1 .and. isub == 0 .and. wbc >= 4) then
      if (use_gpv) then
        !$omp do schedule(runtime) private(j,k)
        do k = 2, nk-1
          do j = 2, nj-2
            w(1,j,k) = w(1,j,k) - wcpx(j,k,1)*(w(2,j,k)-w(1,j,k)) &
              - dmpdt * (w(1,j,k) - (wgpv(1,j,k) + wtd(1,j,k)*tpdt))
          end do
        end do
        !$omp end do
      else
        !$omp do schedule(runtime) private(j,k)
        do k = 2, nk-1
          do j = 2, nj-2
            w(1,j,k) = w(1,j,k) - wcpx(j,k,1)*(w(2,j,k)-w(1,j,k)) - dmpdt*w(1,j,k)
          end do
        end do
        !$omp end do
      end if
    end if

    ! East boundary
    if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4) then
      if (use_gpv) then
        !$omp do schedule(runtime) private(j,k)
        do k = 2, nk-1
          do j = 2, nj-2
            w(nim1,j,k) = w(nim1,j,k) + wcpx(j,k,2)*(w(nim2,j,k)-w(nim1,j,k)) &
              - dmpdt * (w(nim1,j,k) - (wgpv(nim1,j,k) + wtd(nim1,j,k)*tpdt))
          end do
        end do
        !$omp end do
      else
        !$omp do schedule(runtime) private(j,k)
        do k = 2, nk-1
          do j = 2, nj-2
            w(nim1,j,k) = w(nim1,j,k) + wcpx(j,k,2)*(w(nim2,j,k)-w(nim1,j,k)) &
              - dmpdt*w(nim1,j,k)
          end do
        end do
        !$omp end do
      end if
    end if

    ! South boundary
    if (ebs == 1 .and. jsub == 0 .and. sbc >= 4) then
      if (use_gpv) then
        !$omp do schedule(runtime) private(i,k)
        do k = 2, nk-1
          do i = 2, ni-2
            w(i,1,k) = w(i,1,k) - wcpy(i,k,1)*(w(i,2,k)-w(i,1,k)) &
              - dmpdt * (w(i,1,k) - (wgpv(i,1,k) + wtd(i,1,k)*tpdt))
          end do
        end do
        !$omp end do
      else
        !$omp do schedule(runtime) private(i,k)
        do k = 2, nk-1
          do i = 2, ni-2
            w(i,1,k) = w(i,1,k) - wcpy(i,k,1)*(w(i,2,k)-w(i,1,k)) - dmpdt*w(i,1,k)
          end do
        end do
        !$omp end do
      end if
    end if

    ! North boundary
    if (ebn == 1 .and. jsub == njsub-1 .and. nbc >= 4) then
      if (use_gpv) then
        !$omp do schedule(runtime) private(i,k)
        do k = 2, nk-1
          do i = 2, ni-2
            w(i,njm1,k) = w(i,njm1,k) + wcpy(i,k,2)*(w(i,njm2,k)-w(i,njm1,k)) &
              - dmpdt * (w(i,njm1,k) - (wgpv(i,njm1,k) + wtd(i,njm1,k)*tpdt))
          end do
        end do
        !$omp end do
      else
        !$omp do schedule(runtime) private(i,k)
        do k = 2, nk-1
          do i = 2, ni-2
            w(i,njm1,k) = w(i,njm1,k) + wcpy(i,k,2)*(w(i,njm2,k)-w(i,njm1,k)) &
              - dmpdt*w(i,njm1,k)
          end do
        end do
        !$omp end do
      end if
    end if

    !$omp end parallel
  end subroutine kernel_rbcw

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

  subroutine read_params(data_dir, ni, nj, nk, nim1, nim2, njm1, njm2, &
                         wbc, ebc, sbc, nbc, nggopt, lspopt, vspopt, &
                         ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
                         gtinc, dmpdt, tpdt, gpvvar)
    character(len=*), intent(in) :: data_dir
    integer, intent(out) :: ni, nj, nk, nim1, nim2, njm1, njm2
    integer, intent(out) :: wbc, ebc, sbc, nbc, nggopt, lspopt, vspopt
    integer, intent(out) :: ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub
    real, intent(out) :: gtinc, dmpdt, tpdt
    character(len=*), intent(out) :: gpvvar
    character(len=256) :: line, key, sval
    integer :: unit_num, ios, eq_pos
    real :: val

    gpvvar = ''
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
          case('gpvvar'); gpvvar = trim(sval)
          case('lbcvar'); continue  ! Skip string variable
          case default
            read(sval, *, iostat=ios) val
            if (ios /= 0) cycle  ! Skip unparseable values
            select case(trim(key))
              case('ni'); ni = int(val)
              case('nj'); nj = int(val)
              case('nk'); nk = int(val)
              case('nim1'); nim1 = int(val)
              case('nim2'); nim2 = int(val)
              case('njm1'); njm1 = int(val)
              case('njm2'); njm2 = int(val)
              case('wbc'); wbc = int(val)
              case('ebc'); ebc = int(val)
              case('sbc'); sbc = int(val)
              case('nbc'); nbc = int(val)
              case('nggopt'); nggopt = int(val)
              case('lspopt'); lspopt = int(val)
              case('vspopt'); vspopt = int(val)
              case('ebw'); ebw = int(val)
              case('ebe'); ebe = int(val)
              case('ebs'); ebs = int(val)
              case('ebn'); ebn = int(val)
              case('isub'); isub = int(val)
              case('jsub'); jsub = int(val)
              case('nisub'); nisub = int(val)
              case('njsub'); njsub = int(val)
              case('gtinc'); gtinc = val
              case('dmpdt'); dmpdt = val
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

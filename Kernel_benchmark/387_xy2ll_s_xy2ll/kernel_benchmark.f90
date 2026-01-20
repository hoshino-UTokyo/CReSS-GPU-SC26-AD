!***********************************************************************
! Kernel Benchmark: xy2ll (s_xy2ll)
!***********************************************************************
!
! Source: Src/xy2ll.f90
! Description: Convert x,y map coordinates to latitude/longitude using
!              various map projections (lat-lon, Polar Stereographic,
!              Lambert, Mercator, etc.)
!
!***********************************************************************
program kernel_benchmark_xy2ll
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: imin, imax, jmin, jmax

  ! Parameters
  integer :: mpopt, nspol, ncpn
  real :: tlon, x0, y0, r2d
  character(len=8) :: pname
  character(len=2) :: xo

  ! Derived parameters
  real :: rpol, r2d2, r2d3, r2d5, tlonw

  ! Arrays
  real, allocatable :: cpj(:), x(:), y(:)
  real, allocatable :: lat(:,:), lon(:,:)
  real, allocatable :: lat_ref(:,:), lon_ref(:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, tmp_max_err
  real :: tolerance
  integer :: error_count, tmp_err_count
  logical :: validation_passed

  ! Loop variables
  integer :: iter

  ! Constants
  real, parameter :: eps = 1.0e-35
  real, parameter :: d2r = 1.745329251994e-2

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt')

  ! Set derived parameters
  rpol = real(nspol)
  r2d2 = 2.0e0 * r2d
  tlonw = tlon + 180.0e0

  ! Set pname as "solver" for this benchmark
  pname = 'solver'
  ncpn = 6
  xo = 'xx'

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: xy2ll'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6)') ' Grid: imin=', imin, ', imax=', imax
  write(*,'(A,I6,A,I6)') '       jmin=', jmin, ', jmax=', jmax
  write(*,'(A,I6)') ' mpopt=', mpopt
  write(*,'(A,ES12.4)') ' tlon=', tlon
  write(*,'(A,ES12.4,A,ES12.4)') ' x0=', x0, ', y0=', y0
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(cpj(1:7))
  allocate(x(imin:imax))
  allocate(y(jmin:jmax))
  allocate(lat(imin:imax, jmin:jmax))
  allocate(lon(imin:imax, jmin:jmax))
  allocate(lat_ref(imin:imax, jmin:jmax))
  allocate(lon_ref(imin:imax, jmin:jmax))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_1d(trim(data_dir)//'/cpj.bin', cpj, 1, 7)
  call read_array_1d(trim(data_dir)//'/x.bin', x, imin, imax)
  call read_array_1d(trim(data_dir)//'/y.bin', y, jmin, jmax)

  ! Update derived parameters that depend on cpj
  r2d3 = cpj(3) * r2d
  r2d5 = cpj(5) * r2d

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_2d(trim(data_dir)//'/lat_ref.bin', lat_ref, imin, imax, jmin, jmax)
  call read_array_2d(trim(data_dir)//'/lon_ref.bin', lon_ref, imin, imax, jmin, jmax)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    lat = 0.0
    lon = 0.0
    call kernel_xy2ll(mpopt, pname, ncpn, x0, y0, cpj, &
                      imin, imax, jmin, jmax, x, y, lat, lon, &
                      rpol, r2d, r2d2, r2d3, r2d5, tlon, tlonw, eps, d2r)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    lat = 0.0
    lon = 0.0
    t_start = omp_get_wtime()

    call kernel_xy2ll(mpopt, pname, ncpn, x0, y0, cpj, &
                      imin, imax, jmin, jmax, x, y, lat, lon, &
                      rpol, r2d, r2d2, r2d3, r2d5, tlon, tlonw, eps, d2r)

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
  error_count = 0

  call validate_output(lat, lat_ref, imin, imax, jmin, jmax, tolerance, tmp_max_err, tmp_err_count)
  if (tmp_max_err > max_error) max_error = tmp_max_err
  error_count = error_count + tmp_err_count
  if (tmp_err_count > 0) write(*,'(A,I6)') '  lat errors: ', tmp_err_count

  call validate_output(lon, lon_ref, imin, imax, jmin, jmax, tolerance, tmp_max_err, tmp_err_count)
  if (tmp_max_err > max_error) max_error = tmp_max_err
  error_count = error_count + tmp_err_count
  if (tmp_err_count > 0) write(*,'(A,I6)') '  lon errors: ', tmp_err_count

  validation_passed = (error_count == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time:  ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:      ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:      ', t_max * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:    ', t_total * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)')    ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(cpj, x, y, lat, lon, lat_ref, lon_ref, times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: xy2ll
  !-------------------------------------------------------------------
  subroutine kernel_xy2ll(mpopt, pname, ncpn, x0, y0, cpj, &
                          imin, imax, jmin, jmax, x, y, lat, lon, &
                          rpol, r2d, r2d2, r2d3, r2d5, tlon, tlonw, eps, d2r)
    implicit none

    integer, intent(in) :: mpopt, ncpn
    character(len=8), intent(in) :: pname
    real, intent(in) :: x0, y0
    real, intent(in) :: cpj(1:7)
    integer, intent(in) :: imin, imax, jmin, jmax
    real, intent(in) :: x(imin:imax), y(jmin:jmax)
    real, intent(out) :: lat(imin:imax, jmin:jmax)
    real, intent(out) :: lon(imin:imax, jmin:jmax)
    real, intent(in) :: rpol, r2d, r2d2, r2d3, r2d5, tlon, tlonw, eps, d2r

    integer :: i, j
    integer :: istr, iend, jstr, jend
    real :: xx, yy, rr

    ! Set loop bounds (simplified getindx - exclude last index based on dump data)
    istr = imin
    iend = imax - 1
    jstr = jmin
    jend = jmax - 1

    !$omp parallel default(shared)

    ! mpopt == 0 or 10: lat-lon coordinates
    if (mpopt.eq.0.or.mpopt.eq.10) then
      if (pname(1:ncpn).eq.'solver') then
        !$omp do schedule(runtime) private(i,j,xx,yy)
        do j = jstr, jend
          do i = istr, iend
            xx = x(i) + x0
            yy = y(j) + y0

            lat(i,j) = max(min(yy*r2d3, 90.e0), -90.e0)
            lon(i,j) = xx * r2d3

            if (lon(i,j).gt.180.e0) lon(i,j) = lon(i,j) - 360.e0
            if (lon(i,j).lt.-180.e0) lon(i,j) = lon(i,j) + 360.e0
          end do
        end do
        !$omp end do
      else
        !$omp do schedule(runtime) private(i,j)
        do j = jstr, jend
          do i = istr, iend
            lat(i,j) = max(min(y(j)+y0, 90.e0), -90.e0)
            lon(i,j) = x(i) + x0

            if (lon(i,j).gt.180.e0) lon(i,j) = lon(i,j) - 360.e0
            if (lon(i,j).lt.-180.e0) lon(i,j) = lon(i,j) + 360.e0
          end do
        end do
        !$omp end do
      end if

    ! mpopt == 1: Polar Stereographic
    else if (mpopt.eq.1) then
      !$omp do schedule(runtime) private(i,j,xx,yy,rr)
      do j = jstr, jend
        do i = istr, iend
          xx = x(i) + x0
          yy = rpol * y(j) + y0

          rr = sqrt(xx*xx + yy*yy) * cpj(3)
          lat(i,j) = max(min(rpol*(90.e0 - r2d2*atan(rr)), 90.e0), -90.e0)

          if (yy.gt.0.e0) then
            lon(i,j) = tlonw + atan(-xx/(yy+eps)) * r2d
          else
            lon(i,j) = tlon + atan(-xx/(yy-eps)) * r2d
          end if

          if (lon(i,j).gt.180.e0) lon(i,j) = lon(i,j) - 360.e0
          if (lon(i,j).lt.-180.e0) lon(i,j) = lon(i,j) + 360.e0
        end do
      end do
      !$omp end do

    ! mpopt == 2: Lambert Conformal Conic
    else if (mpopt.eq.2) then
      !$omp do schedule(runtime) private(i,j,xx,yy,rr)
      do j = jstr, jend
        do i = istr, iend
          xx = x(i) + x0
          yy = rpol * y(j) + y0

          rr = cpj(2) * exp(cpj(5) * log(sqrt(xx*xx + yy*yy) * cpj(7) + eps))
          lat(i,j) = max(min(rpol*(90.e0 - r2d2*atan(rr)), 90.e0), -90.e0)

          if (yy.gt.0.e0) then
            lon(i,j) = tlonw + atan(-xx/(yy+eps)) * r2d5
          else
            lon(i,j) = tlon + atan(-xx/(yy-eps)) * r2d5
          end if

          if (lon(i,j).gt.180.e0) lon(i,j) = lon(i,j) - 360.e0
          if (lon(i,j).lt.-180.e0) lon(i,j) = lon(i,j) + 360.e0
        end do
      end do
      !$omp end do

    ! mpopt == 3 or 13: Mercator
    else if (mpopt.eq.3.or.mpopt.eq.13) then
      !$omp do schedule(runtime) private(i,j,xx,yy)
      do j = jstr, jend
        do i = istr, iend
          xx = x(i) + x0
          yy = y(j) + y0

          lat(i,j) = max(min(90.e0 - 2.e0*atan(exp(-yy*cpj(3)))*r2d, 90.e0), -90.e0)
          lon(i,j) = xx * r2d3

          if (lon(i,j).gt.180.e0) lon(i,j) = lon(i,j) - 360.e0
          if (lon(i,j).lt.-180.e0) lon(i,j) = lon(i,j) + 360.e0
        end do
      end do
      !$omp end do

    ! mpopt == 4: direct coordinates
    else if (mpopt.eq.4) then
      !$omp do schedule(runtime) private(i,j,xx,yy)
      do j = jstr, jend
        do i = istr, iend
          xx = x(i) + x0
          yy = y(j) + y0

          lat(i,j) = max(min(cpj(2)*yy, 90.e0), -90.e0)
          lon(i,j) = cpj(2)*xx/cos(lat(i,j)*d2r) + tlon

          if (lon(i,j).gt.180.e0) lon(i,j) = lon(i,j) - 360.e0
          if (lon(i,j).lt.-180.e0) lon(i,j) = lon(i,j) + 360.e0
        end do
      end do
      !$omp end do

    ! mpopt == 5: circular cylinder
    else if (mpopt.eq.5) then
      !$omp do schedule(runtime) private(i,j,xx)
      do j = jstr, jend
        do i = istr, iend
          xx = x(i) + x0

          lat(i,j) = max(min(cpj(2)*y0, 90.e0), -90.e0)
          lon(i,j) = cpj(2)*xx/cos(lat(i,j)*d2r)

          if (lon(i,j).gt.180.e0) lon(i,j) = lon(i,j) - 360.e0
          if (lon(i,j).lt.-180.e0) lon(i,j) = lon(i,j) + 360.e0
        end do
      end do
      !$omp end do
    end if

    !$omp end parallel

  end subroutine kernel_xy2ll

  !-------------------------------------------------------------------
  ! Read benchmark configuration
  !-------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios == 0) then
      read(10, '(A)', iostat=ios) data_dir
      read(10, *, iostat=ios) num_iter
      read(10, *, iostat=ios) warmup_iter
      read(10, *, iostat=ios) tol
      close(10)
    end if
  end subroutine read_config

  !-------------------------------------------------------------------
  ! Read parameters
  !-------------------------------------------------------------------
  subroutine read_parameters(filename)
    implicit none
    character(len=*), intent(in) :: filename

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Set defaults
    mpopt = 0
    nspol = 1
    tlon = 0.0
    x0 = 0.0
    y0 = 0.0
    r2d = 57.29578

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
          case ('mpopt')
            read(val, *) mpopt
          case ('nspol')
            read(val, *) nspol
          case ('tlon')
            read(val, *) tlon
          case ('ncpn')
            read(val, *) ncpn
          case ('imin')
            read(val, *) imin
          case ('imax')
            read(val, *) imax
          case ('jmin')
            read(val, *) jmin
          case ('jmax')
            read(val, *) jmax
          case ('x0')
            read(val, *) x0
          case ('y0')
            read(val, *) y0
          case ('r2d')
            read(val, *) r2d
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 1D array
  !-------------------------------------------------------------------
  subroutine read_array_1d(filename, arr, is, ie)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie
    real, intent(out) :: arr(is:ie)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_1d

  !-------------------------------------------------------------------
  ! Read 2D array
  !-------------------------------------------------------------------
  subroutine read_array_2d(filename, arr, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: arr(is:ie, js:je)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_2d

  !-------------------------------------------------------------------
  ! Validate output
  !-------------------------------------------------------------------
  subroutine validate_output(output, reference, is, ie, js, je, tol, max_err, err_count)
    implicit none
    integer, intent(in) :: is, ie, js, je
    real, intent(in) :: output(is:ie, js:je)
    real, intent(in) :: reference(is:ie, js:je)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j

    max_err = 0.0
    err_count = 0

    do j = js, je
      do i = is, ie
        if (abs(reference(i,j)) > 1.0e-30) then
          rel_err = abs(output(i,j) - reference(i,j)) / abs(reference(i,j))
        else
          rel_err = abs(output(i,j) - reference(i,j))
        end if
        if (rel_err > max_err) max_err = rel_err
        if (rel_err > tol) err_count = err_count + 1
      end do
    end do
  end subroutine validate_output

end program kernel_benchmark_xy2ll

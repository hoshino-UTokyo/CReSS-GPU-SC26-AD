!-----------------------------------------------------------------------
! Kernel Benchmark Program: bulksfc - GPU Version
! Source: Src/bulksfc.f90
! Subroutine: s_bulksfc
! Description: Calculate bulk coefficients of surface flux
!-----------------------------------------------------------------------
program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters (dimensions)
  integer :: ni, nj

  ! Physical constants from params.txt
  real :: kappa, wkappa, prnumg, prnumw
  real :: icz0m, icz0h, rmg, rms, rhg, rhs
  real :: oned3, oned9, oned27, tend3, cc

  ! Derived constants from params.txt
  real :: cc05, prgiv, prwiv, kp2, wkp2
  real :: kprg, wkprw, kprg3, wkprw3
  real :: icz0mv, icz0hv, rmg3v, rms3v
  real :: cqg1, cqg2, cpg1, cpg2
  real :: cqs1, cqs2, cps1, cps2

  ! Arrays
  real, allocatable :: za(:,:)
  integer, allocatable :: land(:,:)
  real, allocatable :: kai(:,:)
  real, allocatable :: z0m(:,:), z0h(:,:)
  real, allocatable :: rch(:,:)
  real, allocatable :: cm(:,:), ch(:,:)
  real, allocatable :: cm_ref(:,:), ch_ref(:,:)

  ! Benchmark variables
  character(len=256) :: data_dir
  integer :: num_iterations, warmup_iterations
  real :: tolerance
  double precision :: start_time, end_time, total_time, avg_time
  double precision :: min_time, max_time
  double precision, allocatable :: times(:)
  integer :: iter, errors
  real :: max_error_cm, max_error_ch

  ! Read benchmark configuration
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  ! Read parameters from dump
  call read_parameters(trim(data_dir)//'/params.txt')

  ! Allocate arrays
  allocate(za(0:ni+1, 0:nj+1))
  allocate(land(0:ni+1, 0:nj+1))
  allocate(kai(0:ni+1, 0:nj+1))
  allocate(z0m(0:ni+1, 0:nj+1))
  allocate(z0h(0:ni+1, 0:nj+1))
  allocate(rch(0:ni+1, 0:nj+1))
  allocate(cm(0:ni+1, 0:nj+1))
  allocate(ch(0:ni+1, 0:nj+1))
  allocate(cm_ref(0:ni+1, 0:nj+1))
  allocate(ch_ref(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  ! Read input arrays
  call read_array_2d(trim(data_dir)//'/za.bin', za, 0, ni+1, 0, nj+1)
  call read_array_2d_int(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/kai.bin', kai, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/z0m.bin', z0m, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/z0h.bin', z0h, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/rch.bin', rch, 0, ni+1, 0, nj+1)

  ! Read reference output for validation
  call read_array_2d(trim(data_dir)//'/cm_ref.bin', cm_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/ch_ref.bin', ch_ref, 0, ni+1, 0, nj+1)

  print '(A)',        '========================================'
  print '(A)',        'Kernel: bulksfc (s_bulksfc) - GPU'
  print '(A)',        '========================================'
  print '(A,I0,A,I0)', 'Grid size: ', ni, ' x ', nj
  print '(A,I0)',     'Warmup iterations: ', warmup_iterations
  print '(A,I0)',     'Benchmark iterations: ', num_iterations
  print '(A)',        '========================================'

  ! Warmup iterations
  print '(A)',        'Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_bulksfc(ni, nj, za, land, kai, z0m, z0h, rch, cm, ch)
    !$acc wait
  end do

  ! Benchmark iterations
  print '(A)',        'Running benchmark iterations...'
  total_time = 0.0d0
  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_bulksfc(ni, nj, za, land, kai, z0m, z0h, rch, cm, ch)
    !$acc wait
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / dble(num_iterations)
  min_time = minval(times)
  max_time = maxval(times)

  ! Validate results
  call validate_results(cm, cm_ref, ch, ch_ref, ni, nj, tolerance, errors, max_error_cm, max_error_ch)

  ! Output results
  print '(A)',        ''
  print '(A)',        '========================================'
  print '(A)',        'Results'
  print '(A)',        '========================================'
  print '(A,F12.6,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', 'Total time:   ', total_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', 'Min time:     ', min_time * 1000.0d0, ' ms'
  print '(A,F12.6,A)', 'Max time:     ', max_time * 1000.0d0, ' ms'
  print '(A)',        '----------------------------------------'
  print '(A,ES12.5)', 'Max error (cm): ', max_error_cm
  print '(A,ES12.5)', 'Max error (ch): ', max_error_ch
  print '(A,ES12.5)', 'Tolerance: ', tolerance
  print '(A,I0)',     'Validation errors: ', errors
  if (errors == 0) then
    print '(A)',      'Validation: PASSED'
  else
    print '(A)',      'Validation: FAILED'
  end if
  print '(A)',        '========================================'

  ! Cleanup
  deallocate(za, land, kai, z0m, z0h, rch, cm, ch, cm_ref, ch_ref)
  deallocate(times)

  if (errors /= 0) stop 1

contains

  !---------------------------------------------------------------------
  ! Kernel subroutine: bulksfc - GPU version
  !---------------------------------------------------------------------
  subroutine kernel_bulksfc(ni, nj, za, land, kai, z0m, z0h, rch, cm, ch)
    implicit none

    integer, intent(in) :: ni, nj
    real, intent(in) :: za(0:ni+1, 0:nj+1)
    integer, intent(in) :: land(0:ni+1, 0:nj+1)
    real, intent(in) :: kai(0:ni+1, 0:nj+1)
    real, intent(in) :: z0m(0:ni+1, 0:nj+1)
    real, intent(in) :: z0h(0:ni+1, 0:nj+1)
    real, intent(in) :: rch(0:ni+1, 0:nj+1)
    real, intent(out) :: cm(0:ni+1, 0:nj+1)
    real, intent(out) :: ch(0:ni+1, 0:nj+1)

    integer :: i, j
    real :: dz0m, dz0h, cmice, chice
    real :: a, b, c, d, e, f

    !$acc data copyin(za, land, kai, z0m, z0h, rch) copyout(cm, ch)
    !$acc kernels
    !$acc loop independent collapse(2) private(dz0m, dz0h, cmice, chice, a, b, c, d, e, f)
    do j = 1, nj-1
      do i = 1, ni-1

        ! Set common used variables
        dz0m = za(i,j) - z0m(i,j)
        dz0h = za(i,j) - z0h(i,j)

        a = max(za(i,j) / z0m(i,j), 1.01e0)
        b = max(za(i,j) / z0h(i,j), 1.01e0)

        ! For the unstable case
        if (rch(i,j) .lt. 0.e0) then

          if (land(i,j) .lt. 3) then
            a = log(a)
            b = log(b)

            c = rch(i,j) * prwiv
            c = c * c

            d = cqs1 + cqs2 * c
            e = cps1 + cps2 * c

            f = e * e - d * d * d

            c = (za(i,j) * dz0h * a * a) / (dz0m * dz0m * b)

            if (f .gt. 0.e0) then
              f = exp(oned3 * log(sqrt(f) + abs(e)))
              f = rms * c * (rms3v - (f + d/f))
            else
              f = sqrt(d)
              e = max(min(e / (d * f), 1.e0), -1.e0)
              f = rms * c * (rms3v - 2.e0 * f * cos(oned3 * acos(e)))
            end if

            f = sqrt(1.e0 - min(f, 1.e0))
            e = sqrt(f)

            c = 2.e0 * log(.5e0 * (1.e0 + e)) &
                + log(.5e0 * (1.e0 + f)) - 2.e0 * atan(e) + cc05
            d = 2.e0 * log(.5e0 * (1.e0 + f))

            if (c .lt. .5e0 * a) then
              cm(i,j) = wkappa / (a - c)
            else
              cm(i,j) = wkp2 / a
            end if

            if (d .lt. .7e0 * b) then
              ch(i,j) = wkprw / (b - d)
            else
              ch(i,j) = wkprw3 / b
            end if

          else
            a = log(a)
            b = log(b)

            c = rch(i,j) * prgiv
            c = c * c

            d = cqg1 + cqg2 * c
            e = cpg1 + cpg2 * c

            f = e * e - d * d * d

            c = (za(i,j) * dz0h * a * a) / (dz0m * dz0m * b)

            if (f .gt. 0.e0) then
              f = exp(oned3 * log(sqrt(f) + abs(e)))
              f = rmg * c * (rmg3v - (f + d/f))
            else
              f = sqrt(d)
              e = max(min(e / (d * f), 1.e0), -1.e0)
              f = rmg * c * (rmg3v - 2.e0 * f * cos(oned3 * acos(e)))
            end if

            f = min(f, 1.e0)
            e = sqrt(sqrt(1.e0 - f))

            c = 2.e0 * log(.5e0 * (1.e0 + e)) &
                + log(.5e0 * (1.e0 + e * e)) - 2.e0 * atan(e) + cc05
            f = sqrt(1.e0 - .6e0 * f)
            d = 2.e0 * log(.5e0 * (1.e0 + f))

            if (c .lt. .5e0 * a) then
              cm(i,j) = kappa / (a - c)
            else
              cm(i,j) = kp2 / a
            end if

            if (d .lt. .7e0 * b) then
              ch(i,j) = kprg / (b - d)
            else
              ch(i,j) = kprg3 / b
            end if
          end if

        ! For the stable case
        else
          if (land(i,j) .lt. 3) then
            a = wkappa / log(a)
            b = wkappa / log(b)

            c = sqrt(1.e0 + 5.e0 * rch(i,j))
            d = sqrt(1.e0 + 10.e0 * rch(i,j) * c)

            cm(i,j) = a / d
            ch(i,j) = b * d / (prnumw * (1.e0 + 15.e0 * rch(i,j) * c))
          else
            a = kappa / log(a)
            b = kappa / log(b)

            c = sqrt(1.e0 + 5.e0 * rch(i,j))
            d = sqrt(1.e0 + 10.e0 * rch(i,j) * c)

            cm(i,j) = a / d
            ch(i,j) = b * d / (prnumg * (1.e0 + 15.e0 * rch(i,j) * c))
          end if
        end if

        ! Mix bulk coefficients for weighted average ice surface
        if (land(i,j) .eq. 1) then
          dz0m = za(i,j) - icz0m
          dz0h = za(i,j) - icz0h

          a = max(za(i,j) * icz0mv, 1.01e0)
          b = max(za(i,j) * icz0hv, 1.01e0)

          if (rch(i,j) .lt. 0.e0) then
            a = log(a)
            b = log(b)

            c = rch(i,j) * prgiv
            c = c * c

            d = cqg1 + cqg2 * c
            e = cpg1 + cpg2 * c

            f = e * e - d * d * d

            c = (za(i,j) * dz0h * a * a) / (dz0m * dz0m * b)

            if (f .gt. 0.e0) then
              f = exp(oned3 * log(sqrt(f) + abs(e)))
              f = rmg * c * (rmg3v - (f + d/f))
            else
              f = sqrt(d)
              e = max(min(e / (d * f), 1.e0), -1.e0)
              f = rmg * c * (rmg3v - 2.e0 * f * cos(oned3 * acos(e)))
            end if

            f = min(f, 1.e0)
            e = sqrt(sqrt(1.e0 - f))

            c = 2.e0 * log(.5e0 * (1.e0 + e)) &
                + log(.5e0 * (1.e0 + e * e)) - 2.e0 * atan(e) + cc05
            f = sqrt(1.e0 - .6e0 * f)
            d = 2.e0 * log(.5e0 * (1.e0 + f))

            if (c .lt. .5e0 * a) then
              cmice = kappa / (a - c)
            else
              cmice = kp2 / a
            end if

            if (d .lt. .7e0 * b) then
              chice = kprg / (b - d)
            else
              chice = kprg3 / b
            end if
          else
            a = kappa / log(a)
            b = kappa / log(b)

            c = sqrt(1.e0 + 5.e0 * rch(i,j))
            d = sqrt(1.e0 + 10.e0 * rch(i,j) * c)

            cmice = a / d
            chice = b * d / (prnumg * (1.e0 + 15.e0 * rch(i,j) * c))
          end if

          a = 1.e0 - kai(i,j)
          cm(i,j) = kai(i,j) * cmice + a * cm(i,j)
          ch(i,j) = kai(i,j) * chice + a * ch(i,j)
        end if

      end do
    end do
    !$acc end kernels
    !$acc end data

  end subroutine kernel_bulksfc

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: unit_num, ios

    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open benchmark.conf'
      stop 1
    end if

    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iter
    read(unit_num, *) warmup_iter
    read(unit_num, *) tol
    close(unit_num)
  end subroutine read_config

  !---------------------------------------------------------------------
  ! Read parameters from dump file
  !---------------------------------------------------------------------
  subroutine read_parameters(filename)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=256) :: line, key, value_str
    integer :: unit_num, ios, eq_pos

    unit_num = 11
    open(unit=unit_num, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    do
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit

      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        value_str = adjustl(line(eq_pos+1:))

        select case (trim(key))
          case ('ni')
            read(value_str, *) ni
          case ('nj')
            read(value_str, *) nj
          case ('kappa')
            read(value_str, *) kappa
          case ('wkappa')
            read(value_str, *) wkappa
          case ('prnumg')
            read(value_str, *) prnumg
          case ('prnumw')
            read(value_str, *) prnumw
          case ('icz0m')
            read(value_str, *) icz0m
          case ('icz0h')
            read(value_str, *) icz0h
          case ('rmg')
            read(value_str, *) rmg
          case ('rms')
            read(value_str, *) rms
          case ('rhg')
            read(value_str, *) rhg
          case ('rhs')
            read(value_str, *) rhs
          case ('oned3')
            read(value_str, *) oned3
          case ('oned9')
            read(value_str, *) oned9
          case ('oned27')
            read(value_str, *) oned27
          case ('tend3')
            read(value_str, *) tend3
          case ('cc')
            read(value_str, *) cc
          case ('cc05')
            read(value_str, *) cc05
          case ('prgiv')
            read(value_str, *) prgiv
          case ('prwiv')
            read(value_str, *) prwiv
          case ('kp2')
            read(value_str, *) kp2
          case ('wkp2')
            read(value_str, *) wkp2
          case ('kprg')
            read(value_str, *) kprg
          case ('wkprw')
            read(value_str, *) wkprw
          case ('kprg3')
            read(value_str, *) kprg3
          case ('wkprw3')
            read(value_str, *) wkprw3
          case ('icz0mv')
            read(value_str, *) icz0mv
          case ('icz0hv')
            read(value_str, *) icz0hv
          case ('rmg3v')
            read(value_str, *) rmg3v
          case ('rms3v')
            read(value_str, *) rms3v
          case ('cqg1')
            read(value_str, *) cqg1
          case ('cqg2')
            read(value_str, *) cqg2
          case ('cpg1')
            read(value_str, *) cpg1
          case ('cpg2')
            read(value_str, *) cpg2
          case ('cqs1')
            read(value_str, *) cqs1
          case ('cqs2')
            read(value_str, *) cqs2
          case ('cps1')
            read(value_str, *) cps1
          case ('cps2')
            read(value_str, *) cps2
        end select
      end if
    end do

    close(unit_num)
  end subroutine read_parameters

  !---------------------------------------------------------------------
  ! Read 2D real array from binary file
  !---------------------------------------------------------------------
  subroutine read_array_2d(filename, array, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: array(is:ie, js:je)
    integer :: unit_num, ios

    unit_num = 12
    open(unit=unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(unit_num) array
    close(unit_num)
  end subroutine read_array_2d

  !---------------------------------------------------------------------
  ! Read 2D integer array from binary file
  !---------------------------------------------------------------------
  subroutine read_array_2d_int(filename, array, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    integer, intent(out) :: array(is:ie, js:je)
    integer :: unit_num, ios

    unit_num = 12
    open(unit=unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      print *, 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(unit_num) array
    close(unit_num)
  end subroutine read_array_2d_int

  !---------------------------------------------------------------------
  ! Validate results against reference
  !---------------------------------------------------------------------
  subroutine validate_results(cm, cm_ref, ch, ch_ref, ni, nj, tolerance, errors, max_err_cm, max_err_ch)
    implicit none
    integer, intent(in) :: ni, nj
    real, intent(in) :: cm(0:ni+1, 0:nj+1), cm_ref(0:ni+1, 0:nj+1)
    real, intent(in) :: ch(0:ni+1, 0:nj+1), ch_ref(0:ni+1, 0:nj+1)
    real, intent(in) :: tolerance
    integer, intent(out) :: errors
    real, intent(out) :: max_err_cm, max_err_ch

    integer :: i, j, max_i_cm, max_j_cm, max_i_ch, max_j_ch
    real :: rel_err, denom

    errors = 0
    max_err_cm = 0.0
    max_err_ch = 0.0
    max_i_cm = 0
    max_j_cm = 0
    max_i_ch = 0
    max_j_ch = 0

    do j = 1, nj-1
      do i = 1, ni-1
        ! Check cm
        denom = max(abs(cm_ref(i,j)), 1.0e-20)
        rel_err = abs(cm(i,j) - cm_ref(i,j)) / denom
        if (rel_err > max_err_cm) then
          max_err_cm = rel_err
          max_i_cm = i
          max_j_cm = j
        end if
        if (rel_err > tolerance) errors = errors + 1

        ! Check ch
        denom = max(abs(ch_ref(i,j)), 1.0e-20)
        rel_err = abs(ch(i,j) - ch_ref(i,j)) / denom
        if (rel_err > max_err_ch) then
          max_err_ch = rel_err
          max_i_ch = i
          max_j_ch = j
        end if
        if (rel_err > tolerance) errors = errors + 1
      end do
    end do

    print '(A,I0,A,I0)', 'Max cm error at (i,j) = (', max_i_cm, ',', max_j_cm, ')'
    print '(A,ES15.8,A,ES15.8)', 'cm GPU=', cm(max_i_cm, max_j_cm), '  ref=', cm_ref(max_i_cm, max_j_cm)
    print '(A,I0,A,I0)', 'Max ch error at (i,j) = (', max_i_ch, ',', max_j_ch, ')'
    print '(A,ES15.8,A,ES15.8)', 'ch GPU=', ch(max_i_ch, max_j_ch), '  ref=', ch_ref(max_i_ch, max_j_ch)

    ! Count zeros
    print '(A,I0)', 'cm zeros: ', count(cm(1:ni-1,1:nj-1) == 0.0)
    print '(A,I0)', 'ch zeros: ', count(ch(1:ni-1,1:nj-1) == 0.0)
    print '(A,I0)', 'cm_ref zeros: ', count(cm_ref(1:ni-1,1:nj-1) == 0.0)
    print '(A,ES15.8,A,ES15.8)', 'cm(10,10)=', cm(10,10), '  ref=', cm_ref(10,10)
    print '(A,ES15.8,A,ES15.8)', 'cm(100,100)=', cm(100,100), '  ref=', cm_ref(100,100)
  end subroutine validate_results

end program kernel_benchmark

program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  integer :: wbc, ebc, sbc, nbc
  integer :: npe, siz, nigrp, njgrp, nigrp1, njgrp1, nsrl

  ! Arrays
  integer, allocatable :: idxbuf(:,:), idxbuf_ref(:,:)
  real, allocatable :: mxnbuf(:), mxnbuf_ref(:)
  real, allocatable :: sbuf(:), sbuf_ref(:)
  real, allocatable :: rbuf(:), rbuf_ref(:)
  integer, allocatable :: grpxy(:,:), grpxy_ref(:,:)
  integer, allocatable :: xgrp(:), xgrp_ref(:)
  integer, allocatable :: ygrp(:), ygrp_ref(:)

  ! Benchmark variables
  integer :: num_iterations, warmup_iterations
  real(8) :: total_time, avg_time
  real(8) :: start_time, end_time
  character(len=256) :: data_dir
  real :: tolerance
  integer :: error_count
  integer :: iter, ijpe, igc_sub, jgc_sub, ijsc

  ! Read benchmark configuration
  open(10, file='benchmark.conf', status='old')
  read(10, '(A)') data_dir
  read(10, *) num_iterations
  read(10, *) warmup_iterations
  read(10, *) tolerance
  close(10)

  ! Read parameters
  call read_params(trim(data_dir)//'/params.txt')

  ! Allocate arrays
  allocate(idxbuf(1:3, 0:npe-1), idxbuf_ref(1:3, 0:npe-1))
  allocate(mxnbuf(0:npe-1), mxnbuf_ref(0:npe-1))
  allocate(sbuf(1:siz), sbuf_ref(1:siz))
  allocate(rbuf(1:siz), rbuf_ref(1:siz))
  allocate(grpxy(0:nigrp+1, 0:njgrp+1), grpxy_ref(0:nigrp+1, 0:njgrp+1))
  allocate(xgrp(0:nsrl-1), xgrp_ref(0:nsrl-1))
  allocate(ygrp(0:nsrl-1), ygrp_ref(0:nsrl-1))

  ! Read reference output
  call read_binary_2d_i(trim(data_dir)//'/idxbuf_ref.bin', idxbuf_ref, 1, 3, 0, npe-1)
  call read_binary_1d(trim(data_dir)//'/mxnbuf_ref.bin', mxnbuf_ref, 0, npe-1)
  call read_binary_1d(trim(data_dir)//'/sbuf_ref.bin', sbuf_ref, 1, siz)
  call read_binary_1d(trim(data_dir)//'/rbuf_ref.bin', rbuf_ref, 1, siz)
  call read_binary_2d_i(trim(data_dir)//'/grpxy_ref.bin', grpxy_ref, 0, nigrp+1, 0, njgrp+1)
  call read_binary_1d_i(trim(data_dir)//'/xgrp_ref.bin', xgrp_ref, 0, nsrl-1)
  call read_binary_1d_i(trim(data_dir)//'/ygrp_ref.bin', ygrp_ref, 0, nsrl-1)

  ! Warmup
  do iter = 1, warmup_iterations
    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(ijpe)
    do ijpe = 0, npe-1
      idxbuf(1,ijpe) = 0
      idxbuf(2,ijpe) = 0
      idxbuf(3,ijpe) = 0
      mxnbuf(ijpe) = 0.e0
    end do
    !$omp end do
    !$omp do schedule(runtime) private(ijpe)
    do ijpe = 1, siz
      sbuf(ijpe) = 0.e0
      rbuf(ijpe) = 0.e0
    end do
    !$omp end do
    !$omp do schedule(runtime) private(igc_sub,jgc_sub)
    do jgc_sub = 0, njgrp+1
      do igc_sub = 0, nigrp+1
        grpxy(igc_sub,jgc_sub) = -1
      end do
    end do
    !$omp end do
    if (abs(wbc).eq.1 .and. abs(ebc).eq.1) then
      !$omp do schedule(runtime) private(jgc_sub)
      do jgc_sub = 0, njgrp+1
        grpxy(0,jgc_sub) = nsrl
        grpxy(nigrp1,jgc_sub) = nsrl
      end do
      !$omp end do
    end if
    if (abs(sbc).eq.1 .and. abs(nbc).eq.1) then
      !$omp do schedule(runtime) private(igc_sub)
      do igc_sub = 0, nigrp+1
        grpxy(igc_sub,0) = nsrl
        grpxy(igc_sub,njgrp1) = nsrl
      end do
      !$omp end do
    end if
    !$omp do schedule(runtime) private(ijsc)
    do ijsc = 0, nsrl-1
      xgrp(ijsc) = -1
      ygrp(ijsc) = -1
    end do
    !$omp end do
    !$omp end parallel
  end do

  ! Benchmark
  total_time = 0.0d0
  do iter = 1, num_iterations
    start_time = omp_get_wtime()

    !$omp parallel default(shared)
    !$omp do schedule(runtime) private(ijpe)
    do ijpe = 0, npe-1
      idxbuf(1,ijpe) = 0
      idxbuf(2,ijpe) = 0
      idxbuf(3,ijpe) = 0
      mxnbuf(ijpe) = 0.e0
    end do
    !$omp end do
    !$omp do schedule(runtime) private(ijpe)
    do ijpe = 1, siz
      sbuf(ijpe) = 0.e0
      rbuf(ijpe) = 0.e0
    end do
    !$omp end do
    !$omp do schedule(runtime) private(igc_sub,jgc_sub)
    do jgc_sub = 0, njgrp+1
      do igc_sub = 0, nigrp+1
        grpxy(igc_sub,jgc_sub) = -1
      end do
    end do
    !$omp end do
    if (abs(wbc).eq.1 .and. abs(ebc).eq.1) then
      !$omp do schedule(runtime) private(jgc_sub)
      do jgc_sub = 0, njgrp+1
        grpxy(0,jgc_sub) = nsrl
        grpxy(nigrp1,jgc_sub) = nsrl
      end do
      !$omp end do
    end if
    if (abs(sbc).eq.1 .and. abs(nbc).eq.1) then
      !$omp do schedule(runtime) private(igc_sub)
      do igc_sub = 0, nigrp+1
        grpxy(igc_sub,0) = nsrl
        grpxy(igc_sub,njgrp1) = nsrl
      end do
      !$omp end do
    end if
    !$omp do schedule(runtime) private(ijsc)
    do ijsc = 0, nsrl-1
      xgrp(ijsc) = -1
      ygrp(ijsc) = -1
    end do
    !$omp end do
    !$omp end parallel

    end_time = omp_get_wtime()
    total_time = total_time + (end_time - start_time)
  end do

  avg_time = total_time / dble(num_iterations)

  ! Validate
  error_count = 0
  do ijpe = 0, npe-1
    if (idxbuf(1,ijpe) /= idxbuf_ref(1,ijpe)) error_count = error_count + 1
    if (idxbuf(2,ijpe) /= idxbuf_ref(2,ijpe)) error_count = error_count + 1
    if (idxbuf(3,ijpe) /= idxbuf_ref(3,ijpe)) error_count = error_count + 1
    if (abs(mxnbuf(ijpe) - mxnbuf_ref(ijpe)) > tolerance) error_count = error_count + 1
  end do
  do ijpe = 1, siz
    if (abs(sbuf(ijpe) - sbuf_ref(ijpe)) > tolerance) error_count = error_count + 1
    if (abs(rbuf(ijpe) - rbuf_ref(ijpe)) > tolerance) error_count = error_count + 1
  end do
  do jgc_sub = 0, njgrp+1
    do igc_sub = 0, nigrp+1
      if (grpxy(igc_sub,jgc_sub) /= grpxy_ref(igc_sub,jgc_sub)) error_count = error_count + 1
    end do
  end do
  do ijsc = 0, nsrl-1
    if (xgrp(ijsc) /= xgrp_ref(ijsc)) error_count = error_count + 1
    if (ygrp(ijsc) /= ygrp_ref(ijsc)) error_count = error_count + 1
  end do

  print '(A)', '=== Benchmark Results ==='
  print '(A,I0)', 'Kernel: allocbuf_s_allocbuf'
  print '(A,I0)', 'Iterations: ', num_iterations
  print '(A,F12.6,A)', 'Average time: ', avg_time * 1000.0d0, ' ms'
  print '(A,I0)', 'Errors: ', error_count
  if (error_count == 0) then
    print '(A)', 'PASSED'
  else
    print '(A)', 'FAILED'
  end if

  deallocate(idxbuf, idxbuf_ref, mxnbuf, mxnbuf_ref, sbuf, sbuf_ref, rbuf, rbuf_ref)
  deallocate(grpxy, grpxy_ref, xgrp, xgrp_ref, ygrp, ygrp_ref)

contains

  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    character(len=256) :: line
    character(len=64) :: key
    integer :: ios, eq_pos

    open(10, file=filename, status='old', iostat=ios)
    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        select case (trim(key))
        case ('wbc'); read(line(eq_pos+1:), *) wbc
        case ('ebc'); read(line(eq_pos+1:), *) ebc
        case ('sbc'); read(line(eq_pos+1:), *) sbc
        case ('nbc'); read(line(eq_pos+1:), *) nbc
        case ('npe'); read(line(eq_pos+1:), *) npe
        case ('siz'); read(line(eq_pos+1:), *) siz
        case ('nigrp'); read(line(eq_pos+1:), *) nigrp
        case ('njgrp'); read(line(eq_pos+1:), *) njgrp
        case ('nigrp1'); read(line(eq_pos+1:), *) nigrp1
        case ('njgrp1'); read(line(eq_pos+1:), *) njgrp1
        case ('nsrl'); read(line(eq_pos+1:), *) nsrl
        end select
      end if
    end do
    close(10)
  end subroutine

  subroutine read_binary_1d(filename, array, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: k1, k2
    real, intent(out) :: array(k1:k2)
    open(10, file=filename, form='unformatted', access='stream', status='old')
    read(10) array
    close(10)
  end subroutine

  subroutine read_binary_1d_i(filename, array, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: k1, k2
    integer, intent(out) :: array(k1:k2)
    open(10, file=filename, form='unformatted', access='stream', status='old')
    read(10) array
    close(10)
  end subroutine

  subroutine read_binary_2d_i(filename, array, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    integer, intent(out) :: array(i1:i2, j1:j2)
    open(10, file=filename, form='unformatted', access='stream', status='old')
    read(10) array
    close(10)
  end subroutine

end program kernel_benchmark

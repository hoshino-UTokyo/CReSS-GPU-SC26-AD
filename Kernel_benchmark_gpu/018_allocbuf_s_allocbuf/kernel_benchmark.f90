program kernel_benchmark
  use omp_lib
  implicit none

  ! Kind parameters
  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  ! Parameters
  integer :: wbc, ebc, sbc, nbc
  integer :: npe, siz, nigrp, njgrp, nigrp1, njgrp1, nsrl

  ! Arrays
  integer, allocatable :: idxbuf(:,:)
  real(sp), allocatable :: mxnbuf(:)
  real(sp), allocatable :: sbuf(:)
  real(sp), allocatable :: rbuf(:)
  integer, allocatable :: grpxy(:,:)
  integer, allocatable :: xgrp(:)
  integer, allocatable :: ygrp(:)

  ! Benchmark parameters
  character(len=256) :: data_dir, line
  integer :: num_iterations, warmup_iterations
  real(dp) :: tolerance
  integer :: io_unit, ios

  ! Timing variables
  real(dp) :: start_time, end_time
  real(dp) :: total_time, avg_time
  real(dp), allocatable :: times(:)

  ! Loop variables
  integer :: iter, ijpe, igc_sub, jgc_sub, ijsc

  ! Read benchmark configuration
  open(newunit=io_unit, file='benchmark.conf', status='old', action='read')
  read(io_unit, '(A)') data_dir
  read(io_unit, *) num_iterations
  read(io_unit, *) warmup_iterations
  read(io_unit, *) tolerance
  close(io_unit)

  data_dir = trim(adjustl(data_dir))

  ! Read parameters
  call read_params(trim(data_dir)//'/params.txt')

  write(*,'(A)') '======================================'
  write(*,'(A)') 'GPU Kernel Benchmark: allocbuf'
  write(*,'(A)') '  Buffer allocation & initialization'
  write(*,'(A)') '======================================'
  write(*,'(A,I6)') 'npe = ', npe
  write(*,'(A,I6)') 'siz = ', siz
  write(*,'(A,I6)') 'nigrp = ', nigrp
  write(*,'(A,I6)') 'njgrp = ', njgrp
  write(*,'(A,I6)') 'nsrl = ', nsrl
  write(*,'(A,I6)') 'Iterations = ', num_iterations
  write(*,'(A,I6)') 'Warmup = ', warmup_iterations
  write(*,'(A)') '======================================'

  ! Allocate arrays
  allocate(idxbuf(1:3, 0:npe-1))
  allocate(mxnbuf(0:npe-1))
  allocate(sbuf(1:siz))
  allocate(rbuf(1:siz))
  allocate(grpxy(0:nigrp+1, 0:njgrp+1))
  allocate(xgrp(0:nsrl-1))
  allocate(ygrp(0:nsrl-1))
  allocate(times(num_iterations))

  ! Warmup iterations
  write(*,'(A)') 'Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_allocbuf()
    !$acc wait
  end do

  ! Benchmark iterations
  write(*,'(A)') 'Running benchmark iterations...'
  total_time = 0.0_dp

  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()

    call kernel_allocbuf()

    !$acc wait
    end_time = omp_get_wtime()
    times(iter) = end_time - start_time
    total_time = total_time + times(iter)
  end do

  avg_time = total_time / num_iterations

  ! Report results
  write(*,'(A)') '======================================'
  write(*,'(A)') 'Results:'
  write(*,'(A)') '======================================'
  write(*,'(A)') ''
  write(*,'(A)') 'Timing Results:'
  write(*,'(A,F12.6,A)') 'Total time: ', total_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Min time: ', minval(times) * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Max time: ', maxval(times) * 1000.0_dp, ' ms'

  ! Validation - check expected values
  write(*,'(A)') ''
  write(*,'(A)') 'Validation:'
  call validate_results()

  ! Cleanup
  deallocate(idxbuf, mxnbuf, sbuf, rbuf, grpxy, xgrp, ygrp, times)

contains

  subroutine kernel_allocbuf()
    ! Loop 1: Initialize idxbuf and mxnbuf
    !$acc kernels
    !$acc loop independent
    do ijpe = 0, npe-1
      idxbuf(1,ijpe) = 0
      idxbuf(2,ijpe) = 0
      idxbuf(3,ijpe) = 0
      mxnbuf(ijpe) = 0.e0
    end do
    !$acc end kernels

    ! Loop 2: Initialize sbuf and rbuf
    !$acc kernels
    !$acc loop independent
    do ijpe = 1, siz
      sbuf(ijpe) = 0.e0
      rbuf(ijpe) = 0.e0
    end do
    !$acc end kernels

    ! Loop 3: Initialize grpxy
    !$acc kernels
    !$acc loop independent
    do jgc_sub = 0, njgrp+1
      !$acc loop independent
      do igc_sub = 0, nigrp+1
        grpxy(igc_sub,jgc_sub) = -1
      end do
    end do
    !$acc end kernels

    ! Loop 4: Set periodic boundary for x-direction
    if (abs(wbc).eq.1 .and. abs(ebc).eq.1) then
      !$acc kernels
      !$acc loop independent
      do jgc_sub = 0, njgrp+1
        grpxy(0,jgc_sub) = nsrl
        grpxy(nigrp1,jgc_sub) = nsrl
      end do
      !$acc end kernels
    end if

    ! Loop 5: Set periodic boundary for y-direction
    if (abs(sbc).eq.1 .and. abs(nbc).eq.1) then
      !$acc kernels
      !$acc loop independent
      do igc_sub = 0, nigrp+1
        grpxy(igc_sub,0) = nsrl
        grpxy(igc_sub,njgrp1) = nsrl
      end do
      !$acc end kernels
    end if

    ! Loop 6: Initialize xgrp and ygrp
    !$acc kernels
    !$acc loop independent
    do ijsc = 0, nsrl-1
      xgrp(ijsc) = -1
      ygrp(ijsc) = -1
    end do
    !$acc end kernels
  end subroutine kernel_allocbuf

  subroutine validate_results()
    integer :: error_count
    error_count = 0

    ! Check idxbuf
    do ijpe = 0, npe-1
      if (idxbuf(1,ijpe) /= 0) error_count = error_count + 1
      if (idxbuf(2,ijpe) /= 0) error_count = error_count + 1
      if (idxbuf(3,ijpe) /= 0) error_count = error_count + 1
    end do

    ! Check mxnbuf
    do ijpe = 0, npe-1
      if (mxnbuf(ijpe) /= 0.e0) error_count = error_count + 1
    end do

    ! Check sbuf and rbuf
    do ijpe = 1, siz
      if (sbuf(ijpe) /= 0.e0) error_count = error_count + 1
      if (rbuf(ijpe) /= 0.e0) error_count = error_count + 1
    end do

    ! Check xgrp and ygrp
    do ijsc = 0, nsrl-1
      if (xgrp(ijsc) /= -1) error_count = error_count + 1
      if (ygrp(ijsc) /= -1) error_count = error_count + 1
    end do

    write(*,'(A,I0)') 'Errors: ', error_count
    if (error_count == 0) then
      write(*,'(A)') '======================================'
      write(*,'(A)') 'Validation: PASSED'
    else
      write(*,'(A)') '======================================'
      write(*,'(A)') 'Validation: FAILED'
    end if
    write(*,'(A)') '======================================'
  end subroutine validate_results

  subroutine read_params(filename)
    character(len=*), intent(in) :: filename
    character(len=256) :: pline
    character(len=64) :: key
    integer :: pios, eq_pos

    open(newunit=io_unit, file=filename, status='old', iostat=pios)
    if (pios /= 0) then
      print *, 'Error opening params file: ', trim(filename)
      stop 1
    end if

    do
      read(io_unit, '(A)', iostat=pios) pline
      if (pios /= 0) exit
      eq_pos = index(pline, '=')
      if (eq_pos > 0) then
        key = adjustl(pline(1:eq_pos-1))
        select case (trim(key))
        case ('wbc'); read(pline(eq_pos+1:), *) wbc
        case ('ebc'); read(pline(eq_pos+1:), *) ebc
        case ('sbc'); read(pline(eq_pos+1:), *) sbc
        case ('nbc'); read(pline(eq_pos+1:), *) nbc
        case ('npe'); read(pline(eq_pos+1:), *) npe
        case ('siz'); read(pline(eq_pos+1:), *) siz
        case ('nigrp'); read(pline(eq_pos+1:), *) nigrp
        case ('njgrp'); read(pline(eq_pos+1:), *) njgrp
        case ('nigrp1'); read(pline(eq_pos+1:), *) nigrp1
        case ('njgrp1'); read(pline(eq_pos+1:), *) njgrp1
        case ('nsrl'); read(pline(eq_pos+1:), *) nsrl
        end select
      end if
    end do
    close(io_unit)
  end subroutine read_params

end program kernel_benchmark

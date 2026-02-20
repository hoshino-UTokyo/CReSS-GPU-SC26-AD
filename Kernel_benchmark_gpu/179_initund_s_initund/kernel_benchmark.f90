!***********************************************************************
! GPU Kernel Benchmark: initund (s_initund)
!***********************************************************************
!
! Source: Src/initund.f90
! Description: Initialize soil and sea temperature arrays (tund, tundp) based on
!              surface type (land/sea), SST data, and atmospheric conditions
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_initund
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk, nund

  ! Parameters
  character(len=108) :: sfcdat
  integer :: sfcopt, advopt

  ! Physical constants
  real, parameter :: rd = 287.04e0    ! Gas constant for dry air
  real, parameter :: cp = 1004.0e0    ! Specific heat at constant pressure
  real, parameter :: p0 = 100000.0e0  ! Reference pressure
  real :: t0                          ! Freezing point temperature

  ! Derived constants
  real :: rddvcp, p0iv, dzgrd, tgdeep, sstcst
  real :: enk, enkm1v

  ! Arrays
  integer, allocatable :: land(:,:)
  real, allocatable :: pbr(:,:,:), ptbr(:,:,:), pp(:,:,:), ptp(:,:,:)
  real, allocatable :: sst(:,:), sst_init(:,:)
  real, allocatable :: tund(:,:,:), tundp(:,:,:)
  real, allocatable :: tund_ref(:,:,:), tundp_ref(:,:,:), sst_ref(:,:)
  real, allocatable :: ek(:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error
  integer :: error_count
  real :: tolerance
  logical :: validation_passed

  ! Loop variables
  integer :: iter

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt')

  ! Compute derived constants
  rddvcp = rd / cp
  p0iv = 1.0e0 / p0
  enk = exp(real(1-nund) * dzgrd)
  enkm1v = 1.0e0 / (exp(real(1-nund) * dzgrd) - 1.0e0)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: initund'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk, ', nund=', nund
  write(*,'(A,A)') ' sfcdat=', trim(sfcdat)
  write(*,'(A,I3,A,I3)') ' sfcopt=', sfcopt, ', advopt=', advopt
  write(*,'(A,F10.4)') ' dzgrd=', dzgrd
  write(*,'(A,F10.2,A,F10.2)') ' tgdeep=', tgdeep, ', sstcst=', sstcst
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(land(0:ni+1, 0:nj+1))
  allocate(pbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(pp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(sst(0:ni+1, 0:nj+1))
  allocate(sst_init(0:ni+1, 0:nj+1))
  allocate(tund(0:ni+1, 0:nj+1, 1:nund))
  allocate(tundp(0:ni+1, 0:nj+1, 1:nund))
  allocate(tund_ref(0:ni+1, 0:nj+1, 1:nund))
  allocate(tundp_ref(0:ni+1, 0:nj+1, 1:nund))
  allocate(sst_ref(0:ni+1, 0:nj+1))
  allocate(ek(1:nund))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_2d_int(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/pbr.bin', pbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pp.bin', pp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/sst_in.bin', sst_init, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/tund_ref.bin', tund_ref, 0, ni+1, 0, nj+1, 1, nund)
  call read_array_3d(trim(data_dir)//'/tundp_ref.bin', tundp_ref, 0, ni+1, 0, nj+1, 1, nund)
  call read_array_2d(trim(data_dir)//'/sst_ref.bin', sst_ref, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    tund = 0.0e0
    tundp = 0.0e0
    sst = sst_init
    ek = 0.0e0
    call kernel_initund(sfcdat, sfcopt, advopt, ni, nj, nk, nund, &
         rddvcp, p0iv, dzgrd, tgdeep, sstcst, enk, enkm1v, t0, &
         land, pbr, ptbr, pp, ptp, sst, tund, tundp, ek)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    tund = 0.0e0
    tundp = 0.0e0
    sst = sst_init
    ek = 0.0e0

    !$acc wait
    t_start = omp_get_wtime()

    call kernel_initund(sfcdat, sfcopt, advopt, ni, nj, nk, nund, &
         rddvcp, p0iv, dzgrd, tgdeep, sstcst, enk, enkm1v, t0, &
         land, pbr, ptbr, pp, ptp, sst, tund, tundp, ek)

    !$acc wait
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
  call validate_output_3d(tundp, tundp_ref, ni, nj, nund, tolerance, max_error, error_count)

  validation_passed = (error_count == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:  ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:      ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:      ', t_max * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Total time:    ', t_total * 1000.0d0, ' ms'
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
  deallocate(land, pbr, ptbr, pp, ptp, sst, sst_init)
  deallocate(tund, tundp, tund_ref, tundp_ref, sst_ref, ek, times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: initund (OpenACC version)
  !-------------------------------------------------------------------
  subroutine kernel_initund(sfcdat, sfcopt, advopt, ni, nj, nk, nund, &
       rddvcp, p0iv, dzgrd, tgdeep, sstcst, enk, enkm1v, t0, &
       land, pbr, ptbr, pp, ptp, sst, tund, tundp, ek)
    implicit none

    character(len=*), intent(in) :: sfcdat
    integer, intent(in) :: sfcopt, advopt
    integer, intent(in) :: ni, nj, nk, nund
    real, intent(in) :: rddvcp, p0iv, dzgrd, tgdeep, sstcst
    real, intent(in) :: enk, enkm1v, t0
    integer, intent(in) :: land(0:ni+1, 0:nj+1)
    real, intent(in) :: pbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: pp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ptp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: sst(0:ni+1, 0:nj+1)
    real, intent(out) :: tund(0:ni+1, 0:nj+1, 1:nund)
    real, intent(out) :: tundp(0:ni+1, 0:nj+1, 1:nund)
    real, intent(inout) :: ek(1:nund)

    integer :: i, j, k

    ! Initialized by diagnostic value
    if (sfcopt == 1 .or. sfcopt == 2 .or. sfcopt == 3) then

      ! Set the surface temperature
      if (sfcdat(2:2) == 'o') then

        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            if (land(i,j) < 3) then
              tundp(i,j,1) = sst(i,j)
            else
              tundp(i,j,1) = (ptbr(i,j,1) + ptp(i,j,1)) &
                   * exp(rddvcp * log(p0iv * (pbr(i,j,1) + pp(i,j,1))))
              tundp(i,j,1) = 0.5e0 * (tundp(i,j,1) + (ptbr(i,j,2) + ptp(i,j,2)) &
                   * exp(rddvcp * log(p0iv * (pbr(i,j,2) + pp(i,j,2)))))
              if (land(i,j) < 10) then
                tundp(i,j,1) = min(tundp(i,j,1), t0)
              end if
            end if
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            if (land(i,j) < 3) then
              tundp(i,j,1) = sstcst
            else
              tundp(i,j,1) = (ptbr(i,j,1) + ptp(i,j,1)) &
                   * exp(rddvcp * log(p0iv * (pbr(i,j,1) + pp(i,j,1))))
              tundp(i,j,1) = 0.5e0 * (tundp(i,j,1) + (ptbr(i,j,2) + ptp(i,j,2)) &
                   * exp(rddvcp * log(p0iv * (pbr(i,j,2) + pp(i,j,2)))))
              if (land(i,j) < 10) then
                tundp(i,j,1) = min(tundp(i,j,1), t0)
              end if
            end if
          end do
        end do
        !$acc end kernels

      end if

      ! Set the soil temperature - compute ek array
      !$acc kernels
      !$acc loop independent
      do k = 2, nund
        ek(k) = exp(real(1-k) * dzgrd)
      end do
      !$acc end kernels

      do k = 2, nund
        !$acc kernels
        !$acc loop independent
        do j = 1, nj-1
          !$acc loop independent
          do i = 1, ni-1
            if (land(i,j) < 10) then
              tundp(i,j,k) = tundp(i,j,1)
            else
              tundp(i,j,k) = ((tgdeep - tundp(i,j,1)) * enkm1v) * ek(k) &
                   + (tundp(i,j,1) * enk - tgdeep) * enkm1v
            end if
          end do
        end do
        !$acc end kernels
      end do

      ! Copy the past value to the present
      if (advopt <= 3) then
        do k = 1, nund
          !$acc kernels
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni-1
              tund(i,j,k) = tundp(i,j,k)
            end do
          end do
          !$acc end kernels
        end do
      end if

    else if (sfcopt > 10) then

      ! Reset the sea temperature
      if (advopt <= 3) then
        if (sfcdat(2:2) == 'o') then
          do k = 1, nund
            !$acc kernels
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                if (land(i,j) < 3) then
                  tund(i,j,k) = sst(i,j)
                  tundp(i,j,k) = sst(i,j)
                end if
              end do
            end do
            !$acc end kernels
          end do
        else
          do k = 1, nund
            !$acc kernels
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                if (land(i,j) < 3) then
                  tund(i,j,k) = sstcst
                  tundp(i,j,k) = sstcst
                end if
              end do
            end do
            !$acc end kernels
          end do
        end if
      else
        if (sfcdat(2:2) == 'o') then
          do k = 1, nund
            !$acc kernels
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                if (land(i,j) < 3) then
                  tundp(i,j,k) = sst(i,j)
                end if
              end do
            end do
            !$acc end kernels
          end do
        else
          do k = 1, nund
            !$acc kernels
            !$acc loop independent
            do j = 1, nj-1
              !$acc loop independent
              do i = 1, ni-1
                if (land(i,j) < 3) then
                  tundp(i,j,k) = sstcst
                end if
              end do
            end do
            !$acc end kernels
          end do
        end if
      end if

    end if

  end subroutine kernel_initund

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
          case ('sfcdat')
            read(val, '(A)') sfcdat
          case ('sfcopt')
            read(val, *) sfcopt
          case ('advopt')
            read(val, *) advopt
          case ('dzgrd')
            read(val, *) dzgrd
          case ('tgdeep')
            read(val, *) tgdeep
          case ('sstcst')
            read(val, *) sstcst
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('nund')
            read(val, *) nund
          case ('t0')
            read(val, *) t0
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 2D integer array
  !-------------------------------------------------------------------
  subroutine read_array_2d_int(filename, arr, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    integer, intent(out) :: arr(is:ie, js:je)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_2d_int

  !-------------------------------------------------------------------
  ! Read 2D real array
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
  ! Read 3D array
  !-------------------------------------------------------------------
  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_3d

  !-------------------------------------------------------------------
  ! Validate 3D output
  !-------------------------------------------------------------------
  subroutine validate_output_3d(output, reference, ni, nj, nund, tol, max_err, err_count)
    implicit none
    integer, intent(in) :: ni, nj, nund
    real, intent(in) :: output(0:ni+1, 0:nj+1, 1:nund)
    real, intent(in) :: reference(0:ni+1, 0:nj+1, 1:nund)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = 1, nund
      do j = 0, nj+1
        do i = 0, ni+1
          if (abs(reference(i,j,k)) > 1.0e-30) then
            rel_err = abs(output(i,j,k) - reference(i,j,k)) / abs(reference(i,j,k))
          else
            rel_err = abs(output(i,j,k) - reference(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do
  end subroutine validate_output_3d

end program kernel_benchmark_gpu_initund

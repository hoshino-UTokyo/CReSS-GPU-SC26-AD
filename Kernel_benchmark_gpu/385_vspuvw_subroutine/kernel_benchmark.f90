!***********************************************************************
! GPU Kernel Benchmark: vspuvw (s_vspuvw)
!***********************************************************************
!
! Source: Src/vspuvw.f90
! Description: Apply vertical sponge damping for velocity components
!              (u, v, w), relaxing toward GPV data or base state value.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_vspuvw_gpu
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: ksp0(1:2)
  integer :: vspopt
  real :: gtinc
  character(len=108) :: gpvvar, vspvar

  ! Arrays
  real, allocatable :: ubr(:,:,:), vbr(:,:,:)
  real, allocatable :: rst8u(:,:,:), rst8v(:,:,:), rst8w(:,:,:)
  real, allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:)
  real, allocatable :: rbct(:,:,:,:)
  real, allocatable :: ugpv(:,:,:), utd(:,:,:)
  real, allocatable :: vgpv(:,:,:), vtd(:,:,:)
  real, allocatable :: wgpv(:,:,:), wtd(:,:,:)
  real, allocatable :: ufrc(:,:,:), vfrc(:,:,:), wfrc(:,:,:)
  real, allocatable :: ufrc_init(:,:,:), vfrc_init(:,:,:), wfrc_init(:,:,:)
  real, allocatable :: ufrc_ref(:,:,:), vfrc_ref(:,:,:), wfrc_ref(:,:,:)
  real, allocatable :: rbct8s(:,:,:), rbct8s_init(:,:,:), rbct8s_ref(:,:,:)

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

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt')

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: vspuvw'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' vspopt=', vspopt
  write(*,'(A,ES12.4)') ' gtinc=', gtinc
  write(*,'(A,A)') ' gpvvar=', trim(gpvvar)
  write(*,'(A,A)') ' vspvar=', trim(vspvar)
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(ubr(0:ni+1, 0:nj+1, 1:nk))
  allocate(vbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(up(0:ni+1, 0:nj+1, 1:nk))
  allocate(vp(0:ni+1, 0:nj+1, 1:nk))
  allocate(wp(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbct(1:ni, 1:nj, 1:nk, 1:2))
  allocate(ugpv(0:ni+1, 0:nj+1, 1:nk))
  allocate(utd(0:ni+1, 0:nj+1, 1:nk))
  allocate(vgpv(0:ni+1, 0:nj+1, 1:nk))
  allocate(vtd(0:ni+1, 0:nj+1, 1:nk))
  allocate(wgpv(0:ni+1, 0:nj+1, 1:nk))
  allocate(wtd(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(wfrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbct8s(1:ni, 1:nj, 1:nk))
  allocate(rbct8s_init(1:ni, 1:nj, 1:nk))
  allocate(rbct8s_ref(1:ni, 1:nj, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call try_read_array_1d_int(trim(data_dir)//'/ksp0.bin', ksp0, 1, 2)
  call read_array_3d(trim(data_dir)//'/ubr.bin', ubr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vbr.bin', vbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8u.bin', rst8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8v.bin', rst8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8w.bin', rst8w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/up.bin', up, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vp.bin', vp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wp.bin', wp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_4d(trim(data_dir)//'/rbct.bin', rbct, 1, ni, 1, nj, 1, nk, 1, 2)
  call read_array_3d(trim(data_dir)//'/ugpv.bin', ugpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/utd.bin', utd, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vgpv.bin', vgpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vtd.bin', vtd, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wgpv.bin', wgpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wtd.bin', wtd, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ufrc_in.bin', ufrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_in.bin', vfrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wfrc_in.bin', wfrc_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d_v2(trim(data_dir)//'/rbct8s_in.bin', rbct8s_init, 1, ni, 1, nj, 1, nk)

  write(*,'(A,I6,A,I6)') ' ksp0(1)=', ksp0(1), ', ksp0(2)=', ksp0(2)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ufrc_ref.bin', ufrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vfrc_ref.bin', vfrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wfrc_ref.bin', wfrc_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d_v2(trim(data_dir)//'/rbct8s_ref.bin', rbct8s_ref, 1, ni, 1, nj, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ufrc = ufrc_init
    vfrc = vfrc_init
    wfrc = wfrc_init
    rbct8s = rbct8s_init
    call kernel_vspuvw(gpvvar, vspvar, vspopt, ksp0, gtinc, ni, nj, nk, &
                       ubr, vbr, rst8u, rst8v, rst8w, up, vp, wp, &
                       rbct, ugpv, utd, vgpv, vtd, wgpv, wtd, &
                       ufrc, vfrc, wfrc, rbct8s)
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
    ufrc = ufrc_init
    vfrc = vfrc_init
    wfrc = wfrc_init
    rbct8s = rbct8s_init

    !$acc wait
    t_start = omp_get_wtime()

    call kernel_vspuvw(gpvvar, vspvar, vspopt, ksp0, gtinc, ni, nj, nk, &
                       ubr, vbr, rst8u, rst8v, rst8w, up, vp, wp, &
                       rbct, ugpv, utd, vgpv, vtd, wgpv, wtd, &
                       ufrc, vfrc, wfrc, rbct8s)

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
  max_error = 0.0
  error_count = 0

  call validate_output(ufrc, ufrc_ref, ni, nj, nk, tolerance, tmp_max_err, tmp_err_count)
  if (tmp_max_err > max_error) max_error = tmp_max_err
  error_count = error_count + tmp_err_count
  if (tmp_err_count > 0) write(*,'(A,I6)') '  ufrc errors: ', tmp_err_count

  call validate_output(vfrc, vfrc_ref, ni, nj, nk, tolerance, tmp_max_err, tmp_err_count)
  if (tmp_max_err > max_error) max_error = tmp_max_err
  error_count = error_count + tmp_err_count
  if (tmp_err_count > 0) write(*,'(A,I6)') '  vfrc errors: ', tmp_err_count

  call validate_output(wfrc, wfrc_ref, ni, nj, nk, tolerance, tmp_max_err, tmp_err_count)
  if (tmp_max_err > max_error) max_error = tmp_max_err
  error_count = error_count + tmp_err_count
  if (tmp_err_count > 0) write(*,'(A,I6)') '  wfrc errors: ', tmp_err_count

  call validate_output_v2(rbct8s, rbct8s_ref, ni, nj, nk, tolerance, tmp_max_err, tmp_err_count)
  if (tmp_max_err > max_error) max_error = tmp_max_err
  error_count = error_count + tmp_err_count
  if (tmp_err_count > 0) write(*,'(A,I6)') '  rbct8s errors: ', tmp_err_count

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
  deallocate(ubr, vbr, rst8u, rst8v, rst8w, up, vp, wp, rbct)
  deallocate(ugpv, utd, vgpv, vtd, wgpv, wtd)
  deallocate(ufrc, vfrc, wfrc, ufrc_init, vfrc_init, wfrc_init)
  deallocate(ufrc_ref, vfrc_ref, wfrc_ref, rbct8s, rbct8s_init, rbct8s_ref, times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: vspuvw (GPU version with OpenACC)
  !-------------------------------------------------------------------
  subroutine kernel_vspuvw(gpvvar, vspvar, vspopt, ksp0, gtinc, ni, nj, nk, &
                           ubr, vbr, rst8u, rst8v, rst8w, up, vp, wp, &
                           rbct, ugpv, utd, vgpv, vtd, wgpv, wtd, &
                           ufrc, vfrc, wfrc, rbct8s)
    implicit none

    character(len=*), intent(in) :: gpvvar, vspvar
    integer, intent(in) :: vspopt
    integer, intent(in) :: ksp0(1:2)
    real, intent(in) :: gtinc
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: ubr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vbr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: up(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: wp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rbct(1:ni, 1:nj, 1:nk, 1:2)
    real, intent(in) :: ugpv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: utd(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vgpv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vtd(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: wgpv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: wtd(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: ufrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: wfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: rbct8s(1:ni, 1:nj, 1:nk)

    integer :: i, j, k

    ! Set the common used variable.
    if (vspvar(1:1).eq.'o'.or.vspvar(2:2).eq.'o') then
      if (vspopt.eq.1) then
        do k = ksp0(1)-1, nk-2
          !$acc kernels
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni-1
              rbct8s(i,j,k) = 0.25e0 * (rbct(i,j,k,1) + rbct(i,j,k+1,1))
            end do
          end do
          !$acc end kernels
        end do
      else
        do k = ksp0(2)-1, nk-2
          !$acc kernels
          !$acc loop independent
          do j = 1, nj-1
            !$acc loop independent
            do i = 1, ni-1
              rbct8s(i,j,k) = 0.25e0 * (rbct(i,j,k,2) + rbct(i,j,k+1,2))
            end do
          end do
          !$acc end kernels
        end do
      end if
    end if

    ! For the x components of velocity.
    if (vspvar(1:1).eq.'o') then
      if (vspopt.eq.1) then
        do k = ksp0(1)-1, nk-2
          !$acc kernels
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) &
                   - rst8u(i,j,k) * (rbct8s(i-1,j,k) + rbct8s(i,j,k)) &
                   * (up(i,j,k) - (ugpv(i,j,k) + utd(i,j,k) * gtinc))
            end do
          end do
          !$acc end kernels
        end do
      else
        do k = ksp0(2)-1, nk-2
          !$acc kernels
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) - rst8u(i,j,k) &
                   * (rbct8s(i-1,j,k) + rbct8s(i,j,k)) * (up(i,j,k) - ubr(i,j,k))
            end do
          end do
          !$acc end kernels
        end do
      end if
    end if

    ! For the y components of velocity.
    if (vspvar(2:2).eq.'o') then
      if (vspopt.eq.1) then
        do k = ksp0(1)-1, nk-2
          !$acc kernels
          !$acc loop independent
          do j = 2, nj-1
            !$acc loop independent
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) &
                   - rst8v(i,j,k) * (rbct8s(i,j-1,k) + rbct8s(i,j,k)) &
                   * (vp(i,j,k) - (vgpv(i,j,k) + vtd(i,j,k) * gtinc))
            end do
          end do
          !$acc end kernels
        end do
      else
        do k = ksp0(2)-1, nk-2
          !$acc kernels
          !$acc loop independent
          do j = 2, nj-1
            !$acc loop independent
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) - rst8v(i,j,k) &
                   * (rbct8s(i,j-1,k) + rbct8s(i,j,k)) * (vp(i,j,k) - vbr(i,j,k))
            end do
          end do
          !$acc end kernels
        end do
      end if
    end if

    ! For the z components of velocity.
    if (vspvar(3:3).eq.'o') then
      if (vspopt.eq.1.and.gpvvar(1:1).eq.'o') then
        do k = ksp0(1), nk-1
          !$acc kernels
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-2
              wfrc(i,j,k) = wfrc(i,j,k) - rbct(i,j,k,1) * rst8w(i,j,k) &
                   * (wp(i,j,k) - (wgpv(i,j,k) + wtd(i,j,k) * gtinc))
            end do
          end do
          !$acc end kernels
        end do
      else
        do k = ksp0(2), nk-1
          !$acc kernels
          !$acc loop independent
          do j = 2, nj-2
            !$acc loop independent
            do i = 2, ni-2
              wfrc(i,j,k) = wfrc(i,j,k) - rbct(i,j,k,2) * rst8w(i,j,k) * wp(i,j,k)
            end do
          end do
          !$acc end kernels
        end do
      end if
    end if

  end subroutine kernel_vspuvw

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

    gpvvar = ''
    vspvar = ''
    vspopt = 0
    gtinc = 0.0
    ! Set default values for ksp0 (matching vspdmp.f90 initialization)
    ksp0(1) = 2
    ksp0(2) = 2

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
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('vspopt')
            read(val, *) vspopt
          case ('gtinc')
            read(val, *) gtinc
          case ('gpvvar')
            gpvvar = trim(val)
          case ('vspvar')
            vspvar = trim(val)
          case ('ksp0_1')
            read(val, *) ksp0(1)
          case ('ksp0_2')
            read(val, *) ksp0(2)
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Try to read 1D integer array (optional - uses default if file not found)
  !-------------------------------------------------------------------
  subroutine try_read_array_1d_int(filename, arr, is, ie)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie
    integer, intent(inout) :: arr(is:ie)
    integer :: ios
    logical :: file_exists

    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
      write(*,'(A,A)') '  Note: File not found, using default values: ', trim(filename)
      return
    end if

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') '  Warning: Cannot open file, using default values: ', trim(filename)
      return
    end if
    read(10) arr
    close(10)
  end subroutine try_read_array_1d_int

  !-------------------------------------------------------------------
  ! Read 3D array (0-indexed in first two dims)
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
  ! Read 3D array (1-indexed)
  !-------------------------------------------------------------------
  subroutine read_array_3d_v2(filename, arr, is, ie, js, je, ks, ke)
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
  end subroutine read_array_3d_v2

  !-------------------------------------------------------------------
  ! Read 4D array
  !-------------------------------------------------------------------
  subroutine read_array_4d(filename, arr, is, ie, js, je, ks, ke, ls, le)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke, ls, le
    real, intent(out) :: arr(is:ie, js:je, ks:ke, ls:le)
    integer :: ios

    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)
  end subroutine read_array_4d

  !-------------------------------------------------------------------
  ! Validate output (0-indexed in first two dims)
  !-------------------------------------------------------------------
  subroutine validate_output(output, reference, ni, nj, nk, tol, max_err, err_count)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: output(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: reference(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = 1, nk
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
  end subroutine validate_output

  !-------------------------------------------------------------------
  ! Validate output (1-indexed)
  !-------------------------------------------------------------------
  subroutine validate_output_v2(output, reference, ni, nj, nk, tol, max_err, err_count)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: output(1:ni, 1:nj, 1:nk)
    real, intent(in) :: reference(1:ni, 1:nj, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count

    real :: rel_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
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
  end subroutine validate_output_v2

end program kernel_benchmark_vspuvw_gpu

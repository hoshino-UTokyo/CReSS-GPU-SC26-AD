!***********************************************************************
! GPU Kernel Benchmark: bc4news (s_bc4news)
!***********************************************************************
!
! Source: Src/bc4news.f90
! Description: Sets boundary conditions at the four corners (SW, SE, NW, NE)
!              by averaging adjacent boundary values.
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_gpu_bc4news
  use omp_lib
  implicit none

  ! Grid dimensions
  integer :: ni, nj, kmax

  ! Boundary indices
  integer :: isw, ise, jss, jsn
  integer :: iswp1, isem1, jssp1, jsnm1

  ! Boundary condition options
  integer :: wbc, ebc, sbc, nbc

  ! Domain decomposition flags (simulated for benchmark)
  integer :: ebsw, ebse, ebnw, ebne
  integer :: isub, jsub, nisub, njsub

  ! Arrays
  real, allocatable :: var(:,:,:)
  real, allocatable :: var_ref(:,:,:)
  real, allocatable :: var_input(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
  integer :: error_count
  logical :: validation_passed

  ! Loop variables
  integer :: iter, i, j, k

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt', &
       ni, nj, kmax, isw, ise, jss, jsn, &
       wbc, ebc, sbc, nbc, &
       ebsw, ebse, ebnw, ebne, isub, jsub, nisub, njsub)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: bc4news'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', kmax=', kmax
  write(*,'(A,I6,A,I6,A,I6,A,I6)') ' Corners: isw=', isw, ', ise=', ise, &
       ', jss=', jss, ', jsn=', jsn
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  ! Derived indices
  iswp1 = isw + 1
  isem1 = ise - 1
  jssp1 = jss + 1
  jsnm1 = jsn - 1

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(var(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var_ref(0:ni+1, 0:nj+1, 1:kmax))
  allocate(var_input(0:ni+1, 0:nj+1, 1:kmax))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/var_in.bin', var_input, 0, ni+1, 0, nj+1, 1, kmax)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/var_ref.bin', var_ref, 0, ni+1, 0, nj+1, 1, kmax)

  !---------------------------------------------------------------------
  ! Warmup iterations (includes GPU JIT compilation)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    var = var_input
    call kernel_bc4news(wbc, ebc, sbc, nbc, isw, ise, jss, jsn, &
         iswp1, isem1, jssp1, jsnm1, ni, nj, kmax, var, &
         ebsw, ebse, ebnw, ebne, isub, jsub, nisub, njsub)
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    var = var_input

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_bc4news(wbc, ebc, sbc, nbc, isw, ise, jss, jsn, &
         iswp1, isem1, jssp1, jsnm1, ni, nj, kmax, var, &
         ebsw, ebse, ebnw, ebne, isub, jsub, nisub, njsub)
    !$acc wait
    t_end = omp_get_wtime()

    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  max_error = 0.0
  error_count = 0

  do k = 1, kmax
    do j = 0, nj+1
      do i = 0, ni+1
        rel_error = abs(var(i,j,k) - var_ref(i,j,k))
        if (abs(var_ref(i,j,k)) > 1.0e-20) then
          rel_error = rel_error / abs(var_ref(i,j,k))
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  validation_passed = (error_count == 0)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.6,A)') ' Average time: ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Total time:   ', t_total * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Min time:     ', minval(times) * 1000.0d0, ' ms'
  write(*,'(A,F12.6,A)') ' Max time:     ', maxval(times) * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,ES12.4)') ' Tolerance:          ', tolerance
  write(*,'(A,I12)') ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(var, var_ref, var_input, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: bc4news (OpenACC version)
  !=====================================================================
  subroutine kernel_bc4news(wbc, ebc, sbc, nbc, isw, ise, jss, jsn, &
       iswp1, isem1, jssp1, jsnm1, ni, nj, kmax, var, &
       ebsw, ebse, ebnw, ebne, isub, jsub, nisub, njsub)
    implicit none

    integer, intent(in) :: wbc, ebc, sbc, nbc
    integer, intent(in) :: isw, ise, jss, jsn
    integer, intent(in) :: iswp1, isem1, jssp1, jsnm1
    integer, intent(in) :: ni, nj, kmax
    real, intent(inout) :: var(0:ni+1, 0:nj+1, 1:kmax)
    integer, intent(in) :: ebsw, ebse, ebnw, ebne
    integer, intent(in) :: isub, jsub, nisub, njsub

    integer :: k

    if(abs(wbc).ne.1.or.abs(ebc).ne.1 &
         .or.abs(sbc).ne.1.or.abs(nbc).ne.1) then

      if(ebsw.eq.1.and.isub.eq.0.and.jsub.eq.0) then
        !$acc kernels
        !$acc loop independent
        do k=1,kmax
          var(isw,jss,k)=.5e0*(var(iswp1,jss,k)+var(isw,jssp1,k))
        end do
        !$acc end kernels
      end if

      if(ebse.eq.1.and.isub.eq.nisub-1.and.jsub.eq.0) then
        !$acc kernels
        !$acc loop independent
        do k=1,kmax
          var(ise,jss,k)=.5e0*(var(isem1,jss,k)+var(ise,jssp1,k))
        end do
        !$acc end kernels
      end if

      if(ebnw.eq.1.and.isub.eq.0.and.jsub.eq.njsub-1) then
        !$acc kernels
        !$acc loop independent
        do k=1,kmax
          var(isw,jsn,k)=.5e0*(var(iswp1,jsn,k)+var(isw,jsnm1,k))
        end do
        !$acc end kernels
      end if

      if(ebne.eq.1.and.isub.eq.nisub-1.and.jsub.eq.njsub-1) then
        !$acc kernels
        !$acc loop independent
        do k=1,kmax
          var(ise,jsn,k)=.5e0*(var(isem1,jsn,k)+var(ise,jsnm1,k))
        end do
        !$acc end kernels
      end if

    end if

  end subroutine kernel_bc4news

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    character(len=256) :: config_file
    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    config_file = 'benchmark.conf'
    inquire(file=config_file, exist=exists)

    if (exists) then
      open(unit=10, file=config_file, status='old', iostat=ios)
      if (ios == 0) then
        read(10, '(A)', iostat=ios) data_dir
        read(10, *, iostat=ios) num_iter
        read(10, *, iostat=ios) warmup_iter
        read(10, *, iostat=ios) tol
        close(10)
      end if
    end if

  end subroutine read_config

  !=====================================================================
  ! Parameter reader
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, kmax, isw, ise, jss, jsn, &
       wbc, ebc, sbc, nbc, ebsw, ebse, ebnw, ebne, isub, jsub, nisub, njsub)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, kmax
    integer, intent(out) :: isw, ise, jss, jsn
    integer, intent(out) :: wbc, ebc, sbc, nbc
    integer, intent(out) :: ebsw, ebse, ebnw, ebne
    integer, intent(out) :: isub, jsub, nisub, njsub

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
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('kmax')
            read(val, *) kmax
          case ('isw')
            read(val, *) isw
          case ('ise')
            read(val, *) ise
          case ('jss')
            read(val, *) jss
          case ('jsn')
            read(val, *) jsn
          case ('wbc')
            read(val, *) wbc
          case ('ebc')
            read(val, *) ebc
          case ('sbc')
            read(val, *) sbc
          case ('nbc')
            read(val, *) nbc
          case ('ebsw')
            read(val, *) ebsw
          case ('ebse')
            read(val, *) ebse
          case ('ebnw')
            read(val, *) ebnw
          case ('ebne')
            read(val, *) ebne
          case ('isub')
            read(val, *) isub
          case ('jsub')
            read(val, *) jsub
          case ('nisub')
            read(val, *) nisub
          case ('njsub')
            read(val, *) njsub
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary 3D array reader
  !=====================================================================
  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_3d

end program kernel_benchmark_gpu_bc4news

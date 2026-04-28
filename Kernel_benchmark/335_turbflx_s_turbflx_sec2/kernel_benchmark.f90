!***********************************************************************
! Kernel Benchmark: turbflx section 2 (s_turbflx_sec2)
!***********************************************************************
!
! Source: Src/turbflx.f90
! Description: Calculate the x, y, z components of turbulent fluxes
!              for optional scalar variable. CPU OpenMP version.
!
!***********************************************************************
program kernel_benchmark_turbflx_sec2
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Extra parameters
  integer :: trnopt, sfcopt
  real :: dxv05, dyv05, dziv, dzv125

  ! Input arrays
  real, allocatable :: jcb(:,:,:)     ! Jacobian
  real, allocatable :: s(:,:,:)       ! Optional scalar variable
  real, allocatable :: rkh8u(:,:,:)   ! Horizontal eddy diffusivity at u points
  real, allocatable :: rkh8v(:,:,:)   ! Horizontal eddy diffusivity at v points
  real, allocatable :: rkv8w(:,:,:)   ! Vertical eddy diffusivity at w points
  real, allocatable :: j31s(:,:,:)    ! 4.0 x j31 x optional scalar variable
  real, allocatable :: j32s(:,:,:)    ! 4.0 x j32 x optional scalar variable
  real, allocatable :: sfrc(:,:,:)    ! Optional scalar forcing term

  ! Output arrays
  real, allocatable :: h1(:,:,:)      ! x components of turbulent fluxes
  real, allocatable :: h2(:,:,:)      ! y components of turbulent fluxes
  real, allocatable :: h3(:,:,:)      ! z components of turbulent fluxes
  real, allocatable :: jcbs(:,:,:)    ! jcb x optional scalar variable

  ! Reference outputs for validation
  real, allocatable :: h1_ref(:,:,:)
  real, allocatable :: h2_ref(:,:,:)
  real, allocatable :: h3_ref(:,:,:)

  ! Backup arrays for iteration
  real, allocatable :: h1_in(:,:,:)
  real, allocatable :: h2_in(:,:,:)
  real, allocatable :: h3_in(:,:,:)
  real, allocatable :: jcbs_in(:,:,:)
  real, allocatable :: j31s_in(:,:,:)
  real, allocatable :: j32s_in(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error_h1, max_error_h2, max_error_h3
  real :: rel_error
  real :: tolerance
  integer :: error_count_h1, error_count_h2, error_count_h3
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
  call read_parameters(trim(data_dir)//'/params.txt', ni, nj, nk, &
                       trnopt, sfcopt, dxv05, dyv05, dziv, dzv125)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: turbflx section 2'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' trnopt: ', trnopt
  write(*,'(A,I6)') ' sfcopt: ', sfcopt
  write(*,'(A,ES12.4)') ' dxv05:  ', dxv05
  write(*,'(A,ES12.4)') ' dyv05:  ', dyv05
  write(*,'(A,ES12.4)') ' dziv:   ', dziv
  write(*,'(A,ES12.4)') ' dzv125: ', dzv125
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(s(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkh8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkh8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(rkv8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(sfrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(j31s(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32s(0:ni+1, 0:nj+1, 1:nk))
  allocate(h1(0:ni+1, 0:nj+1, 1:nk))
  allocate(h2(0:ni+1, 0:nj+1, 1:nk))
  allocate(h3(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcbs(0:ni+1, 0:nj+1, 1:nk))
  allocate(h1_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(h2_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(h3_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(h1_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(h2_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(h3_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcbs_in(0:ni+1, 0:nj+1, 1:nk))
  if (trnopt >= 1) then
    allocate(j31s_in(0:ni+1, 0:nj+1, 1:nk))
    allocate(j32s_in(0:ni+1, 0:nj+1, 1:nk))
  end if
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s.bin', s, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rkh8u.bin', rkh8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rkh8v.bin', rkh8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rkv8w.bin', rkv8w, 0, ni+1, 0, nj+1, 1, nk)
  if (trnopt >= 1) then
    call read_array_3d(trim(data_dir)//'/j31s_in.bin', j31s, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/j32s_in.bin', j32s, 0, ni+1, 0, nj+1, 1, nk)
    j31s_in = j31s
    j32s_in = j32s
  end if
  if (sfcopt >= 1) then
    call read_array_3d(trim(data_dir)//'/sfrc.bin', sfrc, 0, ni+1, 0, nj+1, 1, nk)
  else
    sfrc = 0.0
  end if
  call read_array_3d(trim(data_dir)//'/h1_in.bin', h1_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/h2_in.bin', h2_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/h3_in.bin', h3_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcbs_in.bin', jcbs_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/h1_ref.bin', h1_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/h2_ref.bin', h2_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/h3_ref.bin', h3_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    h1 = h1_in
    h2 = h2_in
    h3 = h3_in
    jcbs = jcbs_in
    if (trnopt >= 1) then
      j31s = j31s_in
      j32s = j32s_in
    end if
    call kernel_turbflx_sec2(ni, nj, nk, trnopt, sfcopt, &
         dxv05, dyv05, dziv, dzv125, &
         jcb, s, rkh8u, rkh8v, rkv8w, j31s, j32s, sfrc, &
         h1, h2, h3, jcbs)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    h1 = h1_in
    h2 = h2_in
    h3 = h3_in
    jcbs = jcbs_in
    if (trnopt >= 1) then
      j31s = j31s_in
      j32s = j32s_in
    end if

    t_start = omp_get_wtime()
    call kernel_turbflx_sec2(ni, nj, nk, trnopt, sfcopt, &
         dxv05, dyv05, dziv, dzv125, &
         jcb, s, rkh8u, rkh8v, rkv8w, j31s, j32s, sfrc, &
         h1, h2, h3, jcbs)
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'

  ! Validate h1
  max_error_h1 = 0.0
  error_count_h1 = 0
  do k = 1, nk-1
    do j = 2, nj-2
      do i = 2, ni-1
        rel_error = abs(h1(i,j,k) - h1_ref(i,j,k))
        if (abs(h1_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(h1_ref(i,j,k))
        end if
        if (rel_error > max_error_h1) max_error_h1 = rel_error
        if (rel_error > tolerance) error_count_h1 = error_count_h1 + 1
      end do
    end do
  end do

  ! Validate h2
  max_error_h2 = 0.0
  error_count_h2 = 0
  do k = 1, nk-1
    do j = 2, nj-1
      do i = 2, ni-2
        rel_error = abs(h2(i,j,k) - h2_ref(i,j,k))
        if (abs(h2_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(h2_ref(i,j,k))
        end if
        if (rel_error > max_error_h2) max_error_h2 = rel_error
        if (rel_error > tolerance) error_count_h2 = error_count_h2 + 1
      end do
    end do
  end do

  ! Validate h3
  max_error_h3 = 0.0
  error_count_h3 = 0
  do k = 2, nk-1
    do j = 2, nj-2
      do i = 2, ni-2
        rel_error = abs(h3(i,j,k) - h3_ref(i,j,k))
        if (abs(h3_ref(i,j,k)) > 1.0e-10) then
          rel_error = rel_error / abs(h3_ref(i,j,k))
        end if
        if (rel_error > max_error_h3) max_error_h3 = rel_error
        if (rel_error > tolerance) error_count_h3 = error_count_h3 + 1
      end do
    end do
  end do

  validation_passed = (error_count_h1 == 0 .and. error_count_h2 == 0 &
                       .and. error_count_h3 == 0)

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
  write(*,'(A,ES12.4)') ' h1 max error:   ', max_error_h1
  write(*,'(A,I12)')    ' h1 error count: ', error_count_h1
  write(*,'(A,ES12.4)') ' h2 max error:   ', max_error_h2
  write(*,'(A,I12)')    ' h2 error count: ', error_count_h2
  write(*,'(A,ES12.4)') ' h3 max error:   ', max_error_h3
  write(*,'(A,I12)')    ' h3 error count: ', error_count_h3
  write(*,'(A,ES12.4)') ' Tolerance:      ', tolerance
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(jcb, s, rkh8u, rkh8v, rkv8w, sfrc)
  deallocate(j31s, j32s)
  deallocate(h1, h2, h3, jcbs)
  deallocate(h1_ref, h2_ref, h3_ref)
  deallocate(h1_in, h2_in, h3_in, jcbs_in)
  if (allocated(j31s_in)) deallocate(j31s_in)
  if (allocated(j32s_in)) deallocate(j32s_in)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: turbflx section 2 - OpenMP version
  !=====================================================================
  subroutine kernel_turbflx_sec2(ni, nj, nk, trnopt, sfcopt, &
       dxv05, dyv05, dziv, dzv125, &
       jcb, s, rkh8u, rkh8v, rkv8w, j31s, j32s, sfrc, &
       h1, h2, h3, jcbs)
    implicit none

    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: trnopt, sfcopt
    real, intent(in) :: dxv05, dyv05, dziv, dzv125
    real, intent(in) :: jcb(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: s(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rkh8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rkh8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rkv8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j31s(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32s(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: sfrc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: h1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: h2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: h3(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: jcbs(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! Compute jcbs = s * jcb
    do k = 1, nk-1
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          jcbs(i,j,k) = s(i,j,k) * jcb(i,j,k)
        end do
      end do
      !$omp end do
    end do

    ! Calculate the x components of the turbulent fluxes
    if (trnopt == 0) then

      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-1
            h1(i,j,k) = rkh8u(i,j,k) * (jcbs(i,j,k) - jcbs(i-1,j,k)) * dxv05
          end do
        end do
        !$omp end do
      end do

    else if (trnopt >= 1) then

      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
          do i = 2, ni-1
            h1(i,j,k) = rkh8u(i,j,k) * ((jcbs(i,j,k) - jcbs(i-1,j,k)) * dxv05 &
                 + (j31s(i,j,k+1) - j31s(i,j,k)) * dzv125)
          end do
        end do
        !$omp end do
      end do

    end if

    ! Calculate the y components of the turbulent fluxes
    if (trnopt == 0) then

      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-1
          do i = 2, ni-2
            h2(i,j,k) = rkh8v(i,j,k) * (jcbs(i,j,k) - jcbs(i,j-1,k)) * dyv05
          end do
        end do
        !$omp end do
      end do

    else if (trnopt >= 1) then

      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-1
          do i = 2, ni-2
            h2(i,j,k) = rkh8v(i,j,k) * ((jcbs(i,j,k) - jcbs(i,j-1,k)) * dyv05 &
                 + (j32s(i,j,k+1) - j32s(i,j,k)) * dzv125)
          end do
        end do
        !$omp end do
      end do

    end if

    ! Get the z components of the turbulent fluxes
    do k = 2, nk-1
      !$omp do schedule(runtime) private(i,j)
      do j = 2, nj-2
        do i = 2, ni-2
          h3(i,j,k) = rkv8w(i,j,k) * (s(i,j,k) - s(i,j,k-1)) * dziv
        end do
      end do
      !$omp end do
    end do

    !$omp end parallel

    if (sfcopt >= 1) then

      !$omp parallel do schedule(runtime) private(i,j)
      do j = 2, nj-2
        do i = 2, ni-2
          h3(i,j,2) = sfrc(i,j,1)
        end do
      end do
      !$omp end parallel do

    end if

  end subroutine kernel_turbflx_sec2

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios
    logical :: exists
    data_dir = './data'; num_iter = 10; warmup_iter = 2; tol = 1.0e-5
    inquire(file='benchmark.conf', exist=exists)
    if (exists) then
      open(unit=10, file='benchmark.conf', status='old', iostat=ios)
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
  subroutine read_parameters(filename, ni, nj, nk, &
       trnopt, sfcopt, dxv05, dyv05, dziv, dzv125)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    integer, intent(out) :: trnopt, sfcopt
    real, intent(out) :: dxv05, dyv05, dziv, dzv125
    character(len=256) :: line, key, val
    integer :: ios, eq_pos
    ni = 1; nj = 1; nk = 1
    trnopt = 0; sfcopt = 0
    dxv05 = 0.0; dyv05 = 0.0; dziv = 0.0; dzv125 = 0.0
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
          case ('trnopt')
            read(val, *) trnopt
          case ('sfcopt')
            read(val, *) sfcopt
          case ('dxv05')
            read(val, *) dxv05
          case ('dyv05')
            read(val, *) dyv05
          case ('dziv')
            read(val, *) dziv
          case ('dzv125')
            read(val, *) dzv125
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !=====================================================================
  ! Binary array reader
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

end program kernel_benchmark_turbflx_sec2

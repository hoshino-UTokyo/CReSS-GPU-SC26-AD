!***********************************************************************
! Kernel Benchmark: advp (s_advp)
!***********************************************************************
! Description:
!   Calculates pressure advection using 2nd or 4th order finite
!   difference methods with Jacobian weighting for terrain-following
!   coordinates. Multiple code paths based on advopt, mpopt, mfcopt,
!   diaopt options.
!
! Original location: Src/advp.f90
! Profiling info: 360 calls, 18.87s total, 52.42ms avg
!***********************************************************************

program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup
  integer :: num_iterations
  real :: tolerance

  ! Array dimensions and parameters
  integer :: ni, nj, nk
  integer :: advopt, mpopt, mfcopt, diaopt
  integer :: iwest, ieast, jsouth, jnorth
  real :: dxiv, dyiv, dziv

  ! Arrays - 2D
  real, allocatable :: mf8u(:,:), mf8v(:,:)

  ! Arrays - 3D inputs
  real, allocatable :: jcb8u(:,:,:), jcb8v(:,:,:), jcb8w(:,:,:)
  real, allocatable :: u(:,:,:), v(:,:,:), wc(:,:,:), pp(:,:,:)

  ! Arrays - 3D input/output (initial values saved for reset)
  real, allocatable :: pfrc(:,:,:), pfrc_init(:,:,:)
  real, allocatable :: jcbxu(:,:,:), jcbxv(:,:,:), jcbxwc(:,:,:)
  real, allocatable :: jcbxu_init(:,:,:), jcbxv_init(:,:,:), jcbxwc_init(:,:,:)
  real, allocatable :: hadv(:,:,:), vadv(:,:,:)
  real, allocatable :: hadv_init(:,:,:), vadv_init(:,:,:)
  real, allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:)
  real, allocatable :: tmp1_init(:,:,:), tmp2_init(:,:,:), tmp3_init(:,:,:)

  ! Reference output arrays
  real, allocatable :: pfrc_ref(:,:,:)
  real, allocatable :: jcbxu_ref(:,:,:), jcbxv_ref(:,:,:), jcbxwc_ref(:,:,:)
  real, allocatable :: hadv_ref(:,:,:), vadv_ref(:,:,:)
  real, allocatable :: tmp1_ref(:,:,:), tmp2_ref(:,:,:), tmp3_ref(:,:,:)

  ! Timing variables
  real(8) :: t_start, t_end, t_total, t_min, t_max, t_avg
  real(8), allocatable :: times(:)

  ! Loop counters and validation
  integer :: iter, ios
  integer :: i, j, k
  real :: max_err, rel_err, ref_max
  logical :: validation_passed

  ! Read configuration
  call read_config()

  ! Read parameters
  call read_parameters()

  ! Allocate arrays
  call allocate_arrays()

  ! Read input data
  call read_input_data()

  ! Read reference output
  call read_reference_data()

  ! Warmup iterations
  write(*,'(A)') '=== Kernel Benchmark: advp (s_advp) ==='
  write(*,'(A,I3,A)') 'Performing ', num_warmup, ' warmup iterations...'

  do iter = 1, num_warmup
    call reset_arrays()
    call kernel_advp()
  end do

  ! Benchmark iterations
  write(*,'(A,I4,A)') 'Performing ', num_iterations, ' benchmark iterations...'
  allocate(times(num_iterations))

  t_total = 0.0d0
  t_min = huge(t_min)
  t_max = 0.0d0

  do iter = 1, num_iterations
    call reset_arrays()

    t_start = omp_get_wtime()
    call kernel_advp()
    t_end = omp_get_wtime()

    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
    t_min = min(t_min, times(iter))
    t_max = max(t_max, times(iter))
  end do

  t_avg = t_total / dble(num_iterations)

  ! Validate results
  call validate_results(validation_passed, max_err)

  ! Report results
  write(*,'(A)') ''
  write(*,'(A)') '=== Results ==='
  write(*,'(A,ES12.5,A)') 'Total time:   ', t_total, ' seconds'
  write(*,'(A,ES12.5,A)') 'Average time: ', t_avg, ' seconds'
  write(*,'(A,ES12.5,A)') 'Min time:     ', t_min, ' seconds'
  write(*,'(A,ES12.5,A)') 'Max time:     ', t_max, ' seconds'
  write(*,'(A)') ''
  write(*,'(A,L1)') 'Validation passed: ', validation_passed
  write(*,'(A,ES12.5)') 'Maximum relative error: ', max_err
  write(*,'(A)') ''
  write(*,'(A,I6,A,I4,A,I4)') 'Grid size: ni=', ni, ', nj=', nj, ', nk=', nk

  ! Cleanup
  deallocate(times)
  call deallocate_arrays()

contains

  !---------------------------------------------------------------------
  subroutine read_config()
    integer :: unit_num
    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'Error: Cannot open benchmark.conf'
      stop 1
    end if
    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iterations
    read(unit_num, *) num_warmup
    read(unit_num, *) tolerance
    close(unit_num)
  end subroutine read_config

  !---------------------------------------------------------------------
  subroutine read_parameters()
    integer :: unit_num, tmp_int
    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'Error: Cannot open params.txt'
      stop 1
    end if
    ! Read integer parameters (fp* values are index references, use as actual values)
    read(unit_num, *) advopt   ! fpadvopt -> actual advopt value
    read(unit_num, *) mpopt    ! fpmpopt -> actual mpopt value
    read(unit_num, *) mfcopt   ! fpmfcopt -> actual mfcopt value
    read(unit_num, *) diaopt   ! fpdiaopt -> actual diaopt value
    read(unit_num, *) iwest    ! fpiwest -> actual iwest value
    read(unit_num, *) ieast    ! fpieast -> actual ieast value
    read(unit_num, *) jsouth   ! fpjsouth -> actual jsouth value
    read(unit_num, *) jnorth   ! fpjnorth -> actual jnorth value
    read(unit_num, *) dxiv     ! fpdxiv -> actual dxiv value
    read(unit_num, *) dyiv     ! fpdyiv -> actual dyiv value
    read(unit_num, *) dziv     ! fpdziv -> actual dziv value
    read(unit_num, *) ni
    read(unit_num, *) nj
    read(unit_num, *) nk
    close(unit_num)
  end subroutine read_parameters

  !---------------------------------------------------------------------
  subroutine allocate_arrays()
    ! 2D arrays
    allocate(mf8u(0:ni+1, 0:nj+1))
    allocate(mf8v(0:ni+1, 0:nj+1))

    ! 3D input arrays
    allocate(jcb8u(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcb8v(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcb8w(0:ni+1, 0:nj+1, 1:nk))
    allocate(u(0:ni+1, 0:nj+1, 1:nk))
    allocate(v(0:ni+1, 0:nj+1, 1:nk))
    allocate(wc(0:ni+1, 0:nj+1, 1:nk))
    allocate(pp(0:ni+1, 0:nj+1, 1:nk))

    ! Input/output arrays and their initial copies
    allocate(pfrc(0:ni+1, 0:nj+1, 1:nk))
    allocate(pfrc_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcbxu(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcbxu_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcbxv(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcbxv_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcbxwc(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcbxwc_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(hadv(0:ni+1, 0:nj+1, 1:nk))
    allocate(hadv_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(vadv(0:ni+1, 0:nj+1, 1:nk))
    allocate(vadv_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp1_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp2_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp3_init(0:ni+1, 0:nj+1, 1:nk))

    ! Reference arrays
    allocate(pfrc_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcbxu_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcbxv_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(jcbxwc_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(hadv_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(vadv_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp1_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp2_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp3_ref(0:ni+1, 0:nj+1, 1:nk))
  end subroutine allocate_arrays

  !---------------------------------------------------------------------
  subroutine read_input_data()
    integer :: unit_num
    unit_num = 12

    ! Read 2D arrays
    open(unit=unit_num, file=trim(data_dir)//'/mf8u.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading mf8u.bin'
    read(unit_num) mf8u
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/mf8v.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading mf8v.bin'
    read(unit_num) mf8v
    close(unit_num)

    ! Read 3D input arrays
    open(unit=unit_num, file=trim(data_dir)//'/jcb8u.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcb8u.bin'
    read(unit_num) jcb8u
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/jcb8v.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcb8v.bin'
    read(unit_num) jcb8v
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/jcb8w.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcb8w.bin'
    read(unit_num) jcb8w
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/u.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading u.bin'
    read(unit_num) u
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/v.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading v.bin'
    read(unit_num) v
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/wc.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading wc.bin'
    read(unit_num) wc
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/pp.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading pp.bin'
    read(unit_num) pp
    close(unit_num)

    ! Read input/output initial values
    open(unit=unit_num, file=trim(data_dir)//'/pfrc_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading pfrc_in.bin'
    read(unit_num) pfrc_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/jcbxu_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcbxu_in.bin'
    read(unit_num) jcbxu_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/jcbxv_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcbxv_in.bin'
    read(unit_num) jcbxv_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/jcbxwc_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcbxwc_in.bin'
    read(unit_num) jcbxwc_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/hadv_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading hadv_in.bin'
    read(unit_num) hadv_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/vadv_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading vadv_in.bin'
    read(unit_num) vadv_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp1_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp1_in.bin'
    read(unit_num) tmp1_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp2_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp2_in.bin'
    read(unit_num) tmp2_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp3_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp3_in.bin'
    read(unit_num) tmp3_init
    close(unit_num)
  end subroutine read_input_data

  !---------------------------------------------------------------------
  subroutine read_reference_data()
    integer :: unit_num
    unit_num = 13

    open(unit=unit_num, file=trim(data_dir)//'/pfrc_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading pfrc_ref.bin'
    read(unit_num) pfrc_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/jcbxu_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcbxu_ref.bin'
    read(unit_num) jcbxu_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/jcbxv_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcbxv_ref.bin'
    read(unit_num) jcbxv_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/jcbxwc_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcbxwc_ref.bin'
    read(unit_num) jcbxwc_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/hadv_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading hadv_ref.bin'
    read(unit_num) hadv_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/vadv_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading vadv_ref.bin'
    read(unit_num) vadv_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp1_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp1_ref.bin'
    read(unit_num) tmp1_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp2_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp2_ref.bin'
    read(unit_num) tmp2_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp3_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp3_ref.bin'
    read(unit_num) tmp3_ref
    close(unit_num)
  end subroutine read_reference_data

  !---------------------------------------------------------------------
  subroutine reset_arrays()
    pfrc = pfrc_init
    jcbxu = jcbxu_init
    jcbxv = jcbxv_init
    jcbxwc = jcbxwc_init
    hadv = hadv_init
    vadv = vadv_init
    tmp1 = tmp1_init
    tmp2 = tmp2_init
    tmp3 = tmp3_init
  end subroutine reset_arrays

  !---------------------------------------------------------------------
  subroutine validate_results(passed, max_error)
    logical, intent(out) :: passed
    real, intent(out) :: max_error
    real :: err, ref_val

    passed = .true.
    max_error = 0.0

    ! Validate pfrc
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          ref_val = abs(pfrc_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(pfrc(i,j,k) - pfrc_ref(i,j,k)) / ref_val
          else
            err = abs(pfrc(i,j,k) - pfrc_ref(i,j,k))
          end if
          max_error = max(max_error, err)
          if (err > tolerance) passed = .false.
        end do
      end do
    end do
  end subroutine validate_results

  !---------------------------------------------------------------------
  subroutine deallocate_arrays()
    deallocate(mf8u, mf8v)
    deallocate(jcb8u, jcb8v, jcb8w)
    deallocate(u, v, wc, pp)
    deallocate(pfrc, pfrc_init, pfrc_ref)
    deallocate(jcbxu, jcbxu_init, jcbxu_ref)
    deallocate(jcbxv, jcbxv_init, jcbxv_ref)
    deallocate(jcbxwc, jcbxwc_init, jcbxwc_ref)
    deallocate(hadv, hadv_init, hadv_ref)
    deallocate(vadv, vadv_init, vadv_ref)
    deallocate(tmp1, tmp1_init, tmp1_ref)
    deallocate(tmp2, tmp2_init, tmp2_ref)
    deallocate(tmp3, tmp3_init, tmp3_ref)
  end subroutine deallocate_arrays

  !---------------------------------------------------------------------
  subroutine kernel_advp()
    ! Module constants from commath
    real, parameter :: oned24 = 1.0e0 / 24.0e0
    real, parameter :: fourd3 = 4.0e0 / 3.0e0

    ! Local variables
    real :: dxv05n, dyv05n, dzv05n
    real :: dxv24, dyv24, dzv24

    ! Set the common used variables
    dxv05n = -0.5e0 * dxiv
    dyv05n = -0.5e0 * dyiv
    dzv05n = -0.5e0 * dziv

    dxv24 = oned24 * dxiv
    dyv24 = oned24 * dyiv
    dzv24 = oned24 * dziv

    !$omp parallel default(shared) private(k)

    ! Perform the centered fdm scheme
    if (advopt <= 3) then

      ! The variables at the u, v and w points multiplied by u, v and w
      if (mfcopt == 0) then

        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
          do i = 1, ni
            jcbxu(i,j,k) = jcb8u(i,j,k) * u(i,j,k)
          end do
          end do
          !$omp end do

          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj
          do i = 1, ni-1
            jcbxv(i,j,k) = jcb8v(i,j,k) * v(i,j,k)
          end do
          end do
          !$omp end do
        end do

        do k = 1, nk
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
          do i = 1, ni-1
            jcbxwc(i,j,k) = jcb8w(i,j,k) * wc(i,j,k)
          end do
          end do
          !$omp end do
        end do

      else

        if (mpopt == 0 .or. mpopt == 10) then
          do k = 1, nk-1
            !$omp do schedule(runtime) private(i,j)
            do j = 1, nj-1
            do i = 1, ni
              jcbxu(i,j,k) = mf8u(i,j) * jcb8u(i,j,k) * u(i,j,k)
            end do
            end do
            !$omp end do

            !$omp do schedule(runtime) private(i,j)
            do j = 1, nj
            do i = 1, ni-1
              jcbxv(i,j,k) = jcb8v(i,j,k) * v(i,j,k)
            end do
            end do
            !$omp end do
          end do

        else if (mpopt == 5) then
          do k = 1, nk-1
            !$omp do schedule(runtime) private(i,j)
            do j = 1, nj-1
            do i = 1, ni
              jcbxu(i,j,k) = jcb8u(i,j,k) * u(i,j,k)
            end do
            end do
            !$omp end do

            !$omp do schedule(runtime) private(i,j)
            do j = 1, nj
            do i = 1, ni-1
              jcbxv(i,j,k) = mf8v(i,j) * jcb8v(i,j,k) * v(i,j,k)
            end do
            end do
            !$omp end do
          end do

        else
          do k = 1, nk-1
            !$omp do schedule(runtime) private(i,j)
            do j = 1, nj-1
            do i = 1, ni
              jcbxu(i,j,k) = mf8u(i,j) * jcb8u(i,j,k) * u(i,j,k)
            end do
            end do
            !$omp end do

            !$omp do schedule(runtime) private(i,j)
            do j = 1, nj
            do i = 1, ni-1
              jcbxv(i,j,k) = mf8v(i,j) * jcb8v(i,j,k) * v(i,j,k)
            end do
            end do
            !$omp end do
          end do
        end if

        do k = 1, nk
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
          do i = 1, ni-1
            jcbxwc(i,j,k) = jcb8w(i,j,k) * wc(i,j,k)
          end do
          end do
          !$omp end do
        end do

      end if

      ! Calculate the 2nd order pressure advection
      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
        do i = 2, ni-1
          tmp1(i,j,k) = jcbxu(i,j,k) * (pp(i,j,k) - pp(i-1,j,k)) * dxv05n
        end do
        end do
        !$omp end do

        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-1
        do i = 2, ni-2
          tmp2(i,j,k) = jcbxv(i,j,k) * (pp(i,j,k) - pp(i,j-1,k)) * dyv05n
        end do
        end do
        !$omp end do
      end do

      do k = 2, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
        do i = 2, ni-2
          tmp3(i,j,k) = jcbxwc(i,j,k) * (pp(i,j,k) - pp(i,j,k-1)) * dzv05n
        end do
        end do
        !$omp end do
      end do

      if (diaopt == 0) then
        if (advopt == 1) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
            do i = 2, ni-2
              pfrc(i,j,k) = (tmp3(i,j,k) + tmp3(i,j,k+1)) &
                + ((tmp1(i,j,k) + tmp1(i+1,j,k)) &
                + (tmp2(i,j,k) + tmp2(i,j+1,k)))
            end do
            end do
            !$omp end do
          end do

        else
          if (advopt == 2) then
            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
              do i = 2, ni-2
                pfrc(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) &
                  + (tmp2(i,j,k) + tmp2(i,j+1,k))
              end do
              end do
              !$omp end do
            end do

          else if (advopt == 3) then
            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
              do i = 2, ni-2
                pfrc(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) &
                  + (tmp2(i,j,k) + tmp2(i,j+1,k))
                vadv(i,j,k) = tmp3(i,j,k) + tmp3(i,j,k+1)
              end do
              end do
              !$omp end do
            end do
          end if
        end if

      else
        if (advopt == 1) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2, nj-2
            do i = 2, ni-2
              pfrc(i,j,k) = pfrc(i,j,k) + ((tmp3(i,j,k) + tmp3(i,j,k+1)) &
                + ((tmp1(i,j,k) + tmp1(i+1,j,k)) &
                + (tmp2(i,j,k) + tmp2(i,j+1,k))))
            end do
            end do
            !$omp end do
          end do

        else
          if (advopt == 2) then
            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
              do i = 2, ni-2
                hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) &
                  + (tmp2(i,j,k) + tmp2(i,j+1,k))
              end do
              end do
              !$omp end do
            end do

          else if (advopt == 3) then
            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
              do i = 2, ni-2
                hadv(i,j,k) = (tmp1(i,j,k) + tmp1(i+1,j,k)) &
                  + (tmp2(i,j,k) + tmp2(i,j+1,k))
                vadv(i,j,k) = tmp3(i,j,k) + tmp3(i,j,k+1)
              end do
              end do
              !$omp end do
            end do
          end if
        end if
      end if

      ! Calculate the 4th order pressure advection
      if (advopt == 2 .or. advopt == 3) then
        do k = 2, nk-2
          !$omp do schedule(runtime) private(i,j)
          do j = 2+jsouth, nj-2-jnorth
          do i = 1+iwest, ni-1-ieast
            tmp1(i,j,k) = (jcbxu(i,j,k) + jcbxu(i+1,j,k)) &
              * (pp(i+1,j,k) - pp(i-1,j,k)) * dxv24
          end do
          end do
          !$omp end do

          !$omp do schedule(runtime) private(i,j)
          do j = 1+jsouth, nj-1-jnorth
          do i = 2+iwest, ni-2-ieast
            tmp2(i,j,k) = (jcbxv(i,j,k) + jcbxv(i,j+1,k)) &
              * (pp(i,j+1,k) - pp(i,j-1,k)) * dyv24
          end do
          end do
          !$omp end do
        end do

        if (diaopt == 0) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-2-ieast
              pfrc(i,j,k) = fourd3 * pfrc(i,j,k) &
                + ((tmp1(i-1,j,k) + tmp1(i+1,j,k)) &
                + (tmp2(i,j-1,k) + tmp2(i,j+1,k)))
            end do
            end do
            !$omp end do
          end do

          if (advopt == 2) then
            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
              do i = 2, ni-2
                pfrc(i,j,k) = pfrc(i,j,k) + (tmp3(i,j,k) + tmp3(i,j,k+1))
              end do
              end do
              !$omp end do
            end do

          else if (advopt == 3) then
            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2+jsouth, nj-2-jnorth
              do i = 2+iwest, ni-2-ieast
                tmp3(i,j,k) = (jcbxwc(i,j,k) + jcbxwc(i,j,k+1)) &
                  * (pp(i,j,k+1) - pp(i,j,k-1)) * dzv24
              end do
              end do
              !$omp end do
            end do

            do k = 3, nk-3
              !$omp do schedule(runtime) private(i,j)
              do j = 2+jsouth, nj-2-jnorth
              do i = 2+iwest, ni-2-ieast
                vadv(i,j,k) = fourd3 * vadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k+1))
              end do
              end do
              !$omp end do
            end do

            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
              do i = 2, ni-2
                pfrc(i,j,k) = pfrc(i,j,k) + vadv(i,j,k)
              end do
              end do
              !$omp end do
            end do
          end if

        else
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i,j)
            do j = 2+jsouth, nj-2-jnorth
            do i = 2+iwest, ni-2-ieast
              hadv(i,j,k) = fourd3 * hadv(i,j,k) &
                + ((tmp1(i-1,j,k) + tmp1(i+1,j,k)) &
                + (tmp2(i,j-1,k) + tmp2(i,j+1,k)))
            end do
            end do
            !$omp end do
          end do

          if (advopt == 2) then
            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
              do i = 2, ni-2
                pfrc(i,j,k) = pfrc(i,j,k) &
                  + (hadv(i,j,k) + (tmp3(i,j,k) + tmp3(i,j,k+1)))
              end do
              end do
              !$omp end do
            end do

          else if (advopt == 3) then
            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2+jsouth, nj-2-jnorth
              do i = 2+iwest, ni-2-ieast
                tmp3(i,j,k) = (jcbxwc(i,j,k) + jcbxwc(i,j,k+1)) &
                  * (pp(i,j,k+1) - pp(i,j,k-1)) * dzv24
              end do
              end do
              !$omp end do
            end do

            do k = 3, nk-3
              !$omp do schedule(runtime) private(i,j)
              do j = 2+jsouth, nj-2-jnorth
              do i = 2+iwest, ni-2-ieast
                vadv(i,j,k) = fourd3 * vadv(i,j,k) + (tmp3(i,j,k-1) + tmp3(i,j,k+1))
              end do
              end do
              !$omp end do
            end do

            do k = 2, nk-2
              !$omp do schedule(runtime) private(i,j)
              do j = 2, nj-2
              do i = 2, ni-2
                pfrc(i,j,k) = pfrc(i,j,k) + (hadv(i,j,k) + vadv(i,j,k))
              end do
              end do
              !$omp end do
            end do
          end if
        end if
      end if

    else
      ! Set the forcing term with 0 in the case the Cubic Lagrange scheme is performed
      do k = 2, nk-2
        !$omp do schedule(runtime) private(i,j)
        do j = 2, nj-2
        do i = 2, ni-2
          pfrc(i,j,k) = 0.0e0
        end do
        end do
        !$omp end do
      end do
    end if

    !$omp end parallel

  end subroutine kernel_advp

end program kernel_benchmark

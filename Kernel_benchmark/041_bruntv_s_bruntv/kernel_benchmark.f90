!***********************************************************************
! Kernel Benchmark: bruntv (s_bruntv)
!***********************************************************************
! Description:
!   Calculates Brunt-Vaisala frequency squared for atmospheric
!   stability using potential temperature and moisture fields.
!
! Original location: Src/bruntv.f90
! Profiling info: 360 calls, 11.98s total, 33.28ms avg
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
  integer :: cphopt
  real :: dziv, thresq
  character(len=5) :: fmois

  ! Physical constants (read from params.txt)
  real :: g, rd, cp, cw, ci, p0, lv0, lf0, t0, tlow, epsav, epsva

  ! Arrays - 3D inputs
  real, allocatable :: jcb8w(:,:,:), pbr(:,:,:), ptbr(:,:,:)
  real, allocatable :: pp(:,:,:), ptp(:,:,:), qv(:,:,:), qall(:,:,:)

  ! Arrays - output
  real, allocatable :: nsq8w(:,:,:)

  ! Arrays - working arrays (input/output)
  real, allocatable :: pt(:,:,:), ptv(:,:,:), a(:,:,:)
  real, allocatable :: pt_init(:,:,:), ptv_init(:,:,:), a_init(:,:,:)
  real, allocatable :: t(:,:), t_init(:,:)

  ! Reference output arrays
  real, allocatable :: nsq8w_ref(:,:,:)
  real, allocatable :: pt_ref(:,:,:), ptv_ref(:,:,:), a_ref(:,:,:)
  real, allocatable :: t_ref(:,:)

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
  write(*,'(A)') '=== Kernel Benchmark: bruntv (s_bruntv) ==='
  write(*,'(A,I3,A)') 'Performing ', num_warmup, ' warmup iterations...'

  do iter = 1, num_warmup
    call reset_arrays()
    call kernel_bruntv()
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
    call kernel_bruntv()
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
  write(*,'(A,A)') 'fmois: ', trim(fmois)
  write(*,'(A,I2)') 'cphopt: ', cphopt

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
    integer :: unit_num, eq_pos
    character(len=256) :: line, name, val
    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'Error: Cannot open params.txt'
      stop 1
    end if
    ! Initialize fmois to default
    fmois = 'moist'
    do while (.true.)
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        name = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))
        select case (trim(name))
          case ('cphopt')
            read(val, *) cphopt
          case ('dziv')
            read(val, *) dziv
          case ('thresq')
            read(val, *) thresq
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('fmois')
            fmois = trim(val)
          case ('g')
            read(val, *) g
          case ('rd')
            read(val, *) rd
          case ('cp')
            read(val, *) cp
          case ('cw')
            read(val, *) cw
          case ('ci')
            read(val, *) ci
          case ('p0')
            read(val, *) p0
          case ('lv0')
            read(val, *) lv0
          case ('lf0')
            read(val, *) lf0
          case ('t0')
            read(val, *) t0
          case ('tlow')
            read(val, *) tlow
          case ('epsav')
            read(val, *) epsav
          case ('epsva')
            read(val, *) epsva
        end select
      end if
    end do
    close(unit_num)
  end subroutine read_parameters

  !---------------------------------------------------------------------
  subroutine allocate_arrays()
    ! 3D input arrays
    allocate(jcb8w(0:ni+1, 0:nj+1, 1:nk))
    allocate(pbr(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
    allocate(pp(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
    allocate(qv(0:ni+1, 0:nj+1, 1:nk))
    allocate(qall(0:ni+1, 0:nj+1, 1:nk))

    ! Output arrays
    allocate(nsq8w(0:ni+1, 0:nj+1, 1:nk))
    allocate(nsq8w_ref(0:ni+1, 0:nj+1, 1:nk))

    ! Working arrays
    allocate(pt(0:ni+1, 0:nj+1, 1:nk))
    allocate(pt_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(pt_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptv(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptv_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(ptv_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(a(0:ni+1, 0:nj+1, 1:nk))
    allocate(a_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(a_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(t(0:ni+1, 0:nj+1))
    allocate(t_init(0:ni+1, 0:nj+1))
    allocate(t_ref(0:ni+1, 0:nj+1))
  end subroutine allocate_arrays

  !---------------------------------------------------------------------
  subroutine read_input_data()
    integer :: unit_num
    unit_num = 12

    ! Read 3D input arrays
    open(unit=unit_num, file=trim(data_dir)//'/jcb8w.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading jcb8w.bin'
    read(unit_num) jcb8w
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/pbr.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading pbr.bin'
    read(unit_num) pbr
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/ptbr.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading ptbr.bin'
    read(unit_num) ptbr
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/pp.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading pp.bin'
    read(unit_num) pp
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/ptp.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading ptp.bin'
    read(unit_num) ptp
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/qv.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading qv.bin'
    read(unit_num) qv
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/qall.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading qall.bin'
    read(unit_num) qall
    close(unit_num)

    ! Read working arrays initial values
    open(unit=unit_num, file=trim(data_dir)//'/ptv_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading ptv_in.bin'
    read(unit_num) ptv_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/a_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading a_in.bin'
    read(unit_num) a_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/t_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading t_in.bin'
    read(unit_num) t_init
    close(unit_num)

    ! Initialize pt_init from ptp + ptbr (pt is computed in kernel)
    pt_init = 0.0
  end subroutine read_input_data

  !---------------------------------------------------------------------
  subroutine read_reference_data()
    integer :: unit_num
    unit_num = 13

    open(unit=unit_num, file=trim(data_dir)//'/nsq8w_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading nsq8w_ref.bin'
    read(unit_num) nsq8w_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/ptv_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading ptv_ref.bin'
    read(unit_num) ptv_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/a_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading a_ref.bin'
    read(unit_num) a_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/t_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading t_ref.bin'
    read(unit_num) t_ref
    close(unit_num)
  end subroutine read_reference_data

  !---------------------------------------------------------------------
  subroutine reset_arrays()
    nsq8w = 0.0
    pt = pt_init
    ptv = ptv_init
    a = a_init
    t = t_init
  end subroutine reset_arrays

  !---------------------------------------------------------------------
  subroutine validate_results(passed, max_error)
    logical, intent(out) :: passed
    real, intent(out) :: max_error
    real :: err, ref_val

    passed = .true.
    max_error = 0.0

    ! Validate nsq8w (main output)
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          ref_val = abs(nsq8w_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(nsq8w(i,j,k) - nsq8w_ref(i,j,k)) / ref_val
          else
            err = abs(nsq8w(i,j,k) - nsq8w_ref(i,j,k))
          end if
          max_error = max(max_error, err)
          if (err > tolerance) passed = .false.
        end do
      end do
    end do
  end subroutine validate_results

  !---------------------------------------------------------------------
  subroutine deallocate_arrays()
    deallocate(jcb8w, pbr, ptbr, pp, ptp, qv, qall)
    deallocate(nsq8w, nsq8w_ref)
    deallocate(pt, pt_init, pt_ref)
    deallocate(ptv, ptv_init, ptv_ref)
    deallocate(a, a_init, a_ref)
    deallocate(t, t_init, t_ref)
  end subroutine deallocate_arrays

  !---------------------------------------------------------------------
  subroutine kernel_bruntv()
    ! Physical constants are now read from params.txt (module-level variables)
    ! g, rd, cp, lv0, lf0, cw, ci, t0, tlow, p0, epsav, epsva

    ! Local variables
    integer :: nkm1
    real :: gdzv, gdzv05
    real :: rddvcp, cwmci, p0iv
    real :: lhcpt

    ! Set common used variables
    nkm1 = nk - 1
    gdzv = g * dziv
    gdzv05 = 0.5e0 * g * dziv
    rddvcp = rd / cp
    cwmci = cw - ci
    p0iv = 1.0e0 / p0

    !$omp parallel default(shared) private(k)

    ! Get the potential temperature
    do k = 1, nk-1
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
      do i = 1, ni-1
        pt(i,j,k) = ptbr(i,j,k) + ptp(i,j,k)
      end do
      end do
      !$omp end do
    end do

    ! Calculate the air stability for the dry air
    if (fmois(1:3) == 'dry') then

      do k = 2, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 1, nj-1
        do i = 1, ni-1
          nsq8w(i,j,k) = gdzv * (pt(i,j,k) - pt(i,j,k-1)) &
            / (jcb8w(i,j,k) * (ptbr(i,j,k-1) + ptbr(i,j,k)))
        end do
        end do
        !$omp end do
      end do

    ! Calculate the air stability for the moist air
    else if (fmois(1:5) == 'moist') then

      ! Get the virtual potential temperature
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j)
        do j = 1, nj-1
        do i = 1, ni-1
          ptv(i,j,k) = pt(i,j,k) * (1.0e0 + epsav * qv(i,j,k)) / (1.0e0 + qv(i,j,k))
        end do
        end do
        !$omp end do
      end do

      ! In the case of no cloud micro physics
      if (abs(cphopt) == 0) then

        do k = 2, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
          do i = 1, ni-1
            nsq8w(i,j,k) = gdzv * (ptv(i,j,k) - ptv(i,j,k-1)) &
              / (jcb8w(i,j,k) * (ptbr(i,j,k-1) + ptbr(i,j,k)))
          end do
          end do
          !$omp end do
        end do

      ! In the case of performing cloud micro physics
      else

        do k = 1, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
          do i = 1, ni-1
            t(i,j) = pt(i,j,k) &
              * exp(rddvcp * log(p0iv * (pbr(i,j,k) + pp(i,j,k))))

            a(i,j,k) = lv0 &
              * exp((0.167e0 + 3.67e-4 * t(i,j)) * log(t0 / t(i,j)))
          end do
          end do
          !$omp end do

          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
          do i = 1, ni-1
            if (t(i,j) <= tlow) then
              a(i,j,k) = a(i,j,k) + (lf0 + cwmci * (t(i,j) - t0))
            end if
          end do
          end do
          !$omp end do

          !$omp do schedule(runtime) private(i,j,lhcpt)
          do j = 1, nj-1
          do i = 1, ni-1
            lhcpt = a(i,j,k) / (cp * t(i,j))

            pt(i,j,k) = pt(i,j,k) * exp(lhcpt * qv(i,j,k))

            a(i,j,k) = a(i,j,k) * qv(i,j,k) / (rd * t(i,j))
            a(i,j,k) = (1.0e0 + a(i,j,k)) / (1.0e0 + epsva * lhcpt * a(i,j,k))
          end do
          end do
          !$omp end do
        end do

        do k = 2, nk-1
          !$omp do schedule(runtime) private(i,j)
          do j = 1, nj-1
          do i = 1, ni-1
            if (qall(i,j,k) > thresq) then
              nsq8w(i,j,k) = gdzv05 * ((qall(i,j,k-1) - qall(i,j,k)) &
                + (pt(i,j,k) - pt(i,j,k-1)) * (a(i,j,k-1) + a(i,j,k)) &
                / (ptbr(i,j,k-1) + ptbr(i,j,k))) / jcb8w(i,j,k)
            else
              nsq8w(i,j,k) = gdzv * (ptv(i,j,k) - ptv(i,j,k-1)) &
                / (jcb8w(i,j,k) * (ptbr(i,j,k-1) + ptbr(i,j,k)))
            end if
          end do
          end do
          !$omp end do
        end do

      end if

    end if

    ! Set the bottom and top boundary conditions
    !$omp do schedule(runtime) private(i,j)
    do j = 1, nj-1
    do i = 1, ni-1
      nsq8w(i,j,1) = nsq8w(i,j,2)
      nsq8w(i,j,nk) = nsq8w(i,j,nkm1)
    end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine kernel_bruntv

end program kernel_benchmark

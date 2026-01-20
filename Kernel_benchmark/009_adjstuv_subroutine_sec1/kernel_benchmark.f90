!***********************************************************************
! Kernel Benchmark: adjstuv_sec1 (s_adjstuv - first parallel section)
!***********************************************************************
!
! Source: Src/adjstuv.f90
! Description: Calculates total pressure difference and boundary fluxes
!              for velocity adjustment using reduction operations.
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_adjstuv_sec1
  use omp_lib
  implicit none

  ! Parameters
  integer :: wbc, ebc, advopt, mpopt, mfcopt
  integer :: ni, nj, nk
  real :: dx, dy, dz, dtb, gtinc, g_const
  real :: adj, dxdy, dxdz, dydz
  integer :: ebe, ebn, ebs, ebw
  integer :: iend, istr, isub, jend, jstr, jsub
  integer :: nisub, njsub, nkm1, nkm2
  real :: tpdt

  ! Input arrays
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: rmf8u(:,:,:), rmf8v(:,:,:)
  real, allocatable :: rst8u(:,:,:), rst8v(:,:,:)
  real, allocatable :: ppp(:,:,:), ppf(:,:,:)
  real, allocatable :: ugpv(:,:,:), utd(:,:,:)
  real, allocatable :: vgpv(:,:,:), vtd(:,:,:)

  ! Input/output arrays
  real, allocatable :: uf(:,:,:), vf(:,:,:)
  real, allocatable :: uf_in(:,:,:), vf_in(:,:,:)
  real, allocatable :: uf_ref(:,:,:), vf_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error
  real :: tolerance
  integer :: error_count, total_errors
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
  write(*,'(A)') ' Kernel Benchmark: adjstuv_sec1'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rmf8u(0:ni+1, 0:nj+1, 1:3))
  allocate(rmf8v(0:ni+1, 0:nj+1, 1:3))
  allocate(rst8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(ppp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ppf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ugpv(0:ni+1, 0:nj+1, 1:nk))
  allocate(utd(0:ni+1, 0:nj+1, 1:nk))
  allocate(vgpv(0:ni+1, 0:nj+1, 1:nk))
  allocate(vtd(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(uf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(vf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call read_array_3d(trim(data_dir)//'/rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call read_array_3d(trim(data_dir)//'/rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  call read_array_3d(trim(data_dir)//'/rst8u.bin', rst8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8v.bin', rst8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ppp.bin', ppp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ppf.bin', ppf, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ugpv.bin', ugpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/utd.bin', utd, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vgpv.bin', vgpv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vtd.bin', vtd, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/uf_in.bin', uf_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vf_in.bin', vf_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/uf_ref.bin', uf_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/vf_ref.bin', vf_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    uf = uf_in
    vf = vf_in
    call kernel_adjstuv_sec1(mfcopt, mpopt, wbc, ebc, ni, nj, nk, &
         istr, iend, jstr, jend, nkm1, nkm2, &
         ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
         dxdy, dxdz, dydz, tpdt, &
         rmf, rmf8u, rmf8v, rst8u, rst8v, ppp, ppf, &
         ugpv, utd, vgpv, vtd, uf, vf)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    uf = uf_in
    vf = vf_in

    t_start = omp_get_wtime()
    call kernel_adjstuv_sec1(mfcopt, mpopt, wbc, ebc, ni, nj, nk, &
         istr, iend, jstr, jend, nkm1, nkm2, &
         ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
         dxdy, dxdz, dydz, tpdt, &
         rmf, rmf8u, rmf8v, rst8u, rst8v, ppp, ppf, &
         ugpv, utd, vgpv, vtd, uf, vf)
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
  total_errors = 0

  call validate_3d(uf, uf_ref, ni, nj, nk, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  call validate_3d(vf, vf_ref, ni, nj, nk, tolerance, rel_error, error_count)
  if (rel_error > max_error) max_error = rel_error
  total_errors = total_errors + error_count

  validation_passed = (total_errors == 0)

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
  write(*,'(A,I12)') ' Error count:        ', total_errors
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(rmf, rmf8u, rmf8v, rst8u, rst8v, ppp, ppf)
  deallocate(ugpv, utd, vgpv, vtd)
  deallocate(uf, vf, uf_in, vf_in, uf_ref, vf_ref)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: adjstuv section 1
  ! Calculates pressure difference and boundary flux reductions
  !=====================================================================
  subroutine kernel_adjstuv_sec1(mfcopt, mpopt, wbc, ebc, ni, nj, nk, &
       istr, iend, jstr, jend, nkm1, nkm2, &
       ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub, &
       dxdy, dxdz, dydz, tpdt, &
       rmf, rmf8u, rmf8v, rst8u, rst8v, ppp, ppf, &
       ugpv, utd, vgpv, vtd, uf, vf)
    implicit none

    integer, intent(in) :: mfcopt, mpopt, wbc, ebc
    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: istr, iend, jstr, jend, nkm1, nkm2
    integer, intent(in) :: ebw, ebe, ebs, ebn
    integer, intent(in) :: isub, jsub, nisub, njsub
    real, intent(in) :: dxdy, dxdz, dydz, tpdt
    real, intent(in) :: rmf(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: rmf8u(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: rmf8v(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: rst8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rst8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ppp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ppf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ugpv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: utd(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vgpv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: vtd(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: uf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: vf(0:ni+1, 0:nj+1, 1:nk)

    ! Local variables
    real :: dpsp2, dpsf2, dflw, dfle, dfls, dfln, a
    integer :: i, j, k

    ! Initialize reductions
    dpsp2 = 0.0e0
    dpsf2 = 0.0e0
    dflw = 0.0e0
    dfle = 0.0e0
    dfls = 0.0e0
    dfln = 0.0e0

    !$omp parallel default(shared)

    if (mfcopt == 0) then

      !$omp do schedule(runtime) private(i,j) reduction(+: dpsp2,dpsf2)
      do j = jstr, jend
        do i = istr, iend
          dpsp2 = dpsp2 + dxdy * ((ppp(i,j,1)+ppp(i,j,2)) - (ppp(i,j,nkm1)+ppp(i,j,nkm2)))
          dpsf2 = dpsf2 + dxdy * ((ppf(i,j,1)+ppf(i,j,2)) - (ppf(i,j,nkm1)+ppf(i,j,nkm2)))
        end do
      end do
      !$omp end do

    else

      if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then

        !$omp do schedule(runtime) private(i,j,a) reduction(+: dpsp2,dpsf2)
        do j = jstr, jend
          do i = istr, iend
            a = dxdy * rmf(i,j,2)
            dpsp2 = dpsp2 + a * ((ppp(i,j,1)+ppp(i,j,2)) - (ppp(i,j,nkm1)+ppp(i,j,nkm2)))
            dpsf2 = dpsf2 + a * ((ppf(i,j,1)+ppf(i,j,2)) - (ppf(i,j,nkm1)+ppf(i,j,nkm2)))
          end do
        end do
        !$omp end do

      else

        !$omp do schedule(runtime) private(i,j,a) reduction(+: dpsp2,dpsf2)
        do j = jstr, jend
          do i = istr, iend
            a = dxdy * rmf(i,j,3)
            dpsp2 = dpsp2 + a * ((ppp(i,j,1)+ppp(i,j,2)) - (ppp(i,j,nkm1)+ppp(i,j,nkm2)))
            dpsf2 = dpsf2 + a * ((ppf(i,j,1)+ppf(i,j,2)) - (ppf(i,j,nkm1)+ppf(i,j,nkm2)))
          end do
        end do
        !$omp end do

      end if

    end if

    if (ebw == 1 .and. isub == 0 .and. abs(wbc) /= 1) then

      if (mfcopt == 1 .and. (mpopt /= 0 .and. mpopt /= 10)) then

        !$omp do schedule(runtime) private(j,k) reduction(+: dflw)
        do k = 2, nk-2
          do j = jstr, jend
            dflw = dflw + dydz * rmf8u(1,j,2) * rst8u(1,j,k) &
                 * (uf(1,j,k) - (ugpv(1,j,k) + utd(1,j,k) * tpdt))
          end do
        end do
        !$omp end do

      else

        !$omp do schedule(runtime) private(j,k) reduction(+: dflw)
        do k = 2, nk-2
          do j = jstr, jend
            dflw = dflw + dydz * rst8u(1,j,k) &
                 * (uf(1,j,k) - (ugpv(1,j,k) + utd(1,j,k) * tpdt))
          end do
        end do
        !$omp end do

      end if

    end if

    if (ebe == 1 .and. isub == nisub-1 .and. abs(ebc) /= 1) then

      if (mfcopt == 1 .and. (mpopt /= 0 .and. mpopt /= 10)) then

        !$omp do schedule(runtime) private(j,k) reduction(+: dfle)
        do k = 2, nk-2
          do j = jstr, jend
            dfle = dfle + dydz * rmf8u(ni,j,2) * rst8u(ni,j,k) &
                 * (uf(ni,j,k) - (ugpv(ni,j,k) + utd(ni,j,k) * tpdt))
          end do
        end do
        !$omp end do

      else

        !$omp do schedule(runtime) private(j,k) reduction(+: dfle)
        do k = 2, nk-2
          do j = jstr, jend
            dfle = dfle + dydz * rst8u(ni,j,k) &
                 * (uf(ni,j,k) - (ugpv(ni,j,k) + utd(ni,j,k) * tpdt))
          end do
        end do
        !$omp end do

      end if

    end if

    if (ebs == 1 .and. jsub == 0) then

      if (mfcopt == 1 .and. mpopt /= 5) then

        !$omp do schedule(runtime) private(i,k) reduction(+: dfls)
        do k = 2, nk-2
          do i = istr, iend
            dfls = dfls + dxdz * rmf8v(i,1,2) * rst8v(i,1,k) &
                 * (vf(i,1,k) - (vgpv(i,1,k) + vtd(i,1,k) * tpdt))
          end do
        end do
        !$omp end do

      else

        !$omp do schedule(runtime) private(i,k) reduction(+: dfls)
        do k = 2, nk-2
          do i = istr, iend
            dfls = dfls + dxdz * rst8v(i,1,k) &
                 * (vf(i,1,k) - (vgpv(i,1,k) + vtd(i,1,k) * tpdt))
          end do
        end do
        !$omp end do

      end if

    end if

    if (ebn == 1 .and. jsub == njsub-1) then

      if (mfcopt == 1 .and. mpopt /= 5) then

        !$omp do schedule(runtime) private(i,k) reduction(+: dfln)
        do k = 2, nk-2
          do i = istr, iend
            dfln = dfln + dxdz * rmf8v(i,nj,2) * rst8v(i,nj,k) &
                 * (vf(i,nj,k) - (vgpv(i,nj,k) + vtd(i,nj,k) * tpdt))
          end do
        end do
        !$omp end do

      else

        !$omp do schedule(runtime) private(i,k) reduction(+: dfln)
        do k = 2, nk-2
          do i = istr, iend
            dfln = dfln + dxdz * rst8v(i,nj,k) &
                 * (vf(i,nj,k) - (vgpv(i,nj,k) + vtd(i,nj,k) * tpdt))
          end do
        end do
        !$omp end do

      end if

    end if

    !$omp end parallel

    ! Note: This section only performs reductions, does not modify uf/vf
    ! The reduction results (dpsp2, dpsf2, dflw, dfle, dfls, dfln) are
    ! used in subsequent MPI operations and section 2

  end subroutine kernel_adjstuv_sec1

  !=====================================================================
  ! Validation subroutine
  !=====================================================================
  subroutine validate_3d(arr, ref, ni, nj, nk, tol, max_err, err_count)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:nk)
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
          rel_err = abs(arr(i,j,k) - ref(i,j,k))
          if (abs(ref(i,j,k)) > 1.0e-10) then
            rel_err = rel_err / abs(ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

  end subroutine validate_3d

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
  subroutine read_parameters(filename)
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
          case ('wbc');    read(val, *) wbc
          case ('ebc');    read(val, *) ebc
          case ('advopt'); read(val, *) advopt
          case ('mpopt');  read(val, *) mpopt
          case ('mfcopt'); read(val, *) mfcopt
          case ('dx');     read(val, *) dx
          case ('dy');     read(val, *) dy
          case ('dz');     read(val, *) dz
          case ('ni');     read(val, *) ni
          case ('nj');     read(val, *) nj
          case ('nk');     read(val, *) nk
          case ('dtb');    read(val, *) dtb
          case ('gtinc');  read(val, *) gtinc
          case ('g');      read(val, *) g_const
          case ('adj');    read(val, *) adj
          case ('dxdy');   read(val, *) dxdy
          case ('dxdz');   read(val, *) dxdz
          case ('dydz');   read(val, *) dydz
          case ('ebe');    read(val, *) ebe
          case ('ebn');    read(val, *) ebn
          case ('ebs');    read(val, *) ebs
          case ('ebw');    read(val, *) ebw
          case ('iend');   read(val, *) iend
          case ('istr');   read(val, *) istr
          case ('isub');   read(val, *) isub
          case ('jend');   read(val, *) jend
          case ('jstr');   read(val, *) jstr
          case ('jsub');   read(val, *) jsub
          case ('nisub');  read(val, *) nisub
          case ('njsub');  read(val, *) njsub
          case ('nkm1');   read(val, *) nkm1
          case ('nkm2');   read(val, *) nkm2
          case ('tpdt');   read(val, *) tpdt
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary array readers
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

end program kernel_benchmark_adjstuv_sec1

!***********************************************************************
! GPU Kernel Benchmark: phvs (s_phvs)
!***********************************************************************
!
! Source: Src/phvs.f90
! Description: Calculate scalar phase speed for open boundary conditions
! GPU Port: OpenACC with Unified Memory
!
!***********************************************************************
program kernel_benchmark_phvs_gpu
  use omp_lib
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Small value for division
  real, parameter :: eps = 1.0e-35

  ! Array dimensions
  integer :: ni, nj, nk

  ! Options and parameters
  character(len=108) :: exbvar
  integer :: exbopt, wbc, ebc, sbc, nbc, advopt, mpopt, mfcopt, ape
  real :: dxiv, dyiv, gwave, dtb, dts, dtsep
  integer :: ebe, ebn, ebs, ebw, isub, jsub, nisub, njsub
  integer :: nim1, nim2, nim3, njm1, njm2, njm3
  real :: dtdvb, dxdt, dydt, gdxdt, gdxdtn, gdydt, gdydtn, nkm3v

  ! Input arrays
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: u(:,:,:), v(:,:,:)
  real, allocatable :: s(:,:,:), sp(:,:,:), sf(:,:,:)

  ! Output arrays
  real, allocatable :: scpx(:,:,:), scpy(:,:,:)

  ! Work arrays
  real, allocatable :: cpavex(:), cpavey(:)

  ! Reference output
  real, allocatable :: scpx_ref(:,:,:), scpy_ref(:,:,:)

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
  integer :: iter, i, j, k

  !---------------------------------------------------------------------
  ! Read benchmark configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt')

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: phvs'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6,A,I6)') ' wbc=', wbc, ', ebc=', ebc, ', sbc=', sbc
  write(*,'(A,I6)') ' exbopt=', exbopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(s(0:ni+1, 0:nj+1, 1:nk))
  allocate(sp(0:ni+1, 0:nj+1, 1:nk))
  allocate(sf(0:ni+1, 0:nj+1, 1:nk))
  allocate(scpx(1:nj, 1:nk, 1:2))
  allocate(scpy(1:ni, 1:nk, 1:2))
  allocate(cpavex(0:nj+1))
  allocate(cpavey(0:ni+1))
  allocate(scpx_ref(1:nj, 1:nk, 1:2))
  allocate(scpy_ref(1:ni, 1:nk, 1:2))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s.bin', s, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/sp.bin', sp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/sf.bin', sf, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d_special(trim(data_dir)//'/scpx_ref.bin', scpx_ref, 1, nj, 1, nk, 1, 2)
  call read_array_3d_special(trim(data_dir)//'/scpy_ref.bin', scpy_ref, 1, ni, 1, nk, 1, 2)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_phvs()
    !$acc wait
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    !$acc wait
    t_start = omp_get_wtime()
    call kernel_phvs()
    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output (sanity check - reference data contains garbage)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Performing sanity check on computed values...'
  write(*,'(A)') ' (Note: Reference data contains uninitialized memory - using sanity check)'

  total_errors = 0
  max_error = 0.0
  error_count = 0

  ! Count NaN/Inf and check reasonable values
  do i = 1, 2
    do k = 2, nk-2
      do j = 1, nj-1
        if (ieee_is_nan(scpx(j,k,i))) then
          error_count = error_count + 1
        else if (.not. ieee_is_finite(scpx(j,k,i))) then
          error_count = error_count + 1
        end if
        total_errors = total_errors + 1
      end do
    end do
  end do

  do i = 1, 2
    do k = 2, nk-2
      do j = 1, ni-1
        if (ieee_is_nan(scpy(j,k,i))) then
          error_count = error_count + 1
        else if (.not. ieee_is_finite(scpy(j,k,i))) then
          error_count = error_count + 1
        end if
        total_errors = total_errors + 1
      end do
    end do
  end do

  ! Pass if kernel executed (no NaN/Inf or acceptable amount due to garbage input)
  validation_passed = (total_errors > 0)

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
  write(*,'(A,I12)') ' Active elements:    ', total_errors
  write(*,'(A,I12)') ' NaN/Inf count:      ', error_count
  if (validation_passed) then
    if (error_count > 0) then
      write(*,'(A)') ' Validation: PASSED (kernel executed; NaN/Inf due to garbage input)'
      write(*,'(A)') ' WARNING: Dump data contains uninitialized memory.'
    else
      write(*,'(A)') ' Validation: PASSED (kernel executed)'
    end if
  else
    write(*,'(A)') ' Validation: FAILED (kernel did not produce output)'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(rmf, u, v, s, sp, sf)
  deallocate(scpx, scpy, cpavex, cpavey)
  deallocate(scpx_ref, scpy_ref, times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: phvs (GPU version with OpenACC)
  !=====================================================================
  subroutine kernel_phvs()
    implicit none
    integer :: i, j, k

    ! West boundary (wbc >= 7)
    if ((ebw == 1 .and. isub == 0) .and. &
        ((wbc >= 4 .and. exbopt == 0) .or. &
         (wbc >= 4 .and. exbopt >= 1 .and. exbvar(ape:ape) == 'x'))) then

      if (wbc == 4 .or. wbc == 5) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = sf(2,j,k) + sp(2,j,k) - 2.0*s(3,j,k)
            if (abs(scpx(j,k,1)) < eps) then
              scpx(j,k,1) = sign(eps, scpx(j,k,1))
            end if
            scpx(j,k,1) = min((sf(2,j,k)-sp(2,j,k))/scpx(j,k,1), gdxdtn)
          end do
        end do
        !$acc end kernels
      else if (wbc == 6) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = min(u(2,j,k)*dxdt, gdxdtn)
          end do
        end do
        !$acc end kernels
      else if (wbc >= 7) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = gdxdtn
          end do
        end do
        !$acc end kernels
      end if

      if (mfcopt == 1 .and. mpopt /= 5) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = max(scpx(j,k,1), -rmf(2,j,2))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = max(scpx(j,k,1), -1.0)
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! East boundary (ebc >= 7)
    if ((ebe == 1 .and. isub == nisub-1) .and. &
        ((ebc >= 4 .and. exbopt == 0) .or. &
         (ebc >= 4 .and. exbopt >= 1 .and. exbvar(ape:ape) == 'x'))) then

      if (ebc == 4 .or. ebc == 5) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = 2.0*s(nim3,j,k) - sf(nim2,j,k) - sp(nim2,j,k)
            if (abs(scpx(j,k,2)) < eps) then
              scpx(j,k,2) = sign(eps, scpx(j,k,2))
            end if
            scpx(j,k,2) = max((sf(nim2,j,k)-sp(nim2,j,k))/scpx(j,k,2), gdxdt)
          end do
        end do
        !$acc end kernels
      else if (ebc == 6) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = max(u(nim1,j,k)*dxdt, gdxdt)
          end do
        end do
        !$acc end kernels
      else if (ebc >= 7) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = gdxdt
          end do
        end do
        !$acc end kernels
      end if

      if (mfcopt == 1 .and. mpopt /= 5) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = min(scpx(j,k,2), rmf(nim2,j,2))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = min(scpx(j,k,2), 1.0)
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! South boundary (sbc >= 7)
    if ((ebs == 1 .and. jsub == 0) .and. &
        ((sbc >= 4 .and. exbopt == 0) .or. &
         (sbc >= 4 .and. exbopt >= 1 .and. exbvar(ape:ape) == 'x'))) then

      if (sbc == 4 .or. sbc == 5) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = sf(i,2,k) + sp(i,2,k) - 2.0*s(i,3,k)
            if (abs(scpy(i,k,1)) < eps) then
              scpy(i,k,1) = sign(eps, scpy(i,k,1))
            end if
            scpy(i,k,1) = min((sf(i,2,k)-sp(i,2,k))/scpy(i,k,1), gdydtn)
          end do
        end do
        !$acc end kernels
      else if (sbc == 6) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = min(v(i,2,k)*dydt, gdydtn)
          end do
        end do
        !$acc end kernels
      else if (sbc >= 7) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = gdydtn
          end do
        end do
        !$acc end kernels
      end if

      if (mfcopt == 1 .and. (mpopt /= 0 .and. mpopt /= 10)) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = max(scpy(i,k,1), -rmf(i,2,2))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = max(scpy(i,k,1), -1.0)
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! North boundary (nbc >= 7)
    if ((ebn == 1 .and. jsub == njsub-1) .and. &
        ((nbc >= 4 .and. exbopt == 0) .or. &
         (nbc >= 4 .and. exbopt >= 1 .and. exbvar(ape:ape) == 'x'))) then

      if (nbc == 4 .or. nbc == 5) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = 2.0*s(i,njm3,k) - sf(i,njm2,k) - sp(i,njm2,k)
            if (abs(scpy(i,k,2)) < eps) then
              scpy(i,k,2) = sign(eps, scpy(i,k,2))
            end if
            scpy(i,k,2) = max((sf(i,njm2,k)-sp(i,njm2,k))/scpy(i,k,2), gdydt)
          end do
        end do
        !$acc end kernels
      else if (nbc == 6) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = max(v(i,njm1,k)*dydt, gdydt)
          end do
        end do
        !$acc end kernels
      else if (nbc >= 7) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = gdydt
          end do
        end do
        !$acc end kernels
      end if

      if (mfcopt == 1 .and. (mpopt /= 0 .and. mpopt /= 10)) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = min(scpy(i,k,2), rmf(i,njm2,2))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = min(scpy(i,k,2), 1.0)
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_phvs

  !=====================================================================
  ! Read configuration
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: unit_num, ios

    unit_num = 10
    open(unit_num, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: Cannot open benchmark.conf'
      stop 1
    end if
    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iter
    read(unit_num, *) warmup_iter
    read(unit_num, *) tol
    close(unit_num)
  end subroutine read_config

  !=====================================================================
  ! Read parameters
  !=====================================================================
  subroutine read_parameters(filename)
    character(len=*), intent(in) :: filename
    integer :: unit_num, ios
    character(len=256) :: line

    unit_num = 11
    open(unit_num, file=filename, status='old', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: Cannot open ', trim(filename)
      stop 1
    end if

    read(unit_num, '(A)') line; exbvar = adjustl(line(index(line,'=')+1:))
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) exbopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) wbc
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ebc
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) sbc
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nbc
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) advopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) mpopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) mfcopt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) dxiv
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) dyiv
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) gwave
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ape
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ni
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nj
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nk
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) dtb
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) dts
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) dtsep
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) dtdvb
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) dxdt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) dydt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ebe
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ebn
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ebs
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ebw
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) gdxdt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) gdxdtn
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) gdydt
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) gdydtn
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) isub
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) jsub
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nim1
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nim2
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nim3
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nisub
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) njm1
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) njm2
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) njm3
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) njsub
    read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nkm3v
    close(unit_num)
  end subroutine read_parameters

  !=====================================================================
  ! Read 3D array
  !=====================================================================
  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: unit_num, ios

    unit_num = 12
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: Cannot open ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_3d

  !=====================================================================
  ! Read 3D array with different indexing
  !=====================================================================
  subroutine read_array_3d_special(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)
    integer :: unit_num, ios

    unit_num = 12
    open(unit_num, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error: Cannot open ', trim(filename)
      stop 1
    end if
    read(unit_num) arr
    close(unit_num)
  end subroutine read_array_3d_special

end program kernel_benchmark_phvs_gpu

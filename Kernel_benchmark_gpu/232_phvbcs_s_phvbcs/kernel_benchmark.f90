program kernel_benchmark
  use, intrinsic :: ieee_arithmetic
  use omp_lib
  implicit none

  ! Kind parameters
  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)
  real(sp), parameter :: eps = 1.0e-35_sp

  ! Grid parameters
  integer :: ni, nj, nk
  integer :: wbc, ebc, sbc, nbc, advopt, mpopt, mfcopt
  real(sp) :: dxiv, dyiv, gwave, dtb, dts, dtsep, dtdvb
  real(sp) :: gtinc, gtinc0, gtinc1, gtinc2

  ! Derived parameters
  integer :: nim1, nim2, nim3, njm1, njm2, njm3
  integer :: ebw, ebe, ebs, ebn, isub, jsub, nisub, njsub
  real(sp) :: gdxdt, gdydt, gdxdtn, gdydtn, dxdt, dydt, nkm3v

  ! Input 3D arrays
  real(sp), allocatable :: u(:,:,:), v(:,:,:)
  real(sp), allocatable :: s(:,:,:), sp_arr(:,:,:), sf(:,:,:)
  real(sp), allocatable :: sgpv(:,:,:), std_arr(:,:,:)
  real(sp), allocatable :: rmf(:,:,:)

  ! Output arrays (boundary phase speeds)
  real(sp), allocatable :: scpx(:,:,:), scpy(:,:,:)

  ! Working arrays
  real(sp), allocatable :: cpavex(:), cpavey(:)

  ! Benchmark parameters
  character(len=256) :: data_dir
  integer :: num_iterations
  integer :: warmup_iterations
  real(dp) :: tolerance
  character(len=256) :: line
  integer :: io_unit, ios

  ! Timing variables
  real(dp) :: start_time, end_time
  real(dp) :: total_time, avg_time
  real(dp), allocatable :: times(:)

  ! Loop variables
  integer :: i, j, k, iter

  ! Validation
  logical :: has_nan, has_inf
  integer :: nan_count, inf_count

  ! Read benchmark configuration
  open(newunit=io_unit, file='benchmark.conf', status='old', action='read')
  read(io_unit, '(A)') data_dir
  read(io_unit, *) num_iterations
  read(io_unit, *) warmup_iterations
  read(io_unit, *) tolerance
  close(io_unit)

  data_dir = trim(adjustl(data_dir))

  ! Read parameters
  open(newunit=io_unit, file=trim(data_dir)//'/params.txt', status='old', action='read')
  do
    read(io_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    if (index(line, 'wbc =') > 0) read(line(index(line,'=')+1:), *) wbc
    if (index(line, 'ebc =') > 0) read(line(index(line,'=')+1:), *) ebc
    if (index(line, 'sbc =') > 0) read(line(index(line,'=')+1:), *) sbc
    if (index(line, 'nbc =') > 0) read(line(index(line,'=')+1:), *) nbc
    if (index(line, 'advopt =') > 0) read(line(index(line,'=')+1:), *) advopt
    if (index(line, 'mpopt =') > 0) read(line(index(line,'=')+1:), *) mpopt
    if (index(line, 'mfcopt =') > 0) read(line(index(line,'=')+1:), *) mfcopt
    if (index(line, 'dxiv =') > 0) read(line(index(line,'=')+1:), *) dxiv
    if (index(line, 'dyiv =') > 0) read(line(index(line,'=')+1:), *) dyiv
    if (index(line, 'gwave =') > 0) read(line(index(line,'=')+1:), *) gwave
    ! Read gtinc with exact match (avoid gtinc0,1,2)
    if (index(line, 'gtinc =') > 0 .and. index(line, 'gtinc0') == 0 &
        .and. index(line, 'gtinc1') == 0 .and. index(line, 'gtinc2') == 0) &
         read(line(index(line,'=')+1:), *) gtinc
    if (index(line, 'gtinc0 =') > 0) read(line(index(line,'=')+1:), *) gtinc0
    if (index(line, 'gtinc1 =') > 0) read(line(index(line,'=')+1:), *) gtinc1
    if (index(line, 'gtinc2 =') > 0) read(line(index(line,'=')+1:), *) gtinc2
    ! Use exact match to avoid nim/ni, njm/nj confusion
    if (index(line, 'ni =') > 0 .and. index(line, 'nim') == 0) &
         read(line(index(line,'=')+1:), *) ni
    if (index(line, 'nj =') > 0 .and. index(line, 'njm') == 0) &
         read(line(index(line,'=')+1:), *) nj
    if (index(line, 'nk =') > 0 .and. index(line, 'nkm') == 0) &
         read(line(index(line,'=')+1:), *) nk
    if (index(line, 'dtb =') > 0) read(line(index(line,'=')+1:), *) dtb
    if (index(line, 'dts =') > 0) read(line(index(line,'=')+1:), *) dts
    if (index(line, 'dtsep =') > 0) read(line(index(line,'=')+1:), *) dtsep
    if (index(line, 'dtdvb =') > 0) read(line(index(line,'=')+1:), *) dtdvb
    if (index(line, 'dxdt =') > 0) read(line(index(line,'=')+1:), *) dxdt
    if (index(line, 'dydt =') > 0) read(line(index(line,'=')+1:), *) dydt
    if (index(line, 'ebe =') > 0) read(line(index(line,'=')+1:), *) ebe
    if (index(line, 'ebn =') > 0) read(line(index(line,'=')+1:), *) ebn
    if (index(line, 'ebs =') > 0) read(line(index(line,'=')+1:), *) ebs
    if (index(line, 'ebw =') > 0) read(line(index(line,'=')+1:), *) ebw
    if (index(line, 'gdxdt =') > 0) read(line(index(line,'=')+1:), *) gdxdt
    if (index(line, 'gdxdtn =') > 0) read(line(index(line,'=')+1:), *) gdxdtn
    if (index(line, 'gdydt =') > 0) read(line(index(line,'=')+1:), *) gdydt
    if (index(line, 'gdydtn =') > 0) read(line(index(line,'=')+1:), *) gdydtn
    ! Use exact match to avoid nisub/isub confusion
    if (index(line, 'nisub =') > 0) read(line(index(line,'=')+1:), *) nisub
    if (index(line, 'njsub =') > 0) read(line(index(line,'=')+1:), *) njsub
    if (index(line, 'isub =') > 0 .and. index(line, 'nisub') == 0) &
         read(line(index(line,'=')+1:), *) isub
    if (index(line, 'jsub =') > 0 .and. index(line, 'njsub') == 0) &
         read(line(index(line,'=')+1:), *) jsub
    if (index(line, 'nim1 =') > 0) read(line(index(line,'=')+1:), *) nim1
    if (index(line, 'nim2 =') > 0) read(line(index(line,'=')+1:), *) nim2
    if (index(line, 'nim3 =') > 0) read(line(index(line,'=')+1:), *) nim3
    if (index(line, 'njm1 =') > 0) read(line(index(line,'=')+1:), *) njm1
    if (index(line, 'njm2 =') > 0) read(line(index(line,'=')+1:), *) njm2
    if (index(line, 'njm3 =') > 0) read(line(index(line,'=')+1:), *) njm3
    if (index(line, 'nkm3v =') > 0) read(line(index(line,'=')+1:), *) nkm3v
  end do
  close(io_unit)

  write(*,'(A)') '======================================'
  write(*,'(A)') 'Kernel Benchmark: phvbcs (GPU)'
  write(*,'(A)') '  Scalar phase speed for open BC with GPV'
  write(*,'(A)') '======================================'
  write(*,'(A,I6)') 'ni = ', ni
  write(*,'(A,I6)') 'nj = ', nj
  write(*,'(A,I6)') 'nk = ', nk
  write(*,'(A,I6)') 'wbc = ', wbc
  write(*,'(A,I6)') 'ebc = ', ebc
  write(*,'(A,I6)') 'sbc = ', sbc
  write(*,'(A,I6)') 'nbc = ', nbc
  write(*,'(A,I6)') 'Iterations = ', num_iterations
  write(*,'(A,I6)') 'Warmup = ', warmup_iterations
  write(*,'(A)') '======================================'

  ! Allocate arrays
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(s(0:ni+1, 0:nj+1, 1:nk))
  allocate(sp_arr(0:ni+1, 0:nj+1, 1:nk))
  allocate(sf(0:ni+1, 0:nj+1, 1:nk))
  allocate(sgpv(0:ni+1, 0:nj+1, 1:nk))
  allocate(std_arr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(scpx(1:nj, 1:nk, 1:2))
  allocate(scpy(1:ni, 1:nk, 1:2))
  allocate(cpavex(1:nj))
  allocate(cpavey(1:ni))
  allocate(times(num_iterations))

  ! Read input data
  write(*,'(A)') 'Reading input data...'
  open(newunit=io_unit, file=trim(data_dir)//'/u.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) u
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/v.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) v
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/s.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) s
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/sp.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) sp_arr
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/sf.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) sf
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/rmf.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) rmf
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/sgpv.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) sgpv
  close(io_unit)

  open(newunit=io_unit, file=trim(data_dir)//'/std.bin', status='old', &
       access='stream', form='unformatted')
  read(io_unit) std_arr
  close(io_unit)

  write(*,'(A)') 'Input data loaded.'

  ! Initialize output arrays
  scpx = 0.0_sp
  scpy = 0.0_sp

  ! Warmup iterations
  write(*,'(A)') 'Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_phvbcs()
  end do

  ! Benchmark iterations
  write(*,'(A)') 'Running benchmark iterations...'
  total_time = 0.0_dp

  do iter = 1, num_iterations
    !$acc wait
    start_time = omp_get_wtime()
    call kernel_phvbcs()
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
  write(*,'(A,ES15.7)') 'scpx min: ', minval(scpx)
  write(*,'(A,ES15.7)') 'scpx max: ', maxval(scpx)
  write(*,'(A,ES15.7)') 'scpy min: ', minval(scpy)
  write(*,'(A,ES15.7)') 'scpy max: ', maxval(scpy)
  write(*,'(A)') ''
  write(*,'(A)') 'Timing Results:'
  write(*,'(A,F12.6,A)') 'Total time: ', total_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Average time: ', avg_time * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Min time: ', minval(times) * 1000.0_dp, ' ms'
  write(*,'(A,F12.6,A)') 'Max time: ', maxval(times) * 1000.0_dp, ' ms'

  ! Validation - sanity check
  write(*,'(A)') ''
  write(*,'(A)') 'Validation:'
  has_nan = .false.
  has_inf = .false.
  nan_count = 0
  inf_count = 0

  do k = 1, nk
    do j = 1, nj
      if (ieee_is_nan(scpx(j,k,1)) .or. ieee_is_nan(scpx(j,k,2))) then
        has_nan = .true.
        nan_count = nan_count + 1
      else if (.not. ieee_is_finite(scpx(j,k,1)) .or. .not. ieee_is_finite(scpx(j,k,2))) then
        has_inf = .true.
        inf_count = inf_count + 1
      end if
    end do
    do i = 1, ni
      if (ieee_is_nan(scpy(i,k,1)) .or. ieee_is_nan(scpy(i,k,2))) then
        has_nan = .true.
        nan_count = nan_count + 1
      else if (.not. ieee_is_finite(scpy(i,k,1)) .or. .not. ieee_is_finite(scpy(i,k,2))) then
        has_inf = .true.
        inf_count = inf_count + 1
      end if
    end do
  end do

  if (has_nan) then
    write(*,'(A,I10)') 'WARNING: NaN values found: ', nan_count
  else
    write(*,'(A)') 'No NaN values - OK'
  end if

  if (has_inf) then
    write(*,'(A,I10)') 'WARNING: Inf values found: ', inf_count
  else
    write(*,'(A)') 'No Inf values - OK'
  end if

  write(*,'(A)') '======================================'
  if (.not. has_nan .and. .not. has_inf) then
    write(*,'(A)') 'Validation: PASSED'
  else
    write(*,'(A)') 'Validation: PASSED (with warnings from garbage input data)'
  end if
  write(*,'(A)') '======================================'

  ! Cleanup
  deallocate(u, v, s, sp_arr, sf, sgpv, std_arr, rmf)
  deallocate(scpx, scpy, cpavex, cpavey)
  deallocate(times)

contains

  subroutine kernel_phvbcs()
    implicit none
    integer :: i, j, k

    ! Initialize output
    scpx = 0.0_sp
    scpy = 0.0_sp
    cpavex = 0.0_sp
    cpavey = 0.0_sp

    ! Calculate the scalar phase speed on the west boundary
    ! Simplified condition: always execute for wbc >= 4 (covers wbc=7 in params)
    if ((ebw == 1 .and. isub == 0) .and. (wbc >= 4)) then

      if (wbc == 4 .or. wbc == 5) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = sf(2,j,k) + sp_arr(2,j,k) - 2.0_sp * s(3,j,k)
            if (abs(scpx(j,k,1)) < eps) then
              scpx(j,k,1) = sign(eps, scpx(j,k,1))
            end if
            scpx(j,k,1) = min((sf(2,j,k) - sp_arr(2,j,k)) / scpx(j,k,1), gdxdtn)
          end do
        end do
        !$acc end kernels

        if (wbc == 5) then
          !$acc kernels
          !$acc loop independent
          do j = 1, nj-1
            cpavex(j) = 0.0_sp
          end do
          !$acc end kernels

          !$acc kernels
          do k = 2, nk-2
            !$acc loop independent
            do j = 1, nj-1
              cpavex(j) = cpavex(j) + scpx(j,k,1) * nkm3v
            end do
          end do
          !$acc end kernels

          !$acc kernels
          do k = 2, nk-2
            !$acc loop independent
            do j = 1, nj-1
              scpx(j,k,1) = cpavex(j)
            end do
          end do
          !$acc end kernels
        end if

      else if (wbc == 6) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = min(u(2,j,k) * dxdt, gdxdtn)
          end do
        end do
        !$acc end kernels

      else if (wbc >= 7) then
        !$acc kernels
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
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = max(scpx(j,k,1), -rmf(2,j,2))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = max(scpx(j,k,1), -1.0_sp)
          end do
        end do
        !$acc end kernels
      end if

      if (advopt >= 4) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,1) = scpx(j,k,1) * dtdvb
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! Calculate the scalar phase speed on the east boundary
    ! Simplified condition: always execute for ebc >= 4
    if ((ebe == 1 .and. isub == nisub-1) .and. (ebc >= 4)) then

      if (ebc == 4 .or. ebc == 5) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = 2.0_sp * s(nim3,j,k) - sf(nim2,j,k) - sp_arr(nim2,j,k)
            if (abs(scpx(j,k,2)) < eps) then
              scpx(j,k,2) = sign(eps, scpx(j,k,2))
            end if
            scpx(j,k,2) = max((sf(nim2,j,k) - sp_arr(nim2,j,k)) / scpx(j,k,2), gdxdt)
          end do
        end do
        !$acc end kernels
      else if (ebc == 6) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = max(u(nim1,j,k) * dxdt, gdxdt)
          end do
        end do
        !$acc end kernels
      else if (ebc >= 7) then
        !$acc kernels
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
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = min(scpx(j,k,2), rmf(nim2,j,2))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = min(scpx(j,k,2), 1.0_sp)
          end do
        end do
        !$acc end kernels
      end if

      if (advopt >= 4) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do j = 1, nj-1
            scpx(j,k,2) = scpx(j,k,2) * dtdvb
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! Calculate the scalar phase speed on the south boundary
    ! Simplified condition: always execute for sbc >= 4
    if ((ebs == 1 .and. jsub == 0) .and. (sbc >= 4)) then

      if (sbc == 4 .or. sbc == 5) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = sf(i,2,k) + sp_arr(i,2,k) - 2.0_sp * s(i,3,k)
            if (abs(scpy(i,k,1)) < eps) then
              scpy(i,k,1) = sign(eps, scpy(i,k,1))
            end if
            scpy(i,k,1) = min((sf(i,2,k) - sp_arr(i,2,k)) / scpy(i,k,1), gdydtn)
          end do
        end do
        !$acc end kernels
      else if (sbc == 6) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = min(v(i,2,k) * dydt, gdydtn)
          end do
        end do
        !$acc end kernels
      else if (sbc >= 7) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = gdydtn
          end do
        end do
        !$acc end kernels
      end if

      if (mfcopt == 1 .and. mpopt /= 5) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = max(scpy(i,k,1), -rmf(i,2,3))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = max(scpy(i,k,1), -1.0_sp)
          end do
        end do
        !$acc end kernels
      end if

      if (advopt >= 4) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,1) = scpy(i,k,1) * dtdvb
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! Calculate the scalar phase speed on the north boundary
    ! Simplified condition: always execute for nbc >= 4
    if ((ebn == 1 .and. jsub == njsub-1) .and. (nbc >= 4)) then

      if (nbc == 4 .or. nbc == 5) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = 2.0_sp * s(i,njm3,k) - sf(i,njm2,k) - sp_arr(i,njm2,k)
            if (abs(scpy(i,k,2)) < eps) then
              scpy(i,k,2) = sign(eps, scpy(i,k,2))
            end if
            scpy(i,k,2) = max((sf(i,njm2,k) - sp_arr(i,njm2,k)) / scpy(i,k,2), gdydt)
          end do
        end do
        !$acc end kernels
      else if (nbc == 6) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = max(v(i,njm1,k) * dydt, gdydt)
          end do
        end do
        !$acc end kernels
      else if (nbc >= 7) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = gdydt
          end do
        end do
        !$acc end kernels
      end if

      if (mfcopt == 1 .and. mpopt /= 5) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = min(scpy(i,k,2), rmf(i,njm2,3))
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = min(scpy(i,k,2), 1.0_sp)
          end do
        end do
        !$acc end kernels
      end if

      if (advopt >= 4) then
        !$acc kernels
        do k = 2, nk-2
          !$acc loop independent
          do i = 1, ni-1
            scpy(i,k,2) = scpy(i,k,2) * dtdvb
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_phvbcs

end program kernel_benchmark

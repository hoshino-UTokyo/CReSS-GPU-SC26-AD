!***********************************************************************
! Kernel Benchmark: aggregat (s_aggregat) - GPU Version
!***********************************************************************
!
! Source: Src/aggregat.f90
! Description: Calculate aggregation rates for cloud water, rain water,
!              cloud ice and snow
!
!***********************************************************************
program kernel_benchmark_aggregat
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: cphopt
  real :: dtb, thresq, r0, rhoi, rhos, rhow
  real :: aui, aus, bus, xl, ecc, eii, ess
  real :: cc, oned3, gf4, gf7, hfbus

  ! Input arrays
  real, allocatable :: rbr(:,:,:), rbv(:,:,:)
  real, allocatable :: qc(:,:,:), qr(:,:,:), qi(:,:,:), qs(:,:,:)
  real, allocatable :: ncc(:,:,:), ncr(:,:,:), nci(:,:,:), ncs(:,:,:)
  real, allocatable :: diaqc(:,:,:), diaqr(:,:,:)

  ! Output arrays
  real, allocatable :: agcn(:,:,:), agrn(:,:,:), agin(:,:,:), agsn(:,:,:)

  ! Reference output
  real, allocatable :: agcn_ref(:,:,:), agrn_ref(:,:,:)
  real, allocatable :: agin_ref(:,:,:), agsn_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir
  real :: tolerance

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max

  ! Validation
  real :: max_error
  integer :: error_count
  logical :: validation_passed

  ! Loop variables
  integer :: iter

  !---------------------------------------------------------------------
  ! Read configuration
  !---------------------------------------------------------------------
  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)

  !---------------------------------------------------------------------
  ! Read parameters
  !---------------------------------------------------------------------
  call read_parameters(trim(data_dir)//'/params.txt')

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: aggregat (GPU)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' cphopt: ', cphopt
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qc(0:ni+1, 0:nj+1, 1:nk))
  allocate(qr(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncc(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncr(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncs(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqc(0:ni+1, 0:nj+1, 1:nk))
  allocate(diaqr(0:ni+1, 0:nj+1, 1:nk))
  allocate(agcn(0:ni+1, 0:nj+1, 1:nk))
  allocate(agrn(0:ni+1, 0:nj+1, 1:nk))
  allocate(agin(0:ni+1, 0:nj+1, 1:nk))
  allocate(agsn(0:ni+1, 0:nj+1, 1:nk))
  allocate(agcn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(agrn_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(agin_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(agsn_ref(0:ni+1, 0:nj+1, 1:nk))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rbv.bin', rbv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qc.bin', qc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qr.bin', qr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qi.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  if (abs(cphopt) == 4) then
    call read_array_3d(trim(data_dir)//'/ncc.bin', ncc, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/ncr.bin', ncr, 0, ni+1, 0, nj+1, 1, nk)
  else
    ncc = 0.0
    ncr = 0.0
  end if
  call read_array_3d(trim(data_dir)//'/nci.bin', nci, 0, ni+1, 0, nj+1, 1, nk)
  if (abs(cphopt) >= 3) then
    call read_array_3d(trim(data_dir)//'/ncs.bin', ncs, 0, ni+1, 0, nj+1, 1, nk)
  else
    ncs = 0.0
  end if
  call read_array_3d(trim(data_dir)//'/diaqc.bin', diaqc, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/diaqr.bin', diaqr, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  if (abs(cphopt) == 4) then
    call read_array_3d(trim(data_dir)//'/agcn_ref.bin', agcn_ref, 0, ni+1, 0, nj+1, 1, nk)
    call read_array_3d(trim(data_dir)//'/agrn_ref.bin', agrn_ref, 0, ni+1, 0, nj+1, 1, nk)
  else
    agcn_ref = 0.0
    agrn_ref = 0.0
  end if
  call read_array_3d(trim(data_dir)//'/agin_ref.bin', agin_ref, 0, ni+1, 0, nj+1, 1, nk)
  if (abs(cphopt) >= 3) then
    call read_array_3d(trim(data_dir)//'/agsn_ref.bin', agsn_ref, 0, ni+1, 0, nj+1, 1, nk)
  else
    agsn_ref = 0.0
  end if

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_aggregat(cphopt, dtb, thresq, ni, nj, nk, &
         rbr, rbv, qc, qr, qi, qs, ncc, ncr, nci, ncs, diaqc, diaqr, &
         agcn, agrn, agin, agsn)
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
    !$acc wait
    t_start = omp_get_wtime()
    call kernel_aggregat(cphopt, dtb, thresq, ni, nj, nk, &
         rbr, rbv, qc, qr, qi, qs, ncc, ncr, nci, ncs, diaqc, diaqr, &
         agcn, agrn, agin, agsn)
    !$acc wait
    t_end = omp_get_wtime()
    t_total = t_total + (t_end - t_start)
    if (t_end - t_start < t_min) t_min = t_end - t_start
    if (t_end - t_start > t_max) t_max = t_end - t_start
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'
  call validate_output(agin, agin_ref, agsn, agsn_ref, ni, nj, nk, &
       tolerance, max_error, error_count, validation_passed)

  !---------------------------------------------------------------------
  ! Report results
  !---------------------------------------------------------------------
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Benchmark Results'
  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:  ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:      ', t_min * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:      ', t_max * 1000.0d0, ' ms'
  write(*,'(A,ES12.4)') ' Max error:     ', max_error
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(rbr, rbv, qc, qr, qi, qs)
  deallocate(ncc, ncr, nci, ncs, diaqc, diaqr)
  deallocate(agcn, agrn, agin, agsn)
  deallocate(agcn_ref, agrn_ref, agin_ref, agsn_ref)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: aggregat (GPU version)
  !-------------------------------------------------------------------
  subroutine kernel_aggregat(cphopt, dtb, thresq, ni, nj, nk, &
       rbr, rbv, qc, qr, qi, qs, ncc, ncr, nci, ncs, diaqc, diaqr, &
       agcn, agrn, agin, agsn)
    implicit none

    integer, intent(in) :: cphopt, ni, nj, nk
    real, intent(in) :: dtb, thresq
    real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncc(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: nci(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ncs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: diaqc(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: diaqr(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: agcn(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: agrn(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: agin(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: agsn(0:ni+1,0:nj+1,1:nk)

    ! Local variables
    real :: expo1, expo2
    real :: cagcn, cagrn1, cagrn2, cagin, cagsn
    real :: diaqc3, ercol
    integer :: i, j, k

    ! Set coefficients (using module-level constants)
    expo1 = oned3 * (2.0e0 + bus)
    expo2 = oned3 * (4.0e0 - bus)
    cagcn = 2.59e15 * gf7 * dtb
    cagrn1 = 2.59e15 * gf7 * dtb
    cagrn2 = 3.03e3 * gf4 * dtb
    cagin = r0 * exp(3.0e0 * log(0.5e0 * aui * eii * xl / rhoi * dtb))
    cagsn = 3.47222e-4 * aus * ess * hfbus * dtb &
         * exp(oned3 * (1.0e0 - bus) * log(cc)) &
         * exp(-oned3 * (2.0e0 + bus) * log(rhos))

    if (nk == 1) then
      ! Case nk = 1
      if (abs(cphopt) == 2) then
        !$acc kernels
        !$acc loop independent collapse(2)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qi(i,j,1) > thresq) then
              agin(i,j,1) = rbr(i,j,1) * nci(i,j,1) &
                   * qi(i,j,1) * exp(oned3 * log(cagin * rbv(i,j,1)))
            else
              agin(i,j,1) = 0.0e0
            end if
          end do
        end do
        !$acc end kernels

      else if (abs(cphopt) == 3) then
        !$acc kernels
        !$acc loop independent collapse(2)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qi(i,j,1) > thresq) then
              agin(i,j,1) = rbr(i,j,1) * nci(i,j,1) &
                   * qi(i,j,1) * exp(oned3 * log(cagin * rbv(i,j,1)))
            else
              agin(i,j,1) = 0.0e0
            end if
            if (qs(i,j,1) > thresq) then
              agsn(i,j,1) = cagsn * rbr(i,j,1) &
                   * exp(expo1 * log(qs(i,j,1))) * exp(expo2 * log(ncs(i,j,1)))
            else
              agsn(i,j,1) = 0.0e0
            end if
          end do
        end do
        !$acc end kernels

      else if (abs(cphopt) == 4) then
        !$acc kernels
        !$acc loop independent collapse(2) private(diaqc3, ercol)
        do j = 1, nj-1
          do i = 1, ni-1
            if (qc(i,j,1) > thresq) then
              diaqc3 = diaqc(i,j,1) * diaqc(i,j,1) * diaqc(i,j,1)
              agcn(i,j,1) = cagcn * diaqc3 * diaqc3 * ncc(i,j,1) * ncc(i,j,1) * rbr(i,j,1)
            else
              agcn(i,j,1) = 0.0e0
            end if
            if (qr(i,j,1) > thresq) then
              if (diaqr(i,j,1) < 6.0e-4) then
                ercol = 1.0e0
              else if (diaqr(i,j,1) >= 6.0e-4 .and. diaqr(i,j,1) < 2.0e-3) then
                ercol = exp(-2.5e3 * (diaqr(i,j,1) - 6.0e-4))
              else
                ercol = 0.0e0
              end if
              agrn(i,j,1) = cagrn2 * diaqr(i,j,1)**3 * ercol * ncr(i,j,1)**2 * rbr(i,j,1)
            else
              agrn(i,j,1) = 0.0e0
            end if
            if (qi(i,j,1) > thresq) then
              agin(i,j,1) = rbr(i,j,1) * nci(i,j,1) &
                   * qi(i,j,1) * exp(oned3 * log(cagin * rbv(i,j,1)))
            else
              agin(i,j,1) = 0.0e0
            end if
            if (qs(i,j,1) > thresq) then
              agsn(i,j,1) = cagsn * rbr(i,j,1) &
                   * exp(expo1 * log(qs(i,j,1))) * exp(expo2 * log(ncs(i,j,1)))
            else
              agsn(i,j,1) = 0.0e0
            end if
          end do
        end do
        !$acc end kernels
      end if

    else
      ! Case nk > 1
      if (abs(cphopt) == 2) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni-1
              if (qi(i,j,k) > thresq) then
                agin(i,j,k) = rbr(i,j,k) * nci(i,j,k) &
                     * qi(i,j,k) * exp(oned3 * log(cagin * rbv(i,j,k)))
              else
                agin(i,j,k) = 0.0e0
              end if
            end do
          end do
        end do
        !$acc end kernels

      else if (abs(cphopt) == 3) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni-1
              if (qi(i,j,k) > thresq) then
                agin(i,j,k) = rbr(i,j,k) * nci(i,j,k) &
                     * qi(i,j,k) * exp(oned3 * log(cagin * rbv(i,j,k)))
              else
                agin(i,j,k) = 0.0e0
              end if
              if (qs(i,j,k) > thresq) then
                agsn(i,j,k) = cagsn * rbr(i,j,k) &
                     * exp(expo1 * log(qs(i,j,k))) * exp(expo2 * log(ncs(i,j,k)))
              else
                agsn(i,j,k) = 0.0e0
              end if
            end do
          end do
        end do
        !$acc end kernels

      else if (abs(cphopt) == 4) then
        !$acc kernels
        !$acc loop independent collapse(3) private(diaqc3, ercol)
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni-1
              if (qc(i,j,k) > thresq) then
                diaqc3 = diaqc(i,j,k) * diaqc(i,j,k) * diaqc(i,j,k)
                agcn(i,j,k) = cagcn * diaqc3 * diaqc3 * ncc(i,j,k) * ncc(i,j,k) * rbr(i,j,k)
              else
                agcn(i,j,k) = 0.0e0
              end if
              if (qr(i,j,k) > thresq) then
                if (diaqr(i,j,k) < 6.0e-4) then
                  ercol = 1.0e0
                else if (diaqr(i,j,k) >= 6.0e-4 .and. diaqr(i,j,k) < 2.0e-3) then
                  ercol = exp(-2.5e3 * (diaqr(i,j,k) - 6.0e-4))
                else
                  ercol = 0.0e0
                end if
                agrn(i,j,k) = cagrn2 * diaqr(i,j,k)**3 * ercol * ncr(i,j,k)**2 * rbr(i,j,k)
              else
                agrn(i,j,k) = 0.0e0
              end if
              if (qi(i,j,k) > thresq) then
                agin(i,j,k) = rbr(i,j,k) * nci(i,j,k) &
                     * qi(i,j,k) * exp(oned3 * log(cagin * rbv(i,j,k)))
              else
                agin(i,j,k) = 0.0e0
              end if
              if (qs(i,j,k) > thresq) then
                agsn(i,j,k) = cagsn * rbr(i,j,k) &
                     * exp(expo1 * log(qs(i,j,k))) * exp(expo2 * log(ncs(i,j,k)))
              else
                agsn(i,j,k) = 0.0e0
              end if
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_aggregat

  !-------------------------------------------------------------------
  ! Read configuration
  !-------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios

    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      data_dir = './data'
      num_iter = 10
      warmup_iter = 2
      tol = 1.0e-5
      return
    end if
    read(10, '(A)') data_dir
    read(10, *) num_iter
    read(10, *) warmup_iter
    read(10, *) tol
    close(10)
  end subroutine read_config

  !-------------------------------------------------------------------
  ! Read parameters (key = value format)
  !-------------------------------------------------------------------
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
          case ('cphopt'); read(val, *) cphopt
          case ('ni'); read(val, *) ni
          case ('nj'); read(val, *) nj
          case ('nk'); read(val, *) nk
          case ('dtb'); read(val, *) dtb
          case ('thresq'); read(val, *) thresq
          case ('r0'); read(val, *) r0
          case ('rhoi'); read(val, *) rhoi
          case ('rhos'); read(val, *) rhos
          case ('rhow'); read(val, *) rhow
          case ('aui'); read(val, *) aui
          case ('aus'); read(val, *) aus
          case ('bus'); read(val, *) bus
          case ('xl'); read(val, *) xl
          case ('ecc'); read(val, *) ecc
          case ('eii'); read(val, *) eii
          case ('ess'); read(val, *) ess
          case ('cc'); read(val, *) cc
          case ('oned3'); read(val, *) oned3
          case ('gf4'); read(val, *) gf4
          case ('gf7'); read(val, *) gf7
          case ('hfbus'); read(val, *) hfbus
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 3D array
  !-------------------------------------------------------------------
  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)

    open(unit=10, file=filename, access='stream', form='unformatted', status='old')
    read(10) arr
    close(10)
  end subroutine read_array_3d

  !-------------------------------------------------------------------
  ! Validate output
  !-------------------------------------------------------------------
  subroutine validate_output(agin, agin_ref, agsn, agsn_ref, ni, nj, nk, &
       tol, max_err, err_count, passed)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: agin(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: agin_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: agsn(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: agsn_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count
    logical, intent(out) :: passed

    real :: rel_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = 1, nk
      do j = 1, nj-1
        do i = 1, ni-1
          if (abs(agin_ref(i,j,k)) > 1.0e-30) then
            rel_err = abs(agin(i,j,k) - agin_ref(i,j,k)) / abs(agin_ref(i,j,k))
          else
            rel_err = abs(agin(i,j,k) - agin_ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1

          if (abs(agsn_ref(i,j,k)) > 1.0e-30) then
            rel_err = abs(agsn(i,j,k) - agsn_ref(i,j,k)) / abs(agsn_ref(i,j,k))
          else
            rel_err = abs(agsn(i,j,k) - agsn_ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

    passed = (err_count == 0)
  end subroutine validate_output

end program kernel_benchmark_aggregat

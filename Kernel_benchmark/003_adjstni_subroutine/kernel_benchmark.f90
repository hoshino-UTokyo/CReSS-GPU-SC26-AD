!***********************************************************************
! Kernel Benchmark: adjstni (s_adjstni)
!***********************************************************************
!
! Source: Src/adjstni.f90
! Description: Adjust concentrations of ice hydrometeors (cloud ice,
!              snow, graupel, hail) to be consistent with mixing ratios
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_adjstni
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: haiopt
  real :: rhog, rhos

  ! Physical constants - derived values from params.txt
  real :: mimiv5, mi0iv2, msmiv2, ms0iv2, mgmiv2, mg0iv2, mhmiv2, mh0iv2
  real :: cdiaqs, cdiaqg, cdiaqh

  ! Input arrays
  real, allocatable :: rbr(:,:,:)
  real, allocatable :: rbv(:,:,:)
  real, allocatable :: qi(:,:,:)
  real, allocatable :: qs(:,:,:)
  real, allocatable :: qg(:,:,:)
  real, allocatable :: qh(:,:,:)

  ! Input/Output arrays
  real, allocatable :: nci(:,:,:)
  real, allocatable :: ncs(:,:,:)
  real, allocatable :: ncg(:,:,:)
  real, allocatable :: nch(:,:,:)

  ! Initial values for resetting
  real, allocatable :: nci_init(:,:,:)
  real, allocatable :: ncs_init(:,:,:)
  real, allocatable :: ncg_init(:,:,:)
  real, allocatable :: nch_init(:,:,:)

  ! Reference output for validation
  real, allocatable :: nci_ref(:,:,:)
  real, allocatable :: ncs_ref(:,:,:)
  real, allocatable :: ncg_ref(:,:,:)
  real, allocatable :: nch_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error
  real :: tolerance
  integer :: error_count
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
  call read_parameters(trim(data_dir)//'/params.txt', &
       haiopt, ni, nj, nk, &
       mimiv5, mi0iv2, msmiv2, ms0iv2, mgmiv2, mg0iv2, mhmiv2, mh0iv2, &
       cdiaqs, cdiaqg, cdiaqh)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: adjstni'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' haiopt: ', haiopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rbv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qi(0:ni+1, 0:nj+1, 1:nk))
  allocate(qs(0:ni+1, 0:nj+1, 1:nk))
  allocate(qg(0:ni+1, 0:nj+1, 1:nk))
  allocate(qh(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncs(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncg(0:ni+1, 0:nj+1, 1:nk))
  allocate(nch(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncs_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncg_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(nch_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(nci_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncs_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncg_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(nch_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rbv.bin', rbv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qi.bin', qi, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qs.bin', qs, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qg.bin', qg, 0, ni+1, 0, nj+1, 1, nk)
  if (haiopt /= 0) then
    call read_array_3d(trim(data_dir)//'/qh.bin', qh, 0, ni+1, 0, nj+1, 1, nk)
  else
    qh = 0.0
  end if
  call read_array_3d(trim(data_dir)//'/nci_in.bin', nci_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncs_in.bin', ncs_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncg_in.bin', ncg_init, 0, ni+1, 0, nj+1, 1, nk)
  if (haiopt /= 0) then
    call read_array_3d(trim(data_dir)//'/nch_in.bin', nch_init, 0, ni+1, 0, nj+1, 1, nk)
  else
    nch_init = 0.0
  end if

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/nci_ref.bin', nci_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncs_ref.bin', ncs_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncg_ref.bin', ncg_ref, 0, ni+1, 0, nj+1, 1, nk)
  if (haiopt /= 0) then
    call read_array_3d(trim(data_dir)//'/nch_ref.bin', nch_ref, 0, ni+1, 0, nj+1, 1, nk)
  else
    nch_ref = 0.0
  end if

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    nci = nci_init
    ncs = ncs_init
    ncg = ncg_init
    nch = nch_init
    call kernel_adjstni(haiopt, ni, nj, nk, &
         mimiv5, mi0iv2, msmiv2, ms0iv2, mgmiv2, mg0iv2, mhmiv2, mh0iv2, &
         cdiaqs, cdiaqg, cdiaqh, &
         rbr, rbv, qi, qs, qg, qh, nci, ncs, ncg, nch)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    nci = nci_init
    ncs = ncs_init
    ncg = ncg_init
    nch = nch_init

    t_start = omp_get_wtime()

    call kernel_adjstni(haiopt, ni, nj, nk, &
         mimiv5, mi0iv2, msmiv2, ms0iv2, mgmiv2, mg0iv2, mhmiv2, mh0iv2, &
         cdiaqs, cdiaqg, cdiaqh, &
         rbr, rbv, qi, qs, qg, qh, nci, ncs, ncg, nch)

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
  call validate_output(nci, nci_ref, ncs, ncs_ref, ncg, ncg_ref, nch, nch_ref, &
       ni, nj, nk, tolerance, max_error, error_count, validation_passed)

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
  deallocate(rbr, rbv, qi, qs, qg, qh)
  deallocate(nci, ncs, ncg, nch)
  deallocate(nci_init, ncs_init, ncg_init, nch_init)
  deallocate(nci_ref, ncs_ref, ncg_ref, nch_ref)
  deallocate(times)

contains

  !-------------------------------------------------------------------
  ! Kernel: adjstni
  !-------------------------------------------------------------------
  subroutine kernel_adjstni(haiopt, ni, nj, nk, &
       mimiv5, mi0iv2, msmiv2, ms0iv2, mgmiv2, mg0iv2, mhmiv2, mh0iv2, &
       cdiaqs, cdiaqg, cdiaqh, &
       rbr, rbv, qi, qs, qg, qh, nci, ncs, ncg, nch)
    implicit none

    integer, intent(in) :: haiopt, ni, nj, nk
    real, intent(in) :: mimiv5, mi0iv2, msmiv2, ms0iv2, mgmiv2, mg0iv2, mhmiv2, mh0iv2
    real, intent(in) :: cdiaqs, cdiaqg, cdiaqh
    real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qh(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: nci(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: ncs(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: ncg(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: nch(0:ni+1,0:nj+1,1:nk)

    real :: ndia
    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    if (haiopt == 0) then
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j,ndia)
        do j = 1, nj-1
          do i = 1, ni-1
            ! Adjust cloud ice concentrations
            nci(i,j,k) = min(max(nci(i,j,k), mimiv5*qi(i,j,k)), mi0iv2*qi(i,j,k))

            ! Adjust snow concentrations
            ndia = sqrt(sqrt(cdiaqs*rbr(i,j,k)*qs(i,j,k))) * rbv(i,j,k)
            ncs(i,j,k) = max(ncs(i,j,k), msmiv2*qs(i,j,k), 5.62341e-3*ndia)
            ncs(i,j,k) = min(ncs(i,j,k), ms0iv2*qs(i,j,k), 1.77838e2*ndia)

            ! Adjust graupel concentrations
            ndia = sqrt(sqrt(cdiaqg*rbr(i,j,k)*qg(i,j,k))) * rbv(i,j,k)
            ncg(i,j,k) = max(ncg(i,j,k), mgmiv2*qg(i,j,k), 5.62341e-3*ndia)
            ncg(i,j,k) = min(ncg(i,j,k), mg0iv2*qg(i,j,k), 1.77838e2*ndia)
          end do
        end do
        !$omp end do
      end do
    else
      do k = 1, nk-1
        !$omp do schedule(runtime) private(i,j,ndia)
        do j = 1, nj-1
          do i = 1, ni-1
            ! Adjust cloud ice concentrations
            nci(i,j,k) = min(max(nci(i,j,k), mimiv5*qi(i,j,k)), mi0iv2*qi(i,j,k))

            ! Adjust snow concentrations
            ndia = sqrt(sqrt(cdiaqs*rbr(i,j,k)*qs(i,j,k))) * rbv(i,j,k)
            ncs(i,j,k) = max(ncs(i,j,k), msmiv2*qs(i,j,k), 5.62341e-3*ndia)
            ncs(i,j,k) = min(ncs(i,j,k), ms0iv2*qs(i,j,k), 1.77838e2*ndia)

            ! Adjust graupel concentrations
            ndia = sqrt(sqrt(cdiaqg*rbr(i,j,k)*qg(i,j,k))) * rbv(i,j,k)
            ncg(i,j,k) = max(ncg(i,j,k), mgmiv2*qg(i,j,k), 5.62341e-3*ndia)
            ncg(i,j,k) = min(ncg(i,j,k), mg0iv2*qg(i,j,k), 1.77838e2*ndia)

            ! Adjust hail concentrations
            ndia = sqrt(sqrt(cdiaqh*rbr(i,j,k)*qh(i,j,k))) * rbv(i,j,k)
            nch(i,j,k) = max(nch(i,j,k), mhmiv2*qh(i,j,k), 5.62341e-3*ndia)
            nch(i,j,k) = min(nch(i,j,k), mh0iv2*qh(i,j,k), 1.77838e2*ndia)
          end do
        end do
        !$omp end do
      end do
    end if

    !$omp end parallel

  end subroutine kernel_adjstni

  !-------------------------------------------------------------------
  ! Read benchmark configuration
  !-------------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
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
  ! Read parameters
  !-------------------------------------------------------------------
  subroutine read_parameters(filename, haiopt, ni, nj, nk, &
       mimiv5, mi0iv2, msmiv2, ms0iv2, mgmiv2, mg0iv2, mhmiv2, mh0iv2, &
       cdiaqs, cdiaqg, cdiaqh)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: haiopt, ni, nj, nk
    real, intent(out) :: mimiv5, mi0iv2, msmiv2, ms0iv2
    real, intent(out) :: mgmiv2, mg0iv2, mhmiv2, mh0iv2
    real, intent(out) :: cdiaqs, cdiaqg, cdiaqh

    character(len=256) :: line
    character(len=64) :: name, val_str
    integer :: ios

    ! Set defaults
    haiopt = 0
    ni = 1
    nj = 1
    nk = 1
    mimiv5 = 0.0
    mi0iv2 = 0.0
    msmiv2 = 0.0
    ms0iv2 = 0.0
    mgmiv2 = 0.0
    mg0iv2 = 0.0
    mhmiv2 = 0.0
    mh0iv2 = 0.0
    cdiaqs = 0.0
    cdiaqg = 0.0
    cdiaqh = 0.0

    open(unit=10, file=filename, status='old')
    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, '=') > 0) then
        read(line, *) name
        val_str = line(index(line, '=')+1:)
        select case (trim(name))
        case ('haiopt')
          read(val_str, *) haiopt
        case ('ni')
          read(val_str, *) ni
        case ('nj')
          read(val_str, *) nj
        case ('nk')
          read(val_str, *) nk
        case ('mimiv5')
          read(val_str, *) mimiv5
        case ('mi0iv2')
          read(val_str, *) mi0iv2
        case ('msmiv2')
          read(val_str, *) msmiv2
        case ('ms0iv2')
          read(val_str, *) ms0iv2
        case ('mgmiv2')
          read(val_str, *) mgmiv2
        case ('mg0iv2')
          read(val_str, *) mg0iv2
        case ('mhmiv2')
          read(val_str, *) mhmiv2
        case ('mh0iv2')
          read(val_str, *) mh0iv2
        case ('cdiaqs')
          read(val_str, *) cdiaqs
        case ('cdiaqg')
          read(val_str, *) cdiaqg
        case ('cdiaqh')
          read(val_str, *) cdiaqh
        end select
      end if
    end do
    close(10)

    ! If ni is not in params, assume square grid
    if (ni == 1 .and. nj > 1) ni = nj
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 3D array
  !-------------------------------------------------------------------
  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
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
  subroutine validate_output(nci, nci_ref, ncs, ncs_ref, ncg, ncg_ref, nch, nch_ref, &
       ni, nj, nk, tol, max_err, err_count, passed)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: nci(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: nci_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ncs(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ncs_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ncg(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ncg_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: nch(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: nch_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count
    logical, intent(out) :: passed

    real :: rel_err
    integer :: i, j, k

    max_err = 0.0
    err_count = 0

    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          ! Check nci
          if (abs(nci_ref(i,j,k)) > 1.0e-30) then
            rel_err = abs(nci(i,j,k) - nci_ref(i,j,k)) / abs(nci_ref(i,j,k))
          else
            rel_err = abs(nci(i,j,k) - nci_ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1

          ! Check ncs
          if (abs(ncs_ref(i,j,k)) > 1.0e-30) then
            rel_err = abs(ncs(i,j,k) - ncs_ref(i,j,k)) / abs(ncs_ref(i,j,k))
          else
            rel_err = abs(ncs(i,j,k) - ncs_ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1

          ! Check ncg
          if (abs(ncg_ref(i,j,k)) > 1.0e-30) then
            rel_err = abs(ncg(i,j,k) - ncg_ref(i,j,k)) / abs(ncg_ref(i,j,k))
          else
            rel_err = abs(ncg(i,j,k) - ncg_ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1

          ! Check nch
          if (abs(nch_ref(i,j,k)) > 1.0e-30) then
            rel_err = abs(nch(i,j,k) - nch_ref(i,j,k)) / abs(nch_ref(i,j,k))
          else
            rel_err = abs(nch(i,j,k) - nch_ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do

    passed = (err_count == 0)
  end subroutine validate_output

end program kernel_benchmark_adjstni

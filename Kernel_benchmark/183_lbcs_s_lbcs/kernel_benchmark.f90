!***********************************************************************
! Kernel Benchmark: lbcs (s_lbcs)
!***********************************************************************
!
! Source: Src/lbcs.f90
! Description: Set lateral boundary conditions for optional scalar variable
!              by copying values from interior to boundary points.
!
!***********************************************************************
program kernel_benchmark_lbcs
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  integer :: wbc, ebc, sbc, nbc, advopt
  integer :: nim1, nim2, njm1, njm2

  ! MPI-related parameters (for single process)
  integer :: ebw, ebe, ebs, ebn
  integer :: isub, jsub, nisub, njsub

  ! Arrays
  real, allocatable :: sf(:,:,:)
  real, allocatable :: sf_init(:,:,:)
  real, allocatable :: sf_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
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
  call read_parameters(trim(data_dir)//'/params.txt')

  ! Set derived parameters
  nim1 = ni - 1
  nim2 = ni - 2
  njm1 = nj - 1
  njm2 = nj - 2

  ! Set MPI parameters for single process execution
  ebw = 1
  ebe = 1
  ebs = 1
  ebn = 1
  isub = 0
  jsub = 0
  nisub = 1
  njsub = 1

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: lbcs'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I3,A,I3,A,I3,A,I3)') ' wbc=', wbc, ', ebc=', ebc, ', sbc=', sbc, ', nbc=', nbc
  write(*,'(A,I3)') ' advopt=', advopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(sf(0:ni+1, 0:nj+1, 1:nk))
  allocate(sf_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(sf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/sf_in.bin', sf_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/sf_ref.bin', sf_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    sf = sf_init
    call kernel_lbcs(wbc, ebc, sbc, nbc, advopt, ni, nj, nk, &
         nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
         isub, jsub, nisub, njsub, sf)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    sf = sf_init
    t_start = omp_get_wtime()

    call kernel_lbcs(wbc, ebc, sbc, nbc, advopt, ni, nj, nk, &
         nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
         isub, jsub, nisub, njsub, sf)

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
  call validate_output(sf, sf_ref, ni, nj, nk, tolerance, max_error, error_count)

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
  deallocate(sf, sf_init, sf_ref, times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: lbcs
  !-------------------------------------------------------------------
  subroutine kernel_lbcs(wbc, ebc, sbc, nbc, advopt, ni, nj, nk, &
       nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
       isub, jsub, nisub, njsub, sf)
    implicit none

    integer, intent(in) :: wbc, ebc, sbc, nbc, advopt
    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: nim1, nim2, njm1, njm2
    integer, intent(in) :: ebw, ebe, ebs, ebn
    integer, intent(in) :: isub, jsub, nisub, njsub
    real, intent(inout) :: sf(0:ni+1,0:nj+1,1:nk)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! West boundary
    if (ebw == 1 .and. isub == 0) then
      if (advopt <= 3) then
        if (wbc == 2 .or. wbc == 3) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(j)
            do j = 1, nj-1
              sf(1,j,k) = sf(2,j,k)
            end do
            !$omp end do
          end do
        end if
      else
        if (wbc >= 2) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(j)
            do j = 1, nj-1
              sf(1,j,k) = sf(2,j,k)
            end do
            !$omp end do
          end do
        end if
      end if
    end if

    ! East boundary
    if (ebe == 1 .and. isub == nisub-1) then
      if (advopt <= 3) then
        if (ebc == 2 .or. ebc == 3) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(j)
            do j = 1, nj-1
              sf(nim1,j,k) = sf(nim2,j,k)
            end do
            !$omp end do
          end do
        end if
      else
        if (ebc >= 2) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(j)
            do j = 1, nj-1
              sf(nim1,j,k) = sf(nim2,j,k)
            end do
            !$omp end do
          end do
        end if
      end if
    end if

    ! South boundary
    if (ebs == 1 .and. jsub == 0) then
      if (advopt <= 3) then
        if (sbc == 2 .or. sbc == 3) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i)
            do i = 1, ni-1
              sf(i,1,k) = sf(i,2,k)
            end do
            !$omp end do
          end do
        end if
      else
        if (sbc >= 2) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i)
            do i = 1, ni-1
              sf(i,1,k) = sf(i,2,k)
            end do
            !$omp end do
          end do
        end if
      end if
    end if

    ! North boundary
    if (ebn == 1 .and. jsub == njsub-1) then
      if (advopt <= 3) then
        if (nbc == 2 .or. nbc == 3) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i)
            do i = 1, ni-1
              sf(i,njm1,k) = sf(i,njm2,k)
            end do
            !$omp end do
          end do
        end if
      else
        if (nbc >= 2) then
          do k = 2, nk-2
            !$omp do schedule(runtime) private(i)
            do i = 1, ni-1
              sf(i,njm1,k) = sf(i,njm2,k)
            end do
            !$omp end do
          end do
        end if
      end if
    end if

    !$omp end parallel

  end subroutine kernel_lbcs

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
          case ('wbc')
            read(val, *) wbc
          case ('ebc')
            read(val, *) ebc
          case ('sbc')
            read(val, *) sbc
          case ('nbc')
            read(val, *) nbc
          case ('advopt')
            read(val, *) advopt
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  !-------------------------------------------------------------------
  ! Read 3D array
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
  ! Validate output
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

end program kernel_benchmark_lbcs

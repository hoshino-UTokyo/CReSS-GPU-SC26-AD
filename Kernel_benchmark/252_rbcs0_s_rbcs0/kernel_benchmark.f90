!***********************************************************************
! Kernel Benchmark: rbcs0 (s_rbcs0)
!***********************************************************************
!
! Source: Src/rbcs0.f90
! Description: Set radiative lateral boundary conditions for optional
!              scalar variable at domain corners and edges.
!
!***********************************************************************
program kernel_benchmark_rbcs0
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Parameters
  character(len=108) :: lbcvar
  integer :: wbc, ebc, sbc, nbc, advopt, apl
  integer :: nim1, nim2, njm1, njm2
  real :: lbnews, dt, dmpdt

  ! MPI-related parameters (for single process)
  integer :: ebw, ebe, ebs, ebn
  integer :: isub, jsub, nisub, njsub

  ! Arrays
  real, allocatable :: s(:,:,:), sp(:,:,:)
  real, allocatable :: scpx(:,:,:), scpy(:,:,:)
  real, allocatable :: sf(:,:,:), sf_init(:,:,:), sf_ref(:,:,:)

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

  ! Compute damping coefficient
  if (advopt <= 3) then
    if (lbcvar(apl:apl) == 'o') then
      dmpdt = 2.0e0 * lbnews * dt
    else
      dmpdt = 0.0e0
    end if
  else
    if (lbcvar(apl:apl) == 'o') then
      dmpdt = lbnews * dt
    else
      dmpdt = 0.0e0
    end if
  end if

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: rbcs0'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I3,A,I3,A,I3,A,I3)') ' wbc=', wbc, ', ebc=', ebc, ', sbc=', sbc, ', nbc=', nbc
  write(*,'(A,I3)') ' advopt=', advopt
  write(*,'(A,I3)') ' apl=', apl
  write(*,'(A,ES12.4)') ' dmpdt=', dmpdt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(s(0:ni+1, 0:nj+1, 1:nk))
  allocate(sp(0:ni+1, 0:nj+1, 1:nk))
  allocate(scpx(1:nj, 1:nk, 1:2))
  allocate(scpy(1:ni, 1:nk, 1:2))
  allocate(sf(0:ni+1, 0:nj+1, 1:nk))
  allocate(sf_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(sf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/s.bin', s, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/sp.bin', sp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d_alt(trim(data_dir)//'/scpx.bin', scpx, 1, nj, 1, nk, 1, 2)
  call read_array_3d_alt(trim(data_dir)//'/scpy.bin', scpy, 1, ni, 1, nk, 1, 2)
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
    call kernel_rbcs0(wbc, ebc, sbc, nbc, advopt, ni, nj, nk, &
         nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
         isub, jsub, nisub, njsub, dmpdt, s, sp, scpx, scpy, sf)
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

    call kernel_rbcs0(wbc, ebc, sbc, nbc, advopt, ni, nj, nk, &
         nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
         isub, jsub, nisub, njsub, dmpdt, s, sp, scpx, scpy, sf)

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
  deallocate(s, sp, scpx, scpy, sf, sf_init, sf_ref, times)

  if (.not. validation_passed) stop 1

contains

  !-------------------------------------------------------------------
  ! Kernel: rbcs0
  !-------------------------------------------------------------------
  subroutine kernel_rbcs0(wbc, ebc, sbc, nbc, advopt, ni, nj, nk, &
       nim1, nim2, njm1, njm2, ebw, ebe, ebs, ebn, &
       isub, jsub, nisub, njsub, dmpdt, s, sp, scpx, scpy, sf)
    implicit none

    integer, intent(in) :: wbc, ebc, sbc, nbc, advopt
    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: nim1, nim2, njm1, njm2
    integer, intent(in) :: ebw, ebe, ebs, ebn
    integer, intent(in) :: isub, jsub, nisub, njsub
    real, intent(in) :: dmpdt
    real, intent(in) :: s(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: sp(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: scpx(1:nj, 1:nk, 1:2)
    real, intent(in) :: scpy(1:ni, 1:nk, 1:2)
    real, intent(inout) :: sf(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    real :: gamma, radwe, radsn

    !$omp parallel default(shared)

    ! Set the boundary conditions at the four corners
    if (ebs == 1 .and. jsub == 0) then
      if (advopt <= 3) then
        if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. sbc >= 4) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-2
            radwe = (s(2,1,k) - sp(1,1,k)) * scpx(1,k,1) / (1.0e0 - scpx(1,k,1))
            radsn = (s(1,2,k) - sp(1,1,k)) * scpy(1,k,1) / (1.0e0 - scpy(1,k,1))
            sf(1,1,k) = sp(1,1,k) - 2.0e0 * (radwe + radsn) - dmpdt * sp(1,1,k)
          end do
          !$omp end do
        end if

        if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. sbc >= 4) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-2
            radwe = (s(nim2,1,k) - sp(nim1,1,k)) * scpx(1,k,2) / (1.0e0 + scpx(1,k,2))
            radsn = (s(nim1,2,k) - sp(nim1,1,k)) * scpy(nim1,k,1) / (1.0e0 - scpy(nim1,k,1))
            sf(nim1,1,k) = sp(nim1,1,k) + 2.0e0 * (radwe - radsn) - dmpdt * sp(nim1,1,k)
          end do
          !$omp end do
        end if
      else
        if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. sbc >= 4) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-2
            radwe = (sp(2,1,k) - sp(1,1,k)) * scpx(1,k,1)
            radsn = (sp(1,2,k) - sp(1,1,k)) * scpy(1,k,1)
            sf(1,1,k) = sp(1,1,k) - (radwe + radsn) - dmpdt * sp(1,1,k)
          end do
          !$omp end do
        end if

        if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. sbc >= 4) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-2
            radwe = (sp(nim2,1,k) - sp(nim1,1,k)) * scpx(1,k,2)
            radsn = (sp(nim1,2,k) - sp(nim1,1,k)) * scpy(nim1,k,1)
            sf(nim1,1,k) = sp(nim1,1,k) + (radwe - radsn) - dmpdt * sp(nim1,1,k)
          end do
          !$omp end do
        end if
      end if
    end if

    if (ebn == 1 .and. jsub == njsub-1) then
      if (advopt <= 3) then
        if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. nbc >= 4) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-2
            radwe = (s(2,njm1,k) - sp(1,njm1,k)) * scpx(njm1,k,1) / (1.0e0 - scpx(njm1,k,1))
            radsn = (s(1,njm2,k) - sp(1,njm1,k)) * scpy(1,k,2) / (1.0e0 + scpy(1,k,2))
            sf(1,njm1,k) = sp(1,njm1,k) - 2.0e0 * (radwe - radsn) - dmpdt * sp(1,njm1,k)
          end do
          !$omp end do
        end if

        if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. nbc >= 4) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-2
            radwe = (s(nim2,njm1,k) - sp(nim1,njm1,k)) * scpx(njm1,k,2) / (1.0e0 + scpx(njm1,k,2))
            radsn = (s(nim1,njm2,k) - sp(nim1,njm1,k)) * scpy(nim1,k,2) / (1.0e0 + scpy(nim1,k,2))
            sf(nim1,njm1,k) = sp(nim1,njm1,k) + 2.0e0 * (radwe + radsn) - dmpdt * sp(nim1,njm1,k)
          end do
          !$omp end do
        end if
      else
        if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. nbc >= 4) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-2
            radwe = (sp(2,njm1,k) - sp(1,njm1,k)) * scpx(njm1,k,1)
            radsn = (sp(1,njm2,k) - sp(1,njm1,k)) * scpy(1,k,2)
            sf(1,njm1,k) = sp(1,njm1,k) - (radwe - radsn) - dmpdt * sp(1,njm1,k)
          end do
          !$omp end do
        end if

        if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. nbc >= 4) then
          !$omp do schedule(runtime) private(k,radwe,radsn)
          do k = 2, nk-2
            radwe = (sp(nim2,njm1,k) - sp(nim1,njm1,k)) * scpx(njm1,k,2)
            radsn = (sp(nim1,njm2,k) - sp(nim1,njm1,k)) * scpy(nim1,k,2)
            sf(nim1,njm1,k) = sp(nim1,njm1,k) + (radwe + radsn) - dmpdt * sp(nim1,njm1,k)
          end do
          !$omp end do
        end if
      end if
    end if

    ! Set the west boundary conditions
    if (ebw == 1 .and. isub == 0) then
      if (wbc >= 4) then
        if (advopt <= 3) then
          !$omp do schedule(runtime) private(j,k,gamma)
          do k = 2, nk-2
            do j = 2, nj-2
              gamma = 2.0e0 * scpx(j,k,1) / (1.0e0 - scpx(j,k,1))
              sf(1,j,k) = sp(1,j,k) - gamma * (s(2,j,k) - sp(1,j,k)) - dmpdt * sp(1,j,k)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(j,k)
          do k = 2, nk-2
            do j = 2, nj-2
              sf(1,j,k) = sp(1,j,k) - scpx(j,k,1) * (sp(2,j,k) - sp(1,j,k)) - dmpdt * sp(1,j,k)
            end do
          end do
          !$omp end do
        end if
      end if
    end if

    ! Set the east boundary conditions
    if (ebe == 1 .and. isub == nisub-1) then
      if (ebc >= 4) then
        if (advopt <= 3) then
          !$omp do schedule(runtime) private(j,k,gamma)
          do k = 2, nk-2
            do j = 2, nj-2
              gamma = 2.0e0 * scpx(j,k,2) / (1.0e0 + scpx(j,k,2))
              sf(nim1,j,k) = sp(nim1,j,k) + gamma * (s(nim2,j,k) - sp(nim1,j,k)) - dmpdt * sp(nim1,j,k)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(j,k)
          do k = 2, nk-2
            do j = 2, nj-2
              sf(nim1,j,k) = sp(nim1,j,k) + scpx(j,k,2) * (sp(nim2,j,k) - sp(nim1,j,k)) - dmpdt * sp(nim1,j,k)
            end do
          end do
          !$omp end do
        end if
      end if
    end if

    ! Set the south boundary conditions
    if (ebs == 1 .and. jsub == 0) then
      if (sbc >= 4) then
        if (advopt <= 3) then
          !$omp do schedule(runtime) private(i,k,gamma)
          do k = 2, nk-2
            do i = 2, ni-2
              gamma = 2.0e0 * scpy(i,k,1) / (1.0e0 - scpy(i,k,1))
              sf(i,1,k) = sp(i,1,k) - gamma * (s(i,2,k) - sp(i,1,k)) - dmpdt * sp(i,1,k)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(i,k)
          do k = 2, nk-2
            do i = 2, ni-2
              sf(i,1,k) = sp(i,1,k) - scpy(i,k,1) * (sp(i,2,k) - sp(i,1,k)) - dmpdt * sp(i,1,k)
            end do
          end do
          !$omp end do
        end if
      end if
    end if

    ! Set the north boundary conditions
    if (ebn == 1 .and. jsub == njsub-1) then
      if (nbc >= 4) then
        if (advopt <= 3) then
          !$omp do schedule(runtime) private(i,k,gamma)
          do k = 2, nk-2
            do i = 2, ni-2
              gamma = 2.0e0 * scpy(i,k,2) / (1.0e0 + scpy(i,k,2))
              sf(i,njm1,k) = sp(i,njm1,k) + gamma * (s(i,njm2,k) - sp(i,njm1,k)) - dmpdt * sp(i,njm1,k)
            end do
          end do
          !$omp end do
        else
          !$omp do schedule(runtime) private(i,k)
          do k = 2, nk-2
            do i = 2, ni-2
              sf(i,njm1,k) = sp(i,njm1,k) + scpy(i,k,2) * (sp(i,njm2,k) - sp(i,njm1,k)) - dmpdt * sp(i,njm1,k)
            end do
          end do
          !$omp end do
        end if
      end if
    end if

    !$omp end parallel

  end subroutine kernel_rbcs0

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
          case ('lbcvar')
            read(val, '(A)') lbcvar
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
          case ('lbnews')
            read(val, *) lbnews
          case ('apl')
            read(val, *) apl
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('dt')
            read(val, *) dt
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
  ! Read 3D array with alternative indexing
  !-------------------------------------------------------------------
  subroutine read_array_3d_alt(filename, arr, is, ie, js, je, ks, ke)
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
  end subroutine read_array_3d_alt

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

end program kernel_benchmark_rbcs0

!***********************************************************************
! GPU Kernel Benchmark: rbcs (s_rbcs)
!***********************************************************************
! Source: Src/rbcs.f90
! Description: Set radiative lateral boundary conditions for scalar
! GPU Port: OpenACC with Unified Memory
!***********************************************************************
program kernel_benchmark_rbcs_gpu
  use omp_lib
  implicit none

  integer :: ni, nj, nk
  character(len=108) :: lbcvar
  integer :: wbc, ebc, sbc, nbc, advopt, apl
  integer :: nim1, nim2, njm1, njm2
  real :: lbnews, dt, dmpdt
  integer :: ebw, ebe, ebs, ebn
  integer :: isub, jsub, nisub, njsub

  real, allocatable :: s(:,:,:), sp(:,:,:)
  real, allocatable :: scpx(:,:,:), scpy(:,:,:)
  real, allocatable :: sf(:,:,:), sf_init(:,:,:), sf_ref(:,:,:)

  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)
  real :: max_error, tolerance
  integer :: error_count, iter
  logical :: validation_passed

  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)
  call read_parameters(trim(data_dir)//'/params.txt')

  nim1 = ni - 1; nim2 = ni - 2; njm1 = nj - 1; njm2 = nj - 2
  ebw = 1; ebe = 1; ebs = 1; ebn = 1
  isub = 0; jsub = 0; nisub = 1; njsub = 1

  if (advopt <= 3) then
    dmpdt = merge(2.0e0 * lbnews * dt, 0.0e0, lbcvar(apl:apl) == 'o')
  else
    dmpdt = merge(lbnews * dt, 0.0e0, lbcvar(apl:apl) == 'o')
  end if

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: rbcs'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk

  allocate(s(0:ni+1, 0:nj+1, 1:nk))
  allocate(sp(0:ni+1, 0:nj+1, 1:nk))
  allocate(scpx(1:nj, 1:nk, 1:2))
  allocate(scpy(1:ni, 1:nk, 1:2))
  allocate(sf(0:ni+1, 0:nj+1, 1:nk))
  allocate(sf_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(sf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/s.bin', s, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/sp.bin', sp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d_alt(trim(data_dir)//'/scpx.bin', scpx, 1, nj, 1, nk, 1, 2)
  call read_array_3d_alt(trim(data_dir)//'/scpy.bin', scpy, 1, ni, 1, nk, 1, 2)
  call read_array_3d(trim(data_dir)//'/sf_in.bin', sf_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/sf_ref.bin', sf_ref, 0, ni+1, 0, nj+1, 1, nk)

  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    sf = sf_init
    call kernel_rbcs()
    !$acc wait
  end do

  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  do iter = 1, num_iterations
    sf = sf_init
    !$acc wait
    t_start = omp_get_wtime()
    call kernel_rbcs()
    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do
  t_avg = t_total / dble(num_iterations)

  write(*,'(A)') ' Validating output...'
  call validate_output(sf, sf_ref, ni, nj, nk, tolerance, max_error, error_count)
  validation_passed = (error_count == 0)

  write(*,'(A)') '=================================================='
  write(*,'(A,F12.3,A)') ' Average time:  ', t_avg * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Min time:      ', minval(times) * 1000.0d0, ' ms'
  write(*,'(A,F12.3,A)') ' Max time:      ', maxval(times) * 1000.0d0, ' ms'
  write(*,'(A)') '--------------------------------------------------'
  write(*,'(A,ES12.4)') ' Max relative error: ', max_error
  write(*,'(A,I12)')    ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  deallocate(s, sp, scpx, scpy, sf, sf_init, sf_ref, times)
  if (.not. validation_passed) stop 1

contains

  subroutine kernel_rbcs()
    integer :: i, j, k
    real :: gamma, radwe, radsn

    ! Corners
    if (ebs == 1 .and. jsub == 0) then
      if (advopt <= 3) then
        if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. sbc >= 4) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            radwe = (s(2,1,k) - sp(1,1,k)) * scpx(1,k,1) / (1.0e0 - scpx(1,k,1))
            radsn = (s(1,2,k) - sp(1,1,k)) * scpy(1,k,1) / (1.0e0 - scpy(1,k,1))
            sf(1,1,k) = sp(1,1,k) - 2.0e0 * (radwe + radsn) - dmpdt * sp(1,1,k)
          end do
          !$acc end kernels
        end if
        if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. sbc >= 4) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            radwe = (s(nim2,1,k) - sp(nim1,1,k)) * scpx(1,k,2) / (1.0e0 + scpx(1,k,2))
            radsn = (s(nim1,2,k) - sp(nim1,1,k)) * scpy(nim1,k,1) / (1.0e0 - scpy(nim1,k,1))
            sf(nim1,1,k) = sp(nim1,1,k) + 2.0e0 * (radwe - radsn) - dmpdt * sp(nim1,1,k)
          end do
          !$acc end kernels
        end if
      else
        if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. sbc >= 4) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            radwe = (sp(2,1,k) - sp(1,1,k)) * scpx(1,k,1)
            radsn = (sp(1,2,k) - sp(1,1,k)) * scpy(1,k,1)
            sf(1,1,k) = sp(1,1,k) - (radwe + radsn) - dmpdt * sp(1,1,k)
          end do
          !$acc end kernels
        end if
        if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. sbc >= 4) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            radwe = (sp(nim2,1,k) - sp(nim1,1,k)) * scpx(1,k,2)
            radsn = (sp(nim1,2,k) - sp(nim1,1,k)) * scpy(nim1,k,1)
            sf(nim1,1,k) = sp(nim1,1,k) + (radwe - radsn) - dmpdt * sp(nim1,1,k)
          end do
          !$acc end kernels
        end if
      end if
    end if

    if (ebn == 1 .and. jsub == njsub-1) then
      if (advopt <= 3) then
        if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. nbc >= 4) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            radwe = (s(2,njm1,k) - sp(1,njm1,k)) * scpx(njm1,k,1) / (1.0e0 - scpx(njm1,k,1))
            radsn = (s(1,njm2,k) - sp(1,njm1,k)) * scpy(1,k,2) / (1.0e0 + scpy(1,k,2))
            sf(1,njm1,k) = sp(1,njm1,k) - 2.0e0 * (radwe - radsn) - dmpdt * sp(1,njm1,k)
          end do
          !$acc end kernels
        end if
        if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. nbc >= 4) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            radwe = (s(nim2,njm1,k) - sp(nim1,njm1,k)) * scpx(njm1,k,2) / (1.0e0 + scpx(njm1,k,2))
            radsn = (s(nim1,njm2,k) - sp(nim1,njm1,k)) * scpy(nim1,k,2) / (1.0e0 + scpy(nim1,k,2))
            sf(nim1,njm1,k) = sp(nim1,njm1,k) + 2.0e0 * (radwe + radsn) - dmpdt * sp(nim1,njm1,k)
          end do
          !$acc end kernels
        end if
      else
        if (ebw == 1 .and. isub == 0 .and. wbc >= 4 .and. nbc >= 4) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            radwe = (sp(2,njm1,k) - sp(1,njm1,k)) * scpx(njm1,k,1)
            radsn = (sp(1,njm2,k) - sp(1,njm1,k)) * scpy(1,k,2)
            sf(1,njm1,k) = sp(1,njm1,k) - (radwe - radsn) - dmpdt * sp(1,njm1,k)
          end do
          !$acc end kernels
        end if
        if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4 .and. nbc >= 4) then
          !$acc kernels
          !$acc loop independent
          do k = 2, nk-2
            radwe = (sp(nim2,njm1,k) - sp(nim1,njm1,k)) * scpx(njm1,k,2)
            radsn = (sp(nim1,njm2,k) - sp(nim1,njm1,k)) * scpy(nim1,k,2)
            sf(nim1,njm1,k) = sp(nim1,njm1,k) + (radwe + radsn) - dmpdt * sp(nim1,njm1,k)
          end do
          !$acc end kernels
        end if
      end if
    end if

    ! West boundary
    if (ebw == 1 .and. isub == 0 .and. wbc >= 4) then
      if (advopt <= 3) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-2
            gamma = 2.0e0 * scpx(j,k,1) / (1.0e0 - scpx(j,k,1))
            sf(1,j,k) = sp(1,j,k) - gamma * (s(2,j,k) - sp(1,j,k)) - dmpdt * sp(1,j,k)
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-2
            sf(1,j,k) = sp(1,j,k) - scpx(j,k,1) * (sp(2,j,k) - sp(1,j,k)) - dmpdt * sp(1,j,k)
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! East boundary
    if (ebe == 1 .and. isub == nisub-1 .and. ebc >= 4) then
      if (advopt <= 3) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-2
            gamma = 2.0e0 * scpx(j,k,2) / (1.0e0 + scpx(j,k,2))
            sf(nim1,j,k) = sp(nim1,j,k) + gamma * (s(nim2,j,k) - sp(nim1,j,k)) - dmpdt * sp(nim1,j,k)
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do j = 2, nj-2
            sf(nim1,j,k) = sp(nim1,j,k) + scpx(j,k,2) * (sp(nim2,j,k) - sp(nim1,j,k)) - dmpdt * sp(nim1,j,k)
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! South boundary
    if (ebs == 1 .and. jsub == 0 .and. sbc >= 4) then
      if (advopt <= 3) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 2, ni-2
            gamma = 2.0e0 * scpy(i,k,1) / (1.0e0 - scpy(i,k,1))
            sf(i,1,k) = sp(i,1,k) - gamma * (s(i,2,k) - sp(i,1,k)) - dmpdt * sp(i,1,k)
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 2, ni-2
            sf(i,1,k) = sp(i,1,k) - scpy(i,k,1) * (sp(i,2,k) - sp(i,1,k)) - dmpdt * sp(i,1,k)
          end do
        end do
        !$acc end kernels
      end if
    end if

    ! North boundary
    if (ebn == 1 .and. jsub == njsub-1 .and. nbc >= 4) then
      if (advopt <= 3) then
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 2, ni-2
            gamma = 2.0e0 * scpy(i,k,2) / (1.0e0 + scpy(i,k,2))
            sf(i,njm1,k) = sp(i,njm1,k) + gamma * (s(i,njm2,k) - sp(i,njm1,k)) - dmpdt * sp(i,njm1,k)
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent
        do k = 2, nk-2
          !$acc loop independent
          do i = 2, ni-2
            sf(i,njm1,k) = sp(i,njm1,k) + scpy(i,k,2) * (sp(i,njm2,k) - sp(i,njm1,k)) - dmpdt * sp(i,njm1,k)
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine kernel_rbcs

  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios
    data_dir = './data'; num_iter = 10; warmup_iter = 2; tol = 1.0e-5
    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios == 0) then
      read(10, '(A)', iostat=ios) data_dir
      read(10, *, iostat=ios) num_iter
      read(10, *, iostat=ios) warmup_iter
      read(10, *, iostat=ios) tol
      close(10)
    end if
  end subroutine read_config

  subroutine read_parameters(filename)
    character(len=*), intent(in) :: filename
    character(len=256) :: line, key, val
    integer :: ios, eq_pos
    open(unit=10, file=filename, status='old', iostat=ios)
    if (ios /= 0) stop 1
    do while (.true.)
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        key = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))
        select case (trim(key))
          case ('lbcvar'); read(val, '(A)') lbcvar
          case ('wbc'); read(val, *) wbc
          case ('ebc'); read(val, *) ebc
          case ('sbc'); read(val, *) sbc
          case ('nbc'); read(val, *) nbc
          case ('advopt'); read(val, *) advopt
          case ('lbnews'); read(val, *) lbnews
          case ('apl'); read(val, *) apl
          case ('ni'); read(val, *) ni
          case ('nj'); read(val, *) nj
          case ('nk'); read(val, *) nk
          case ('dt'); read(val, *) dt
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: ios
    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) stop 1
    read(10) arr
    close(10)
  end subroutine read_array_3d

  subroutine read_array_3d_alt(filename, arr, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    integer :: ios
    open(unit=10, file=filename, access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) stop 1
    read(10) arr
    close(10)
  end subroutine read_array_3d_alt

  subroutine validate_output(output, reference, ni, nj, nk, tol, max_err, err_count)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: output(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: reference(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_count
    real :: rel_err
    integer :: i, j, k
    max_err = 0.0; err_count = 0
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

end program kernel_benchmark_rbcs_gpu

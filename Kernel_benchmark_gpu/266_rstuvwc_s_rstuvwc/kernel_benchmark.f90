!***********************************************************************
! Kernel Benchmark: rstuvwc (s_rstuvwc) - OpenACC GPU Version
!***********************************************************************
!
! Source: Src/rstuvwc.f90
! Description: Multiply base state density x Jacobian by velocity
!              components u, v, and wc
!
!***********************************************************************
program kernel_benchmark_rstuvwc
  use omp_lib
  implicit none

  integer :: ni, nj, nk
  integer :: mpopt, mfcopt, iwest, ieast, jsouth, jnorth

  real, allocatable :: mf8u(:,:), mf8v(:,:)
  real, allocatable :: rst8u(:,:,:), rst8v(:,:,:), rst8w(:,:,:)
  real, allocatable :: u(:,:,:), v(:,:,:), wc(:,:,:)
  real, allocatable :: rstxu(:,:,:), rstxv(:,:,:), rstxwc(:,:,:)
  real, allocatable :: rstxu_ref(:,:,:), rstxv_ref(:,:,:), rstxwc_ref(:,:,:)

  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)
  real :: max_error, tolerance
  integer :: error_count, iter
  logical :: validation_passed

  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)
  call read_parameters(trim(data_dir)//'/params.txt', &
       mpopt, mfcopt, iwest, ieast, jsouth, jnorth, ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: rstuvwc (OpenACC GPU)'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  allocate(mf8u(0:ni+1, 0:nj+1), mf8v(0:ni+1, 0:nj+1))
  allocate(rst8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxu(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxv(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxwc(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxu_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxv_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(rstxwc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  write(*,'(A)') ' Loading input data...'
  call read_array_2d(trim(data_dir)//'/mf8u.bin', mf8u, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/mf8v.bin', mf8v, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/rst8u.bin', rst8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8v.bin', rst8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst8w.bin', rst8w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/wc.bin', wc, 0, ni+1, 0, nj+1, 1, nk)

  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/rstxu_ref.bin', rstxu_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rstxv_ref.bin', rstxv_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rstxwc_ref.bin', rstxwc_ref, 0, ni+1, 0, nj+1, 1, nk)

  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    call kernel_rstuvwc(mpopt, mfcopt, iwest, ieast, jsouth, jnorth, &
         ni, nj, nk, mf8u, mf8v, rst8u, rst8v, rst8w, u, v, wc, &
         rstxu, rstxv, rstxwc)
  end do

  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    !$acc wait
    t_start = omp_get_wtime()
    call kernel_rstuvwc(mpopt, mfcopt, iwest, ieast, jsouth, jnorth, &
         ni, nj, nk, mf8u, mf8v, rst8u, rst8v, rst8w, u, v, wc, &
         rstxu, rstxv, rstxwc)
    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
    if (times(iter) < t_min) t_min = times(iter)
    if (times(iter) > t_max) t_max = times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  write(*,'(A)') ' Validating output...'
  call validate_output(rstxu, rstxu_ref, rstxv, rstxv_ref, rstxwc, rstxwc_ref, &
       ni, nj, nk, tolerance, max_error, error_count, validation_passed)

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

  deallocate(mf8u, mf8v, rst8u, rst8v, rst8w, u, v, wc)
  deallocate(rstxu, rstxv, rstxwc, rstxu_ref, rstxv_ref, rstxwc_ref, times)

contains

  subroutine kernel_rstuvwc(mpopt, mfcopt, iwest, ieast, jsouth, jnorth, &
       ni, nj, nk, mf8u, mf8v, rst8u, rst8v, rst8w, u, v, wc, &
       rstxu, rstxv, rstxwc)
    implicit none
    integer, intent(in) :: mpopt, mfcopt, iwest, ieast, jsouth, jnorth
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: mf8u(0:ni+1,0:nj+1), mf8v(0:ni+1,0:nj+1)
    real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: wc(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: rstxu(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: rstxv(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: rstxwc(0:ni+1,0:nj+1,1:nk)

    integer :: i, j, k

    if (mfcopt == 0) then
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 1, nk-1
        do j = jsouth, nj-jnorth
          do i = iwest, ni+1-ieast
            rstxu(i,j,k) = rst8u(i,j,k) * u(i,j,k)
          end do
        end do
      end do
      !$acc end kernels
      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 1, nk-1
        do j = jsouth, nj+1-jnorth
          do i = iwest, ni-ieast
            rstxv(i,j,k) = rst8v(i,j,k) * v(i,j,k)
          end do
        end do
      end do
      !$acc end kernels
    else
      if (mpopt == 0 .or. mpopt == 10) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 1, nk-1
          do j = jsouth, nj-jnorth
            do i = iwest, ni+1-ieast
              rstxu(i,j,k) = mf8u(i,j) * rst8u(i,j,k) * u(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 1, nk-1
          do j = jsouth, nj+1-jnorth
            do i = iwest, ni-ieast
              rstxv(i,j,k) = rst8v(i,j,k) * v(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      else if (mpopt == 5) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 1, nk-1
          do j = jsouth, nj-jnorth
            do i = iwest, ni+1-ieast
              rstxu(i,j,k) = rst8u(i,j,k) * u(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 1, nk-1
          do j = jsouth, nj+1-jnorth
            do i = iwest, ni-ieast
              rstxv(i,j,k) = mf8v(i,j) * rst8v(i,j,k) * v(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      else
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 1, nk-1
          do j = jsouth, nj-jnorth
            do i = iwest, ni+1-ieast
              rstxu(i,j,k) = mf8u(i,j) * rst8u(i,j,k) * u(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 1, nk-1
          do j = jsouth, nj+1-jnorth
            do i = iwest, ni-ieast
              rstxv(i,j,k) = mf8v(i,j) * rst8v(i,j,k) * v(i,j,k)
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk
      do j = jsouth, nj-jnorth
        do i = iwest, ni-ieast
          rstxwc(i,j,k) = rst8w(i,j,k) * wc(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_rstuvwc

  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol
    integer :: ios
    open(unit=10, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      data_dir = './data'; num_iter = 10; warmup_iter = 2; tol = 1.0e-5
      return
    end if
    read(10, '(A)') data_dir
    read(10, *) num_iter
    read(10, *) warmup_iter
    read(10, *) tol
    close(10)
  end subroutine read_config

  subroutine read_parameters(filename, mpopt, mfcopt, iwest, ieast, jsouth, jnorth, ni, nj, nk)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: mpopt, mfcopt, iwest, ieast, jsouth, jnorth, ni, nj, nk
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
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
          case ('iwest')
            read(val, *) iwest
          case ('ieast')
            read(val, *) ieast
          case ('jsouth')
            read(val, *) jsouth
          case ('jnorth')
            read(val, *) jnorth
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

  subroutine read_array_2d(filename, arr, is, ie, js, je)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: arr(is:ie, js:je)
    open(unit=10, file=filename, access='stream', form='unformatted', status='old')
    read(10) arr
    close(10)
  end subroutine read_array_2d

  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    open(unit=10, file=filename, access='stream', form='unformatted', status='old')
    read(10) arr
    close(10)
  end subroutine read_array_3d

  subroutine validate_output(rstxu, rstxu_ref, rstxv, rstxv_ref, rstxwc, rstxwc_ref, &
       ni, nj, nk, tol, max_err, err_count, passed)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: rstxu(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rstxu_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rstxv(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rstxv_ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rstxwc(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: rstxwc_ref(0:ni+1, 0:nj+1, 1:nk)
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
        do i = 1, ni
          if (abs(rstxu_ref(i,j,k)) > 1.0e-30) then
            rel_err = abs(rstxu(i,j,k) - rstxu_ref(i,j,k)) / abs(rstxu_ref(i,j,k))
          else
            rel_err = abs(rstxu(i,j,k) - rstxu_ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do
    passed = (err_count == 0)
  end subroutine validate_output

end program kernel_benchmark_rstuvwc

!***********************************************************************
! Kernel Benchmark: upwnp (s_upwnp)
!***********************************************************************
!
! Source: Src/upwnp.f90
! Description: Calculate sedimentation for optional precipitation
!              concentrations using upwind flux divergence scheme
!
!***********************************************************************
program kernel_benchmark_upwnp
  use omp_lib
  implicit none

  integer :: ni, nj, nk
  real :: dziv, dtp

  real, allocatable :: rbr(:,:,:), rst(:,:,:), un(:,:,:)
  real, allocatable :: ncf(:,:,:), ncflx(:,:,:)
  real, allocatable :: ncf_init(:,:,:)
  real, allocatable :: ncf_ref(:,:,:)

  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir
  real(8) :: t_start, t_end, t_total, t_avg, t_min, t_max
  real(8), allocatable :: times(:)
  real :: max_error, tolerance
  integer :: error_count, iter
  logical :: validation_passed

  call read_config(data_dir, num_iterations, warmup_iterations, tolerance)
  call read_parameters(trim(data_dir)//'/params.txt', dziv, dtp, ni, nj, nk)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: upwnp'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  allocate(rbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(un(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncf(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncflx(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncf_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(ncf_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  write(*,'(A)') ' Loading input data...'
  call read_array_3d(trim(data_dir)//'/rbr.bin', rbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/rst.bin', rst, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/un.bin', un, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ncf_in.bin', ncf_init, 0, ni+1, 0, nj+1, 1, nk)

  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/ncf_ref.bin', ncf_ref, 0, ni+1, 0, nj+1, 1, nk)

  ncflx = 0.0

  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ncf = ncf_init
    call kernel_upwnp(dziv, dtp, ni, nj, nk, rbr, rst, un, ncf, ncflx)
  end do

  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0
  t_min = 1.0d30
  t_max = 0.0d0

  do iter = 1, num_iterations
    ncf = ncf_init
    t_start = omp_get_wtime()
    call kernel_upwnp(dziv, dtp, ni, nj, nk, rbr, rst, un, ncf, ncflx)
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
    if (times(iter) < t_min) t_min = times(iter)
    if (times(iter) > t_max) t_max = times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  write(*,'(A)') ' Validating output...'
  call validate_output(ncf, ncf_ref, ni, nj, nk, tolerance, max_error, error_count, validation_passed)

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

  deallocate(rbr, rst, un, ncf, ncflx, ncf_init, ncf_ref, times)

contains

  subroutine kernel_upwnp(dziv, dtp, ni, nj, nk, rbr, rst, un, ncf, ncflx)
    implicit none
    real, intent(in) :: dziv, dtp
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: un(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: ncf(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: ncflx(0:ni+1,0:nj+1,1:nk)

    integer :: nkm1, nkm2
    real :: dzvdt
    integer :: i, j, k

    nkm1 = nk - 1
    nkm2 = nk - 2
    dzvdt = dziv * dtp

    !$omp parallel default(shared) private(k)

    ! Calculate flux
    do k = 2, nkm1
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          ncflx(i,j,k) = rbr(i,j,k) * un(i,j,k) * ncf(i,j,k)
        end do
      end do
      !$omp end do
    end do

    ! Set boundary flux
    !$omp do schedule(runtime) private(i,j)
    do j = 1, nj-1
      do i = 1, ni-1
        ncflx(i,j,1) = 0.0e0
        ncflx(i,j,nk) = ncflx(i,j,nkm1)
      end do
    end do
    !$omp end do

    ! Update concentration
    do k = 1, nkm2
      !$omp do schedule(runtime) private(i,j)
      do j = 1, nj-1
        do i = 1, ni-1
          ncf(i,j,k) = max(ncf(i,j,k) &
               + (ncflx(i,j,k+1) - ncflx(i,j,k)) * dzvdt / rst(i,j,k), 0.0e0)
        end do
      end do
      !$omp end do
    end do

    ! Copy boundary
    !$omp do schedule(runtime) private(i,j)
    do j = 1, nj-1
      do i = 1, ni-1
        ncf(i,j,nkm1) = ncf(i,j,nkm2)
      end do
    end do
    !$omp end do

    !$omp end parallel
  end subroutine kernel_upwnp

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

  subroutine read_parameters(filename, dziv, dtp, ni, nj, nk)
    character(len=*), intent(in) :: filename
    real, intent(out) :: dziv, dtp
    integer, intent(out) :: ni, nj, nk
    character(len=256) :: line
    integer :: ios, eq_pos

    ! Read parameters in name=value format
    open(unit=10, file=filename, status='old')
    do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        select case(trim(adjustl(line(1:eq_pos-1))))
          case('dziv'); read(line(eq_pos+1:), *) dziv
          case('dtp');  read(line(eq_pos+1:), *) dtp
          case('ni');   read(line(eq_pos+1:), *) ni
          case('nj');   read(line(eq_pos+1:), *) nj
          case('nk');   read(line(eq_pos+1:), *) nk
        end select
      end if
    end do
    close(10)
  end subroutine read_parameters

  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)
    open(unit=10, file=filename, access='stream', form='unformatted', status='old')
    read(10) arr
    close(10)
  end subroutine read_array_3d

  subroutine validate_output(ncf, ncf_ref, ni, nj, nk, tol, max_err, err_count, passed)
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: ncf(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ncf_ref(0:ni+1, 0:nj+1, 1:nk)
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
          if (abs(ncf_ref(i,j,k)) > 1.0e-30) then
            rel_err = abs(ncf(i,j,k) - ncf_ref(i,j,k)) / abs(ncf_ref(i,j,k))
          else
            rel_err = abs(ncf(i,j,k) - ncf_ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_count = err_count + 1
        end do
      end do
    end do
    passed = (err_count == 0)
  end subroutine validate_output

end program kernel_benchmark_upwnp

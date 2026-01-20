!***********************************************************************
! Kernel Benchmark: steptund_sec1 (s_steptund section 1)
!***********************************************************************
!
! Source: Src/steptund.f90
! Description: Set boundary conditions and coefficient matrix for
!              soil/sea temperature tridiagonal equation solver
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_steptund_sec1
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nund, nundm1
  integer :: nk  ! For work arrays (same as nund)

  ! Parameters
  integer :: sfcopt
  real :: ctg1, ctgm1, cts1, ctsm1
  real :: s1g, skg, rkg, tkg
  real :: s1s, sks, rks, tks
  real :: dtsoil, stinc, t0

  ! Input arrays (2D)
  integer, allocatable :: land(:,:)
  real, allocatable :: cap(:,:)
  real, allocatable :: nuu(:,:)
  real, allocatable :: sst(:,:)
  real, allocatable :: sstd(:,:)
  real, allocatable :: hs(:,:)
  real, allocatable :: le(:,:)
  real, allocatable :: rsd(:,:)
  real, allocatable :: rld(:,:)
  real, allocatable :: rlu(:,:)

  ! Input arrays (3D)
  real, allocatable :: t(:,:,:)
  real, allocatable :: ss_in(:,:,:)

  ! Input/output arrays
  real, allocatable :: tundp(:,:,:)
  real, allocatable :: tundp_init(:,:,:)
  real, allocatable :: rr(:,:,:)
  real, allocatable :: rr_init(:,:,:)
  real, allocatable :: tt(:,:,:)
  real, allocatable :: tt_init(:,:,:)
  real, allocatable :: tmp1(:,:,:)

  ! Output array
  real, allocatable :: tundf(:,:,:)

  ! Reference outputs for validation
  real, allocatable :: tundf_ref(:,:,:)
  real, allocatable :: tundp_ref(:,:,:)
  real, allocatable :: rr_ref(:,:,:)
  real, allocatable :: tt_ref(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_error, rel_error, tol_val
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
  call read_parameters(trim(data_dir)//'/params.txt', &
       sfcopt, ni, nj, nund, nundm1, &
       ctg1, ctgm1, cts1, ctsm1, &
       s1g, skg, rkg, tkg, s1s, sks, rks, tks, &
       dtsoil, stinc, t0)

  nk = nund  ! Work arrays use nund size

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: steptund_sec1'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nund=', nund
  write(*,'(A,I6)') ' sfcopt=', sfcopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(land(0:ni+1, 0:nj+1))
  allocate(cap(0:ni+1, 0:nj+1))
  allocate(nuu(0:ni+1, 0:nj+1))
  allocate(sst(0:ni+1, 0:nj+1))
  allocate(sstd(0:ni+1, 0:nj+1))
  allocate(hs(0:ni+1, 0:nj+1))
  allocate(le(0:ni+1, 0:nj+1))
  allocate(rsd(0:ni+1, 0:nj+1))
  allocate(rld(0:ni+1, 0:nj+1))
  allocate(rlu(0:ni+1, 0:nj+1))
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(ss_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tundp(0:ni+1, 0:nj+1, 1:nund))
  allocate(tundp_init(0:ni+1, 0:nj+1, 1:nund))
  allocate(rr(0:ni+1, 0:nj+1, 1:nk))
  allocate(rr_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(tt(0:ni+1, 0:nj+1, 1:nk))
  allocate(tt_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tundf(0:ni+1, 0:nj+1, 1:nund))
  allocate(tundf_ref(0:ni+1, 0:nj+1, 1:nund))
  allocate(tundp_ref(0:ni+1, 0:nj+1, 1:nund))
  allocate(rr_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tt_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_2d_int(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/cap.bin', cap, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/nuu.bin', nuu, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/sst.bin', sst, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/sstd.bin', sstd, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/hs.bin', hs, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/le.bin', le, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/rsd.bin', rsd, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/rld.bin', rld, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/rlu.bin', rlu, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/t.bin', t, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ss.bin', ss_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tundp_in.bin', tundp_init, 0, ni+1, 0, nj+1, 1, nund)
  call read_array_3d(trim(data_dir)//'/rr_in.bin', rr_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tt_in.bin', tt_init, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/tundf_ref.bin', tundf_ref, 0, ni+1, 0, nj+1, 1, nund)
  call read_array_3d(trim(data_dir)//'/tundp_ref.bin', tundp_ref, 0, ni+1, 0, nj+1, 1, nund)
  call read_array_3d(trim(data_dir)//'/rr_ref.bin', rr_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tt_ref.bin', tt_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    ! Reset input arrays
    tundp = tundp_init
    rr = rr_init
    tt = tt_init

    call kernel_steptund_sec1(sfcopt, ni, nj, nund, nundm1, &
         ctg1, ctgm1, cts1, ctsm1, &
         s1g, skg, rkg, tkg, s1s, sks, rks, tks, &
         t0, stinc, &
         t, land, cap, nuu, sst, sstd, hs, le, rsd, rld, rlu, &
         tundp, tundf, rr, ss_in, tt)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    ! Reset input arrays
    tundp = tundp_init
    rr = rr_init
    tt = tt_init

    t_start = omp_get_wtime()

    call kernel_steptund_sec1(sfcopt, ni, nj, nund, nundm1, &
         ctg1, ctgm1, cts1, ctsm1, &
         s1g, skg, rkg, tkg, s1s, sks, rks, tks, &
         t0, stinc, &
         t, land, cap, nuu, sst, sstd, hs, le, rsd, rld, rlu, &
         tundp, tundf, rr, ss_in, tt)

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
  error_count = 0

  ! Validate tundf
  do k = 1, nund
    do j = 1, nj-1
      do i = 1, ni-1
        rel_error = abs(tundf(i,j,k) - tundf_ref(i,j,k))
        tol_val = abs(tundf_ref(i,j,k))
        if (tol_val > 1.0e-10) then
          rel_error = rel_error / tol_val
        end if
        if (rel_error > max_error) max_error = rel_error
        if (rel_error > tolerance) error_count = error_count + 1
      end do
    end do
  end do

  validation_passed = (error_count == 0)

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
  write(*,'(A,I12)') ' Error count:        ', error_count
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(land, cap, nuu, sst, sstd, hs, le, rsd, rld, rlu)
  deallocate(t, ss_in, tundp, tundp_init, tundf)
  deallocate(rr, rr_init, tt, tt_init, tmp1)
  deallocate(tundf_ref, tundp_ref, rr_ref, tt_ref)
  deallocate(times)

contains

  !-----------------------------------------------------------------
  ! Kernel: steptund section 1
  !-----------------------------------------------------------------
  subroutine kernel_steptund_sec1(sfcopt, ni, nj, nund, nundm1, &
       ctg1, ctgm1, cts1, ctsm1, &
       s1g, skg, rkg, tkg, s1s, sks, rks, tks, &
       t0, stinc, &
       t, land, cap, nuu, sst, sstd, hs, le, rsd, rld, rlu, &
       tundp, tundf, rr, ss, tt)
    implicit none

    integer, intent(in) :: sfcopt, ni, nj, nund, nundm1
    real, intent(in) :: ctg1, ctgm1, cts1, ctsm1
    real, intent(in) :: s1g, skg, rkg, tkg, s1s, sks, rks, tks
    real, intent(in) :: t0, stinc
    real, intent(in) :: t(0:ni+1,0:nj+1,1:nund)
    integer, intent(in) :: land(0:ni+1,0:nj+1)
    real, intent(in) :: cap(0:ni+1,0:nj+1)
    real, intent(in) :: nuu(0:ni+1,0:nj+1)
    real, intent(in) :: sst(0:ni+1,0:nj+1)
    real, intent(in) :: sstd(0:ni+1,0:nj+1)
    real, intent(in) :: hs(0:ni+1,0:nj+1)
    real, intent(in) :: le(0:ni+1,0:nj+1)
    real, intent(in) :: rsd(0:ni+1,0:nj+1)
    real, intent(in) :: rld(0:ni+1,0:nj+1)
    real, intent(in) :: rlu(0:ni+1,0:nj+1)
    real, intent(inout) :: tundp(0:ni+1,0:nj+1,1:nund)
    real, intent(out) :: tundf(0:ni+1,0:nj+1,1:nund)
    real, intent(inout) :: rr(0:ni+1,0:nj+1,1:nund)
    real, intent(in) :: ss(0:ni+1,0:nj+1,1:nund)
    real, intent(inout) :: tt(0:ni+1,0:nj+1,1:nund)

    integer :: i, j, k

    !$omp parallel default(shared) private(k)

    ! Set the ice and snow surface temperature
    !$omp do schedule(runtime) private(i,j)
    do j=1,nj-1
      do i=1,ni-1
        if(land(i,j).ge.3.and.land(i,j).lt.10) then
          tundp(i,j,1)=min(.5e0*(t(i,j,1)+t(i,j,2)),t0)
        end if
      end do
    end do
    !$omp end do

    do k=2,nund
      !$omp do schedule(runtime) private(i,j)
      do j=1,nj-1
        do i=1,ni-1
          if(land(i,j).ge.3.and.land(i,j).lt.10) then
            tundp(i,j,k)=tundp(i,j,1)
          end if
        end do
      end do
      !$omp end do
    end do

    ! Copy the past value to future and convert Kelvin to Celsius
    do k=1,nund
      !$omp do schedule(runtime) private(i,j)
      do j=1,nj-1
        do i=1,ni-1
          tundf(i,j,k)=tundp(i,j,k)-t0
        end do
      end do
      !$omp end do
    end do

    ! Set the top and bottom boundary conditions
    if(sfcopt.eq.1.or.sfcopt.eq.11) then
      !$omp do schedule(runtime) private(i,j)
      do j=1,nj-1
        do i=1,ni-1
          if(land(i,j).lt.3) then
            tundf(i,j,1)=tundf(i,j,1) &
                 +cts1*(rsd(i,j)+rld(i,j)-rlu(i,j)-hs(i,j)-le(i,j))/cap(i,j)
            tundf(i,j,nundm1)=(1.e0+ctsm1*nuu(i,j))*tundf(i,j,nundm1)
          end if
          if(land(i,j).ge.10) then
            tundf(i,j,1)=tundf(i,j,1) &
                 +ctg1*(rsd(i,j)+rld(i,j)-rlu(i,j)-hs(i,j)-le(i,j))/cap(i,j)
            tundf(i,j,nundm1)=tundf(i,j,nundm1)+ctgm1*nuu(i,j)*tundf(i,j,nund)
          end if
        end do
      end do
      !$omp end do

    else if(sfcopt.eq.2.or.sfcopt.eq.12) then
      !$omp do schedule(runtime) private(i,j)
      do j=1,nj-1
        do i=1,ni-1
          if(land(i,j).ge.10) then
            tundf(i,j,1)=tundf(i,j,1) &
                 +ctg1*(rsd(i,j)+rld(i,j)-rlu(i,j)-hs(i,j)-le(i,j))/cap(i,j)
            tundf(i,j,nundm1)=tundf(i,j,nundm1)+ctgm1*nuu(i,j)*tundf(i,j,nund)
          end if
        end do
      end do
      !$omp end do

    else if(sfcopt.eq.3.or.sfcopt.eq.13) then
      !$omp do schedule(runtime) private(i,j)
      do j=1,nj-1
        do i=1,ni-1
          if(land(i,j).lt.3) then
            tundf(i,j,1)=(sst(i,j)+sstd(i,j)*stinc)-t0
          end if
          if(land(i,j).ge.10) then
            tundf(i,j,1)=tundf(i,j,1) &
                 +ctg1*(rsd(i,j)+rld(i,j)-rlu(i,j)-hs(i,j)-le(i,j))/cap(i,j)
            tundf(i,j,nundm1)=tundf(i,j,nundm1)+ctgm1*nuu(i,j)*tundf(i,j,nund)
          end if
        end do
      end do
      !$omp end do
    end if

    ! Set the constant sea temperature
    if(sfcopt.eq.3.or.sfcopt.eq.13) then
      do k=2,nund
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            if(land(i,j).lt.3) then
              tundf(i,j,k)=tundf(i,j,1)
            end if
          end do
        end do
        !$omp end do
      end do
    end if

    ! Set the coefficient matrix
    if(sfcopt.eq.1.or.sfcopt.eq.11) then
      do k=1,nund-1
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            if(land(i,j).lt.3) then
              rr(i,j,k)=rks*nuu(i,j)
              tt(i,j,k)=tks*nuu(i,j)
            else if(land(i,j).ge.10) then
              rr(i,j,k)=rkg*nuu(i,j)
              tt(i,j,k)=tkg*nuu(i,j)
            else
              rr(i,j,k)=0.e0
              tt(i,j,k)=0.e0
            end if
          end do
        end do
        !$omp end do
      end do

    else if(sfcopt.eq.2.or.sfcopt.eq.3.or.sfcopt.ge.12) then
      do k=1,nund-1
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            if(land(i,j).ge.10) then
              rr(i,j,k)=rkg*nuu(i,j)
              tt(i,j,k)=tkg*nuu(i,j)
            else
              rr(i,j,k)=0.e0
              tt(i,j,k)=0.e0
            end if
          end do
        end do
        !$omp end do
      end do
    end if

    !$omp end parallel

  end subroutine kernel_steptund_sec1

  !-----------------------------------------------------------------
  ! Read benchmark configuration
  !-----------------------------------------------------------------
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    implicit none
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    open(unit=10, file='benchmark.conf', status='old', action='read')
    read(10,'(A)') data_dir
    read(10,*) num_iter
    read(10,*) warmup_iter
    read(10,*) tol
    close(10)
  end subroutine read_config

  !-----------------------------------------------------------------
  ! Read parameters from file
  !-----------------------------------------------------------------
  subroutine read_parameters(filename, sfcopt, ni, nj, nund, nundm1, &
       ctg1, ctgm1, cts1, ctsm1, &
       s1g, skg, rkg, tkg, s1s, sks, rks, tks, &
       dtsoil, stinc, t0)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: sfcopt, ni, nj, nund, nundm1
    real, intent(out) :: ctg1, ctgm1, cts1, ctsm1
    real, intent(out) :: s1g, skg, rkg, tkg, s1s, sks, rks, tks
    real, intent(out) :: dtsoil, stinc, t0

    character(len=256) :: line
    character(len=32) :: name
    real :: dummy_r
    integer :: dummy_i

    ! Read parameters in the order they appear in params.txt
    open(unit=10, file=filename, status='old', action='read')
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) sfcopt
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) dummy_r  ! dzgrd (skip)
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) dummy_r  ! dzsea (skip)
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) ni
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) nj
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) dummy_i  ! nk (skip)
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) nund
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) dtsoil
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) stinc
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) t0
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) ctg1
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) ctgm1
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) cts1
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) ctsm1
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) nundm1
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) rkg
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) rks
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) s1g
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) s1s
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) skg
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) sks
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) tkg
    read(10,'(A)') line; read(line(index(line,'=')+1:),*) tks
    close(10)
  end subroutine read_parameters

  !-----------------------------------------------------------------
  ! Read 2D real array from binary file
  !-----------------------------------------------------------------
  subroutine read_array_2d(filename, arr, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    real, intent(out) :: arr(is:ie, js:je)

    open(unit=10, file=filename, status='old', access='stream', form='unformatted')
    read(10) arr
    close(10)
  end subroutine read_array_2d

  !-----------------------------------------------------------------
  ! Read 2D integer array from binary file
  !-----------------------------------------------------------------
  subroutine read_array_2d_int(filename, arr, is, ie, js, je)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je
    integer, intent(out) :: arr(is:ie, js:je)

    open(unit=10, file=filename, status='old', access='stream', form='unformatted')
    read(10) arr
    close(10)
  end subroutine read_array_2d_int

  !-----------------------------------------------------------------
  ! Read 3D real array from binary file
  !-----------------------------------------------------------------
  subroutine read_array_3d(filename, arr, is, ie, js, je, ks, ke)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: is, ie, js, je, ks, ke
    real, intent(out) :: arr(is:ie, js:je, ks:ke)

    open(unit=10, file=filename, status='old', access='stream', form='unformatted')
    read(10) arr
    close(10)
  end subroutine read_array_3d

end program kernel_benchmark_steptund_sec1

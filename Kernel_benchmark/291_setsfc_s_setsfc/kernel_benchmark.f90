!***********************************************************************
! Kernel Benchmark: setsfc (s_setsfc)
!***********************************************************************
!
! Source: Src/setsfc.f90
! Description: Calculate surface parameters including pressure, temperature,
!              virtual potential temperature, saturation mixing ratio,
!              and velocity magnitude at lowest levels.
!
! This benchmark:
!   1. Reads parameters and input arrays from files
!   2. Executes the kernel specified number of times
!   3. Validates output against reference data
!   4. Reports timing statistics
!
!***********************************************************************
program kernel_benchmark_setsfc
  use omp_lib
  implicit none

  ! Constants (from m_comphy and m_commath)
  real, parameter :: es0 = 610.78
  real, parameter :: t0 = 273.16
  real, parameter :: epsva = 0.622
  real, parameter :: epsav = 0.608
  real, parameter :: tlow = 236.15
  real, parameter :: lim35n = 238.15
  real, parameter :: p0 = 100000.0
  real, parameter :: vamin = 0.4e0
  real, parameter :: icbeta = 1.0

  ! Array dimensions
  integer :: ni, nj, nk, nund

  ! Parameters
  character(len=5) :: fmois
  integer :: levpbl, tubopt, cphopt
  real :: rddvcp, p0iv

  ! Input arrays (2D)
  integer, allocatable :: land(:,:)
  real, allocatable :: beta(:,:)
  real, allocatable :: kai(:,:)
  real, allocatable :: fall(:,:)

  ! Input arrays (3D)
  real, allocatable :: pbr(:,:,:)
  real, allocatable :: ptbr(:,:,:)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)
  real, allocatable :: w(:,:,:)
  real, allocatable :: pp(:,:,:)
  real, allocatable :: ptp(:,:,:)
  real, allocatable :: qv(:,:,:)
  real, allocatable :: tund(:,:,:)

  ! Input/output arrays
  real, allocatable :: pt(:,:,:)
  real, allocatable :: pt_init(:,:,:)
  real, allocatable :: ps(:,:)
  real, allocatable :: ps_init(:,:)
  real, allocatable :: qvsts(:,:)
  real, allocatable :: qvsts_init(:,:)
  real, allocatable :: qvsice(:,:)
  real, allocatable :: qvsice_init(:,:)

  ! Output arrays
  real, allocatable :: p(:,:,:)
  real, allocatable :: t(:,:,:)
  real, allocatable :: ptv(:,:,:)
  real, allocatable :: qvsfc(:,:)
  real, allocatable :: tice(:,:)
  real, allocatable :: va(:,:)

  ! Reference outputs for validation
  real, allocatable :: p_ref(:,:,:)
  real, allocatable :: t_ref(:,:,:)
  real, allocatable :: ptv_ref(:,:,:)
  real, allocatable :: qvsfc_ref(:,:)
  real, allocatable :: tice_ref(:,:)
  real, allocatable :: va_ref(:,:)

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
       fmois, levpbl, tubopt, cphopt, ni, nj, nk, nund, rddvcp, p0iv)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel Benchmark: setsfc'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,A)') ' fmois=', trim(fmois)
  write(*,'(A,I6,A,I6,A,I6)') ' levpbl=', levpbl, ', tubopt=', tubopt, ', cphopt=', cphopt
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A,I6)') ' OpenMP threads: ', omp_get_max_threads()
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(land(0:ni+1, 0:nj+1))
  allocate(beta(0:ni+1, 0:nj+1))
  allocate(kai(0:ni+1, 0:nj+1))
  allocate(fall(0:ni+1, 0:nj+1))
  allocate(pbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptbr(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(w(0:ni+1, 0:nj+1, 1:nk))
  allocate(pp(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptp(0:ni+1, 0:nj+1, 1:nk))
  allocate(qv(0:ni+1, 0:nj+1, 1:nk))
  allocate(tund(0:ni+1, 0:nj+1, 1:nund))
  allocate(pt(0:ni+1, 0:nj+1, 1:nk))
  allocate(pt_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(ps(0:ni+1, 0:nj+1))
  allocate(ps_init(0:ni+1, 0:nj+1))
  allocate(qvsts(0:ni+1, 0:nj+1))
  allocate(qvsts_init(0:ni+1, 0:nj+1))
  allocate(qvsice(0:ni+1, 0:nj+1))
  allocate(qvsice_init(0:ni+1, 0:nj+1))
  allocate(p(0:ni+1, 0:nj+1, 1:nk))
  allocate(t(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptv(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsfc(0:ni+1, 0:nj+1))
  allocate(tice(0:ni+1, 0:nj+1))
  allocate(va(0:ni+1, 0:nj+1))
  allocate(p_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(t_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(ptv_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(qvsfc_ref(0:ni+1, 0:nj+1))
  allocate(tice_ref(0:ni+1, 0:nj+1))
  allocate(va_ref(0:ni+1, 0:nj+1))
  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'
  call read_array_2d_int(trim(data_dir)//'/land.bin', land, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/beta.bin', beta, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/kai.bin', kai, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/fall.bin', fall, 0, ni+1, 0, nj+1)
  call read_array_3d(trim(data_dir)//'/pbr.bin', pbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptbr.bin', ptbr, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/w.bin', w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/pp.bin', pp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptp.bin', ptp, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/qv.bin', qv, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tund.bin', tund, 0, ni+1, 0, nj+1, 1, nund)
  call read_array_3d(trim(data_dir)//'/pt.bin', pt_init, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/ps_in.bin', ps_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/qvsts_in.bin', qvsts_init, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/qvsice_in.bin', qvsice_init, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/p.bin', p_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/t.bin', t_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/ptv_ref.bin', ptv_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_2d(trim(data_dir)//'/qvsfc_ref.bin', qvsfc_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/tice_ref.bin', tice_ref, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/va_ref.bin', va_ref, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    pt = pt_init
    ps = ps_init
    qvsts = qvsts_init
    qvsice = qvsice_init

    call kernel_setsfc(fmois, levpbl, tubopt, cphopt, &
         ni, nj, nk, nund, rddvcp, p0iv, &
         pbr, ptbr, u, v, w, pp, ptp, qv, &
         land, beta, kai, tund, fall, &
         p, t, ptv, qvsfc, tice, va, pt, ps, qvsts, qvsice)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    pt = pt_init
    ps = ps_init
    qvsts = qvsts_init
    qvsice = qvsice_init

    t_start = omp_get_wtime()

    call kernel_setsfc(fmois, levpbl, tubopt, cphopt, &
         ni, nj, nk, nund, rddvcp, p0iv, &
         pbr, ptbr, u, v, w, pp, ptp, qv, &
         land, beta, kai, tund, fall, &
         p, t, ptv, qvsfc, tice, va, pt, ps, qvsts, qvsice)

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

  ! Validate va (surface velocity magnitude)
  do j = 1, nj-1
    do i = 1, ni-1
      rel_error = abs(va(i,j) - va_ref(i,j))
      tol_val = abs(va_ref(i,j))
      if (tol_val > 1.0e-10) then
        rel_error = rel_error / tol_val
      end if
      if (rel_error > max_error) max_error = rel_error
      if (rel_error > tolerance) error_count = error_count + 1
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
  deallocate(land, beta, kai, fall)
  deallocate(pbr, ptbr, u, v, w, pp, ptp, qv, tund)
  deallocate(pt, pt_init, ps, ps_init, qvsts, qvsts_init, qvsice, qvsice_init)
  deallocate(p, t, ptv, qvsfc, tice, va)
  deallocate(p_ref, t_ref, ptv_ref, qvsfc_ref, tice_ref, va_ref)
  deallocate(times)

contains

  !-----------------------------------------------------------------
  ! Kernel: setsfc
  !-----------------------------------------------------------------
  subroutine kernel_setsfc(fmois, levpbl, tubopt, cphopt, &
       ni, nj, nk, nund, rddvcp, p0iv, &
       pbr, ptbr, u, v, w, pp, ptp, qv, &
       land, beta, kai, tund, fall, &
       p, t, ptv, qvsfc, tice, va, pt, ps, qvsts, qvsice)
    implicit none

    character(len=5), intent(in) :: fmois
    integer, intent(in) :: levpbl, tubopt, cphopt
    integer, intent(in) :: ni, nj, nk, nund
    real, intent(in) :: rddvcp, p0iv
    real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
    real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
    integer, intent(in) :: land(0:ni+1,0:nj+1)
    real, intent(in) :: beta(0:ni+1,0:nj+1)
    real, intent(in) :: kai(0:ni+1,0:nj+1)
    real, intent(in) :: tund(0:ni+1,0:nj+1,1:nund)
    real, intent(in) :: fall(0:ni+1,0:nj+1)
    real, intent(out) :: p(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: t(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: ptv(0:ni+1,0:nj+1,1:nk)
    real, intent(out) :: qvsfc(0:ni+1,0:nj+1)
    real, intent(out) :: tice(0:ni+1,0:nj+1)
    real, intent(out) :: va(0:ni+1,0:nj+1)
    real, intent(inout) :: pt(0:ni+1,0:nj+1,1:nk)
    real, intent(inout) :: ps(0:ni+1,0:nj+1)
    real, intent(inout) :: qvsts(0:ni+1,0:nj+1)
    real, intent(inout) :: qvsice(0:ni+1,0:nj+1)

    integer :: i, j, k
    real :: ests, esice, u8s, v8s, w8s

    !$omp parallel default(shared) private(k)

    ! Get the pressure, potential temperature and air temperature
    do k=1,nk-1
      !$omp do schedule(runtime) private(i,j)
      do j=1,nj-1
        do i=1,ni-1
          p(i,j,k)=pbr(i,j,k)+pp(i,j,k)
          pt(i,j,k)=ptbr(i,j,k)+ptp(i,j,k)
          t(i,j,k)=pt(i,j,k)*exp(rddvcp*log(p0iv*p(i,j,k)))
        end do
      end do
      !$omp end do
    end do

    ! Get the surface pressure and ice surface temperature
    !$omp do schedule(runtime) private(i,j)
    do j=1,nj-1
      do i=1,ni-1
        ps(i,j)=.5e0*(p(i,j,1)+p(i,j,2))
        if(land(i,j).eq.1) then
          tice(i,j)=min(.5e0*(t(i,j,1)+t(i,j,2)),t0)
        else
          tice(i,j)=lim35n
        end if
      end do
    end do
    !$omp end do

    ! Calculate water vapor mixing ratio on the surface
    if(fmois(1:3).eq.'dry') then
      !$omp do schedule(runtime) private(i,j)
      do j=1,nj-1
        do i=1,ni-1
          qvsfc(i,j)=0.e0
        end do
      end do
      !$omp end do
    else if(fmois(1:5).eq.'moist') then
      !$omp do schedule(runtime) private(i,j,ests,esice)
      do j=1,nj-1
        do i=1,ni-1
          if(land(i,j).lt.0) then
            ests=es0*exp(17.269e0*(tund(i,j,1)-t0)/(tund(i,j,1)-35.86e0))
            qvsts(i,j)=epsva*ests/(ps(i,j)-ests)
          else if(land(i,j).eq.1) then
            if(tice(i,j).gt.tlow) then
              ests=es0*exp(17.269e0*(tund(i,j,1)-t0)/(tund(i,j,1)-35.86e0))
              esice=es0*exp(17.269e0*(tice(i,j)-t0)/(tice(i,j)-35.86e0))
              qvsts(i,j)=epsva*ests/(ps(i,j)-ests)
              qvsice(i,j)=epsva*esice/(ps(i,j)-esice)
            else
              ests=es0*exp(17.269e0*(tund(i,j,1)-t0)/(tund(i,j,1)-35.86e0))
              esice=es0*exp(21.875e0*(tice(i,j)-t0)/(tice(i,j)-7.66e0))
              qvsts(i,j)=epsva*ests/(ps(i,j)-ests)
              qvsice(i,j)=epsva*esice/(ps(i,j)-esice)
            end if
          else
            if(tund(i,j,1).gt.tlow) then
              ests=es0*exp(17.269e0*(tund(i,j,1)-t0)/(tund(i,j,1)-35.86e0))
              qvsts(i,j)=epsva*ests/(ps(i,j)-ests)
            else
              ests=es0*exp(21.875e0*(tund(i,j,1)-t0)/(tund(i,j,1)-7.66e0))
              qvsts(i,j)=epsva*ests/(ps(i,j)-ests)
            end if
          end if
        end do
      end do
      !$omp end do

      if(abs(cphopt).eq.0) then
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            if(land(i,j).eq.1) then
              qvsfc(i,j)=qv(i,j,2)+icbeta*kai(i,j)*(qvsice(i,j)-qv(i,j,2)) &
                   +beta(i,j)*(1.e0-kai(i,j))*(qvsts(i,j)-qv(i,j,2))
            else
              qvsfc(i,j)=beta(i,j)*(qvsts(i,j)-qv(i,j,2))+qv(i,j,2)
            end if
          end do
        end do
        !$omp end do
      else
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            if(land(i,j).eq.1) then
              if(fall(i,j).gt.0.e0) then
                qvsfc(i,j)=kai(i,j)*qvsice(i,j)+(1.e0-kai(i,j))*qvsts(i,j)
              else
                qvsfc(i,j)=qv(i,j,2)+icbeta*kai(i,j)*(qvsice(i,j)-qv(i,j,2)) &
                     +beta(i,j)*(1.e0-kai(i,j))*(qvsts(i,j)-qv(i,j,2))
              end if
            else
              if(fall(i,j).gt.0.e0) then
                qvsfc(i,j)=beta(i,j)*(qvsts(i,j)-qv(i,j,2))+qv(i,j,2)
              else
                qvsfc(i,j)=beta(i,j)*(qvsts(i,j)-qv(i,j,2))+qv(i,j,2)
              end if
            end if
          end do
        end do
        !$omp end do
      end if
    end if

    ! Calculate the virtual potential temperature
    if(fmois(1:3).eq.'dry') then
      !$omp do schedule(runtime) private(i,j)
      do j=1,nj-1
        do i=1,ni-1
          if(land(i,j).eq.1) then
            ptv(i,j,1)=exp(rddvcp*log(p0/ps(i,j))) &
                 *(kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1))
          else
            ptv(i,j,1)=exp(rddvcp*log(p0/ps(i,j)))*tund(i,j,1)
          end if
        end do
      end do
      !$omp end do
      do k=2,levpbl+1
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            ptv(i,j,k)=pt(i,j,k)
          end do
        end do
        !$omp end do
      end do
    else if(fmois(1:5).eq.'moist') then
      !$omp do schedule(runtime) private(i,j)
      do j=1,nj-1
        do i=1,ni-1
          if(land(i,j).eq.1) then
            ptv(i,j,1)=exp(rddvcp*log(p0/ps(i,j))) &
                 *(1.e0+epsav*qvsfc(i,j))/(1.e0+qvsfc(i,j)) &
                 *(kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1))
          else
            ptv(i,j,1)=exp(rddvcp*log(p0/ps(i,j))) &
                 *tund(i,j,1)*(1.e0+epsav*qvsfc(i,j))/(1.e0+qvsfc(i,j))
          end if
        end do
      end do
      !$omp end do
      do k=2,levpbl+1
        !$omp do schedule(runtime) private(i,j)
        do j=1,nj-1
          do i=1,ni-1
            ptv(i,j,k)=pt(i,j,k)*(1.e0+epsav*qv(i,j,k))/(1.e0+qv(i,j,k))
          end do
        end do
        !$omp end do
      end do
    end if

    ! Calculate the magnitude of velocity
    if(tubopt.eq.0) then
      !$omp do schedule(runtime) private(i,j,u8s,v8s)
      do j=1,nj-1
        do i=1,ni-1
          u8s=u(i,j,2)+u(i+1,j,2)
          v8s=v(i,j,2)+v(i,j+1,2)
          va(i,j)=max(.5e0*sqrt(u8s*u8s+v8s*v8s),vamin)
        end do
      end do
      !$omp end do
    else
      !$omp do schedule(runtime) private(i,j,u8s,v8s,w8s)
      do j=1,nj-1
        do i=1,ni-1
          u8s=u(i,j,2)+u(i+1,j,2)
          v8s=v(i,j,2)+v(i,j+1,2)
          w8s=w(i,j,2)+w(i,j,3)
          va(i,j)=max(.5e0*sqrt(u8s*u8s+v8s*v8s+w8s*w8s),vamin)
        end do
      end do
      !$omp end do
    end if

    !$omp end parallel

  end subroutine kernel_setsfc

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
  ! Read parameters from file (by variable name)
  !-----------------------------------------------------------------
  subroutine read_parameters(filename, fmois, levpbl, tubopt, cphopt, &
       ni, nj, nk, nund, rddvcp, p0iv)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=5), intent(out) :: fmois
    integer, intent(out) :: levpbl, tubopt, cphopt, ni, nj, nk, nund
    real, intent(out) :: rddvcp, p0iv

    character(len=256) :: line
    character(len=64) :: varname
    integer :: ios, eq_pos

    ! Initialize defaults
    fmois = 'moist'
    levpbl = 1
    tubopt = 0
    cphopt = 0
    ni = 0
    nj = 0
    nk = 0
    nund = 1
    rddvcp = 0.0
    p0iv = 0.0

    open(unit=10, file=filename, status='old', action='read')
    do
      read(10,'(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        varname = trim(adjustl(line(1:eq_pos-1)))
        select case (trim(varname))
        case ('fmois')
          fmois = trim(adjustl(line(eq_pos+1:)))
        case ('levpbl')
          read(line(eq_pos+1:),*) levpbl
        case ('tubopt')
          read(line(eq_pos+1:),*) tubopt
        case ('cphopt')
          read(line(eq_pos+1:),*) cphopt
        case ('ni')
          read(line(eq_pos+1:),*) ni
        case ('nj')
          read(line(eq_pos+1:),*) nj
        case ('nk')
          read(line(eq_pos+1:),*) nk
        case ('nund')
          read(line(eq_pos+1:),*) nund
        case ('rddvcp')
          read(line(eq_pos+1:),*) rddvcp
        case ('p0iv')
          read(line(eq_pos+1:),*) p0iv
        end select
      end if
    end do
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

end program kernel_benchmark_setsfc

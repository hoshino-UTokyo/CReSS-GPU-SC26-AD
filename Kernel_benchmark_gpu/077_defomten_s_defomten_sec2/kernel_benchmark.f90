!***********************************************************************
! GPU Kernel Benchmark: defomten section 2 (s_defomten_sec2)
!***********************************************************************
!
! Description: Calculate all components of the deformation tensor:
!              s11, s22, s33, s12, s13, s23, s31, s32.
!              This is the second parallel region of defomten.
!              GPU version using OpenACC with Unified Memory.
!
!***********************************************************************
program kernel_benchmark_defomten_sec2
  use omp_lib
  implicit none

  ! Array dimensions
  integer :: ni, nj, nk

  ! Scalar parameters
  integer :: trnopt, mfcopt, mpopt
  integer :: iwest, ieast, jsouth, jnorth
  real :: dxiv, dyiv, dziv, dxv25, dyv25, dzv25, dzv125

  ! Input 3D arrays (0:ni+1, 0:nj+1, 1:nk)
  real, allocatable :: j31(:,:,:), j32(:,:,:)
  real, allocatable :: jcb(:,:,:), jcb8u(:,:,:), jcb8v(:,:,:), jcb8w(:,:,:)
  real, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)

  ! Input 2D arrays (0:ni+1, 0:nj+1)
  real, allocatable :: mf(:,:), mf8u(:,:), mf8v(:,:)

  ! Output/temp 3D arrays (0:ni+1, 0:nj+1, 1:nk)
  real, allocatable :: s11(:,:,:), s22(:,:,:), s33(:,:,:)
  real, allocatable :: s12(:,:,:), s13(:,:,:), s23(:,:,:)
  real, allocatable :: s31(:,:,:), s32(:,:,:)
  real, allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:), tmp4(:,:,:)

  ! Reference outputs for validation
  real, allocatable :: s11_ref(:,:,:), s22_ref(:,:,:), s33_ref(:,:,:)
  real, allocatable :: s12_ref(:,:,:), s13_ref(:,:,:), s23_ref(:,:,:)
  real, allocatable :: s31_ref(:,:,:), s32_ref(:,:,:)

  ! Backup arrays for iteration (input copies)
  real, allocatable :: s11_in(:,:,:), s22_in(:,:,:)
  real, allocatable :: s33_in(:,:,:), s12_in(:,:,:)
  real, allocatable :: s13_in(:,:,:), s23_in(:,:,:)
  real, allocatable :: s31_in(:,:,:), s32_in(:,:,:)
  real, allocatable :: tmp1_in(:,:,:), tmp2_in(:,:,:)
  real, allocatable :: tmp3_in(:,:,:), tmp4_in(:,:,:)

  ! Benchmark control
  integer :: num_iterations, warmup_iterations
  character(len=256) :: data_dir

  ! Timing
  real(8) :: t_start, t_end, t_total, t_avg
  real(8), allocatable :: times(:)

  ! Validation
  real :: max_err, rel_err, tolerance
  real :: max_err_s11, max_err_s22, max_err_s33
  real :: max_err_s12, max_err_s13, max_err_s23
  real :: max_err_s31, max_err_s32
  integer :: err_cnt_s11, err_cnt_s22, err_cnt_s33
  integer :: err_cnt_s12, err_cnt_s13, err_cnt_s23
  integer :: err_cnt_s31, err_cnt_s32
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
       ni, nj, nk, trnopt, mfcopt, mpopt, &
       iwest, ieast, jsouth, jnorth, &
       dxiv, dyiv, dziv, dxv25, dyv25, dzv25, dzv125)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' GPU Kernel Benchmark: defomten_sec2'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6,A,I6)') ' trnopt=', trnopt, ', mfcopt=', mfcopt, ', mpopt=', mpopt
  write(*,'(A,I6,A,I6,A,I6,A,I6)') ' iwest=', iwest, ', ieast=', ieast, &
       ', jsouth=', jsouth, ', jnorth=', jnorth
  write(*,'(A,I6)') ' Warmup iterations: ', warmup_iterations
  write(*,'(A,I6)') ' Benchmark iterations: ', num_iterations
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  ! 3D input arrays
  allocate(j31(0:ni+1, 0:nj+1, 1:nk))
  allocate(j32(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(w(0:ni+1, 0:nj+1, 1:nk))

  ! 2D input arrays
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(mf8u(0:ni+1, 0:nj+1))
  allocate(mf8v(0:ni+1, 0:nj+1))

  ! Output/temp 3D arrays
  allocate(s11(0:ni+1, 0:nj+1, 1:nk))
  allocate(s22(0:ni+1, 0:nj+1, 1:nk))
  allocate(s33(0:ni+1, 0:nj+1, 1:nk))
  allocate(s12(0:ni+1, 0:nj+1, 1:nk))
  allocate(s13(0:ni+1, 0:nj+1, 1:nk))
  allocate(s23(0:ni+1, 0:nj+1, 1:nk))
  allocate(s31(0:ni+1, 0:nj+1, 1:nk))
  allocate(s32(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp4(0:ni+1, 0:nj+1, 1:nk))

  ! Reference outputs
  allocate(s11_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(s22_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(s33_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(s12_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(s13_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(s23_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(s31_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(s32_ref(0:ni+1, 0:nj+1, 1:nk))

  ! Backup arrays for iteration reset
  allocate(s11_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(s22_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(s33_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(s12_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(s13_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(s23_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(s31_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(s32_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp4_in(0:ni+1, 0:nj+1, 1:nk))

  allocate(times(num_iterations))

  !---------------------------------------------------------------------
  ! Read input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading input data...'

  ! 3D input arrays
  call read_array_3d(trim(data_dir)//'/j31.bin', j31, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/j32.bin', j32, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb.bin', jcb, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8u.bin', jcb8u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8v.bin', jcb8v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/jcb8w.bin', jcb8w, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/w.bin', w, 0, ni+1, 0, nj+1, 1, nk)

  ! 2D input arrays
  call read_array_2d(trim(data_dir)//'/mf.bin', mf, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/mf8u.bin', mf8u, 0, ni+1, 0, nj+1)
  call read_array_2d(trim(data_dir)//'/mf8v.bin', mf8v, 0, ni+1, 0, nj+1)

  ! Input state arrays (s11, s22 are pre-computed when trnopt>=1)
  call read_array_3d(trim(data_dir)//'/s11_in.bin', s11_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s22_in.bin', s22_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s33_in.bin', s33_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s12_in.bin', s12_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s13_in.bin', s13_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s23_in.bin', s23_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s31_in.bin', s31_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s32_in.bin', s32_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp1_in.bin', tmp1_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp2_in.bin', tmp2_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp3_in.bin', tmp3_in, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/tmp4_in.bin', tmp4_in, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Read reference outputs
  !---------------------------------------------------------------------
  write(*,'(A)') ' Loading reference output...'
  call read_array_3d(trim(data_dir)//'/s11_ref.bin', s11_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s22_ref.bin', s22_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s33_ref.bin', s33_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s12_ref.bin', s12_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s13_ref.bin', s13_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s23_ref.bin', s23_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s31_ref.bin', s31_ref, 0, ni+1, 0, nj+1, 1, nk)
  call read_array_3d(trim(data_dir)//'/s32_ref.bin', s32_ref, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Warmup iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running warmup iterations...'
  do iter = 1, warmup_iterations
    s11 = s11_in; s22 = s22_in; s33 = s33_in; s12 = s12_in
    s13 = s13_in; s23 = s23_in; s31 = s31_in; s32 = s32_in
    tmp1 = tmp1_in; tmp2 = tmp2_in; tmp3 = tmp3_in; tmp4 = tmp4_in
    call kernel_defomten_sec2(ni, nj, nk, &
         trnopt, mfcopt, mpopt, iwest, ieast, jsouth, jnorth, &
         dxiv, dyiv, dziv, dxv25, dyv25, dzv25, dzv125, &
         j31, j32, jcb, jcb8u, jcb8v, jcb8w, mf, mf8u, mf8v, &
         u, v, w, &
         s11, s22, s33, s12, s13, s23, s31, s32, &
         tmp1, tmp2, tmp3, tmp4)
  end do

  !---------------------------------------------------------------------
  ! Benchmark iterations
  !---------------------------------------------------------------------
  write(*,'(A)') ' Running benchmark iterations...'
  t_total = 0.0d0

  do iter = 1, num_iterations
    s11 = s11_in; s22 = s22_in; s33 = s33_in; s12 = s12_in
    s13 = s13_in; s23 = s23_in; s31 = s31_in; s32 = s32_in
    tmp1 = tmp1_in; tmp2 = tmp2_in; tmp3 = tmp3_in; tmp4 = tmp4_in

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_defomten_sec2(ni, nj, nk, &
         trnopt, mfcopt, mpopt, iwest, ieast, jsouth, jnorth, &
         dxiv, dyiv, dziv, dxv25, dyv25, dzv25, dzv125, &
         j31, j32, jcb, jcb8u, jcb8v, jcb8w, mf, mf8u, mf8v, &
         u, v, w, &
         s11, s22, s33, s12, s13, s23, s31, s32, &
         tmp1, tmp2, tmp3, tmp4)
    !$acc wait
    t_end = omp_get_wtime()
    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
  end do

  t_avg = t_total / dble(num_iterations)

  !---------------------------------------------------------------------
  ! Validate output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Validating output...'

  call validate_array(s11, s11_ref, ni, nj, nk, tolerance, max_err_s11, err_cnt_s11)
  call validate_array(s22, s22_ref, ni, nj, nk, tolerance, max_err_s22, err_cnt_s22)
  call validate_array(s33, s33_ref, ni, nj, nk, tolerance, max_err_s33, err_cnt_s33)
  call validate_array(s12, s12_ref, ni, nj, nk, tolerance, max_err_s12, err_cnt_s12)
  call validate_array(s13, s13_ref, ni, nj, nk, tolerance, max_err_s13, err_cnt_s13)
  call validate_array(s23, s23_ref, ni, nj, nk, tolerance, max_err_s23, err_cnt_s23)
  call validate_array(s31, s31_ref, ni, nj, nk, tolerance, max_err_s31, err_cnt_s31)
  call validate_array(s32, s32_ref, ni, nj, nk, tolerance, max_err_s32, err_cnt_s32)

  validation_passed = (err_cnt_s11 == 0 .and. err_cnt_s22 == 0 .and. &
       err_cnt_s33 == 0 .and. err_cnt_s12 == 0 .and. &
       err_cnt_s13 == 0 .and. err_cnt_s23 == 0 .and. &
       err_cnt_s31 == 0 .and. err_cnt_s32 == 0)

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
  write(*,'(A,ES12.4,A,I8)') ' s11 max_err: ', max_err_s11, '  err_cnt: ', err_cnt_s11
  write(*,'(A,ES12.4,A,I8)') ' s22 max_err: ', max_err_s22, '  err_cnt: ', err_cnt_s22
  write(*,'(A,ES12.4,A,I8)') ' s33 max_err: ', max_err_s33, '  err_cnt: ', err_cnt_s33
  write(*,'(A,ES12.4,A,I8)') ' s12 max_err: ', max_err_s12, '  err_cnt: ', err_cnt_s12
  write(*,'(A,ES12.4,A,I8)') ' s13 max_err: ', max_err_s13, '  err_cnt: ', err_cnt_s13
  write(*,'(A,ES12.4,A,I8)') ' s23 max_err: ', max_err_s23, '  err_cnt: ', err_cnt_s23
  write(*,'(A,ES12.4,A,I8)') ' s31 max_err: ', max_err_s31, '  err_cnt: ', err_cnt_s31
  write(*,'(A,ES12.4,A,I8)') ' s32 max_err: ', max_err_s32, '  err_cnt: ', err_cnt_s32
  write(*,'(A,ES12.4)') ' Tolerance:   ', tolerance
  if (validation_passed) then
    write(*,'(A)') ' Validation: PASSED'
  else
    write(*,'(A)') ' Validation: FAILED'
  end if
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(j31, j32, jcb, jcb8u, jcb8v, jcb8w)
  deallocate(u, v, w, mf, mf8u, mf8v)
  deallocate(s11, s22, s33, s12, s13, s23, s31, s32)
  deallocate(tmp1, tmp2, tmp3, tmp4)
  deallocate(s11_ref, s22_ref, s33_ref, s12_ref)
  deallocate(s13_ref, s23_ref, s31_ref, s32_ref)
  deallocate(s11_in, s22_in, s33_in, s12_in)
  deallocate(s13_in, s23_in, s31_in, s32_in)
  deallocate(tmp1_in, tmp2_in, tmp3_in, tmp4_in)
  deallocate(times)

  if (.not. validation_passed) stop 1

contains

  !=====================================================================
  ! Kernel: defomten_sec2 - GPU version using OpenACC
  !=====================================================================
  subroutine kernel_defomten_sec2(ni, nj, nk, &
       trnopt, mfcopt, mpopt, iwest, ieast, jsouth, jnorth, &
       dxiv, dyiv, dziv, dxv25, dyv25, dzv25, dzv125, &
       j31, j32, jcb, jcb8u, jcb8v, jcb8w, mf, mf8u, mf8v, &
       u, v, w, &
       s11, s22, s33, s12, s13, s23, s31, s32, &
       tmp1, tmp2, tmp3, tmp4)
    implicit none

    integer, intent(in) :: ni, nj, nk
    integer, intent(in) :: trnopt, mfcopt, mpopt
    integer, intent(in) :: iwest, ieast, jsouth, jnorth
    real, intent(in) :: dxiv, dyiv, dziv, dxv25, dyv25, dzv25, dzv125

    real, intent(in) :: j31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: j32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: jcb8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: mf(0:ni+1, 0:nj+1)
    real, intent(in) :: mf8u(0:ni+1, 0:nj+1)
    real, intent(in) :: mf8v(0:ni+1, 0:nj+1)
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: w(0:ni+1, 0:nj+1, 1:nk)

    real, intent(inout) :: s11(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: s22(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: s33(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: s12(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: s13(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: s23(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: s31(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: s32(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp3(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp4(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k
    real :: jcbiv2, mfdvj2

    !-------------------------------------------------------------------
    ! Calculate the diagonal and the x-y components of the deformation
    ! tensor.
    !-------------------------------------------------------------------

    ! tmp1 = u * jcb8u
    !$acc kernels
    !$acc loop independent
    do k=1,nk-1
      !$acc loop independent
      do j=jsouth,nj-jnorth
        !$acc loop independent
        do i=iwest,ni+1-ieast
          tmp1(i,j,k)=u(i,j,k)*jcb8u(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! tmp2 = v * jcb8v
    !$acc kernels
    !$acc loop independent
    do k=1,nk-1
      !$acc loop independent
      do j=jsouth,nj+1-jnorth
        !$acc loop independent
        do i=iwest,ni-ieast
          tmp2(i,j,k)=v(i,j,k)*jcb8v(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! tmp3 (and optionally tmp4) for s12 computation
    if(mfcopt.eq.0) then

      !$acc kernels
      !$acc loop independent
      do k=1,nk-1
        !$acc loop independent
        do j=1+jsouth,nj-jnorth
          !$acc loop independent
          do i=1+iwest,ni-ieast
            tmp3(i,j,k)=4.e0/((jcb(i-1,j-1,k)+jcb(i,j,k))             &
     &        +(jcb(i-1,j,k)+jcb(i,j-1,k)))
          end do
        end do
      end do
      !$acc end kernels

    else

      if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

        !$acc kernels
        !$acc loop independent
        do j=1+jsouth,nj-jnorth
          !$acc loop independent
          do i=1+iwest,ni-ieast
            tmp4(i,j,1)=(mf(i-1,j-1)+mf(i,j))+(mf(i-1,j)+mf(i,j-1))
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              tmp3(i,j,k)=4.e0/((jcb(i-1,j-1,k)+jcb(i,j,k))           &
     &          +(jcb(i-1,j,k)+jcb(i,j-1,k)))
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        !$acc loop independent
        do j=1+jsouth,nj-jnorth
          !$acc loop independent
          do i=1+iwest,ni-ieast
            tmp4(i,j,1)=(mf(i-1,j-1)+mf(i,j))+(mf(i-1,j)+mf(i,j-1))
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              tmp3(i,j,k)=tmp4(i,j,1)/((jcb(i-1,j-1,k)+jcb(i,j,k))    &
     &          +(jcb(i-1,j,k)+jcb(i,j-1,k)))
            end do
          end do
        end do
        !$acc end kernels

      end if

    end if

    !-------------------------------------------------------------------
    ! Compute s11, s22, s12 depending on trnopt/mfcopt/mpopt
    !-------------------------------------------------------------------

    if(trnopt.eq.0) then

      if(mfcopt.eq.1.and.(mpopt.eq.0.or.mpopt.eq.10)) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv &
     &          +tmp4(i,j,1)*(tmp2(i,j,k)-tmp2(i-1,j,k))*dxv25)
            end do
          end do
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1,ni-1
              s11(i,j,k)=(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv
              s22(i,j,k)=(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv
            end do
          end do
        end do
        !$acc end kernels

      else if(mfcopt.eq.1.and.mpopt.eq.5) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((tmp2(i,j,k)-tmp2(i-1,j,k))*dxiv &
     &          +tmp4(i,j,1)*(tmp1(i,j,k)-tmp1(i,j-1,k))*dyv25)
            end do
          end do
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1,ni-1
              s11(i,j,k)=(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv
              s22(i,j,k)=(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv &
     &          +(tmp2(i,j,k)-tmp2(i-1,j,k))*dxiv)
            end do
          end do
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1,ni-1
              s11(i,j,k)=(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv
              s22(i,j,k)=(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv
            end do
          end do
        end do
        !$acc end kernels

      end if

    else if(trnopt.ge.1) then

      if(mfcopt.eq.1.and.(mpopt.eq.0.or.mpopt.eq.10)) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s33(i,j,k)=.25e0*tmp4(i,j,1)                              &
     &          *(s22(i-1,j,k)+s22(i,j,k))*(j31(i,j-1,k)+j31(i,j,k))  &
     &          +(s11(i,j-1,k)+s11(i,j,k))*(j32(i-1,j,k)+j32(i,j,k))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((s33(i,j,k+1)-s33(i,j,k))*dzv125 &
     &          +(tmp4(i,j,1)*(tmp2(i,j,k)-tmp2(i-1,j,k))*dxv25        &
     &          +(tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv))
            end do
          end do
        end do
        !$acc end kernels

      else if(mfcopt.eq.1.and.mpopt.eq.5) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s33(i,j,k)=.25e0*tmp4(i,j,1)                              &
     &          *(s11(i,j-1,k)+s11(i,j,k))*(j32(i-1,j,k)+j32(i,j,k))  &
     &          +(s22(i-1,j,k)+s22(i,j,k))*(j31(i,j-1,k)+j31(i,j,k))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((s33(i,j,k+1)-s33(i,j,k))*dzv125 &
     &          +(tmp4(i,j,1)*(tmp1(i,j,k)-tmp1(i,j-1,k))*dyv25        &
     &          +(tmp2(i,j,k)-tmp2(i-1,j,k))*dxiv))
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        !$acc loop independent
        do k=1,nk
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s33(i,j,k)                                                 &
     &          =(s22(i-1,j,k)+s22(i,j,k))*(j31(i,j-1,k)+j31(i,j,k))  &
     &          +(s11(i,j-1,k)+s11(i,j,k))*(j32(i-1,j,k)+j32(i,j,k))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((s33(i,j,k+1)-s33(i,j,k))*dzv125 &
     &          +((tmp2(i,j,k)-tmp2(i-1,j,k))*dxiv                      &
     &          +(tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv))
            end do
          end do
        end do
        !$acc end kernels

      end if

      ! s11*j31, s22*j32
      !$acc kernels
      !$acc loop independent
      do k=1,nk
        !$acc loop independent
        do j=1,nj-1
          !$acc loop independent
          do i=1,ni
            s11(i,j,k)=s11(i,j,k)*j31(i,j,k)
          end do
        end do
        !$acc loop independent
        do j=1,nj
          !$acc loop independent
          do i=1,ni-1
            s22(i,j,k)=s22(i,j,k)*j32(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      ! s13, s23 intermediate
      !$acc kernels
      !$acc loop independent
      do k=1,nk
        !$acc loop independent
        do j=1,nj-1
          !$acc loop independent
          do i=1,ni-1
            s13(i,j,k)=s11(i,j,k)+s11(i+1,j,k)
            s23(i,j,k)=s22(i,j,k)+s22(i,j+1,k)
          end do
        end do
      end do
      !$acc end kernels

      ! Final s11, s22
      !$acc kernels
      !$acc loop independent
      do k=1,nk-1
        !$acc loop independent
        do j=1,nj-1
          !$acc loop independent
          do i=1,ni-1
            s11(i,j,k)=(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv                &
     &        +(s13(i,j,k+1)-s13(i,j,k))*dzv25

            s22(i,j,k)=(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv                &
     &        +(s23(i,j,k+1)-s23(i,j,k))*dzv25
          end do
        end do
      end do
      !$acc end kernels

    end if

    !-------------------------------------------------------------------
    ! Apply jcb scaling and compute s33
    !-------------------------------------------------------------------

    if(mfcopt.eq.0) then

      !$acc kernels
      !$acc loop independent
      do k=1,nk-1
        !$acc loop independent
        do j=1,nj-1
          !$acc loop independent
          do i=1,ni-1
            jcbiv2=2.e0/jcb(i,j,k)
            s11(i,j,k)=jcbiv2*s11(i,j,k)
            s22(i,j,k)=jcbiv2*s22(i,j,k)
            s33(i,j,k)=jcbiv2*(w(i,j,k+1)-w(i,j,k))*dziv
          end do
        end do
      end do
      !$acc end kernels

    else

      if(mpopt.eq.0.or.mpopt.eq.10) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1,ni-1
              jcbiv2=2.e0/jcb(i,j,k)
              s11(i,j,k)=jcbiv2*mf(i,j)*s11(i,j,k)
              s22(i,j,k)=jcbiv2*s22(i,j,k)
              s33(i,j,k)=jcbiv2*(w(i,j,k+1)-w(i,j,k))*dziv
            end do
          end do
        end do
        !$acc end kernels

      else if(mpopt.eq.5) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1,ni-1
              jcbiv2=2.e0/jcb(i,j,k)
              s11(i,j,k)=jcbiv2*s11(i,j,k)
              s22(i,j,k)=jcbiv2*mf(i,j)*s22(i,j,k)
              s33(i,j,k)=jcbiv2*(w(i,j,k+1)-w(i,j,k))*dziv
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1,ni-1
              jcbiv2=2.e0/jcb(i,j,k)
              mfdvj2=jcbiv2*mf(i,j)
              s11(i,j,k)=mfdvj2*s11(i,j,k)
              s22(i,j,k)=mfdvj2*s22(i,j,k)
              s33(i,j,k)=jcbiv2*(w(i,j,k+1)-w(i,j,k))*dziv
            end do
          end do
        end do
        !$acc end kernels

      end if

    end if

    !-------------------------------------------------------------------
    ! Calculate the x-z, y-z, z-x and z-y components of the deformation
    ! tensor.
    !-------------------------------------------------------------------

    ! Set tmp1 = w * jcb8w
    !$acc kernels
    !$acc loop independent
    do k=1,nk
      !$acc loop independent
      do j=jsouth,nj-jnorth
        !$acc loop independent
        do i=iwest,ni-ieast
          tmp1(i,j,k)=w(i,j,k)*jcb8w(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Calculate the x-z and z-x components

    ! tmp3 for x-z
    !$acc kernels
    !$acc loop independent
    do k=1,nk
      !$acc loop independent
      do j=1,nj-1
        !$acc loop independent
        do i=1+iwest,ni-ieast
          tmp3(i,j,k)=2.e0/(jcb8w(i-1,j,k)+jcb8w(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    if(trnopt.eq.0) then

      if(mfcopt.eq.1.and.mpopt.ne.5) then

        !$acc kernels
        !$acc loop independent
        do k=2,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s13(i,j,k)=(mf8u(i,j)*(tmp1(i,j,k)-tmp1(i-1,j,k))*dxiv  &
     &          +(u(i,j,k)-u(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        !$acc loop independent
        do k=2,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s13(i,j,k)=((tmp1(i,j,k)-tmp1(i-1,j,k))*dxiv             &
     &          +(u(i,j,k)-u(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      end if

    else if(trnopt.ge.1) then

      !$acc kernels
      !$acc loop independent
      do k=1,nk
        !$acc loop independent
        do j=1,nj-1
          !$acc loop independent
          do i=1+iwest,ni-ieast
            tmp2(i,j,k)=(w(i-1,j,k)+w(i,j,k))*j31(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      if(mfcopt.eq.1.and.mpopt.ne.5) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1+iwest,ni-ieast
              tmp4(i,j,k)=tmp2(i,j,k)+tmp2(i,j,k+1)
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k=2,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s13(i,j,k)=((u(i,j,k)-u(i,j,k-1))*dziv                   &
     &          +mf8u(i,j)*((tmp1(i,j,k)-tmp1(i-1,j,k))*dxiv           &
     &          +(tmp4(i,j,k)-tmp4(i,j,k-1))*dzv25))*tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1+iwest,ni-ieast
              tmp4(i,j,k)=u(i,j,k)+.25e0*(tmp2(i,j,k)+tmp2(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k=2,nk-1
          !$acc loop independent
          do j=1,nj-1
            !$acc loop independent
            do i=1+iwest,ni-ieast
              s13(i,j,k)=((tmp1(i,j,k)-tmp1(i-1,j,k))*dxiv             &
     &          +(tmp4(i,j,k)-tmp4(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      end if

    end if

    ! s31 = s13
    !$acc kernels
    !$acc loop independent
    do k=2,nk-1
      !$acc loop independent
      do j=1,nj-1
        !$acc loop independent
        do i=1+iwest,ni-ieast
          s31(i,j,k)=s13(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Calculate the y-z and z-y components

    ! tmp3 for y-z
    !$acc kernels
    !$acc loop independent
    do k=1,nk
      !$acc loop independent
      do j=1+jsouth,nj-jnorth
        !$acc loop independent
        do i=1,ni-1
          tmp3(i,j,k)=2.e0/(jcb8w(i,j-1,k)+jcb8w(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    if(trnopt.eq.0) then

      if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

        !$acc kernels
        !$acc loop independent
        do k=2,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1,ni-1
              s23(i,j,k)=(mf8v(i,j)*(tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv  &
     &          +(v(i,j,k)-v(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        !$acc loop independent
        do k=2,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1,ni-1
              s23(i,j,k)=((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv             &
     &          +(v(i,j,k)-v(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      end if

    else if(trnopt.ge.1) then

      !$acc kernels
      !$acc loop independent
      do k=1,nk
        !$acc loop independent
        do j=1+jsouth,nj-jnorth
          !$acc loop independent
          do i=1,ni-1
            tmp2(i,j,k)=(w(i,j-1,k)+w(i,j,k))*j32(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1,ni-1
              tmp4(i,j,k)=tmp2(i,j,k)+tmp2(i,j,k+1)
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k=2,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1,ni-1
              s23(i,j,k)=((v(i,j,k)-v(i,j,k-1))*dziv                   &
     &          +mf8v(i,j)*((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv           &
     &          +(tmp4(i,j,k)-tmp4(i,j,k-1))*dzv25))*tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else

        !$acc kernels
        !$acc loop independent
        do k=1,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1,ni-1
              tmp4(i,j,k)=v(i,j,k)+.25e0*(tmp2(i,j,k)+tmp2(i,j,k+1))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent
        do k=2,nk-1
          !$acc loop independent
          do j=1+jsouth,nj-jnorth
            !$acc loop independent
            do i=1,ni-1
              s23(i,j,k)=((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv             &
     &          +(tmp4(i,j,k)-tmp4(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      end if

    end if

    ! s32 = s23
    !$acc kernels
    !$acc loop independent
    do k=2,nk-1
      !$acc loop independent
      do j=1+jsouth,nj-jnorth
        !$acc loop independent
        do i=1,ni-1
          s32(i,j,k)=s23(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

  end subroutine kernel_defomten_sec2

  !=====================================================================
  ! Validation helper
  !=====================================================================
  subroutine validate_array(arr, ref, ni, nj, nk, tol, max_err, err_cnt)
    implicit none
    integer, intent(in) :: ni, nj, nk
    real, intent(in) :: arr(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: ref(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: tol
    real, intent(out) :: max_err
    integer, intent(out) :: err_cnt

    integer :: i, j, k
    real :: rel_err

    max_err = 0.0
    err_cnt = 0
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 2, ni-2
          rel_err = abs(arr(i,j,k) - ref(i,j,k))
          if (abs(ref(i,j,k)) > 1.0e-10) then
            rel_err = rel_err / abs(ref(i,j,k))
          end if
          if (rel_err > max_err) max_err = rel_err
          if (rel_err > tol) err_cnt = err_cnt + 1
        end do
      end do
    end do
  end subroutine validate_array

  !=====================================================================
  ! Configuration reader
  !=====================================================================
  subroutine read_config(data_dir, num_iter, warmup_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: num_iter, warmup_iter
    real, intent(out) :: tol

    integer :: ios
    logical :: exists

    data_dir = './data'
    num_iter = 10
    warmup_iter = 2
    tol = 1.0e-5

    inquire(file='benchmark.conf', exist=exists)
    if (exists) then
      open(unit=10, file='benchmark.conf', status='old', iostat=ios)
      if (ios == 0) then
        read(10, '(A)', iostat=ios) data_dir
        read(10, *, iostat=ios) num_iter
        read(10, *, iostat=ios) warmup_iter
        read(10, *, iostat=ios) tol
        close(10)
      end if
    end if

  end subroutine read_config

  !=====================================================================
  ! Parameter reader
  !=====================================================================
  subroutine read_parameters(filename, ni, nj, nk, &
       trnopt, mfcopt, mpopt, iwest, ieast, jsouth, jnorth, &
       dxiv, dyiv, dziv, dxv25, dyv25, dzv25, dzv125)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: ni, nj, nk
    integer, intent(out) :: trnopt, mfcopt, mpopt
    integer, intent(out) :: iwest, ieast, jsouth, jnorth
    real, intent(out) :: dxiv, dyiv, dziv, dxv25, dyv25, dzv25, dzv125

    character(len=256) :: line, key, val
    integer :: ios, eq_pos

    ! Defaults
    ni = 1; nj = 1; nk = 1
    trnopt = 0; mfcopt = 0; mpopt = 0
    iwest = 0; ieast = 0; jsouth = 0; jnorth = 0
    dxiv = 0.0; dyiv = 0.0; dziv = 0.0
    dxv25 = 0.0; dyv25 = 0.0; dzv25 = 0.0; dzv125 = 0.0

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
          case ('ni');      read(val, *) ni
          case ('nj');      read(val, *) nj
          case ('nk');      read(val, *) nk
          case ('trnopt');  read(val, *) trnopt
          case ('mfcopt');  read(val, *) mfcopt
          case ('mpopt');   read(val, *) mpopt
          case ('iwest');   read(val, *) iwest
          case ('ieast');   read(val, *) ieast
          case ('jsouth');  read(val, *) jsouth
          case ('jnorth');  read(val, *) jnorth
          case ('dxiv');    read(val, *) dxiv
          case ('dyiv');    read(val, *) dyiv
          case ('dziv');    read(val, *) dziv
          case ('dxv25');   read(val, *) dxv25
          case ('dyv25');   read(val, *) dyv25
          case ('dzv25');   read(val, *) dzv25
          case ('dzv125');  read(val, *) dzv125
        end select
      end if
    end do
    close(10)

  end subroutine read_parameters

  !=====================================================================
  ! Binary 3D array reader
  !=====================================================================
  subroutine read_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(out) :: arr(i1:i2, j1:j2, k1:k2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_3d

  !=====================================================================
  ! Binary 2D array reader
  !=====================================================================
  subroutine read_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(out) :: arr(i1:i2, j1:j2)

    integer :: ios

    open(unit=10, file=filename, status='old', access='stream', &
         form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    read(10) arr
    close(10)

  end subroutine read_array_2d

end program kernel_benchmark_defomten_sec2

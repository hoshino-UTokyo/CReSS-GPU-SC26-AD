!***********************************************************************
! Kernel Benchmark: curveuvw (s_curveuvw) - GPU Version
!***********************************************************************
! Description:
!   Calculate earth curvature forcing terms for u, v, w velocity
!   equations using temporary arrays and map scale factors.
!
! Original location: Src/curveuvw.f90
!***********************************************************************

program kernel_benchmark
  use omp_lib
  implicit none

  ! Parameters
  character(len=256) :: data_dir
  integer :: num_warmup
  integer :: num_iterations
  real :: tolerance

  ! Array dimensions and parameters
  integer :: ni, nj, nk
  integer :: mpopt, mfcopt
  real :: rearth  ! Earth's radius (read from params.txt)

  ! Arrays - 3D inputs (rmf8u, rmf8v have special k dimension 1:3)
  real, allocatable :: rmf8u(:,:,:), rmf8v(:,:,:)
  real, allocatable :: rst(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:)

  ! Arrays - 3D input/output (initial values saved for reset)
  real, allocatable :: ufrc(:,:,:), vfrc(:,:,:), wfrc(:,:,:)
  real, allocatable :: ufrc_init(:,:,:), vfrc_init(:,:,:), wfrc_init(:,:,:)
  real, allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:)
  real, allocatable :: tmp4(:,:,:), tmp5(:,:,:)
  real, allocatable :: tmp1_init(:,:,:), tmp2_init(:,:,:), tmp3_init(:,:,:)
  real, allocatable :: tmp4_init(:,:,:), tmp5_init(:,:,:)

  ! Reference output arrays
  real, allocatable :: ufrc_ref(:,:,:), vfrc_ref(:,:,:), wfrc_ref(:,:,:)
  real, allocatable :: tmp1_ref(:,:,:), tmp2_ref(:,:,:), tmp3_ref(:,:,:)
  real, allocatable :: tmp4_ref(:,:,:), tmp5_ref(:,:,:)

  ! Timing variables
  real(8) :: t_start, t_end, t_total, t_min, t_max, t_avg
  real(8), allocatable :: times(:)

  ! Loop counters and validation
  integer :: iter, ios
  integer :: i, j, k
  real :: max_err, rel_err, ref_max
  logical :: validation_passed

  ! Read configuration
  call read_config()

  ! Read parameters
  call read_parameters()

  ! Allocate arrays
  call allocate_arrays()

  ! Read input data
  call read_input_data()

  ! Read reference output
  call read_reference_data()

  ! Warmup iterations
  write(*,'(A)') '=== Kernel Benchmark: curveuvw (s_curveuvw) GPU ==='
  write(*,'(A,I3,A)') 'Performing ', num_warmup, ' warmup iterations...'

  do iter = 1, num_warmup
    call reset_arrays()
    call kernel_curveuvw()
  end do

  ! Benchmark iterations
  write(*,'(A,I4,A)') 'Performing ', num_iterations, ' benchmark iterations...'
  allocate(times(num_iterations))

  t_total = 0.0d0
  t_min = huge(t_min)
  t_max = 0.0d0

  do iter = 1, num_iterations
    call reset_arrays()

    !$acc wait
    t_start = omp_get_wtime()
    call kernel_curveuvw()
    !$acc wait
    t_end = omp_get_wtime()

    times(iter) = t_end - t_start
    t_total = t_total + times(iter)
    t_min = min(t_min, times(iter))
    t_max = max(t_max, times(iter))
  end do

  t_avg = t_total / dble(num_iterations)

  ! Validate results
  call validate_results(validation_passed, max_err)

  ! Report results
  write(*,'(A)') ''
  write(*,'(A)') '=== Results ==='
  write(*,'(A,ES12.5,A)') 'Total time:   ', t_total, ' seconds'
  write(*,'(A,ES12.5,A)') 'Average time: ', t_avg, ' seconds'
  write(*,'(A,ES12.5,A)') 'Min time:     ', t_min, ' seconds'
  write(*,'(A,ES12.5,A)') 'Max time:     ', t_max, ' seconds'
  write(*,'(A)') ''
  write(*,'(A,L1)') 'Validation passed: ', validation_passed
  write(*,'(A,ES12.5)') 'Maximum relative error: ', max_err
  write(*,'(A)') ''
  write(*,'(A,I6,A,I4,A,I4)') 'Grid size: ni=', ni, ', nj=', nj, ', nk=', nk

  ! Cleanup
  deallocate(times)
  call deallocate_arrays()

contains

  !---------------------------------------------------------------------
  subroutine read_config()
    integer :: unit_num
    unit_num = 10
    open(unit=unit_num, file='benchmark.conf', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'Error: Cannot open benchmark.conf'
      stop 1
    end if
    read(unit_num, '(A)') data_dir
    read(unit_num, *) num_iterations
    read(unit_num, *) num_warmup
    read(unit_num, *) tolerance
    close(unit_num)
  end subroutine read_config

  !---------------------------------------------------------------------
  subroutine read_parameters()
    integer :: unit_num, eq_pos
    character(len=256) :: line, name, val
    unit_num = 11
    open(unit=unit_num, file=trim(data_dir)//'/params.txt', status='old', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') 'Error: Cannot open params.txt'
      stop 1
    end if
    do while (.true.)
      read(unit_num, '(A)', iostat=ios) line
      if (ios /= 0) exit
      eq_pos = index(line, '=')
      if (eq_pos > 0) then
        name = adjustl(line(1:eq_pos-1))
        val = adjustl(line(eq_pos+1:))
        select case (trim(name))
          case ('mpopt')
            read(val, *) mpopt
          case ('mfcopt')
            read(val, *) mfcopt
          case ('ni')
            read(val, *) ni
          case ('nj')
            read(val, *) nj
          case ('nk')
            read(val, *) nk
          case ('rearth')
            read(val, *) rearth
        end select
      end if
    end do
    close(unit_num)
  end subroutine read_parameters

  !---------------------------------------------------------------------
  subroutine allocate_arrays()
    ! 3D input arrays - rmf8u/rmf8v have special k dimension 1:3
    allocate(rmf8u(0:ni+1, 0:nj+1, 1:3))
    allocate(rmf8v(0:ni+1, 0:nj+1, 1:3))
    allocate(rst(0:ni+1, 0:nj+1, 1:nk))
    allocate(u(0:ni+1, 0:nj+1, 1:nk))
    allocate(v(0:ni+1, 0:nj+1, 1:nk))
    allocate(w(0:ni+1, 0:nj+1, 1:nk))

    ! Input/output arrays and their initial copies
    allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
    allocate(ufrc_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(vfrc(0:ni+1, 0:nj+1, 1:nk))
    allocate(vfrc_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(wfrc(0:ni+1, 0:nj+1, 1:nk))
    allocate(wfrc_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp1_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp2_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp3_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp4(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp4_init(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp5(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp5_init(0:ni+1, 0:nj+1, 1:nk))

    ! Reference arrays
    allocate(ufrc_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(vfrc_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(wfrc_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp1_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp2_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp3_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp4_ref(0:ni+1, 0:nj+1, 1:nk))
    allocate(tmp5_ref(0:ni+1, 0:nj+1, 1:nk))
  end subroutine allocate_arrays

  !---------------------------------------------------------------------
  subroutine read_input_data()
    integer :: unit_num
    unit_num = 12

    ! Read 3D input arrays
    open(unit=unit_num, file=trim(data_dir)//'/rmf8u.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading rmf8u.bin'
    read(unit_num) rmf8u
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/rmf8v.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading rmf8v.bin'
    read(unit_num) rmf8v
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/rst.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading rst.bin'
    read(unit_num) rst
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/u.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading u.bin'
    read(unit_num) u
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/v.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading v.bin'
    read(unit_num) v
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/w.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading w.bin'
    read(unit_num) w
    close(unit_num)

    ! Read input/output initial values
    open(unit=unit_num, file=trim(data_dir)//'/ufrc_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading ufrc_in.bin'
    read(unit_num) ufrc_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/vfrc_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading vfrc_in.bin'
    read(unit_num) vfrc_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/wfrc_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading wfrc_in.bin'
    read(unit_num) wfrc_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp1_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp1_in.bin'
    read(unit_num) tmp1_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp2_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp2_in.bin'
    read(unit_num) tmp2_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp3_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp3_in.bin'
    read(unit_num) tmp3_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp4_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp4_in.bin'
    read(unit_num) tmp4_init
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp5_in.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp5_in.bin'
    read(unit_num) tmp5_init
    close(unit_num)
  end subroutine read_input_data

  !---------------------------------------------------------------------
  subroutine read_reference_data()
    integer :: unit_num
    unit_num = 13

    open(unit=unit_num, file=trim(data_dir)//'/ufrc_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading ufrc_ref.bin'
    read(unit_num) ufrc_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/vfrc_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading vfrc_ref.bin'
    read(unit_num) vfrc_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/wfrc_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading wfrc_ref.bin'
    read(unit_num) wfrc_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp1_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp1_ref.bin'
    read(unit_num) tmp1_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp2_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp2_ref.bin'
    read(unit_num) tmp2_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp3_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp3_ref.bin'
    read(unit_num) tmp3_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp4_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp4_ref.bin'
    read(unit_num) tmp4_ref
    close(unit_num)

    open(unit=unit_num, file=trim(data_dir)//'/tmp5_ref.bin', &
         form='unformatted', access='stream', status='old', iostat=ios)
    if (ios /= 0) stop 'Error reading tmp5_ref.bin'
    read(unit_num) tmp5_ref
    close(unit_num)
  end subroutine read_reference_data

  !---------------------------------------------------------------------
  subroutine reset_arrays()
    ufrc = ufrc_init
    vfrc = vfrc_init
    wfrc = wfrc_init
    tmp1 = tmp1_init
    tmp2 = tmp2_init
    tmp3 = tmp3_init
    tmp4 = tmp4_init
    tmp5 = tmp5_init
  end subroutine reset_arrays

  !---------------------------------------------------------------------
  subroutine validate_results(passed, max_error)
    logical, intent(out) :: passed
    real, intent(out) :: max_error
    real :: err, ref_val

    passed = .true.
    max_error = 0.0

    ! Validate ufrc
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          ref_val = abs(ufrc_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(ufrc(i,j,k) - ufrc_ref(i,j,k)) / ref_val
          else
            err = abs(ufrc(i,j,k) - ufrc_ref(i,j,k))
          end if
          max_error = max(max_error, err)
          if (err > tolerance) passed = .false.
        end do
      end do
    end do

    ! Validate vfrc
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          ref_val = abs(vfrc_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(vfrc(i,j,k) - vfrc_ref(i,j,k)) / ref_val
          else
            err = abs(vfrc(i,j,k) - vfrc_ref(i,j,k))
          end if
          max_error = max(max_error, err)
          if (err > tolerance) passed = .false.
        end do
      end do
    end do

    ! Validate wfrc
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          ref_val = abs(wfrc_ref(i,j,k))
          if (ref_val > 1.0e-30) then
            err = abs(wfrc(i,j,k) - wfrc_ref(i,j,k)) / ref_val
          else
            err = abs(wfrc(i,j,k) - wfrc_ref(i,j,k))
          end if
          max_error = max(max_error, err)
          if (err > tolerance) passed = .false.
        end do
      end do
    end do
  end subroutine validate_results

  !---------------------------------------------------------------------
  subroutine deallocate_arrays()
    deallocate(rmf8u, rmf8v)
    deallocate(rst, u, v, w)
    deallocate(ufrc, ufrc_init, ufrc_ref)
    deallocate(vfrc, vfrc_init, vfrc_ref)
    deallocate(wfrc, wfrc_init, wfrc_ref)
    deallocate(tmp1, tmp1_init, tmp1_ref)
    deallocate(tmp2, tmp2_init, tmp2_ref)
    deallocate(tmp3, tmp3_init, tmp3_ref)
    deallocate(tmp4, tmp4_init, tmp4_ref)
    deallocate(tmp5, tmp5_init, tmp5_ref)
  end subroutine deallocate_arrays

  !---------------------------------------------------------------------
  subroutine kernel_curveuvw()
    ! Local variables
    real :: rev125

    ! Set the common used variable
    rev125 = 0.125e0 / rearth

    ! Set the common used array - tmp1, tmp2
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-1
      do j = 1, nj-1
        do i = 1, ni-1
          tmp1(i,j,k) = u(i,j,k) + u(i+1,j,k)
          tmp2(i,j,k) = v(i,j,k) + v(i,j+1,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Set tmp3
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 1, nj-1
        do i = 1, ni-1
          tmp3(i,j,k) = rst(i,j,k) * (w(i,j,k) + w(i,j,k+1))
        end do
      end do
    end do
    !$acc end kernels

    ! For x and y components of velocity - tmp4, tmp5
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 1, ni-1
          tmp4(i,j,k) = tmp1(i,j,k) * tmp3(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 1, nj-1
        do i = 2, ni-2
          tmp5(i,j,k) = tmp2(i,j,k) * tmp3(i,j,k)
        end do
      end do
    end do
    !$acc end kernels

    ! Update ufrc, vfrc
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-2
        do i = 2, ni-1
          ufrc(i,j,k) = ufrc(i,j,k) - rev125 * (tmp4(i-1,j,k) + tmp4(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-2
      do j = 2, nj-1
        do i = 2, ni-2
          vfrc(i,j,k) = vfrc(i,j,k) - rev125 * (tmp5(i,j-1,k) + tmp5(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    ! Calculate the curvature of earth in the z components of velocity equation
    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 1, nk-1
      do j = 2, nj-2
        do i = 2, ni-2
          tmp3(i,j,k) = rst(i,j,k) &
            * (tmp1(i,j,k) * tmp1(i,j,k) + tmp2(i,j,k) * tmp2(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    !$acc kernels
    !$acc loop independent collapse(3)
    do k = 2, nk-1
      do j = 2, nj-2
        do i = 2, ni-2
          wfrc(i,j,k) = wfrc(i,j,k) + rev125 * (tmp3(i,j,k-1) + tmp3(i,j,k))
        end do
      end do
    end do
    !$acc end kernels

    ! Add the terms in the case of turning on the option for map scale factor
    if (mfcopt == 1) then

      if (mpopt == 0 .or. mpopt == 10) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 1, nj-1
            do i = 1, ni-1
              tmp3(i,j,k) = rmf8v(i,j,3) * rst(i,j,k) * tmp1(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else if (mpopt == 5) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 1, nj-1
            do i = 1, ni-1
              tmp3(i,j,k) = rmf8u(i,j,3) * rst(i,j,k) * tmp2(i,j,k)
            end do
          end do
        end do
        !$acc end kernels

      else
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 1, nj-1
            do i = 1, ni-1
              tmp3(i,j,k) = rst(i,j,k) &
                * (rmf8v(i,j,3) * tmp1(i,j,k) - rmf8u(i,j,3) * tmp2(i,j,k))
            end do
          end do
        end do
        !$acc end kernels
      end if

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 1, nj-1
          do i = 2, ni-2
            tmp1(i,j,k) = tmp1(i,j,k) * tmp3(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels
      !$acc loop independent collapse(3)
      do k = 2, nk-2
        do j = 2, nj-2
          do i = 1, ni-1
            tmp2(i,j,k) = tmp2(i,j,k) * tmp3(i,j,k)
          end do
        end do
      end do
      !$acc end kernels

      if (mpopt == 5) then
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) - (tmp2(i-1,j,k) + tmp2(i,j,k))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) + (tmp1(i,j-1,k) + tmp1(i,j,k))
            end do
          end do
        end do
        !$acc end kernels

      else
        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-2
            do i = 2, ni-1
              ufrc(i,j,k) = ufrc(i,j,k) + (tmp2(i-1,j,k) + tmp2(i,j,k))
            end do
          end do
        end do
        !$acc end kernels

        !$acc kernels
        !$acc loop independent collapse(3)
        do k = 2, nk-2
          do j = 2, nj-1
            do i = 2, ni-2
              vfrc(i,j,k) = vfrc(i,j,k) - (tmp1(i,j-1,k) + tmp1(i,j,k))
            end do
          end do
        end do
        !$acc end kernels
      end if

    end if

  end subroutine kernel_curveuvw

end program kernel_benchmark

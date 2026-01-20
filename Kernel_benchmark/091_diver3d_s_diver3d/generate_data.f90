!***********************************************************************
! Data Generator for diver3d Kernel Benchmark
!***********************************************************************
!
! Generates test input data and reference output for validation.
! Run this program once before running the benchmark.
!
!***********************************************************************
program generate_data_diver3d
  implicit none

  ! Array dimensions (configurable)
  integer :: ni, nj, nk

  ! Parameters
  integer :: mpopt, mfcopt
  real :: dxiv, dyiv, dziv

  ! Arrays
  real, allocatable :: mf(:,:)
  real, allocatable :: rmf(:,:,:)
  real, allocatable :: rmf8u(:,:,:)
  real, allocatable :: rmf8v(:,:,:)
  real, allocatable :: var8u(:,:,:)
  real, allocatable :: var8v(:,:,:)
  real, allocatable :: var8w(:,:,:)
  real, allocatable :: u(:,:,:)
  real, allocatable :: v(:,:,:)
  real, allocatable :: wc(:,:,:)
  real, allocatable :: div3d(:,:,:)
  real, allocatable :: tmp1(:,:,:)
  real, allocatable :: tmp2(:,:,:)
  real, allocatable :: tmp3(:,:,:)

  integer :: i, j, k, seed_size
  integer, allocatable :: seed(:)
  character(len=256) :: data_dir

  !---------------------------------------------------------------------
  ! Configuration - adjust these for your test case
  !---------------------------------------------------------------------
  ! Grid size (similar to real simulation)
  ni = 902
  nj = 902
  nk = 128

  ! Physics options
  mpopt = 10    ! Map projection option
  mfcopt = 1    ! Map scale factor option

  ! Grid spacing inverses (typical values)
  dxiv = 1.0 / 1000.0   ! 1/dx [1/m]
  dyiv = 1.0 / 1000.0   ! 1/dy [1/m]
  dziv = 1.0 / 200.0    ! 1/dz [1/m]

  data_dir = './data'

  !---------------------------------------------------------------------
  ! Initialize random seed
  !---------------------------------------------------------------------
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 12345  ! Fixed seed for reproducibility
  call random_seed(put=seed)

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Data Generator for diver3d Kernel Benchmark'
  write(*,'(A)') '=================================================='
  write(*,'(A,I6,A,I6,A,I6)') ' Grid size: ni=', ni, ', nj=', nj, ', nk=', nk
  write(*,'(A,I6,A,I6)') ' Options: mpopt=', mpopt, ', mfcopt=', mfcopt

  !---------------------------------------------------------------------
  ! Create data directory
  !---------------------------------------------------------------------
  call execute_command_line('mkdir -p '//trim(data_dir))

  !---------------------------------------------------------------------
  ! Allocate arrays
  !---------------------------------------------------------------------
  allocate(mf(0:ni+1, 0:nj+1))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))
  allocate(rmf8u(0:ni+1, 0:nj+1, 1:3))
  allocate(rmf8v(0:ni+1, 0:nj+1, 1:3))
  allocate(var8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8v(0:ni+1, 0:nj+1, 1:nk))
  allocate(var8w(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  allocate(v(0:ni+1, 0:nj+1, 1:nk))
  allocate(wc(0:ni+1, 0:nj+1, 1:nk))
  allocate(div3d(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))

  !---------------------------------------------------------------------
  ! Generate input data
  !---------------------------------------------------------------------
  write(*,'(A)') ' Generating input data...'

  ! Map factors (close to 1.0 for small domains)
  call random_number(mf)
  mf = 0.95 + 0.1 * mf  ! Range: 0.95-1.05

  call random_number(rmf)
  rmf = 0.95 + 0.1 * rmf

  call random_number(rmf8u)
  rmf8u = 0.95 + 0.1 * rmf8u

  call random_number(rmf8v)
  rmf8v = 0.95 + 0.1 * rmf8v

  ! Jacobian-like variables (positive values)
  call random_number(var8u)
  var8u = 0.8 + 0.4 * var8u  ! Range: 0.8-1.2

  call random_number(var8v)
  var8v = 0.8 + 0.4 * var8v

  call random_number(var8w)
  var8w = 0.8 + 0.4 * var8w

  ! Velocity components (typical atmospheric values, m/s)
  call random_number(u)
  u = -20.0 + 40.0 * u  ! Range: -20 to 20 m/s

  call random_number(v)
  v = -20.0 + 40.0 * v

  call random_number(wc)
  wc = -2.0 + 4.0 * wc  ! Range: -2 to 2 m/s (smaller vertical velocity)

  ! Initialize output array
  div3d = 0.0
  tmp1 = 0.0
  tmp2 = 0.0
  tmp3 = 0.0

  !---------------------------------------------------------------------
  ! Write parameters file
  !---------------------------------------------------------------------
  write(*,'(A)') ' Writing parameters...'
  open(unit=10, file=trim(data_dir)//'/params.txt', status='replace')
  write(10, '(I8)') ni
  write(10, '(I8)') nj
  write(10, '(I8)') nk
  write(10, '(I8)') mpopt
  write(10, '(I8)') mfcopt
  write(10, '(ES20.12)') dxiv
  write(10, '(ES20.12)') dyiv
  write(10, '(ES20.12)') dziv
  close(10)

  !---------------------------------------------------------------------
  ! Write input arrays (binary format)
  !---------------------------------------------------------------------
  write(*,'(A)') ' Writing input arrays...'
  call write_array_2d(trim(data_dir)//'/mf.bin', mf, 0, ni+1, 0, nj+1)
  call write_array_3d(trim(data_dir)//'/rmf.bin', rmf, 0, ni+1, 0, nj+1, 1, 4)
  call write_array_3d(trim(data_dir)//'/rmf8u.bin', rmf8u, 0, ni+1, 0, nj+1, 1, 3)
  call write_array_3d(trim(data_dir)//'/rmf8v.bin', rmf8v, 0, ni+1, 0, nj+1, 1, 3)
  call write_array_3d(trim(data_dir)//'/var8u.bin', var8u, 0, ni+1, 0, nj+1, 1, nk)
  call write_array_3d(trim(data_dir)//'/var8v.bin', var8v, 0, ni+1, 0, nj+1, 1, nk)
  call write_array_3d(trim(data_dir)//'/var8w.bin', var8w, 0, ni+1, 0, nj+1, 1, nk)
  call write_array_3d(trim(data_dir)//'/u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
  call write_array_3d(trim(data_dir)//'/v.bin', v, 0, ni+1, 0, nj+1, 1, nk)
  call write_array_3d(trim(data_dir)//'/wc.bin', wc, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Compute reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Computing reference output...'
  call compute_diver3d_ref(mpopt, mfcopt, dxiv, dyiv, dziv, ni, nj, nk, &
       mf, rmf, rmf8u, rmf8v, var8u, var8v, var8w, &
       u, v, wc, div3d, tmp1, tmp2, tmp3)

  !---------------------------------------------------------------------
  ! Write reference output
  !---------------------------------------------------------------------
  write(*,'(A)') ' Writing reference output...'
  call write_array_3d(trim(data_dir)//'/div3d_ref.bin', div3d, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Write benchmark configuration
  !---------------------------------------------------------------------
  write(*,'(A)') ' Writing benchmark configuration...'
  open(unit=10, file='benchmark.conf', status='replace')
  write(10, '(A)') trim(data_dir)
  write(10, '(I8)') 10    ! Number of iterations
  write(10, '(I8)') 2     ! Warmup iterations
  write(10, '(ES12.4)') 1.0e-5  ! Tolerance
  close(10)

  !---------------------------------------------------------------------
  ! Summary
  !---------------------------------------------------------------------
  write(*,'(A)') ''
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Data generation complete!'
  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Files created in: '//trim(data_dir)
  write(*,'(A)') '   - params.txt'
  write(*,'(A)') '   - mf.bin, rmf.bin, rmf8u.bin, rmf8v.bin'
  write(*,'(A)') '   - var8u.bin, var8v.bin, var8w.bin'
  write(*,'(A)') '   - u.bin, v.bin, wc.bin'
  write(*,'(A)') '   - div3d_ref.bin'
  write(*,'(A)') ' Config file: benchmark.conf'
  write(*,'(A)') ''
  write(*,'(A)') ' To run the benchmark:'
  write(*,'(A)') '   ./kernel_benchmark'
  write(*,'(A)') '=================================================='

  !---------------------------------------------------------------------
  ! Cleanup
  !---------------------------------------------------------------------
  deallocate(mf, rmf, rmf8u, rmf8v)
  deallocate(var8u, var8v, var8w)
  deallocate(u, v, wc)
  deallocate(div3d, tmp1, tmp2, tmp3)
  deallocate(seed)

contains

  !=====================================================================
  ! Reference computation (serial, for generating correct output)
  !=====================================================================
  subroutine compute_diver3d_ref(mpopt, mfcopt, dxiv, dyiv, dziv, ni, nj, nk, &
       mf, rmf, rmf8u, rmf8v, var8u, var8v, var8w, &
       u, v, wc, div3d, tmp1, tmp2, tmp3)
    implicit none

    integer, intent(in) :: mpopt, mfcopt
    real, intent(in) :: dxiv, dyiv, dziv
    integer, intent(in) :: ni, nj, nk

    real, intent(in) :: mf(0:ni+1, 0:nj+1)
    real, intent(in) :: rmf(0:ni+1, 0:nj+1, 1:4)
    real, intent(in) :: rmf8u(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: rmf8v(0:ni+1, 0:nj+1, 1:3)
    real, intent(in) :: var8u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: var8v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: var8w(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: u(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: v(0:ni+1, 0:nj+1, 1:nk)
    real, intent(in) :: wc(0:ni+1, 0:nj+1, 1:nk)

    real, intent(out) :: div3d(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp1(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp2(0:ni+1, 0:nj+1, 1:nk)
    real, intent(inout) :: tmp3(0:ni+1, 0:nj+1, 1:nk)

    integer :: i, j, k

    ! Phase 1: Multiply velocities by optional variables and map factors
    if (mfcopt == 0) then
      do k = 1, nk-1
        do j = 1, nj-1
          do i = 1, ni
            tmp1(i,j,k) = var8u(i,j,k) * u(i,j,k)
          end do
        end do
        do j = 1, nj
          do i = 1, ni-1
            tmp2(i,j,k) = var8v(i,j,k) * v(i,j,k)
          end do
        end do
      end do
    else
      if (mpopt == 0 .or. mpopt == 10) then
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni
              tmp1(i,j,k) = var8u(i,j,k) * u(i,j,k)
            end do
          end do
          do j = 1, nj
            do i = 1, ni-1
              tmp2(i,j,k) = rmf8v(i,j,2) * var8v(i,j,k) * v(i,j,k)
            end do
          end do
        end do
      else if (mpopt == 5) then
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni
              tmp1(i,j,k) = rmf8u(i,j,2) * var8u(i,j,k) * u(i,j,k)
            end do
          end do
          do j = 1, nj
            do i = 1, ni-1
              tmp2(i,j,k) = var8v(i,j,k) * v(i,j,k)
            end do
          end do
        end do
      else
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni
              tmp1(i,j,k) = rmf8u(i,j,2) * var8u(i,j,k) * u(i,j,k)
            end do
          end do
          do j = 1, nj
            do i = 1, ni-1
              tmp2(i,j,k) = rmf8v(i,j,2) * var8v(i,j,k) * v(i,j,k)
            end do
          end do
        end do
      end if
    end if

    do k = 1, nk
      do j = 1, nj-1
        do i = 1, ni-1
          tmp3(i,j,k) = var8w(i,j,k) * wc(i,j,k)
        end do
      end do
    end do

    ! Phase 2: Calculate divergence
    if (mfcopt == 0) then
      do k = 1, nk-1
        do j = 1, nj-1
          do i = 1, ni-1
            div3d(i,j,k) = (tmp3(i,j,k) - tmp3(i,j,k+1)) * dziv &
                 + ((tmp1(i,j,k) - tmp1(i+1,j,k)) * dxiv &
                 + (tmp2(i,j,k) - tmp2(i,j+1,k)) * dyiv)
          end do
        end do
      end do
    else
      if (mpopt == 0 .or. mpopt == 5 .or. mpopt == 10) then
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni-1
              div3d(i,j,k) = mf(i,j) * ((tmp1(i,j,k) - tmp1(i+1,j,k)) * dxiv &
                   + (tmp2(i,j,k) - tmp2(i,j+1,k)) * dyiv) &
                   + (tmp3(i,j,k) - tmp3(i,j,k+1)) * dziv
            end do
          end do
        end do
      else
        do k = 1, nk-1
          do j = 1, nj-1
            do i = 1, ni-1
              div3d(i,j,k) = rmf(i,j,1) * ((tmp1(i,j,k) - tmp1(i+1,j,k)) * dxiv &
                   + (tmp2(i,j,k) - tmp2(i,j+1,k)) * dyiv) &
                   + (tmp3(i,j,k) - tmp3(i,j,k+1)) * dziv
            end do
          end do
        end do
      end if
    end if

  end subroutine compute_diver3d_ref

  !=====================================================================
  ! Binary array writers
  !=====================================================================
  subroutine write_array_2d(filename, arr, i1, i2, j1, j2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2
    real, intent(in) :: arr(i1:i2, j1:j2)

    open(unit=10, file=filename, status='replace', access='stream', &
         form='unformatted')
    write(10) arr
    close(10)

  end subroutine write_array_2d

  subroutine write_array_3d(filename, arr, i1, i2, j1, j2, k1, k2)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    real, intent(in) :: arr(i1:i2, j1:j2, k1:k2)

    open(unit=10, file=filename, status='replace', access='stream', &
         form='unformatted')
    write(10) arr
    close(10)

  end subroutine write_array_3d

end program generate_data_diver3d

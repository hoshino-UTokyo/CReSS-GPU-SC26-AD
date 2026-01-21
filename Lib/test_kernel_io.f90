!***********************************************************************
! Test program for Kernel I/O Library
!***********************************************************************
program test_kernel_io
  use m_kernel_io
  implicit none

  ! Test arrays
  real, allocatable :: arr1d(:), arr1d_loaded(:)
  real, allocatable :: arr2d(:,:), arr2d_loaded(:,:)
  real, allocatable :: arr3d(:,:,:), arr3d_loaded(:,:,:)
  real, allocatable :: arr4d(:,:,:,:), arr4d_loaded(:,:,:,:)
  integer, allocatable :: iarr1d(:), iarr1d_loaded(:)

  ! Test parameters
  integer :: ni, nj, nk, param_int
  real :: param_real

  ! Test variables
  integer :: i, j, k, l
  real :: maxdiff
  logical :: all_passed

  all_passed = .true.

  write(*,'(A)') '=================================================='
  write(*,'(A)') ' Kernel I/O Library Test Suite'
  write(*,'(A)') '=================================================='
  write(*,'(A)') ''

  !---------------------------------------------------------------------
  ! Test 1: Size check function
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 1: Size check function ---'

  if (kio_check_size(100000000_8)) then
    write(*,'(A)') '  100M elements (400MB): OK (within limit)'
  else
    write(*,'(A)') '  100M elements (400MB): ERROR'
    all_passed = .false.
  end if

  if (.not. kio_check_size(300000000_8)) then
    write(*,'(A)') '  300M elements (1.2GB): OK (exceeds limit)'
  else
    write(*,'(A)') '  300M elements (1.2GB): ERROR'
    all_passed = .false.
  end if

  write(*,'(A,I0,A)') '  KIO_MAX_ELEMENTS_4BYTE = ', KIO_MAX_ELEMENTS_4BYTE, ' elements'
  write(*,'(A,I0,A)') '  KIO_BUFFER_LIMIT_BYTES = ', KIO_BUFFER_LIMIT_BYTES, ' bytes'
  write(*,'(A)') ''

  !---------------------------------------------------------------------
  ! Test 2: Dump and load scalars
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 2: Dump scalars ---'

  ni = 10
  nj = 20
  nk = 5

  call kio_dump_init('test_kernel')
  call kio_dump_scalar_i('ni', ni)
  call kio_dump_scalar_i('nj', nj)
  call kio_dump_scalar_i('nk', nk)
  call kio_dump_scalar_r('dx', 100.0)
  call kio_dump_scalar_d('pi', 3.14159265358979d0)
  call kio_dump_scalar_c('name', 'test_value')
  write(*,'(A)') ''

  !---------------------------------------------------------------------
  ! Test 3: Dump and load 1D array
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 3: Dump 1D array ---'

  allocate(arr1d(1:100))
  do i = 1, 100
    arr1d(i) = real(i) * 0.1
  end do
  call kio_dump_array_1d('arr1d.bin', arr1d, 1, 100)

  !---------------------------------------------------------------------
  ! Test 4: Dump and load 2D array
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 4: Dump 2D array ---'

  allocate(arr2d(0:ni+1, 0:nj+1))
  do j = 0, nj+1
    do i = 0, ni+1
      arr2d(i,j) = real(i*100 + j)
    end do
  end do
  call kio_dump_array_2d('arr2d.bin', arr2d, 0, ni+1, 0, nj+1)

  !---------------------------------------------------------------------
  ! Test 5: Dump and load 3D array
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 5: Dump 3D array ---'

  allocate(arr3d(0:ni+1, 0:nj+1, 1:nk))
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        arr3d(i,j,k) = real(i*10000 + j*100 + k)
      end do
    end do
  end do
  call kio_dump_array_3d('arr3d.bin', arr3d, 0, ni+1, 0, nj+1, 1, nk)

  !---------------------------------------------------------------------
  ! Test 6: Dump and load 4D array (with 2D slice writing)
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 6: Dump 4D array (2D slices) ---'

  allocate(arr4d(0:ni+1, 0:nj+1, 1:nk, 1:3))
  do l = 1, 3
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          arr4d(i,j,k,l) = real(i*1000000 + j*10000 + k*100 + l)
        end do
      end do
    end do
  end do
  call kio_dump_array_4d('arr4d.bin', arr4d, 0, ni+1, 0, nj+1, 1, nk, 1, 3)

  !---------------------------------------------------------------------
  ! Test 7: Dump integer array
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 7: Dump integer array ---'

  allocate(iarr1d(1:50))
  do i = 1, 50
    iarr1d(i) = i * 10
  end do
  call kio_dump_array_1d_int('iarr1d.bin', iarr1d, 1, 50)

  call kio_dump_finalize()
  write(*,'(A)') ''

  !---------------------------------------------------------------------
  ! Test 8: Load and verify parameters
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 8: Load parameters ---'

  call kio_load_params('./kernel_dump/test_kernel/params.txt')

  if (kio_params_loaded()) then
    write(*,'(A)') '  Parameters loaded successfully'
  else
    write(*,'(A)') '  ERROR: Parameters not loaded'
    all_passed = .false.
  end if

  call kio_get_param_i('ni', param_int)
  write(*,'(A,I0)') '  Loaded ni = ', param_int
  if (param_int /= ni) then
    write(*,'(A)') '  ERROR: ni mismatch!'
    all_passed = .false.
  end if

  call kio_get_param_r('dx', param_real)
  write(*,'(A,F10.2)') '  Loaded dx = ', param_real
  if (abs(param_real - 100.0) > 0.001) then
    write(*,'(A)') '  ERROR: dx mismatch!'
    all_passed = .false.
  end if
  write(*,'(A)') ''

  !---------------------------------------------------------------------
  ! Test 9: Load and verify 1D array
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 9: Load and verify 1D array ---'

  allocate(arr1d_loaded(1:100))
  call kio_load_array_1d('./kernel_dump/test_kernel/arr1d.bin', arr1d_loaded, 1, 100)

  maxdiff = 0.0
  do i = 1, 100
    maxdiff = max(maxdiff, abs(arr1d(i) - arr1d_loaded(i)))
  end do
  write(*,'(A,ES12.4)') '  Max difference: ', maxdiff
  if (maxdiff < 1.0e-6) then
    write(*,'(A)') '  1D array verification: PASSED'
  else
    write(*,'(A)') '  1D array verification: FAILED'
    all_passed = .false.
  end if
  write(*,'(A)') ''

  !---------------------------------------------------------------------
  ! Test 10: Load and verify 3D array
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 10: Load and verify 3D array ---'

  allocate(arr3d_loaded(0:ni+1, 0:nj+1, 1:nk))
  call kio_load_array_3d('./kernel_dump/test_kernel/arr3d.bin', arr3d_loaded, &
                         0, ni+1, 0, nj+1, 1, nk)

  maxdiff = 0.0
  do k = 1, nk
    do j = 0, nj+1
      do i = 0, ni+1
        maxdiff = max(maxdiff, abs(arr3d(i,j,k) - arr3d_loaded(i,j,k)))
      end do
    end do
  end do
  write(*,'(A,ES12.4)') '  Max difference: ', maxdiff
  if (maxdiff < 1.0e-6) then
    write(*,'(A)') '  3D array verification: PASSED'
  else
    write(*,'(A)') '  3D array verification: FAILED'
    all_passed = .false.
  end if
  write(*,'(A)') ''

  !---------------------------------------------------------------------
  ! Test 11: Load and verify 4D array (2D slice reading)
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 11: Load and verify 4D array ---'

  allocate(arr4d_loaded(0:ni+1, 0:nj+1, 1:nk, 1:3))
  call kio_load_array_4d('./kernel_dump/test_kernel/arr4d.bin', arr4d_loaded, &
                         0, ni+1, 0, nj+1, 1, nk, 1, 3)

  maxdiff = 0.0
  do l = 1, 3
    do k = 1, nk
      do j = 0, nj+1
        do i = 0, ni+1
          maxdiff = max(maxdiff, abs(arr4d(i,j,k,l) - arr4d_loaded(i,j,k,l)))
        end do
      end do
    end do
  end do
  write(*,'(A,ES12.4)') '  Max difference: ', maxdiff
  if (maxdiff < 1.0e-6) then
    write(*,'(A)') '  4D array verification: PASSED'
  else
    write(*,'(A)') '  4D array verification: FAILED'
    all_passed = .false.
  end if
  write(*,'(A)') ''

  !---------------------------------------------------------------------
  ! Test 12: Load integer array
  !---------------------------------------------------------------------
  write(*,'(A)') '--- Test 12: Load and verify integer array ---'

  allocate(iarr1d_loaded(1:50))
  call kio_load_array_1d_int('./kernel_dump/test_kernel/iarr1d.bin', iarr1d_loaded, 1, 50)

  maxdiff = 0.0
  do i = 1, 50
    maxdiff = max(maxdiff, real(abs(iarr1d(i) - iarr1d_loaded(i))))
  end do
  write(*,'(A,ES12.4)') '  Max difference: ', maxdiff
  if (maxdiff < 0.5) then
    write(*,'(A)') '  Integer array verification: PASSED'
  else
    write(*,'(A)') '  Integer array verification: FAILED'
    all_passed = .false.
  end if
  write(*,'(A)') ''

  !---------------------------------------------------------------------
  ! Summary
  !---------------------------------------------------------------------
  write(*,'(A)') '=================================================='
  if (all_passed) then
    write(*,'(A)') ' All tests PASSED!'
  else
    write(*,'(A)') ' Some tests FAILED!'
  end if
  write(*,'(A)') '=================================================='

  ! Cleanup
  deallocate(arr1d, arr1d_loaded)
  deallocate(arr2d)
  deallocate(arr3d, arr3d_loaded)
  deallocate(arr4d, arr4d_loaded)
  deallocate(iarr1d, iarr1d_loaded)

  if (.not. all_passed) stop 1

end program test_kernel_io

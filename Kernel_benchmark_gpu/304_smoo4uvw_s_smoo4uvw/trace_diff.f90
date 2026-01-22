program trace_diff
  implicit none
  
  integer :: ni, nj, nk
  integer :: smtopt, iwest, ieast, jsouth, jnorth
  real :: smhcoe, smvcoe
  integer :: nkm1, nkm2
  
  real, allocatable :: ufrc(:,:,:), ufrc_init(:,:,:), ufrc_ref(:,:,:)
  real, allocatable :: tmp1(:,:,:), tmp2(:,:,:), tmp3(:,:,:), tmp4(:,:,:), tmp5(:,:,:)
  real, allocatable :: jcb8u(:,:,:), rst8u(:,:,:), ubr(:,:,:), u(:,:,:)
  
  integer :: i, j, k, ios, eq_pos
  character(len=256) :: line, key, val, data_dir
  real :: rel_err
  
  data_dir = './data'
  
  ! Read params
  open(unit=10, file=trim(data_dir)//'/params.txt', status='old')
  do while (.true.)
    read(10, '(A)', iostat=ios) line
    if (ios /= 0) exit
    eq_pos = index(line, '=')
    if (eq_pos > 0) then
      key = adjustl(line(1:eq_pos-1))
      val = adjustl(line(eq_pos+1:))
      select case (trim(key))
        case ('ni'); read(val, *) ni
        case ('nj'); read(val, *) nj
        case ('nk'); read(val, *) nk
        case ('smtopt'); read(val, *) smtopt
        case ('iwest'); read(val, *) iwest
        case ('ieast'); read(val, *) ieast
        case ('jsouth'); read(val, *) jsouth
        case ('jnorth'); read(val, *) jnorth
        case ('smhcoe'); read(val, *) smhcoe
        case ('smvcoe'); read(val, *) smvcoe
      end select
    end if
  end do
  close(10)
  
  nkm1 = nk - 1
  nkm2 = nk - 2
  
  allocate(ufrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc_init(0:ni+1, 0:nj+1, 1:nk))
  allocate(ufrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp1(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp2(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp3(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp4(0:ni+1, 0:nj+1, 1:nk))
  allocate(tmp5(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst8u(0:ni+1, 0:nj+1, 1:nk))
  allocate(ubr(0:ni+1, 0:nj+1, 1:nk))
  allocate(u(0:ni+1, 0:nj+1, 1:nk))
  
  ! Read all data
  open(unit=10, file=trim(data_dir)//'/ufrc_in.bin', access='stream', form='unformatted', status='old')
  read(10) ufrc_init; close(10)
  open(unit=10, file=trim(data_dir)//'/ufrc_ref.bin', access='stream', form='unformatted', status='old')
  read(10) ufrc_ref; close(10)
  open(unit=10, file=trim(data_dir)//'/jcb8u.bin', access='stream', form='unformatted', status='old')
  read(10) jcb8u; close(10)
  open(unit=10, file=trim(data_dir)//'/rst8u.bin', access='stream', form='unformatted', status='old')
  read(10) rst8u; close(10)
  open(unit=10, file=trim(data_dir)//'/ubr.bin', access='stream', form='unformatted', status='old')
  read(10) ubr; close(10)
  open(unit=10, file=trim(data_dir)//'/u.bin', access='stream', form='unformatted', status='old')
  read(10) u; close(10)
  
  ufrc = ufrc_init
  tmp1 = 0.0; tmp2 = 0.0; tmp3 = 0.0; tmp4 = 0.0; tmp5 = 0.0
  
  ! CPU version (serial) - exact same as original
  do k = 1, nk-1
    do j = jsouth, nj-jnorth
      do i = iwest, ni+1-ieast
        tmp4(i,j,k) = rst8u(i,j,k) * (u(i,j,k) - ubr(i,j,k)) / jcb8u(i,j,k)
        tmp5(i,j,k) = 2.0 * tmp4(i,j,k)
      end do
    end do
  end do
  
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 1+iwest, ni-ieast
        tmp1(i,j,k) = (tmp4(i+1,j,k) + tmp4(i-1,j,k)) - tmp5(i,j,k)
      end do
    end do
    
    do j = 1+jsouth, nj-1-jnorth
      do i = 2, ni-1
        tmp2(i,j,k) = (tmp4(i,j+1,k) + tmp4(i,j-1,k)) - tmp5(i,j,k)
      end do
    end do
    
    do j = 2, nj-2
      do i = 2, ni-1
        tmp3(i,j,k) = (tmp4(i,j,k+1) + tmp4(i,j,k-1)) - tmp5(i,j,k)
      end do
    end do
  end do
  
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-1
        ufrc(i,j,k) = ufrc(i,j,k) + smhcoe * (tmp1(i,j,k) + tmp2(i,j,k)) + smvcoe * tmp3(i,j,k)
      end do
    end do
  end do
  
  if (mod(smtopt, 10) /= 2) then
    do j = 2+jsouth, nj-2-jnorth
      do i = 2+iwest, ni-1-ieast
        tmp3(i,j,1) = tmp3(i,j,2)
        tmp3(i,j,nkm1) = tmp3(i,j,nkm2)
      end do
    end do
    
    do k = 2, nk-2
      do j = 2+jsouth, nj-2-jnorth
        do i = 2+iwest, ni-1-ieast
          ufrc(i,j,k) = ufrc(i,j,k) &
               + (smvcoe * (tmp3(i,j,k) - (tmp3(i,j,k+1) + tmp3(i,j,k-1))) &
               + smhcoe * ((tmp1(i,j,k) - (tmp1(i+1,j,k) + tmp1(i-1,j,k))) &
               + (tmp2(i,j,k) - (tmp2(i,j+1,k) + tmp2(i,j-1,k)))))
        end do
      end do
    end do
  end if
  
  ! Find first error
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-1
        rel_err = abs(ufrc(i,j,k) - ufrc_ref(i,j,k))
        if (abs(ufrc_ref(i,j,k)) > 1.0e-10) rel_err = rel_err / abs(ufrc_ref(i,j,k))
        if (rel_err > 1.0e-5) then
          write(*,*) 'First error at i,j,k = ', i, j, k
          write(*,*) '  ufrc_cpu  = ', ufrc(i,j,k)
          write(*,*) '  ufrc_ref  = ', ufrc_ref(i,j,k)
          write(*,*) '  ufrc_init = ', ufrc_init(i,j,k)
          write(*,*) '  tmp1      = ', tmp1(i,j,k)
          write(*,*) '  tmp2      = ', tmp2(i,j,k)
          write(*,*) '  tmp3      = ', tmp3(i,j,k)
          stop
        end if
      end do
    end do
  end do
  
  write(*,*) 'All values match!'
  
end program

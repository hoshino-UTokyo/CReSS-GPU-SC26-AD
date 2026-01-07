program debug4
  implicit none
  real, parameter :: oned3 = 1.0/3.0, eps = 1.0e-20
  integer :: ni, nj, nk, i, j, k
  real, allocatable :: tkefrc(:,:,:), tkefrc_in(:,:,:), tkefrc_ref(:,:,:)
  real, allocatable :: priv(:,:,:), rst(:,:,:), jcb(:,:,:), rmf(:,:,:), tke(:,:,:)
  real :: ln, ds308, dz05, dz, dx, dy
  real :: max_err, err
  integer :: err_cnt

  ni = 899; nj = 899; nk = 128
  dx = 2292.649414063
  dy = 2499.999267578
  dz = 100.0
  dz05 = 0.5 * dz
  ds308 = 0.125 * dx * dy * dz

  allocate(tkefrc(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(priv(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))

  open(10, file='data/tkefrc_ref.bin', access='stream', form='unformatted')
  read(10) tkefrc_ref; close(10)
  open(10, file='data/tkefrc_in.bin', access='stream', form='unformatted')
  read(10) tkefrc_in; close(10)
  open(10, file='data/tke.bin', access='stream', form='unformatted')
  read(10) tke; close(10)
  open(10, file='data/priv.bin', access='stream', form='unformatted')
  read(10) priv; close(10)
  open(10, file='data/rst.bin', access='stream', form='unformatted')
  read(10) rst; close(10)
  open(10, file='data/jcb.bin', access='stream', form='unformatted')
  read(10) jcb; close(10)
  open(10, file='data/rmf.bin', access='stream', form='unformatted')
  read(10) rmf; close(10)

  ! Copy input
  tkefrc = tkefrc_in

  ! Execute kernel (isoopt=1, mfcopt=1, mpopt=0)
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        ln = (priv(i,j,k) - 1.0) * exp(oned3 * log(ds308 * rmf(i,j,2) * jcb(i,j,k))) + eps
        tkefrc(i,j,k) = tkefrc(i,j,k) - (0.37 * priv(i,j,k) - 0.18) &
             * rst(i,j,k) * tke(i,j,k) * sqrt(tke(i,j,k)) / ln
      end do
    end do
  end do

  ! Validate
  max_err = 0.0
  err_cnt = 0
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        err = abs(tkefrc(i,j,k) - tkefrc_ref(i,j,k))
        if (abs(tkefrc_ref(i,j,k)) > 1.0e-20) err = err / abs(tkefrc_ref(i,j,k))
        if (err > max_err) max_err = err
        if (err > 1.0e-5) err_cnt = err_cnt + 1
      end do
    end do
  end do

  write(*,*) 'Max relative error:', max_err
  write(*,*) 'Error count:', err_cnt
  if (err_cnt == 0) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program

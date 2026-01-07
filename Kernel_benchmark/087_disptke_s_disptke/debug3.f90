program debug3
  implicit none
  integer :: ni, nj, nk
  real, allocatable :: tkefrc_ref(:,:,:), tkefrc_in(:,:,:), tke(:,:,:)
  real, allocatable :: priv(:,:,:), rst(:,:,:), jcb(:,:,:), rmf(:,:,:)
  integer :: i, j, k, cnt
  real :: ln, ds308, dz05, dz, dx, dy
  real, parameter :: oned3 = 1.0/3.0, eps = 1.0e-20
  real :: computed, expected, diff

  ni = 899; nj = 899; nk = 128
  dx = 2292.649414063
  dy = 2499.999267578
  dz = 100.0
  dz05 = 0.5 * dz
  ds308 = 0.125 * dx * dy * dz

  allocate(tkefrc_ref(0:ni+1, 0:nj+1, 1:nk))
  allocate(tkefrc_in(0:ni+1, 0:nj+1, 1:nk))
  allocate(tke(0:ni+1, 0:nj+1, 1:nk))
  allocate(priv(0:ni+1, 0:nj+1, 1:nk))
  allocate(rst(0:ni+1, 0:nj+1, 1:nk))
  allocate(jcb(0:ni+1, 0:nj+1, 1:nk))
  allocate(rmf(0:ni+1, 0:nj+1, 1:4))

  open(10, file='data/tkefrc_ref.bin', access='stream', form='unformatted')
  read(10) tkefrc_ref
  close(10)
  open(10, file='data/tkefrc_in.bin', access='stream', form='unformatted')
  read(10) tkefrc_in
  close(10)
  open(10, file='data/tke.bin', access='stream', form='unformatted')
  read(10) tke
  close(10)
  open(10, file='data/priv.bin', access='stream', form='unformatted')
  read(10) priv
  close(10)
  open(10, file='data/rst.bin', access='stream', form='unformatted')
  read(10) rst
  close(10)
  open(10, file='data/jcb.bin', access='stream', form='unformatted')
  read(10) jcb
  close(10)
  open(10, file='data/rmf.bin', access='stream', form='unformatted')
  read(10) rmf
  close(10)

  ! Find first non-zero tke
  do k = 2, nk-2
    do j = 2, nj-2
      do i = 2, ni-2
        if (tke(i,j,k) > 1.0e-5) then
          write(*,'(A,3I5)') 'First non-zero tke at:', i, j, k
          write(*,'(A,ES15.7)') '  tke      = ', tke(i,j,k)
          write(*,'(A,ES15.7)') '  priv     = ', priv(i,j,k)
          write(*,'(A,ES15.7)') '  rst      = ', rst(i,j,k)
          write(*,'(A,ES15.7)') '  jcb      = ', jcb(i,j,k)
          write(*,'(A,ES15.7)') '  rmf(2)   = ', rmf(i,j,2)
          write(*,'(A,ES15.7)') '  tkefrc_in= ', tkefrc_in(i,j,k)
          write(*,'(A,ES15.7)') '  tkefrc_ref=', tkefrc_ref(i,j,k)

          ! Compute like isoopt=1, mfcopt=1, mpopt=0
          ln = (priv(i,j,k)-1.0)*exp(oned3*log(ds308*rmf(i,j,2)*jcb(i,j,k)))+eps
          computed = tkefrc_in(i,j,k) - (0.37*priv(i,j,k)-0.18) &
               * rst(i,j,k)*tke(i,j,k)*sqrt(tke(i,j,k))/ln
          write(*,'(A,ES15.7)') '  ln       = ', ln
          write(*,'(A,ES15.7)') '  computed = ', computed
          write(*,'(A,ES15.7)') '  diff     = ', computed - tkefrc_ref(i,j,k)
          stop
        end if
      end do
    end do
  end do

end program

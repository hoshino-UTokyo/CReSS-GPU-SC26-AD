function bbplk(wl,tmp,naplk,aplk)
  implicit none
  integer::naplk,j
  real(8)::wl,tmp,aplk(1:naplk),b,xx,bbplk
!--exec      
  xx=1.d0/(wl*tmp)
  b=aplk(1)
  do j=2,naplk
     b=b+aplk(j)*xx**(j-1)
  enddo
  bbplk=1.d0/(exp(b)*wl**3*xx)
end function bbplk

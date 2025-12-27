subroutine gtrh(nln,pl,tl,cgas,rh)
!--- history
! 97.04.24 created
!--- input
! nln       I        number of layers.
! pl      R(nln)     atmospheric pressure at the layer (hPa)
! tl      R(nln)     temperature of layers (K) 
! cgas    R(nln)     gas concentration in the sublayer (ppmv)
!--- output
! rh      r(nln)      relative humidity
!---
  use m_commstrn, only: tstd
  use m_commpi
  real(8),parameter::e0=611.d-2 ,rl=2.5d6,rv=461.d0
  integer,intent(in)::nln
  real(8),intent(in)::cgas(nln),pl(nln),tl(nln)
  real(8),intent(out)::rh(nln)
  integer::l
  real(8)::e(nln),es(nln)
!--exec  
  do l=1,nln
     e(l)  = cgas(l)*1.0d-6*pl(l)
     es(l) = e0*exp( rl/rv*(1.d0/tstd - 1.d0/tl(l)) )
     rh(l) = e(l)/es(l)
  enddo
  return
end subroutine gtrh

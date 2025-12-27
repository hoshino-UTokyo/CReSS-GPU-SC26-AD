subroutine sdp2(nln,p,dp,np)
! derive absorption coefficient from pt table and fitting.
!--history
! 02.02.27  modified knupt in ckdng.f (Zhang Hua)
! 02.11.25  modified from knpt1.      (Miho Sekiguchi)
! 03.09.09  modified from sdp for free source form.(Miho Sekiguchi)
!--
  implicit none
  integer,intent(in)::nln
  real(8),intent(in)::p(nln)
  integer,intent(out)::np(nln)
  real(8),intent(out)::dp(nln)
! works
  integer,parameter::kp=26
  real(8)::prs(1:kp)=(/ &
       1.0130d+3,  0.63096d+3, 0.39811d+3,  0.25119d+3,  0.15849d+3, &
       1.0000d+2,  0.63096d+2, 0.39811d+2,  0.25119d+2,  0.15849d+2, &
       1.0000d+1,  0.63096d+1, 0.39811d+1,  0.25119d+1,  0.15849d+1, &
       1.0000d+0,  0.63096d+0, 0.39811d+0,  0.25119d+0,  0.15849d+0, &
       1.0000d-1,  0.63096d-1, 0.39811d-1,  0.25119d-1,  0.15849d-1,&
       1.0000d-2/)
  integer::l,ip
!--exec
  do l=1,nln
     do ip = 1, kp
        if (p(l) >= prs(ip)) exit
     enddo
     np(l)=ip
     if(np(l) > kp) np(l)=kp
     if(np(l) <= 1) np(l)=2
     dp(l)=log10(p(l)/prs(np(l)-1))/log10(prs(np(l))/prs(np(l)-1))
  enddo
end subroutine sdp2

subroutine tdok2(nln,akt,t,knu)
!----------------------------------------------------------------------c
!                using shi's formula: (t/to)**(a+bt)
!----------------------------------------------------------------------c
  implicit none
  integer,intent(in)::nln
  real(8),intent(inout)::akt(3,nln)
  real(8),intent(in)::t(nln)
  real(8),intent(out)::knu(nln)
  integer::l
  real(8)::a,b,at1,at2,bt1,bt2
!--exec
  at1 = log10(200.d0/260.d0); at2 = log10(320.d0/260.d0)
  
  do l=1,nln
     if(akt(2,l)<=-20.d0.and.akt(1,l)>-20.d0.and.akt(3,l)>-20.d0) &
          akt(2,l)=(akt(1,l)+akt(3,l))*0.5d0
     bt1 = akt(1,l)-akt(2,l)
     bt2 = akt(3,l)-akt(2,l)
     if(bt1==0.d0 .and. bt2==0.d0) then
        knu(l)=10**akt(2,l); cycle
     endif
     b   = (bt2/at2-bt1/at1)/120.d0
     a   = bt2/at2 - b*320.d0
     
     knu(l)=10**(akt(2,l))*(t(l)/260.d0)**(a+b*t(l))
  enddo
  return
end subroutine tdok2

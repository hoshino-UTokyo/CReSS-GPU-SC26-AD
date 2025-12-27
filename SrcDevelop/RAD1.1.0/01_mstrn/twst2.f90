subroutine twst2(am0,nln,thk,omg,gdcfrc,g,galb,fsol,nplk1,cplk,bgnd, &
     nb,flxd,flxu,err)
! solve the radiative transfer in atmosphere system for
! dynamical models.  flux at interfaces. two stream method
! by Teruyuki Nakajima
!--- history
! 92  2.20 created from rtrn1.
!     6.25 bug for if
! 93. 5.29 single scattering approx. for thin layer.
! 94. 12.2 partial cloudiness approx. by semi-random method.(Y.Tsushima)
! 95. 5.30 introduce sumwm to improve ground reflection.(Y.Tsushima)
! 02.12.05 integrated twst12 and twst32 (M.Sekiguchi)
!--- input
! nda        I        # of streams
! am0        R        consine of solar zenith angles .gt. 0.
! nln        I        number of atmospheric sublayers.
! thk   R(knln,2)     optical thickness of sublayers from top to bottom
! omg   R(knln,2)     single scattering albedo
! gdcfrc  R(knln)     cloud cover rate of sublayer.
! g    R(3,knln,2)    legendre moments of phase function
!                      g(i), i=1,nlgn1;  nlgn1=2*nda+1=3
! galb     R(2)       ground albedo (lambertian)
! fsol       R        solar irradiance at the system top.
! nplk1      I        number of order to approximate plank + 1.
!                      if 0 then no thermal.
! cplk R(kplk1,knln,2) thermal emission
!                      = sum(k1=1,nplk1) cplk(k1,l) * tau**(k1-1)
!                      tau is the optical depth measured from the top
!                      of the sublayer (same unit as fsol).
! bgnd     R(2)        thermal emission from the ground
!                      = (1-galb)*plank function
!                      (same unit as fsol). if 0 then no emission
! nb         I        flag of solar or thermal region
!--- output
! flxd     R(knln1)   downward flux at interfaces.
! flxu     R(knln1)   upward   flux at interfaces.
! err      c*64       error indicater. if " " then normal end.
!--
  use m_commstrn, only: kplnk
  implicit none
  integer,intent(in)::nln,nplk1,nb
  real(8),intent(in)::am0,thk(nln,2),omg(nln,2),gdcfrc(nln),g(3,nln,2), &
       galb(2),fsol,cplk(kplnk,nln,2),bgnd(2)
  real(8),intent(out)::flxd(nln+1),flxu(nln+1)
  character::err*64
  real(8),parameter::pi=3.141592654d0,rad=pi/180.0d0

! two-stream
  real(8)::amua,wa,wmp,wmm,eps,sumwm,pla(2),pl0(2),pa(2),p0(2)

! delta-m method
  real(8)::expdm(nln+1),exptm,thkt(nln,2)
  real(8)::re(nln+1),te(nln+1),ser(nln+1),set(nln+1),rl(nln+1),sl(nln+1)
  real(8)::dptt(nln+1),g0,ff(nln),gt(nln),omgt(nln),cplkt(kplnk,nln)

! sublayer
  integer::ipk,iw0
  real(8)::t1,t2,tau,ww,gg,cpl(kplnk),pt,pr,pt0,pr0,sum,r,t,sum1,er,et, &
       x,y,q,qi,zeig,zeig2,xi,dp(kplnk),dm(kplnk),cplk1,sp,sm,g1,gam,trns0, &
       sum2,vp,vm,e0,exp2,exp4,cf,sf,slf,sl1f,sum3,sum4,ap,am,bp,bm,c11,c22, &
       dp1,dm1,dp0,dm0,vp0,vm1,rn(nln+1),rr,tt,r2d,su,bb,sd,rd, &
       rup(nln+2),rdn(nln+2),s2u,r1d,s1d,rt,st,fu1,fd1
  
! work
  integer::l,ic,k1,i,j,nln1,nlt,nlt1,l1
!--exec
  err=''
  if(nb==1) then
! amua=1/1.66
     amua=0.602409638d0; wa=1.0d0
! wmp=sqrt(wa*amua); wmm=sqrt(wa/amua)
     wmp =0.776150525d0; wmm=1.288409873d0
     eps=1.d-10
  elseif(nb==2) then
     amua=0.577350269d0; wa=1.0d0
     wmp =0.759835686d0; wmm=1.316074013d0
     eps=1.d-5
  else
     err='not decided the region'; return
  endif

  sumwm=amua*wa
! set legendre polynomials.
  pla(1)=wmm ; pla(2)=wmm*amua
  pl0(1)=1.d0; pl0(2)=am0
  pa(1)=0.5d0*pla(1)**2; pa(2)=1.5d0*pla(2)**2
  p0(1)=0.5d0*pla(1)*pl0(1); p0(2)=1.5d0*pla(2)*pl0(2)

! delta-m moments
  expdm(1)=1.d0
  do l=1,nln
     thkt(l,1:2)=(1.d0-g(3,l,1:2)/g(1,l,1:2)*omg(l,1:2))*thk(l,1:2)
     exptm=gdcfrc(l)*exp(-thkt(l,1)/am0)+(1.d0-gdcfrc(l))*exp(-thkt(l,2)/am0)
     expdm(l+1)=expdm(l)*exptm
  enddo

  re(1:nln+1)=0.d0; te(1:nln+1)=0.d0
  set(1:nln+1)=0.d0; ser(1:nln+1)=0.d0
  rl(1:nln+1)=0.d0; sl(1:nln+1)=0.d0
  do ic=1,2
     dptt(1)=0.d0
     if(ic==1) rn(1:nln)=gdcfrc(1:nln)
     if(ic==2) rn(1:nln)=1.d0-gdcfrc(1:nln)
     do l=1,nln
        g0=g(1,l,ic)
        ff(l)=g(3,l,ic)/g0
        gt(l)=(g(2,l,ic)/g0-ff(l))/(1-ff(l))
        omgt(l)=(1-ff(l))*omg(l,ic)/(1-ff(l)*omg(l,ic))
        dptt(l+1)=dptt(l)+thkt(l,ic)
     enddo

! scaling coefficients for thermal emission
     if(nplk1>0) then
        do k1=1,nplk1
           cplkt(k1,1:nln)=cplk(k1,1:nln,ic)/(1.d0-omg(1:nln,ic)*ff(1:nln)) &
                **(k1-1)
        enddo
     endif

! solve the eigenvalue problem of atmospheric sublayers.
! loop for sublayers
     do l=1,nln
        t1=dptt(l); t2=dptt(l+1); tau=thkt(l,ic)
        ww=omgt(l); gg=gt(l); ipk=0
        if(1.d0-ww .le. eps) then
           iw0=1
        else
           iw0=0
           if(nplk1>0) ipk=1
        endif
        if(nplk1>0) then
           cpl(1:nplk1)=2*pi*(1.d0-ww)*cplkt(1:nplk1,l)
        endif

!! scattering media
!! pt, pr
        pt =pa(1)+gg*pa(2); pr =pa(1)-gg*pa(2)
        pt0=p0(1)+gg*p0(2); pr0=p0(1)-gg*p0(2)
        if(tau<=1.0d-4) then
           sum=0.d0
           if(nplk1>0) then
              do i=1,nplk1
                 sum=sum+cpl(i)*tau**(i-1)
              enddo
              sum=wmm*(sum+cpl(1))*0.5d0*tau
           endif
           r=ww*pr*tau
           t=1-tau/amua+ww*pt*tau
           sum1=ww*tau*expdm(l)*exp(-thkt(l,ic)*0.5d0/am0)*fsol
           er=pr0*sum1+sum; et=pt0*sum1+sum
        else

!! eigenvalue problem
! x, y matrices
           x=1.d0/amua-ww*(pt-pr); y=1.d0/amua-ww*(pt+pr)
           q=sqrt(x); qi=1/q
           zeig2=x*y; zeig=sqrt(zeig2)
           xi=1.d0/x

! thermal source
           if(ipk==0.d0) then
              if(nplk1>0) then
                 dp(1:nplk1)=0.d0; dm(1:nplk1)=0.d0
              endif
           else
              do j=nplk1,1,-1
                 if(j+2>nplk1) then
                    cplk1=0
                 else
                    cplk1=dp(j+2)
                 endif
                 dp(j)=((j+1)*j*cplk1+q*wmm*cpl(j))/zeig2
                 dm(j)=q*dp(j)
              enddo
              do j=1,nplk1
                 sum=0.d0
                 if(j+1<=nplk1) sum=qi*dp(j+1)
                 dp(j)=dm(j)-j*sum; dm(j)=dm(j)+j*sum
              enddo
           endif
           sp=ww*(pt0+pr0); sm=ww*(pt0-pr0)
           g1=-x*sp-sm/am0; gam=qi*g1/(1.0/am0**2-zeig2)
           trns0=expdm(l)*fsol
           sum1=q*gam; sum2=qi*gam/am0+xi*sm
           vp=(sum1+sum2)/2.0*trns0; vm=(sum1-sum2)/2.0*trns0
           e0=exp(-tau/am0)

! base function  cf(tau) and sf(tau).
           exp2=zeig*tau; exp4=exp(-exp2)
           cf  =(1+exp4)*0.5d0; sf  =(1-exp4)*0.5d0
           slf =zeig*sf
           if(abs(exp2)<=eps) then
              sl1f =tau*0.5d0
           else
              sl1f =sf/zeig
           endif

! a+-, b+-
           sum1=q*cf; sum2=qi*slf; sum3=q*sl1f; sum4=qi*cf
           ap=sum1-sum2; am=sum1+sum2
           bp=sum3-sum4; bm=sum3+sum4

! c11 and c22 -> their inversion.
           c11=1.0d0/(2*am); c22=1.0d0/(2*bm)

! r, t matrices
           sum1=ap*c11; sum2=bp*c22
           r=sum1+sum2; t=sum1-sum2

! er, et matrices
           dp1=0.d0; dm1=0.d0
           if(ipk==1) then
              do j=1,nplk1
                 dp1=dp1+dp(j)*tau**(j-1)
                 dm1=dm1+dm(j)*tau**(j-1)
              enddo
              dp0=dp(1); dm0=dm(1)
           else
              dp0=0.d0; dm0=0.d0
           endif
           vp0=vp   +dp0; vm1=vm*e0+dm1
           sum1=r*vp0+t*vm1; sum2=t*vp0+r*vm1
           er=vm   +dm0-sum1; et=vp*e0+dp1-sum2
        endif
!! set to buffer
        re (l)=re(l)+rn(l)*r; te (l)=te(l)+rn(l)*t
        ser(l)=ser(l)+rn(l)*er; set(l)=set(l)+rn(l)*et
      enddo

! lambert surface
      nln1=nln+1
      if(galb(ic)<=0.d0 .and. bgnd(ic)<=0.d0) then
         nlt=nln
      else
         nlt=nln1
         t=0.d0; r=galb(ic)*wmp**2/sumwm
         x=galb(ic)*am0*expdm(nlt)*fsol/sumwm+2*pi*bgnd(ic)
         et=0; er=wmp*x
         rn(nln1)=rn(nln)
         re (nlt)=re(nlt)+rn(nlt)*r; te (nlt)=te(nlt)+rn(nlt)*t
         ser(nlt)=ser(nlt)+rn(nlt)*er;set(nlt)=set(nlt)+rn(nlt)*et
      endif
   enddo
   nlt1=nlt+1
! upward adding
   rl(nlt)=re(nlt)
   sl(nlt)=ser(nlt)
   if(nlt >= 2) then
      do l=nlt-1,1,-1
         l1=l+1
         rr=re(l); tt=te(l)
         r2d=rl(l1); su=r2d*set(l)+sl(l1)
         bb=1/(1-r2d*rr)
         sd=bb*su; su=tt*sd+ser(l)
         rd=rr+tt*bb*r2d*tt
         rl(l)=rd; sl(l)=su
      enddo
   endif

! field at the top of the system
   rup(1)=sl(1);  rdn(1)=0.d0

! downward adding
   rt=re(1); st=set(1)
   if(nlt>=2) then
      do l=2,nlt
! internal field
         s2u=sl(l); r2d=rl(l)
         sd=st+rt*s2u
         r1d=1.d0/(1.d0-rt*r2d); s1d=r1d*sd
         rdn(l)=s1d; rup(l)=r2d*s1d+s2u

         rr=re(l); tt=te(l)
         bb=1/(1-rt*rr); sd=bb*(rt*ser(l)+st)
         st=tt*sd+set(l); rt=rr+tt*bb*rt*tt
      enddo
   endif
! field at the bottom of the system
   rdn(nlt1)=st; rup(nlt1)=0.d0

! flux
   do l=1,nln1
      fu1=wmp*rup(l); fd1=wmp*rdn(l)
!      if(nb==1) then
      fu1=fu1*(nb-1)+fu1/amua*0.5d0*(2-nb)
      fd1=fd1*(nb-1)+fd1/amua*0.5d0*(2-nb)
!      endif
      flxu(l)=fu1; flxd(l)=fd1+am0*expdm(l)*fsol
   enddo
   return
end subroutine twst2

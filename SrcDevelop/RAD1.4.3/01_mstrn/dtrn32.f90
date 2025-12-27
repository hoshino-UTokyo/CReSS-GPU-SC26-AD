subroutine dtrn3(iug,iup,iuv,ams,nln,pb,pl,tb,tl,gtmp,cpcl,gdcfrc,cgas,ccfc, &
     prg,fd,fu,err)
! flux calculation
!--- history
! 93. 3.18 (quarante) hitran 1992  k-distributi
!     4.14 (teruyuki) registered
! 94. 6.11 change fitting yy to sqrt(sqrt(b))
!          add layer temperature for planck 2-dim. expansion term
! 94. 6.15 use layer and layer boundary p,t and avoid z
!     6.23 revise cplk
!    11.25 change h2o continuum treatment
!    11.29 change planck function treatment
! 94.12. 2 (yoko)     partial cloudiness approximation by semi-random method.
!    12. 5 polish up h2o continuum treatment
!    12.11 kd treatment; pt-interpolation -> function fitting
!    12.14 add cfcs. number of optical properties ifbnd changes 7 -> 8
! 95. 3. 8 revise o2 continuum optical thickness calculation
! 97.04.28  modified from mstr7e kpcl=8 -> kpcl=11
!				 iup -> iupp,iupg,iupv
!				 rmode interpolation(linear) by Y.Tsushima
! 02.11.21 modidfied from dtrn22 by M.Sekiguchi
! 04. 4. 7 Add CFCs for mstrnX by M.Sekiguchi.
! 10.12.30 Bug fixed by M.Sekiguchi (from Hashino-kun). <== add by CTI(2014/06/10)
!--- input
! iug      I           read unit number of gas line absorption para file
! iup      I           read unit number of other extinction para file
! iuv      I           read unit number of source file of iup
! ams       R          cos(solar zenith angle). if ams.le.0 then no solar.
! nln       I        number of layers.
! pb      R(nln1)    atmospheric pressure at the interface (mb) (94. 6.15)
! pl      R(nln)     atmospheric pressure at the layer (mb) (94. 6.15)
! tb      R(nln1)    temperature of layer boundaries(k)
! tl      R(nln)     temperature of layers(k) (94. 6.11)
! gtmp     R         ground temperature (k)
! gdcfrc  R(knln)     cloud cover rate of sublayer.
! cpcl R(nln,npcl,3) parameter packet for particulates in sublayers. (l,i,j)
!                    l: layer number (from top to bottom)
!                       1: water cloud               2: ice cloud
!                       3: dust-like                 4: soot
!                       5: volcanic-ash              6: h2so4
!                       7: rural                     8: sea salt
!                       9: urban                    10: tropo.
!                      11: yellow dust
!                    j=1: concentration in ppmv
!                         = aerosol volume(cm3/cm3)*t*p0/t0/p*1.0e6
!                         with t0=273, p0=1atm.
!                         type 7,8,9,10 = dry concentration(parapc.voldp)
!                                       = total concentration(parapc.volp)
!                    j=2  mode radius(cm)
!                    j=3 undecided for future use.
!
! cgas  R(nln,ngas)   gas concentration in the sublayer (ppmv)
!  1: h2o, 2: co2, 3: o3, 4: n2o, 5: co, 6: ch4, 7: o2
!
! ccfc   R(kcfc)     CFCs concentrations
!   1: CFC-11  ,  2: CFC-12   ,  3: CFC-13   ,  4: CFC-14  ,    5: CFC-113  , 
!   6: CFC-114 ,  7: CFC-115  ,  8: HCFC-21  ,  9: HCFC-22 ,   10: HCFC-123 ,
!  11: HCFC-124, 12: HCFC-141b, 13: HCFC-142b, 14: HCFC-225ca, 15: HCFC-225cb, 
!  16: HFC-32  , 17: HFC-125  , 18: HFC-134  , 19: HFC-134a  , 20: HFC-143a  ,
!  21: HFC-152a, 22: SF6      , 23: ClONO2   , 24: CCl4    ,   25: N2O5     , 
!  26: C2F6    , 27: HNO4     , 28: SF5CF3  
!
! prg     r(5)       parameter packet for surface reflection.
!                      (1)             (2)           (3)
!                      1 ocean         -ams-         -tr- (automatically set
!                      2 wet land      water content      inside this loutine)
!                      3 dry land       -              -
!                      4 low plants     -              -
!                      5 forest         -              -
!                      6 snow           -              -
!                      7 ice            -              -
!--- output
! fd      r(nln1,2)  downward flux at the sublayer interface (wm-2)
!                     (from top to bottom)
!***                    j=2 shorter than 4 micron, j=1 longer than 4 micron
! fu      r(nln1,2)   upward   flux at the sublayer interface.
! err     c*64       error code.  if ' ' then no error.
!--
  use m_commstrn
  use m_commpi
  use m_defmpi
  implicit none
  integer,intent(in)::iug,iup,iuv,nln
  real(8),intent(in)::pb(nln+1),pl(nln),tb(nln+1),tl(nln),gtmp, &
       cpcl(nln,kpcl,kpclc),gdcfrc(nln),cgas(nln,kmol),ccfc(kcfc),prg(kprg)
  real(8),intent(inout)::ams
  real(8),intent(out)::fd(nln+1,2),fu(nln+1,2)
  character::err*64

! gtrh
  real(8)::rh(nln)

! gttbl
  integer,save::nbnd,nda,np,nt,nabs(kbnd),iabs(kabs,kbnd),nflg, &
       iflgb(kflg,kbnd),nch(kbnd),nsfc,nplk,naplk,npcl,nvar(kpcl,2)
  real(8),save::wv(kbnd1),prs(kp),temp(kt),wgt(kch,kbnd), &
       akd(kch,kp,kt,kabs,kbnd),skd(kch,kp,kt,kbnd),acfc(kcfc,kbnd), &
       aplnk(kbnd,kaplk),fsol(kbnd),sr(kbnd,ksfc),ry(kbnd), &
       qmol(kmax2),q(kbnd,kpcl,kvar,kmax2),r0(kpcl,kvar)
! 
  real(8)::aplnk_w(kbnd*kaplk)

! gas concentration for absorption
  real(8)::amtpb(kmol,nln),wbrd(nln)

! p-t table
  integer::mp(nln)
  real(8)::akt(kt,nln),knu(nln),sku(nln),tkd(nln),ddp(nln)

! particle scattering
  integer::mvar
  real(8)::taur,cz,qq(nln,kmax2,2),r(nln),qqq(1),taup(nln,2),taut(2), &
       scap(nln,2),g(3,nln,2)
      
! surface
  integer::is
  real(8)::tr(2),galb(2),wet

! planck function
  integer::nplk1
  real(8)::bgnd(2),wl,bp(nln+1),bl(nln),bgnd0,cplk(kplnk,nln,2),dbp1,dbp2
  
! sublayers
  real(8)::tau(nln,2),omg(nln,2)

! twst
  real(8)::flxd(nln+1),flxu(nln+1)

! work
  integer,save::init=1
  integer::nda2,nln1,iw,l,k,ipcl,ia,ib,ic,ich
  real(8)::dp(nln),rhoair(nln),dz(nln)
  real(8)::c0=rgas/(grav*airm)
  real(8)::ssrfc,bbplk
  integer ierr

  integer i,j,cnt_min,cnt_max
  real(8)::b,ww
!--- exec
! initialization
  if (init>0) then
     init=0

! read optical parameters

     if(mype.eq.root) then

       rewind(iug)
       rewind(iup)
       rewind(iuv)

       call gttbl(iug,iup,iuv,nbnd,nda,np,nt,nabs,iabs,nflg,iflgb,nch, &
         wv,prs,temp,wgt,akd,skd,acfc,nsfc,nplk,naplk,aplnk,fsol,sr,ry,&
         qmol,q,npcl,nvar,r0,err)

     end if

     if(err/='') return

     call mpi_bcast(nbnd,1,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(nda,1,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(np,1,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(nt,1,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(nabs,kbnd,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(iabs,kabs*kbnd,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(nflg,1,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(iflgb,kflg*kbnd,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(nch,kbnd,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(wv,kbnd1,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(prs,kp,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(temp,kt,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(wgt,kch*kbnd,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(akd,kch*kp*kt*kabs*kbnd,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(skd,kch*kp*kt*kbnd,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(acfc,kcfc*kbnd,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(nsfc,1,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(nplk,1,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(naplk,1,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(aplnk,kbnd*kaplk,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(fsol,kbnd,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(sr,kbnd*ksfc,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(ry,kbnd,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(qmol,kmax2,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(q,kbnd*kpcl*kvar*kmax2,mpi_real8,root,mpi_comm_cress,ierr)
     call mpi_bcast(npcl,1,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(nvar,2*kpcl,mpi_integer,root,mpi_comm_cress,ierr)
     call mpi_bcast(r0,kpcl*kvar,mpi_real8,root,mpi_comm_cress,ierr)

  endif

! calculate relative humidity in each layer
  call gtrh(nln,pl,tl,cgas(1:nln,1),rh)

! write(*,'(A,5I12)') 'DUMP-CHECK(nbnd,nda2,npcl,nln,naplk) : ',nbnd,nda2,npcl,nln,naplk

  do l=1,nln
!! get p, t, dp, dz from grid data (94. 6.15)
     dp(l)=abs(pb(l)-pb(l+1))
!!  c0=rgas/grav/airm
     dz(l)=c0*tstd*dp(l)/pstd
     rhoair(l)=pl(l)*tstd/(pstd*tl(l))
     amtpb(1:kmol,l)=cgas(l,1:kmol)*dz(l)*1.d-1
     wbrd(l)=dz(l)*1.d5
  enddo

! for interpolation of layers (p,T)
  call sdp2(nln,pl,ddp,mp)

  ams=max(ams, 1.0d-6); nda2=2*nda; nln1 = nln+1
  fd(1:nln1,1:2) = 0.d0; fu(1:nln1,1:2)=0.d0

! loop for bands
  do iw = 1,nbnd
     taut(1:2)=0.d0

!! loop for sublayers
     do l=1,nln

        taur=ry(iw)*dp(l)/pstd
        cz=0.1d0*dz(l)
        do k=1,nda2+2
           qq(l,k,1:2)=taur*qmol(k)
           do ipcl=1,npcl
              if(cpcl(l,ipcl,1) > 0.d0) then

!! bug fixed
                 mvar=nvar(ipcl,2)
                    IF(cpcl(l,ipcl,2)<r0(ipcl,1)) THEN
                       ww=(r0(ipcl,1)-cpcl(l,ipcl,2))/(r0(ipcl,2)-r0(ipcl,1))
                       qqq(1)=q(iw,ipcl,1,k)-(q(iw,ipcl,2,k)-q(iw,ipcl,1,k))*ww
                    ELSEIF(cpcl(l,ipcl,2)<=r0(ipcl,mvar)) THEN
!cdir novector
                       DO J=1,mvar-1
                          IF(r0(ipcl,J)<=cpcl(l,ipcl,2) .and. r0(ipcl,J+1)>cpcl(l,ipcl,2)) then
                             ww=(cpcl(l,ipcl,2)-r0(ipcl,J))/(r0(ipcl,J+1)-r0(ipcl,J))
                             qqq(1)=q(iw,ipcl,J,k)*(1-ww)+q(iw,ipcl,J+1,k)*ww
                             EXIT
                          ENDIF
                       ENDDO
                    ELSE
                       ww=(cpcl(l,ipcl,2)-r0(ipcl,mvar))/(r0(ipcl,mvar)-r0(ipcl,mvar-1))
                       qqq(1)=q(iw,ipcl,mvar,k)+(q(iw,ipcl,mvar,k)-q(iw,ipcl,mvar-1,k))*ww
                    ENDIF

              endif

              qq(l,k,1)=qq(l,k,1)+qqq(1)*cpcl(l,ipcl,1)*cz
              if(ipcl>=3) qq(l,k,2)=qq(l,k,2)+qqq(1)*cpcl(l,ipcl,1)*cz

           enddo
        enddo
     enddo
     do l=1,nln
        do ic=1,2
           taup(l,ic)=qq(l,1,ic)
           taut(ic)=taut(ic)+taup(l,ic)
           scap(l,ic)=qq(l,1,ic)-qq(l,2,ic)
           g(1,l,ic)=1.d0
           if(scap(l,ic)==0.d0) then
              g(2:nda+1,l,ic)=0.d0
           else
              g(2:nda2+1,l,ic)=qq(l,3:nda2+2,ic)/scap(l,ic)
           endif
        enddo
     enddo
! loop-sublayer end

! ground surface
!! atmospheric transmissivity
     is=int(prg(1))
     if(is==1) then
!! ocean
        if(ams>0.d0) then
           do ic=1,2
              if(taut(ic)>0.d0) then
!! 101230 Bug fixed. (Thanks Hashino-kun)
!ORG             tr(1)=ams/(4.d0*taut(ic));tr(1)=min(tr(1),1.d0)
!ORG             galb(ic)=ssrfc(ams,tr(1))
                 tr(ic)=ams/(4.d0*taut(ic)); tr(ic)=min(tr(ic),1.d0)
                 galb(ic)=ssrfc(ams,tr(ic))
              else
                 galb(ic)=0.1d0
              endif
           enddo
        else
           galb(1:2)=0.1d0
        endif
     elseif(is==2) then
!! land surface
        wet=prg(2)
        galb(1:2)=wet*sr(iw,2)+(1-wet)*sr(iw,3)
     elseif(is<=7) then
!! others
        galb(1:2)=sr(iw,is)
     else
        err = 'no such land surface'
        return
     endif


! planck functions at sublayer interfaces and at layers
! not t table interpolation but function fitting
!

!    aplnk_w(1:naplk) = aplnk(iw,1:naplk)

     if(iflgb(3,iw)==0) then
        nplk1=0; bgnd(1:2)=0.d0; galb(1:2)=0.1d0
     else
        nplk1=3; wl=1.d4/sqrt(wv(iw)*wv(iw+1))
!!! !cdir novector
        do l=1,nln
         ! bp(l)=bbplk(wl,tb(l),naplk,aplnk(iw,1:naplk))
         ! bl(l)=bbplk(wl,tl(l),naplk,aplnk(iw,1:naplk))
         ! bp(l)=bbplk(wl,tb(l),naplk,aplnk_w)
         ! bl(l)=bbplk(wl,tl(l),naplk,aplnk_w)
           b=aplnk(iw,1)
           do j=2,naplk
              b=b+aplnk(iw,j)/(wl*tb(l))**(j-1)
           enddo
           bp(l)=1.d0/(exp(b)*wl**3/(wl*tb(l)))
           b=aplnk(iw,1)
           do j=2,naplk
              b=b+aplnk(iw,j)/(wl*tl(l))**(j-1)
           enddo
           bl(l)=1.d0/(exp(b)*wl**3/(wl*tl(l)))

        enddo
      ! bp(nln1)=bbplk(wl,tb(nln1),naplk,aplnk(iw,1:naplk))
      ! bgnd0=bbplk(wl,gtmp,naplk,aplnk(iw,1:naplk))
      ! bp(nln1)=bbplk(wl,tb(nln1),naplk,aplnk_w)
      ! bgnd0=bbplk(wl,gtmp,naplk,aplnk_w)
        b=aplnk(iw,1)
        do j=2,naplk
           b=b+aplnk(iw,j)/(wl*tb(nln1))**(j-1)
        enddo
        bp(nln1)=1.d0/(exp(b)*wl**3/(wl*tb(nln1)))
        b=aplnk(iw,1)
        do j=2,naplk
           b=b+aplnk(iw,j)/(wl*gtmp)**(j-1)
        enddo
        bgnd0=1.d0/(exp(b)*wl**3/(wl*gtmp))

        galb(1:2)=0.d0
        bgnd(1:2)=(1-galb(1:2))*bgnd0
     endif

! loop for channels
     
     do ich = 1, nch(iw)
        tkd(1:nln)=0.d0
        if(iflgb(1,iw)>0) then
           if(iflgb(7,iw)>0) then
              do l=1,nln
                 tkd(l)=sum(10**acfc(1:kcfc,iw)*ccfc(1:kcfc))*wbrd(l)
              enddo
           endif
           do ia=1,nabs(iw)
              do l=1,nln
                 akt(1:kt,l)=akd(ich,mp(l)-1,1:kt,ia,iw) &
                      +(akd(ich,mp(l),1:kt,ia,iw)-akd(ich,mp(l)-1,1:kt,ia,iw))&
                      *ddp(l)
              enddo
              call tdok2(nln,akt(1:kt,1:nln),tl,knu(1:nln))
              tkd(1:nln)=tkd(1:nln)+knu(1:nln)*amtpb(iabs(ia,iw),1:nln)
              
           enddo

           if(iflgb(5,iw)>0) then
              do l=1,nln
                 akt(1:kt,l)=skd(ich,mp(l)-1,1:kt,iw)+(skd(ich,mp(l),1:kt,iw)&
                      -skd(ich,mp(l)-1,1:kt,iw))*ddp(l)                     
              enddo
              call tdok2(nln,akt(1:3,1:nln),tl,sku(1:nln))
              tkd(1:nln)=tkd(1:nln)+sku(1:nln)*amtpb(1,1:nln)**2/ &
                   (amtpb(1,1:nln)+wbrd(1:nln))
           endif
        endif

! loop for sublayers
        do l=1,nln
! total tau and omega
           tau(l,1:2)=tkd(l)+taup(l,1:2)
           do ic=1,2
              if(tau(l,ic)<=0.d0) then
                 tau(l,ic)=0.d0; omg(l,ic)=1.d0
              else
                 omg(l,ic)=scap(l,ic)/tau(l,ic)
              endif
           enddo
!!         omg(l,1:2)=0.d0 <-- This is bug !!

! get plank expansion coefficients
           cplk(1:3,l,1:2)=0.d0
           if(iflgb(3,iw)>0)then
              dbp1=4*bl(l)-bp(l+1)-3*bp(l)
              dbp2=bp(l+1)+bp(l)-2*bl(l)
              do ic=1,2
                 cplk(1,l,ic)=bp(l)
                 if(tau(l,ic)>0.d0)then
                    cplk(2,l,ic)=dbp1/tau(l,ic)
                    cplk(3,l,ic)=dbp2/(tau(l,ic)**2)*2.d0
                 endif
              enddo
           endif
        enddo

! two-stream transfer
        ib=1
        if(iflgb(4,iw)>0) ib=2
        call twst2(ams,nln,tau,omg,gdcfrc,g,galb,fsol(iw),nplk1,cplk, &
             bgnd,ib,flxd,flxu,err)
        if (err/='') return
        do l=1,nln1
           fd(l,ib)=fd(l,ib)+wgt(ich,iw)*flxd(l)
           fu(l,ib)=fu(l,ib)+wgt(ich,iw)*flxu(l)
        enddo
     enddo
  enddo
  err = ''
  return
end subroutine dtrn3
!

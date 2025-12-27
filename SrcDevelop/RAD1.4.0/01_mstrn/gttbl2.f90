subroutine gttbl(iug,iup,iuv,nbnd,nda,np,nt,nabs,iabs,nflg,iflgb,nch,wv, &
     prs,temp,wgt,akd,skd,acfc,nsfc,nplk,naplk,aplnk,fsol,sr,ry,qmol,q,npcl, &
     nvar,r0,err)
! read optical parameters
!--- history
! 93. 3.28 created by quarante
!     4.26 registered by T. Nakajima
! 94.11.25 change h2o continuum treatment
!    11.29 change planck function treatment
!    12.14 add cfcs. number of optical properties ifbnd changes 7 -> 8
! 97.04.28  modified from mstr7e kpcl=8 -> kpcl=11
!				 iup -> iup,iupg,iuv
!				 rmode interpolation(linear) by Y.Tsushima
! 02.11.18 modified from mstrn8 by M. Sekiguchi 
! 04. 4. 7 Add CFCs for mstrnX by M. Sekiguchi
!--- input
! iug       I        read unit number of parameters of gas line absorption
! iup       I        read unit number of parameters of other extinction
! iuv       I        read unit number of source file of iup
!--- output
! nbnd      I        # of bands
! nda       I        # of streams
! np        I        # of pressures for P-T table
! nt        I        # of temperatures for P-T table
! nabs   I(kbnd)     # of absorbers in each band
! iabs  I(kabs,kbnd) # of absorbers spicies in each band
! nflg  I(kbnd)      # of flags
! iflgb I(kflg,kbnd) optical flag
! 1:line abs, 2:scat, 3:plnk, 4:solar, 5:h2o self, 6:CFC
! wv     R(kbnd1)    boundary wavenumbers [cm-1]
! prs    R(kp)       pressures for P-T table
! temp   R(kt)       temperatures for P-T table
! nch    I(kbnd)     # of subintervals in each band
! wgt   R(kch,kbnd)   weights of subintervals in each band
! akd  R(kch,kp,kt,kabs,kbnd)  absoption coefficients
! skd  R(kch,kp,kt,kbnd)  absoption coefficients for H2O self broadening
! acfc  R(kcfc,kbnd) ansorption coefficients for CFC
! npcl      I        # of particle spicies
! nplk      I        order of planck function fitting
! naplk     I        # of planck function fitting
! aplnk R(kbnd,kaplk) coefficients of planck function fitting
! fsol   R(kbnd)     solar constants in each bands [W/m**2]
! sr    R(kbnd,ksfc) surface condition in each bands
! ry     R(kbnd)     rayleigh scattering in each bands
! qmol   R(kmax2)    
! q   R(kbnd,kpcl,kvar,kmax2) coefficients of particles
! nvar  I(kpcl,2)    1: effective radius change or humidity growth
!                    2: # of values
! r0   R(kpcl,kvar)  values in each particles
!--
  use m_commstrn, only: kbnd,kbnd1,kabs,kp,kt,kch,kcfc,kaplk,ksfc,kmax2, &
       kpcl,kvar,kflg
  implicit none

! file number
  integer,intent(in)::iug,iup,iuv

! gas parameters
  integer,intent(out)::nbnd,nda,np,nt,nabs(kbnd),iabs(kabs,kbnd),nflg, &
       iflgb(kflg,kbnd),nch(kbnd)
  real(8),intent(out)::wv(kbnd1),prs(kp),temp(kt),wgt(kch,kbnd), &
       akd(kch,kp,kt,kabs,kbnd),skd(kch,kp,kt,kbnd),acfc(kcfc,kbnd)
  integer::iw,ip,it,ncfc

! particle parameters 
  integer::mbnd,mda,im,mpcl,ia
  integer,intent(out)::nsfc,nplk,naplk
  real(8),intent(out)::aplnk(kbnd,kaplk),fsol(kbnd),sr(kbnd,ksfc),ry(kbnd), &
       qmol(kmax2),q(kbnd,kpcl,kvar,kmax2)

! vardata      
  integer::ipcl
  integer,intent(out)::npcl,nvar(kpcl,2)
  real(8),intent(out)::r0(kpcl,kvar)

! work
  character::err*64

!--- exec
! get vardata
  read(iuv,*) npcl
  do ipcl=1,npcl
     read(iuv,*) 
     read(iuv,*) nvar(ipcl,1:2)
     read(iuv,*) r0(ipcl,1:nvar(ipcl,2))
  enddo

! read gas parameters
  read(iug,*) nbnd,nda,np,nt,nflg,ncfc

! band boundaries
  read(iug,*)
  read(iug,*) wv(1:nbnd+1)

! read particle parameters
  read(iup,*) mbnd,nsfc,mpcl,mda,nplk,naplk
  read(iup,*)
  read(iup,*) wv(1:mbnd+1)
  if(mbnd/=nbnd.or.mda/=nda) then
     err='Different between gas and particle tables!'
     return
  endif
  if(mpcl/=npcl) then
     err='Defferent between particle and vardata table!'
     return
  endif

! log(pressure) grids
  read(iug,*)
  read(iug,*) prs(1:np)
! temperature grids
  read(iug,*)
  read(iug,*) temp(1:nt)

! quantities for each band
  do iw=1,nbnd

!! optical properties flag
     read(iug,*) 
     read(iug,*) iflgb(1:nflg,iw)

!! number of subintervals
     read(iug,*) 
     read(iug,*) nch(iw)

!! weights for channels
     read(iug,*)
     read(iug,*) wgt(1:nch(iw),iw)

!! molecules
     read(iug,*)
     read(iug,*) nabs(iw)
!!! major gas absorption
     if(nabs(iw)>0) then
        do ia=1,nabs(iw)
           read(iug,*) iabs(ia,iw)
           do it=1,nt; do ip=1,np
              read(iug,*) akd(1:nch(iw),ip,it,ia,iw)
           enddo; enddo
        enddo
     endif
!!! H2O continuum
     if(iflgb(5,iw)>0) then
        read(iug,*) 
        do it=1,nt; do ip=1,np
           read(iug,*) skd(1:nch(iw),ip,it,iw)
        enddo; enddo
     endif
!!! CFC absorption
     if(iflgb(7,iw)>0) then
        read(iug,*) 
        read(iug,*) acfc(1:kcfc,iw)
     endif

!! plank functions
     read(iup,*)
     read(iup,*) aplnk(iw,1:naplk)

!! solar insolation
     read(iup,*)
     read(iup,*) fsol(iw)

!! surface properties
     read(iup,*)
     read(iup,*) sr(iw,1:nsfc)

!! rayleigh scattering
     read(iup,*)
     read(iup,*) ry(iw)
     read(iup,*)

     do im=1,kmax2
!!! moments for rayleigh scattering phase function
        read(iup,*) qmol(im)
!!! moments for particle scattering phase function
        do ipcl=1,npcl
           read(iup,*) q(iw,ipcl,1:nvar(ipcl,2),im)
        enddo
     enddo

   enddo
   return
end subroutine gttbl

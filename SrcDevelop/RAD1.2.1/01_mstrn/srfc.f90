FUNCTION SRFC(WL,PRG)
! Surface Albedo
!--- History
! 93. 3.31(quarante) C KWL -> KWL1
!     4. 4 Set IS
!--- INPUT
! WL      R          WAVELENGTH (MICRON)
! PRG     R(5)   PARAMETERS FOR SURFACE CONDITION
!                    (1)             (2)             (3)
!                    1 OCEAN         COS(SOL.ZEN)    ATMOS.TRANS.
!                    2 BARE LAND     WATER CONTENT
!                    3 DRY LAND
!                    4 LOW PLANTS
!                    5 FOREST
!                    6 SNOW
!                    7 ICE
!--- OUTPUT
! SRFC      R        SURFACE ALBEDO
!                    IF 99 THEN ERROR
!--
  implicit none
  integer,parameter::kwl1=21,ksr=7
  real(8),intent(in)::wl,prg(5)
  real(8)::srfc
!  WORKING AREA
  real(8),save::wlmn=0.2d0, wlmx=200.d0
  real(8),save::r(kwl1,ksr)
  integer::is,iw,iw1
  real(8)::dwl,ddwl,am,tr,w,r1,r2,ssrfc
!--exec
!  OCEAN
  r(1:kwl1,1)=0.06d0
!  WET-LAND
  r(1:kwl1,2)=0.06d0
!  DRY-LAND
  r(1:kwl1,3)=(/0.30,0.30,0.30,0.30,0.30,0.30,0.20,0.10,0.06,0.06 &
       ,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06/)
!  LOW PLANTS
  r(1:kwl1,4)=(/0.20,0.20,0.20,0.20,0.30,0.30,0.20,0.10,0.06,0.06 &
       ,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06/)
!  FOREST
  r(1:kwl1,5)=(/0.06,0.06,0.06,0.10,0.30,0.30,0.20,0.10,0.06,0.06 &
       ,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06/)
!  SNOW
  r(1:kwl1,6)=(/0.80,0.80,0.80,0.80,0.80,0.80,0.60,0.30,0.20,0.10 &
       ,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06/)
!  ICE
  r(1:kwl1,7)=(/0.50,0.50,0.50,0.50,0.50,0.50,0.40,0.20, 0.10,0.06 &
       ,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06/)

  IS=int(PRG(1)+0.1)
  DWL=LOG(WL/WLMN)/LOG(WLMX/WLMN)*dble(KWL1-1)
  IW=int(DWL+1); DDWL=DWL+1.d0-dble(IW)
  IW=MAX(IW,1); IW=MIN(IW,KWL1)
  IW1=IW+1

! OCEAN
  IF(IS==1) THEN
     AM=PRG(2)
     TR=PRG(3)
     SRFC=SSRFC(AM,TR)
! BARE LAND
  elseif(IS==2) THEN
     W=PRG(2)
     R1=W*R(IW ,2)+(1-W)*R(IW ,3)
     R2=W*R(IW1,2)+(1-W)*R(IW1,3)
     SRFC=(1-DDWL)*R1+DDWL*R2
! DRY LAND, LOW PLANTS, FOREST, SNOW, ICE
  else
     SRFC=(1-DDWL)*R(IW,IS)+DDWL*R(IW1,IS)
  ENDIF
  RETURN
END function srfc

FUNCTION SSRFC(AM,TR)
! Sea surface reflectance by Payne
!--- HISTORY 
! 92. 4.28 CREATED
!--- INPUT
! AM       R       COS(SOLAR ZENITH ANGLE) (0-1)
! TR       R       FLUX ATMOSPHERIC TRANSMITTANCE (0-1)
!--- OUTPUT
! SSRFC    RF      ALBEDO OF OCEAN SURFACE
!--
  implicit none
  real(8),intent(in)::am,tr
  real(8)::ssrfc
  integer,parameter::m1=3,m2=5
  real(8),save::c(m2,m1)
  integer::k,l
  real(8)::am1,tr1,s,a,t

  c(1:m2,1)=(/-2.8108d+00, -1.3651d+00,  2.9210d+01, -4.3907d+01,  1.8125d+01/)
  c(1:m2,2)=(/ 6.5626d-01, -8.7253d+00, -2.7749d+01,  4.9486d+01, -1.8345d+01/)
  c(1:m2,3)=(/-6.5423d-01,  9.9967d+00,  2.7769d+00, -1.7620d+01,  7.0838d+00/)
!--exec
  AM1=MAX(AM,0.0349d0); AM1=MIN(AM1,0.961d0); TR1=MAX(TR,0.05d0)
  S=0.d0; a=1.d0
  do k=1,m1
     t=1.d0
     do l=1,m2
        s=s+c(l,k)*a*t
        t=t*tr1
     enddo
     a=a*am1
  enddo
  SSRFC=EXP(S)
  RETURN
END function ssrfc

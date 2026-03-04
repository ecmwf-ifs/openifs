! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE SUDCMIP12_SPEC(YDGEOMETRY,KTESTCASE,PTEMP,PDIV,PVOR,PLNSP,POROG)

!     Purpose.
!     --------
!           Initialize spectral fields for DCMIP12 test cases

!        Explicit arguments :
!        --------------------


!     Author.
!     -------
!        S.Malardel *ECMWF*

!     Modifications.
!     --------------
!        F.Voitus 23.04.18 

USE GEOMETRY_MOD, ONLY : GEOMETRY
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN    , ONLY : NULOUT, NULERR
USE YOMCST    , ONLY : RPI, RA, RD, RG, ROMEGA, RCPD, RATM
USE YOMDYNCORE, ONLY : RT00_DYN,RP00_DYN,RU00_DYN,RST0_DYN,RAMP
USE INTDYN_MOD, ONLY : YYTXYB

IMPLICIT NONE

TYPE(GEOMETRY),     INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN)    :: KTESTCASE
REAL(KIND=JPRB),    INTENT(INOUT) :: PTEMP(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB),    INTENT(INOUT) :: PDIV(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB),    INTENT(INOUT) :: PVOR(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB),    INTENT(INOUT) :: PLNSP(YDGEOMETRY%YRDIM%NSPEC2),POROG(YDGEOMETRY%YRDIM%NSPEC2)

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JLEV, JWORD
REAL(KIND=JPRB) :: ZLONC, ZLATC, ZR, ZLOND, ZLATD, ZLF
! Mountain waves
REAL(KIND=JPRB) :: ZA, ZB, ZLAMBDA, ZKHI, ZAMP
REAL(KIND=JPRB) :: ZZ, ZRAD, ZRADX, ZRADY
REAL(KIND=JPRB) :: ZVORG(YDGEOMETRY%YRGEM%NGPTOT,1:YDGEOMETRY%YRDIMV%NFLEVG)
! baroclinic case
REAL(KIND=JPRB) :: ZCOSLAT,ZSINLAT
REAL(KIND=JPRB) :: ZETA0,ZETAt,ZETAs,ZETAv(YDGEOMETRY%YRDIMV%NFLEVG),ZUP
REAL(KIND=JPRB) :: ZU0, ZT0, ZQ0, ZVP00, ZST,ZH_TROPO,ZTH_TROPO
REAL(KIND=JPRB) :: ZLAPS,ZGAMMA,ZDT0,ZPHIW,ZVPW
REAL(KIND=JPRB) :: ZPHI_F(YDGEOMETRY%YRGEM%NGPTOT,1:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZU(YDGEOMETRY%YRGEM%NGPTOT,1:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZV(YDGEOMETRY%YRGEM%NGPTOT,1:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZT(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZOROGG(YDGEOMETRY%YRGEM%NGPTOT), ZPSP(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZQ(YDGEOMETRY%YRGEM%NGPTOT,1:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPHIVERT(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTVERT(YDGEOMETRY%YRDIMV%NFLEVG)
! Tropical Cyclone
REAL(KIND=JPRB) ::    ZQT,ZZQ1,ZZQ2,ZZT,ZZZ,ZZ2
REAL(KIND=JPRB) ::    ZTV0,ZTVT
REAL(KIND=JPRB) ::    ZP0,ZDELTAP,ZZP,ZRP
REAL(KIND=JPRB) ::    ZD1,ZD2,ZD
REAL(KIND=JPRB) ::    ZWIND,ZFC
INTEGER(KIND=JPIM) :: ILOOP
REAL(KIND=JPRB) ::    ZXYB   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB) ::    ZPRESH(YDGEOMETRY%YRGEM%NGPTOT,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZPRESF (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZPHI(YDGEOMETRY%YRGEM%NGPTOT,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZPHIF(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZRDGAZ(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) ::    ZGELAM(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) ::    ZN1,ZN2
REAL(KIND=JPRB) ::    ZP_TROPO(YDGEOMETRY%YRGEM%NGPTOT),ZT_TROPO(YDGEOMETRY%YRGEM%NGPTOT)
!     ------------------------------------------------------------------

#include "abor1.intfb.h"

#include "reespe.intfb.h"
#include "uvspe.intfb.h"
#include "gphpre.intfb.h"
#include "gpgeo.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUDCMIP12_SPEC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, &
 & YDCVER=>YDGEOMETRY%YRCVER,YDVERT_GEOM=>YDGEOMETRY%YRVERT_GEOM)

ASSOCIATE( NRESOL=>YDDIM%NRESOL, NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, &
 & NGPTOT=>YDGEM%NGPTOT)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Montain waves
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF( KTESTCASE == 20 ) THEN
! Steady state with orography

! Constant static stability
  ZGAMMA=0.0065_JPRB
  ZT0=300._JPRB
  ZVP00=100000._JPRB
! mountain position (centre)
  ZLATC =  0._JPRB
  ZLONC =  3._JPRB/2._JPRB*RPI
  ! mountain height
  ZAMP  =  2000.0_JPRB
  ! mountain radius
  ZR   = 3._JPRB/4._JPRB*RPI
  ! Mountain Oscillation half-width
  ZKHI = RPI/16._JPRB

! orography
DO JWORD=1,NGPTOT
  ZZ   = SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
  &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
  ZRAD = ACOS(ZZ)
  IF(ZRAD<ZR) THEN
  ZOROGG(JWORD) = 0.5_JPRB*RG*ZAMP*(1._JPRB+cos(RPI*ZRAD/ZR)) &
 &  *cos(RPI*ZRAD/ZKHI)**2
  ELSE
  ZOROGG(JWORD) = 0.0_JPRB
  ENDIF
ENDDO
! set spectral orography
CALL REESPE(YDGEOMETRY,1,1,POROG,ZOROGG)

! "balanced" surface pressure
DO JWORD=1,NGPTOT
  ZPRESH(JWORD,NFLEVG) = ZVP00*(1._JPRB-ZGAMMA/ZT0*ZOROGG(JWORD)/RG)**(RG/(RD*ZGAMMA))
  ZPSP(JWORD) = LOG(ZPRESH(JWORD,NFLEVG))
ENDDO

! spectral surface pressure
  CALL REESPE(YDGEOMETRY,1,1,PLNSP,ZPSP)

! Compute pressure
  CALL GPHPRE (NGPTOT,NFLEVG,1,NGPTOT,YDVAB,YDCVER,ZPRESH,PXYB=ZXYB,PRESF=ZPRESF)

! temperature profile
DO JLEV=1,NFLEVG
  DO JWORD=1,NGPTOT
    ZT(JWORD,JLEV)= ZT0*(ZPRESF(JWORD,JLEV)/ZVP00)**((RD*ZGAMMA)/RG)
  ENDDO
ENDDO

! spectral temperature
CALL REESPE(YDGEOMETRY,NFLEVL,NFLEVG,PTEMP,ZT)

! no wind (ZVOR/ZDIV are spectral arrays)
PVOR=0._JPRB
PDIV=0._JPRB

ELSEIF( (KTESTCASE == 211) .OR. &
 &      (KTESTCASE == 212) .OR. &
 &      (KTESTCASE == 213) .OR. &
 &      (KTESTCASE == 214) .OR. &
 &      (KTESTCASE == 215) .OR. &
 &      (KTESTCASE == 216)      &
 &                             ) THEN
! Mountain wave cases, isothermal, no shear

!!!
! KTESTCASE=211 : Schar Mountain (DCMIP12 test case 2.1)
! KTESTCASE=212 : Agnesi (?) Mountain
! KTESTCASE=213 : Mountain as in case 20 (scale in radian, so rescale
!                 automatically with the planet radius
! KTESTCASE=214 : Schar Mountain (DCMIP12 test case 2.1 modified as in Klemp et al, 2015)
! KTESTCASE=215 : Zaengl MWR 2012 steep orography test
!!!
!  ZT0=300._JPRB
!  ZVP00=100000._JPRB
!  ZVP00=100907._JPRB
!  ZU0=20._JPRB
   ZT0=RT00_DYN
   ZU0=RU00_DYN
   ZVP00=RP00_DYN
   ZST=RST0_DYN
   ZH_TROPO=0._JPRB !isothermal case

!===========================
!       Prescribe orography
!===========================
  IF (KTESTCASE==211) THEN
  ! mountain position (centre)
  ZLATC =  0._JPRB
!  ZLONC =  0.25_JPRB*RPI
  ZLONC =  RPI
!  ZLATC =  RPI/6._JPRB
!  ZLONC =  1.5_JPRB*RPI
  ! mountain height
  ZAMP  =  250.0_JPRB
  ! mountain radius
  ZR   = 5000._JPRB
!  ZR   = 5000._JPRB*166.7_JPRB/RUSPLRADI
  ! off centering parameter
  ZLAMBDA=4000._JPRB
!  ZLAMBDA=4000._JPRB*166.7_JPRB/RUSPLRADI

  DO JWORD=1,NGPTOT
    ZZ   = SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
     &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
    ZRAD = RA*ACOS(ZZ)
    ZOROGG(JWORD) = RG*ZAMP*EXP(-ZRAD**2/ZR**2)*(COS(RPI*ZRAD/ZLAMBDA))**2 
  ENDDO

  ELSEIF (KTESTCASE==212) THEN
  ! no shear/shear
  IF (ZST==0._JPRB) THEN
    ZH_TROPO=0._JPRB
  ELSE
    ZH_TROPO=11000._JPRB
  ENDIF

  ! mountain height
  !ZAMP = 175.0_JPRB
  ! ZAMP  =  500.0_JPRB

  ! DBG JOVI
  ZAMP = 100.0   ! LI NH regime
  ! ZAMP = 1000.0  ! NL NH regime

  ! 2 foci points N/S from equator, determine the width of the mountain
  ! centre of ellipse at equator
  ZLATC =  0._JPRB
  ZLONC =  1.5_JPRB*RPI

  ! half width of 2d mountain 
  ! (elliptical shape allowed)
  ! long distance of ellipse, x-direction
  ZA = 1000.0_JPRB
  ! long distance of ellipse, y-direction
  ZB = 3000._JPRB

  DO JWORD=1,NGPTOT
    ! projection on x, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(ZLATC)+&
     &COS(ZLATC)*COS(ZLATC)*COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
    ZRADX = RA*ACOS(ZZ)      
    ! projection on y, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
     &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*COS(ZLONC-ZLONC)
    ZRADY = RA*ACOS(ZZ)
    ZOROGG(JWORD) = RG*ZAMP/(sqrt(1._JPRB+(ZRADX/ZA)**2+(ZRADY/ZB)**2))         
  ENDDO

  ! add a narrow montain upstream of the main one
  ! mountain height
  !ZAMP = 175.0_JPRB
  ZAMP  =  0.0_JPRB
  ! 2 foci points N/S from equator, determine the width of the mountain
  ! centre of ellipse at equator
  ZLATC =  0._JPRB
  ZLONC =  1.5_JPRB*RPI-7500._JPRB/RA
  ! long distance of ellipse, x-direction
  ZA = 1000.0_JPRB
  ! long distance of ellipse, y-direction
  ZB = 3000._JPRB
  DO JWORD=1,NGPTOT
    ! projection on x, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(ZLATC)+&
     &COS(ZLATC)*COS(ZLATC)*COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
    ZRADX = RA*ACOS(ZZ)      
    ! projection on y, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
     &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*COS(ZLONC-ZLONC)
    ZRADY = RA*ACOS(ZZ)
    ZOROGG(JWORD) = ZOROGG(JWORD)+RG*ZAMP/(sqrt(1._JPRB+(ZRADX/ZA)**2+(ZRADY/ZB)**2))         
  ENDDO

! add a narrow montain downstream of the main one
  ! mountain height
  !ZAMP = 175.0_JPRB
  ZAMP  =  0.0_JPRB
  ! 2 foci points N/S from equator, determine the width of the mountain
  ! centre of ellipse at equator
  ZLATC =  0._JPRB
  ZLONC =  1.5_JPRB*RPI+7500._JPRB/RA
  ! long distance of ellipse, x-direction
  ZA = 1500.0_JPRB
  ! long distance of ellipse, y-direction
  ZB    =40000._JPRB
  DO JWORD=1,NGPTOT
    ! projection on x, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(ZLATC)+&
     &COS(ZLATC)*COS(ZLATC)*COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
    ZRADX = RA*ACOS(ZZ)      
    ! projection on y, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
     &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*COS(ZLONC-ZLONC)
    ZRADY = RA*ACOS(ZZ)
    ZOROGG(JWORD) = ZOROGG(JWORD)+RG*ZAMP/(sqrt(1._JPRB+(ZRADX/ZA)**2+(ZRADY/ZB)**2))         
  ENDDO

  ELSEIF (KTESTCASE==213) THEN
! mountain position (centre)
  ZLATC =  0._JPRB
  ZLONC =  3._JPRB/2._JPRB*RPI
  ! mountain height
  ZAMP  =  2000.0_JPRB
  ! mountain radius
  ZR   = 3._JPRB/4._JPRB*RPI
  ! Mountain Oscillation half-width
  ZKHI = RPI/16._JPRB

! orography
  DO JWORD=1,NGPTOT
  ZZ   = SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
  &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
  ZRAD = ACOS(ZZ)
    IF(ZRAD<ZR) THEN
    ZOROGG(JWORD) = 0.5_JPRB*RG*ZAMP*(1._JPRB+cos(RPI*ZRAD/ZR)) &
 &    *cos(RPI*ZRAD/ZKHI)**2
    ELSE
    ZOROGG(JWORD) = 0.0_JPRB
    ENDIF
  ENDDO

  ELSEIF (KTESTCASE==214) THEN
  ! mountain position (centre)
  ZLATC =  0._JPRB
  ZLONC =  0._JPRB
  ! mountain height
  ZAMP  =  250.0_JPRB
  ! mountain radius
!  ZR   = 5000._JPRB*166.7_JPRB/RUSPLRADI
  ZR   = 5000._JPRB
  ! off centering parameter
!  ZLAMBDA=4000._JPRB*166.7_JPRB/RUSPLRADI
  ZLAMBDA=4000._JPRB
  DO JWORD=1,NGPTOT
    ! original GELAM between 0 and 2PI
    ! new GELAM between -PI and PI

    IF (YDGSGEOM_NB%GELAM(JWORD)>RPI) THEN
      ZGELAM(JWORD)=YDGSGEOM_NB%GELAM(JWORD)-2.0_JPRB*RPI
    ELSE
      ZGELAM(JWORD)=YDGSGEOM_NB%GELAM(JWORD)
    ENDIF
    ZRAD = RA*(ZGELAM(JWORD)-ZLONC)
    ZOROGG(JWORD) = RG*ZAMP*EXP(-ZRAD**2/ZR**2)*(COS(RPI*ZRAD/ZLAMBDA))**2 &
 &                * COS(YDGSGEOM_NB%GELAT(JWORD))
  ENDDO

  ELSEIF (KTESTCASE==215) THEN
  ! position at equator
  ZLATC =  0._JPRB
  ZLONC =  0._JPRB

  ! mountain height
!  ZAMP = 1750._JPRB
  ZAMP = RAMP
  ! half-width
  ZR = 2000._JPRB

  DO JWORD=1,NGPTOT
    ZZ   = SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
     &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
    ZRAD = RA*ACOS(ZZ)
   !ZOROGG(JWORD) = RG*ZAMP/(sqrt(1._JPRB+(ZRAD/ZR)**2))**3
    ZOROGG(JWORD) = RG*ZAMP*EXP(-(ZRAD/ZR)**2)
  ENDDO
  ELSEIF (KTESTCASE==216) THEN
  ! no shear/shear
  IF (ZST==0._JPRB) THEN
    ZH_TROPO=0._JPRB
  ELSE
    ZH_TROPO=10500._JPRB
  ENDIF

  ! mountain height
  !ZAMP = 175.0_JPRB
  ZAMP  =  500.0_JPRB
  ! 2 foci points N/S from equator, determine the width of the mountain
  ! centre of ellipse at equator
  ZLATC =  0._JPRB
  ZLONC =  1.5_JPRB*RPI
  ZLATD =  RPI/3._JPRB
  ZLOND =  1.5_JPRB*RPI
  ZZ   = SIN(ZLATD)*SIN(ZLATC)+&
     &COS(ZLATD)*COS(ZLATC)*COS(ZLOND-ZLONC) 
  ZLF = RA*ACOS(ZZ) 
  ! long distance of ellipse, x-direction
  ZA = 2500.0_JPRB
  ! long distance of ellipse, y-direction
!  ZB    =40000._JPRB
  ZB    =ZLF
  DO JWORD=1,NGPTOT
    ! projection on x, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(ZLATC)+&
     &COS(ZLATC)*COS(ZLATC)*COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
    ZRADX = RA*ACOS(ZZ)      
    ! projection on y, distance to the centre
    ZZ   = SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
     &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*COS(ZLONC-ZLONC)
    ZRADY = RA*ACOS(ZZ)
    ZOROGG(JWORD) = RG*ZAMP/(sqrt(1._JPRB+(ZRADX/ZA)**2+(ZRADY/ZB)**2))         
  ENDDO
  ENDIF

  ! set spectral orography
  CALL REESPE(YDGEOMETRY,1,1,POROG,ZOROGG)

  !===========================================
  ! prescribe temperature and surface pressure
  !===========================================


! Compure surface pressure as if uniform T=ZT0 and U=ZU0

  ZZ = 0.5_JPRB*(RA*ROMEGA*2.0_JPRB+ZU0)*ZU0
  ZZ=ZZ/(RD*ZT0)
  DO JWORD=1,NGPTOT
    ZPSP(JWORD) = LOG(ZVP00)-ZZ*YDGSGEOM_NB%GEMU(JWORD)**2-ZOROGG(JWORD)/(RD*ZT0)
  ENDDO

  DO JWORD = 1,NGPTOT
    ZPRESH(JWORD,NFLEVG)=EXP(ZPSP(JWORD))
  ENDDO

  CALL GPHPRE (NGPTOT,NFLEVG,1,NGPTOT,YDVAB,YDCVER,ZPRESH,PXYB=ZXYB,PRESF=ZPRESF)

IF (ZH_TROPO == 0._JPRB) THEN
!!!! isothermal case

  DO JLEV=1,NFLEVG
    DO JWORD=1,NGPTOT
      ZT(JWORD,JLEV)= ZT0
    ENDDO
  ENDDO

ELSE
! Brunt-Vaisala frequency 
! ZN1(in practice here, dtheta/dz)
    ZN1 = (1.E-2_JPRB)**2*ZT0/RG
! ZN1(in practice here, dtheta/dz * T/TH): isothermal strato
    ZN2 = (RG/RCPD)

! Compute first guess geopotential with isothermal atmosphere T=ZT0
! 1st Guess for T, to compute first guess for Z
  ZT = ZT0 ! T at surface
  ZRDGAZ = RD

! Theta at tropo
  ZTH_TROPO= (ZT0+ZN1*ZH_TROPO)

  DO JWORD = 1,NGPTOT
    ZPHI(JWORD,NFLEVG) = ZOROGG(JWORD)
  ENDDO

!!! LOOP for adjustment Z/T
  DO ILOOP=1,10
!!! LOOP for adjustment Z/T

! compute geopotentiel of reference atmosphere
    CALL GPGEO (NGPTOT,1,NGPTOT,NFLEVG,&
              & ZPHI, ZPHIF, &
              & ZT,ZRDGAZ,ZXYB(1,1,YYTXYB%M_LNPR),&
              & ZXYB(1,1,YYTXYB%M_ALPH), YDVERT_GEOM)

! compute pressure at ZH_TROPO
    DO JWORD = 1,NGPTOT
    DO JLEV = NFLEVG-1,0,-1
        ZZZ = ZPHI(JWORD,JLEV)/RG
      IF( ZZZ > ZH_TROPO ) THEN
        ZP_TROPO(JWORD)=  ZPRESH(JWORD,JLEV) +(ZPRESH(JWORD,JLEV+1)-ZPRESH(JWORD,JLEV)) &
 &                        /(ZPHI(JWORD,JLEV+1)-ZPHI(JWORD,JLEV)) &
 &                        *(ZH_TROPO*RG-ZPHI(JWORD,JLEV))
        ZT_TROPO(JWORD)=  ZTH_TROPO*(ZP_TROPO(JWORD)/RATM)**(RD/RCPD)
        EXIT
      ENDIF
    ENDDO
    ENDDO
! Temperature at ZH_TROPO

    DO JWORD = 1,NGPTOT
    DO JLEV = 1, NFLEVG
    ZZZ = ZPHIF(JWORD,JLEV)/RG
      IF( ZZZ >= ZH_TROPO ) THEN
        ZT(JWORD,JLEV)= ZT_TROPO(JWORD)
        ZT(JWORD,JLEV)=MAX(180._JPRB,ZT(JWORD,JLEV))
      ELSE
        ZT(JWORD,JLEV)= (ZT0+ZN1*ZZZ)*(ZPRESF(JWORD,JLEV)/RATM)**(RD/RCPD)
      ENDIF
    ENDDO
    ENDDO
!!! LOOP for adjustment Z/T
  ENDDO 
!!! LOOP for adjustment Z/T
ENDIF

  ! spectral
  CALL REESPE(YDGEOMETRY,NFLEVL,NFLEVG,PTEMP,ZT)
  CALL REESPE(YDGEOMETRY,1,1,PLNSP,ZPSP)

!-------------------------------------------------------------------
!!! VOR/DIV in spectral space

  DO JWORD=1,NSPEC2
    PDIV(:,JWORD)=0.0_JPRB
    PVOR(:,JWORD)=0.0_JPRB
  ENDDO

!!!--- No shear case
IF (RST0_DYN == 0) THEN
  DO JLEV=NFLEVG,1,-1
    DO JWORD=1,NGPTOT
      ZVORG(JWORD,JLEV) = 2.0_JPRB*ZU0/RA*YDGSGEOM_NB%GEMU(JWORD)
    ENDDO
  ENDDO
  ! spectral
  CALL REESPE(YDGEOMETRY,NFLEVL,NFLEVG,PVOR,ZVORG)
ELSE
  ! profile and wind shear following T. Keller (1994) JAS, 51, 1915--1929

  ! Ri = 16     == N^2/(ZST)^2
  ! set ZU0 = 10.0
  ! RA = 128000./(2*pi),  dx - 250 m
  ! use IFS L91
  ! N  = 0.01  - T0 == 957 K
  ! c  = -2.5E-03
  ! ZST = ZU0 * c
  ! upper threshold above which U=const (in meter)

    DO JWORD = 1,NGPTOT
    DO JLEV = 1, NFLEVG
      ZZZ=ZPHIF(JWORD,JLEV )/RG
      IF( (ZZZ >= ZH_TROPO) ) THEN
        ZU(JWORD,JLEV) = ZU0 + ZST*ZH_TROPO
      ELSE
        ZU(JWORD,JLEV) = ZU0 + ZST*ZZZ
      ENDIF
    ENDDO
    ENDDO

  DO JLEV=NFLEVG,1,-1
    DO JWORD=1,NGPTOT
      ZVORG(JWORD,JLEV) = 2.0_JPRB*ZU(JWORD,JLEV)/RA*YDGSGEOM_NB%GEMU(JWORD)
    ENDDO
  ENDDO
  ! spectral
  CALL REESPE(YDGEOMETRY,NFLEVL,NFLEVG,PVOR,ZVORG)

ENDIF

!==================================================================
! Case 4.1 and 4.2:  Baroclinic waves (dry or moist) 
!==================================================================
ELSEIF( (KTESTCASE == 41) .OR. (KTESTCASE == 42) .OR. (KTESTCASE == 43)) THEN

    ZETA0=0.252_JPRB
    ZETAt=0.2_JPRB
    ZETAs=1._JPRB
    ZU0=35._JPRB
    ZT0=288._JPRB
    ZDT0=4.8E5_JPRB
    ZGAMMA=0.005_JPRB

    ZLONC=RPI/9._JPRB
    ZLATC=2._JPRB*RPI/9._JPRB

!!! for zonal jet only (no initial perturbation)
!    ZUP=0._JPRB
!!! Amplitude of initial zonal wind perturbation
    ZUP=1._JPRB

  DO JLEV=1,NFLEVG
    ZETAv(JLEV)=0.5_JPRB*RPI*(YDVETA%VETAF(JLEV)-ZETA0)
  ENDDO

! temperature, vert. profile

  ZLAPS = RD*ZGAMMA/RG

  DO JLEV=1,NFLEVG
    IF (YDVETA%VETAF(JLEV)>=ZETAt)  THEN 
    ZTVERT(JLEV)=ZT0*YDVETA%VETAF(JLEV)**ZLAPS
    ZPHIVERT(JLEV)=ZT0*RG/ZGAMMA*(1._JPRB-YDVETA%VETAF(JLEV)**ZLAPS)
    ELSE
    ZTVERT(JLEV)=ZT0*YDVETA%VETAF(JLEV)**ZLAPS+ &
&                ZDT0*(ZETAt-YDVETA%VETAF(JLEV))**5
    ZPHIVERT(JLEV)=ZT0*RG/ZGAMMA*(1._JPRB-YDVETA%VETAF(JLEV)**ZLAPS) &
&  - RD*ZDT0*( (log(YDVETA%VETAF(JLEV)/ZETAt)+137._JPRB/60._JPRB)*ZETAt**5 &
&              - 5._JPRB*ZETAt**4*YDVETA%VETAF(JLEV) &
&              + 5._JPRB*ZETAt**3*YDVETA%VETAF(JLEV)**2 &
&              - 10._JPRB/3._JPRB*ZETAt**2*YDVETA%VETAF(JLEV)**3 &
&              + 5._JPRB/4._JPRB*ZETAt*YDVETA%VETAF(JLEV)**4 &
&              - 1._JPRB/5._JPRB*YDVETA%VETAF(JLEV)**5 )
    ENDIF
  ENDDO

! Prescribe orography (not flat, in order to balance uniform surf. pres.)

  DO JWORD=1,NGPTOT
    ZR= RA*ACOS(SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+&
   & COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))*&
   & COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC))
    ZCOSLAT=cos(YDGSGEOM_NB%GELAT(JWORD))
    ZSINLAT=sin(YDGSGEOM_NB%GELAT(JWORD))

      ZOROGG(JWORD) = ZU0*cos(0.5_JPRB*RPI*(ZETAs-ZETA0))**(1.5_JPRB) * &
&       ( (-2._JPRB*ZSINLAT**6* &
&       (ZCOSLAT**2+1._JPRB/3._JPRB)+10._JPRB/63._JPRB)* &
&       ZU0*cos(0.5_JPRB*RPI*(ZETAs-ZETA0))**(1.5_JPRB) + &
&       (8._JPRB/5._JPRB*ZCOSLAT**3* &
&       (ZSINLAT**2+2._JPRB/3._JPRB)-RPI/4._JPRB)* &
&       RA*ROMEGA )

!       Prescribe wind and temperature

    DO JLEV=1,NFLEVG
      ZU(JWORD,JLEV) = ZU0* &
&       cos(ZETAv(JLEV))**(1.5_JPRB) * &
&       sin(2._JPRB*YDGSGEOM_NB%GELAT(JWORD))**2 + &
&       ZUP*exp(-(ZR/(RA/10._JPRB))**2)

      ZV(JWORD,JLEV) = 0.0_JPRB

      ZT(JWORD,JLEV)= ZTVERT(JLEV) + &
&     3._JPRB/4._JPRB*YDVETA%VETAF(JLEV)*RPI*ZU0/RD*sin(ZETAv(JLEV))*&
&     SQRT(cos(ZETAv(JLEV)))*&
&       ( (-2._JPRB*ZSINLAT**6* &
&       (ZCOSLAT**2+1._JPRB/3._JPRB)+10._JPRB/63._JPRB)* &
&       2._JPRB*ZU0*cos(ZETAv(JLEV))**(1.5_JPRB) + &
&       (8._JPRB/5._JPRB*ZCOSLAT**3* &
&       (ZSINLAT**2+2._JPRB/3._JPRB)-RPI/4._JPRB)* &
&       RA*ROMEGA )
      ZPHI_F(JWORD,JLEV)= ZPHIVERT(JLEV) + &
&       ZU0*cos(ZETAv(JLEV))**(1.5_JPRB) * &
&       ( (-2._JPRB*ZSINLAT**6* &
&       (ZCOSLAT**2+1._JPRB/3._JPRB)+10._JPRB/63._JPRB)* &
&       ZU0*cos(ZETAv(JLEV))**(1.5_JPRB) + &
&       (8._JPRB/5._JPRB*ZCOSLAT**3* &
&       (ZSINLAT**2+2._JPRB/3._JPRB)-RPI/4._JPRB)* &
&       RA*ROMEGA )
   ENDDO
  ENDDO

! compute virtual temperature

ZPHIW=2._JPRB*RPI/9._JPRB
! if non dry case
!ZQ0=0.021_JPRB
! dry case
ZQ0=0.0_JPRB
ZVP00 = 1.E5_JPRB
ZVPW = 34000._JPRB

DO JWORD = 1,NGPTOT
  DO JLEV = 1, NFLEVG
      ZQ(JWORD,JLEV) = ZQ0 * exp(-(YDGSGEOM_NB%GELAT(JWORD)/ZPHIW)**4) &
 &                     *exp(-((YDVETA%VETAF(JLEV)-1._JPRB)*ZVP00/ZVPW)**2)

      ZT(JWORD,JLEV)=ZT(JWORD,JLEV)/(1._JPRB+0.608_JPRB*ZQ(JWORD,JLEV))
  ENDDO
ENDDO

! Prescribe surface pressure
! constant surface pressure
! ZVP00 = 1.E5_JPRB
  ZVP00 = 101325._JPRB
  DO JWORD=1,NGPTOT
    ZPSP(JWORD) = LOG ( ZVP00 )
  ENDDO

!!!! From GP space to SP space
! set spectral orography
CALL REESPE(YDGEOMETRY,1,1,POROG,ZOROGG)
!!! VOR/DIV in spectral space
PVOR=0._JPRB
PDIV=0._JPRB
CALL UVSPE(YDGEOMETRY,PVOR,PDIV,ZU,ZV,NFLEVL,NFLEVG,1)
! spectral temperature
CALL REESPE(YDGEOMETRY,NFLEVL,NFLEVG,PTEMP,ZT)
! spectral surface pressure
CALL REESPE(YDGEOMETRY,1,1,PLNSP,ZPSP)

  WRITE(NULOUT,'(''END OF GMV computation in suspecg2 for NTESTCASE=41/42'')')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Tropical Cyclone (Reed and Jablonowski, 2011)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSEIF( (KTESTCASE == 51) .OR. (KTESTCASE == 52)) THEN

!!! TC Characteristics

ZLATC=RPI/18._JPRB
ZLONC=RPI

ZQ0=21.E-3_JPRB
ZQT=1.E-11_JPRB
ZZQ1=3000._JPRB
ZZQ2=8000._JPRB
ZZT=15000._JPRB

ZGAMMA=0.007_JPRB
ZT0=302.15_JPRB
ZTV0=302.15_JPRB*(1._JPRB+0.608_JPRB*ZQ0)
ZTVT=ZTV0-ZGAMMA*ZZT

ZP0=101500._JPRB
ZDELTAP=1115._JPRB
!ZDELTAP=0._JPRB

ZZP=7000._JPRB
ZRP=282000._JPRB

ZFC=2._JPRB*ROMEGA*SIN(ZLATC)

!===========================
!       Prescribe flat orography
!===========================
  DO JWORD=1,NGPTOT
    ZOROGG(JWORD) = 0.0_JPRB    
  ENDDO
! set spectral orography
CALL REESPE(YDGEOMETRY,1,1,POROG,ZOROGG)

! 1st Guess for T, to compute first guess for Z
  ZT = ZT0

ZRDGAZ = RD

DO JWORD = 1,NGPTOT
  ZPHI(JWORD,NFLEVG) = ZOROGG(JWORD)
ENDDO

!!! SURFACE PRESSURE
DO JWORD = 1,NGPTOT
ZR   = RA*ACOS( SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+ &
     & COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))* &
     & COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC) )
ZPRESH(JWORD,NFLEVG)=ZP0-ZDELTAP*EXP(-(ZR/ZRP)**(1.5_JPRB))
ENDDO

  CALL GPHPRE (NGPTOT,NFLEVG,1,NGPTOT,YDVAB,YDCVER,ZPRESH,PXYB=ZXYB,PRESF=ZPRESF)

DO JWORD=1,NGPTOT
    ZPSP(JWORD) = LOG(ZPRESH(JWORD,NFLEVG) )
  ENDDO
  ! spectral
  CALL REESPE(YDGEOMETRY,1,1,PLNSP,ZPSP)

!!! LOOP for adjustment Z/T
  DO ILOOP=1,100
!!! LOOP for adjustment Z/T

! compute geopotentiel of reference atmosphere
  CALL GPGEO (NGPTOT,1,NGPTOT,NFLEVG,&
              & ZPHI, ZPHIF, &
              & ZT,ZRDGAZ,ZXYB(1,1,YYTXYB%M_LNPR),&
              & ZXYB(1,1,YYTXYB%M_ALPH), YDVERT_GEOM)

    DO JWORD = 1,NGPTOT
      DO JLEV = 1, NFLEVG

        ZZZ  = ZPHIF(JWORD,JLEV )/RG

        ZR   = RA*ACOS( SIN(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))+ &
             &COS(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))* &
             &COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC) )

        ZD1  = SIN(ZLATC)*COS(YDGSGEOM_NB%GELAT(JWORD))- &
             &COS(ZLATC)*SIN(YDGSGEOM_NB%GELAT(JWORD))*&
             &COS(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
        ZD2  = COS(YDGSGEOM_NB%GELAT(JWORD))*SIN(YDGSGEOM_NB%GELAM(JWORD)-ZLONC)
        ZD   = MAX(1.E-25_JPRB,SQRT(ZD1**2+ZD2**2))

        IF (ZZZ<ZZT) THEN
          ZQ(JWORD,JLEV) = ZQ0*EXP(-(ZZZ/ZZQ1))*EXP(-(ZZZ/ZZQ2)**2)
          ZRDGAZ(JWORD,JLEV)=RD*(1._JPRB+0.608_JPRB*ZQ(JWORD,JLEV))
          IF (ZDELTAP==0._JPRB) THEN
            ZT(JWORD,JLEV) = (ZTV0-ZGAMMA*ZZZ)/(1._JPRB+0.608_JPRB*ZQ(JWORD,JLEV))
          ELSE
            ZT(JWORD,JLEV) = (ZTV0-ZGAMMA*ZZZ)/(1._JPRB+0.608_JPRB*ZQ(JWORD,JLEV))/ &
                           & ( 1._JPRB+ ( 2._JPRB*RD*(ZTV0-ZGAMMA*ZZZ)*ZZZ) / &
                           & ( RG*ZZP**2* (1._JPRB-ZP0/ZDELTAP*EXP((ZR/ZRP)**(1.5_JPRB))* &
                           &  EXP((ZZZ/ZZP)**(2._JPRB))) ) &
                           & )
          ENDIF
          IF (ZDELTAP==0._JPRB) THEN
            ZWIND= 0.0_JPRB
          ELSE
            ZWIND= - ZFC*ZR/2._JPRB + SQRT( (ZFC*ZR/2._JPRB)**2 -  &
                 &  1.5_JPRB*RD*(ZTV0-ZGAMMA*ZZZ)*(ZR/ZRP)**(1.5_JPRB) / &
                 &  ( 1._JPRB+ 2._JPRB*RD*(ZTV0-ZGAMMA*ZZZ)*ZZZ/(RG*ZZP**2) &
                 & -ZP0/ZDELTAP*EXP((ZR/ZRP)**(1.5_JPRB))* &
                 &  EXP((ZZZ/ZZP)**(2._JPRB)) ) &
                 &   )
          ENDIF

          ZU(JWORD,JLEV) =ZWIND*ZD1/ZD
          ZV(JWORD,JLEV) =ZWIND*ZD2/ZD

         ELSE
           ZT(JWORD,JLEV) = ZTVT
           ZU(JWORD,JLEV) = 0._JPRB
           ZV(JWORD,JLEV) = 0._JPRB
        ENDIF
      ENDDO
    ENDDO
!!! LOOP for adjustment Z/T
  ENDDO
!!! LOOP for adjustment Z/T

! T to spectral space
  CALL REESPE(YDGEOMETRY,NFLEVL,NFLEVG,PTEMP,ZT)

! Spectral vorticity/Divergence from GP wind

  CALL UVSPE(YDGEOMETRY,PVOR,PDIV,ZU,ZV,NFLEVL,NFLEVG,1)

ELSE
  WRITE(NULERR,'('' INVALID SETUP NTESTCASE FOR DCMIP12 FAMILY'')')
  CALL ABOR1(' INVALID SETUP NTESTCASE FOR DCMIP12 FAMILY')
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDCMIP12_SPEC',1,ZHOOK_HANDLE)
END SUBROUTINE SUDCMIP12_SPEC

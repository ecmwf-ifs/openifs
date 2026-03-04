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

SUBROUTINE LARCHE2(KPROMA,KSTART,KPROF,&
 & KSTTYP,PDSTRET,PC2M1,PC2P1,PI,PDEPI,&
 & PLOCEN,PMUCEN,YDGSGEOM,YDCSGEOM,&
 & PCOSCO,PSINCO,PSINLA,PCOPHI,&
 & KROT,PLON,PLAT,PQX,PQY)  

!**** *LARCHE2  - semi-LAgrangian scheme:
!                 Research of the Coordinates (of the medium or origin
!                 point). 2D model.

!     Purpose.
!     --------
!       Computes the longitude and latitude of the interpolation
!       point from its cartesian coordinates.
!       Then computes the vector displacement matrix
!                        I po pq I
!                        I       I
!                        I-pq po I
!       from the interpolation point to the grid point.

!**   Interface.
!     ----------
!        *CALL* *LARCHE2(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA   - horizontal dimension.
!          KSTART   - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          KSTTYP   - 1: Not tilted pole;  2: Tilted pole.
!          PDSTRET  - 2*c (where c is the stretching factor).
!          PC2M1    - c*c-1.
!          PC2P1    - c*c+1.
!          PI       - number PI
!          PDEPI    - 2 * PI
!          PLOCEN   - geographic longitude of the stretching pole.
!          PMUCEN   - sinus of the geographic latitude of the stretching pole.
!          YDGSGEOM - grid point geometry.
!          YDCSGEOM - computational sphere geometry.
!          PCOSCO   - cos(Longitude-Longitude(grid-point))*cos(Latitude)
!                     of the interpolation point (geographic lon and lat).
!          PSINCO   - sin(Longitude-Longitude(grid-point))*cos(Latitude)
!                     of the interpolation point (geographic lon and lat).
!          PSINLA   - sinus of the interpolation point geographic latitude.
!          PCOPHI   - cosinus of the geographic angle between the
!                     interpolation point and the grid-point.
!          KROT     - KROT=1: computation of the elements po and pq
!                     of the wind displacement matrix.
!                     KROT=0: no computation.

!        OUTPUT:
!          PLON     - longitude of the interpolation point in the
!                     computational sphere.
!          PLAT     - latitude of the interpolation point in the
!                     computational sphere.
!          PQX      - first element of the wind displacement matrix (p,q).
!          PQY      - second element of the wind displacement matrix (p,q).

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        No external.

!     Reference.
!     ----------

!     Author.
!     -------
!      K. YESSAD, from LARCHE.
!      Original : SEP 2010. 

!     Modifications.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCSGEOM , ONLY : TCSGEOM
USE YOMGSGEOM , ONLY : TGSGEOM
USE YOMJFH    , ONLY : N_VMASS
USE LARCHE_HLP, ONLY : JFH_VMOD, JFH_VIF2

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTTYP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSTRET
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2M1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2P1
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDEPI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLOCEN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMUCEN
TYPE(TGSGEOM)     ,INTENT(IN)    :: YDGSGEOM
TYPE(TCSGEOM)     ,INTENT(IN)    :: YDCSGEOM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOSCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINLA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOPHI(KPROMA) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KROT 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLON(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLAT(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQX(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQY(KPROMA) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZTMP1(KPROMA)
REAL(KIND=JPRB) :: ZTMP2(KPROMA)
REAL(KIND=JPRB) :: ZTMP4(KPROMA)
REAL(KIND=JPRB) :: ZTMP5(KPROMA)
REAL(KIND=JPRB) :: ZTMP6(KPROMA)
REAL(KIND=JPRB) :: ZTMP7(KPROMA)

INTEGER(KIND=JPIM) :: JROF, ILEN

REAL(KIND=JPRB) :: Z1ST, Z1STT, ZCA, ZCOSCO, ZCOSL2, ZCOSLA,&
 & ZCOSLR, ZDECO, ZDI, ZEPS, ZFACT, ZMU, ZP, &
 & ZPCLAP, ZPCLOP, ZPSLAP, ZPSLOP, ZQ, ZSA, &
 & ZSINCO, ZSINLA, ZXX, ZYY, ZZP, ZZW, &
 & ZZX, ZZY, ZZZ  

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARCHE2',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF LAT LON OF THE INTERPOLATION POINT.
!              IF KROT=1 COMPUTATION OF THE WIND DISPLACEMENT MATRIX
!              FROM THE INTERPOLATION POINT TO THE FINAL POINT.
!              ( T FOR LATITUDE THETA, L FOR LONGITUDE LAMBDA).
!              PQX= ( 1 / (1+cos(PHI)) )
!                  *( cos(TG)*cos(T) + (1+sin(TG)*sin(T))*cos(L-LG) )
!              PQY= (-1 / (1+cos(PHI)) )
!                  *( sin(TG)+sin(T) ) * sin(L-LG)
!     ------------------------------------------------------------------

ZEPS=1.E-14_JPRB

IF(KSTTYP == 1) THEN

  !      POLE STRETCHING, NO POLE TILTING.

  IF(N_VMASS <= 0) THEN

    ! ---- Not using MASS library vector routines

    DO JROF=KSTART,KPROF
      ZCOSLR=SQRT(PCOSCO(JROF)*PCOSCO(JROF)&
       & +PSINCO(JROF)*PSINCO(JROF))  
      ZCOSLR=MAX(ZEPS,ZCOSLR)
      ZFACT = PDSTRET/(PC2P1-PC2M1*PSINLA(JROF))
      ZCOSCO= PCOSCO(JROF)*ZFACT
      ZSINCO= PSINCO(JROF)*ZFACT
      ZCOSLA= ZCOSLR*ZFACT
      ZSINLA= (PC2P1*PSINLA(JROF)-PC2M1)/(PC2P1-PC2M1*PSINLA(JROF))
      ZXX   = YDCSGEOM%RCOLON(JROF)*ZCOSCO-YDCSGEOM%RSILON(JROF)*ZSINCO
      ZYY   = YDCSGEOM%RSILON(JROF)*ZCOSCO+YDCSGEOM%RCOLON(JROF)*ZSINCO
      ZMU   = MAX(-1.0_JPRB,MIN(1.0_JPRB,ZSINLA))
      ZZW   = MAX(-1.0_JPRB,MIN(1.0_JPRB,ZXX/ZCOSLA))
      PLAT(JROF)=ASIN(ZMU)
      PLON(JROF)=MOD(PI+SIGN(1.0_JPRB,ZYY)*(ACOS(ZZW)-PI),PDEPI)
      IF (KROT == 1) THEN
        ZDECO =1.0_JPRB/(1.0_JPRB+PCOPHI(JROF))
        Z1ST  =1.0_JPRB/ZCOSLR
        PQX(JROF)=(YDGSGEOM%GSQM2(JROF)*ZCOSLR &
         & +(1.0_JPRB+YDGSGEOM%GEMU(JROF)*PSINLA(JROF))&
         & *PCOSCO(JROF)*Z1ST )*ZDECO  
        PQY(JROF)=-(YDGSGEOM%GEMU(JROF)+PSINLA(JROF))&
         & *(PSINCO(JROF)*Z1ST)*ZDECO  
      ENDIF
    ENDDO

  ELSE

    ! ---- Using MASS library vector routines

    ILEN=KPROF-KSTART+1
    DO JROF=KSTART,KPROF
      ZTMP1(JROF-KSTART+1)=PC2P1-PC2M1*PSINLA(JROF)
    ENDDO

    CALL VREC(ZTMP2,ZTMP1,ILEN)

    DO JROF=KSTART,KPROF
      ZCOSLR=SQRT(PCOSCO(JROF)*PCOSCO(JROF)&
       & +PSINCO(JROF)*PSINCO(JROF))  
      ZCOSLR=MAX(ZEPS,ZCOSLR)
      ZFACT = PDSTRET*ZTMP2(JROF-KSTART+1)
      ZCOSCO= PCOSCO(JROF)*ZFACT
      ZSINCO= PSINCO(JROF)*ZFACT
      ZCOSLA= ZCOSLR*ZFACT
      ZSINLA= (PC2P1*PSINLA(JROF)-PC2M1)*ZTMP2(JROF-KSTART+1)
      ZXX   = YDCSGEOM%RCOLON(JROF)*ZCOSCO-YDCSGEOM%RSILON(JROF)*ZSINCO
      ZYY   = YDCSGEOM%RSILON(JROF)*ZCOSCO+YDCSGEOM%RCOLON(JROF)*ZSINCO
      PLON(JROF)=SIGN(1.0_JPRB,ZYY)
      ZTMP2(JROF-KSTART+1)   =ZSINLA 
      ZTMP7(JROF-KSTART+1)  = ZXX/ZCOSLA
      ! Note: replacing ZXX/ZCOSLA by vdiv not good for adjoint test
      IF(KROT == 1) THEN
        ZTMP4(JROF-KSTART+1) = 1.0_JPRB+PCOPHI(JROF)
        ZTMP5(JROF-KSTART+1) = ZCOSLR
      ENDIF
    ENDDO

    CALL JFH_VIF2(ZTMP2,KPROF-KSTART+1)
    CALL VASIN(ZTMP6,ZTMP2,ILEN)
    CALL JFH_VIF2(ZTMP2,ILEN)
    CALL VASIN(PLAT(KSTART),ZTMP2,ILEN)
    CALL JFH_VIF2(ZTMP7,ILEN)
    CALL VACOS(ZTMP1,ZTMP7,ILEN)

    DO JROF=KSTART,KPROF
      PLON(JROF)=PI+PLON(JROF)*(ZTMP1(JROF-KSTART+1)-PI)
    ENDDO

    CALL JFH_VMOD(PLON(KSTART),PDEPI,KPROF-KSTART+1)

    IF (KROT == 1) THEN

      CALL VREC(ZTMP1,ZTMP4,ILEN)
      CALL VREC(ZTMP2,ZTMP5,ILEN)

      DO JROF=KSTART,KPROF
        ZCOSLR=ZTMP5(JROF-KSTART+1)
        ZDECO =ZTMP1(JROF-KSTART+1)
        Z1ST  =ZTMP2(JROF-KSTART+1)
        PQX(JROF)=(YDGSGEOM%GSQM2(JROF)*ZCOSLR &
         & +(1.0_JPRB+YDGSGEOM%GEMU(JROF)*PSINLA(JROF))&
         & *PCOSCO(JROF)*Z1ST )*ZDECO  
        PQY(JROF)=-(YDGSGEOM%GEMU(JROF)+PSINLA(JROF))&
         & *(PSINCO(JROF)*Z1ST)*ZDECO  

      ENDDO
    ENDIF
  ENDIF !End using MASS library vector routines

ELSEIF(KSTTYP == 2) THEN

  !      POLE STRETCHING AND TILTING.

  IF(N_VMASS > 0) THEN
    CALL ABOR1(' LARCHE2: VMASS lib not yet coded with NSTTYP=2')
  ENDIF

  ZPSLAP = PMUCEN
  ZPCLAP = SQRT(1.0_JPRB-PMUCEN*PMUCEN)
  ZPSLOP = SIN(PLOCEN)
  ZPCLOP = COS(PLOCEN)

  DO JROF=KSTART,KPROF

    ZZY   = YDGSGEOM%GESLO(JROF)*PCOSCO(JROF)+YDGSGEOM%GECLO(JROF)*PSINCO(JROF)
    ZZX   = YDGSGEOM%GECLO(JROF)*PCOSCO(JROF)-YDGSGEOM%GESLO(JROF)*PSINCO(JROF)
    ZCOSL2= ZZY*ZZY+ZZX*ZZX
    ZCOSLR= SQRT(ZCOSL2)
    ZCOSLR= MAX(ZEPS,ZCOSLR)
    ZZZ   = ZZX*ZPCLOP+ZPSLOP*ZZY
    ZZP   = ZPSLAP*PSINLA(JROF)+ZPCLAP*ZZZ
    ZDI= 1.0_JPRB/(PC2P1-PC2M1*ZZP)
    ZYY   = PDSTRET*ZDI*(-ZPCLOP*ZZY+ZPSLOP*ZZX)
    ZXX   = PDSTRET*ZDI*(ZPCLAP*PSINLA(JROF)-ZPSLAP*ZZZ)
    ZCOSLA= SQRT(ZXX*ZXX+ZYY*ZYY)
    ZCOSLA= MAX(ZEPS,ZCOSLA)
    ZSINLA=-(PC2M1-PC2P1*ZZP)*ZDI
    ZMU   = MAX(-1.0_JPRB,MIN(1.0_JPRB,ZSINLA))
    ZZW   = MAX(-1.0_JPRB,MIN(1.0_JPRB,ZXX/ZCOSLA))

    PLAT(JROF)=ASIN(ZMU)
    PLON(JROF)=MOD(PI+SIGN(1.0_JPRB,ZYY)*(ACOS(ZZW)-PI),PDEPI)

    IF (KROT == 1) THEN

      ZDECO =1.0_JPRB/(1.0_JPRB+PCOPHI(JROF))
      Z1ST  =1.0_JPRB/ZCOSLR
      Z1STT =Z1ST/ZCOSLA

      ZP = (YDGSGEOM%GSQM2(JROF)*ZCOSLR &
       & +(1.0_JPRB+YDGSGEOM%GEMU(JROF)*PSINLA(JROF))&
       & *PCOSCO(JROF)*Z1ST )*ZDECO  
      ZQ =-(YDGSGEOM%GEMU(JROF)+PSINLA(JROF))*(PSINCO(JROF)*Z1ST )*ZDECO
      ZCA=PDSTRET*ZDI*(ZPSLAP*ZCOSL2-ZPCLAP*PSINLA(JROF)*ZZZ)*Z1STT
      ZSA=-ZYY*ZPCLAP*Z1STT
      PQX(JROF)=ZP*ZCA+ZQ*ZSA
      PQY(JROF)=ZQ*ZCA-ZP*ZSA

    ENDIF
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARCHE2',1,ZHOOK_HANDLE)

END SUBROUTINE LARCHE2


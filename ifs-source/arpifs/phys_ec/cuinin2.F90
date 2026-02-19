! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUININ2 &
 & (YDECUMF2,YDEPHLI,  KIDIA,    KFDIA,    KLON,    KLEV,&
 & LDRAIN1D,&
 & PTEN,     PQEN,     PQSEN,    PUEN,     PVEN,&
 & PVERVEL,  PGEO,     PAPH,&
 & KLWMIN,   KLAB,&
 & PTENH,    PQENH,    PQSENH,   PGEOH,&
 & PTU,      PQU,      PTD,      PQD,&
 & PUU,      PVU,      PUD,      PVD,&
 & PLU  )  

!          VERSION FOR SIMPLIFIED CONVECTION SCHEME

!          PURPOSE
!          -------

!          THIS ROUTINE INTERPOLATES LARGE-SCALE FIELDS OF T,Q ETC.
!          TO HALF LEVELS (I.E. GRID FOR MASSFLUX SCHEME),
!          DETERMINES LEVEL OF MAXIMUM VERTICAL VELOCITY
!          AND INITIALIZES VALUES FOR UPDRAFTS AND DOWNDRAFTS (1dvar)

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTRN2*.

!          METHOD.
!          --------
!          FOR EXTRAPOLATION TO HALF LEVELS SEE TIEDTKE(1989)

!     PARAMETER     DESCRIPTION                                   UNITS 
!     ---------     -----------                                   ----- 
!     INPUT PARAMETERS (INTEGER): 

!    *KIDIA*        START POINT 
!    *KFDIA*        END POINT 
!    *KLON*         NUMBER OF GRID POINTS PER PACKET 
!    *KLEV*         NUMBER OF LEVELS 

!    INPUT PARAMETERS (REAL): 

!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K 
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG 
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG 
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA

!    OUTPUT PARAMETERS (INTEGER):

!    *KLWMIN*       LEVEL OF MAXIMUM VERTICAL VELOCITY 
!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CONDENSATION LEVEL

!    OUTPUT PARAMETERS (REAL):

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS         K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS    KG/KG
!    *PQSENH*       ENV. SPEC. SATURATION HUMIDITY (T+1)
!                   ON HALF LEVELS                              KG/KG
!    *PTU*          TEMPERATURE IN UPDRAFTS                       K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                  KG/KG
!    *PTD*          TEMPERATURE IN DOWNDRAFTS                     K
!    *PQU*          SPEC. HUMIDITY IN DOWNDRAFTS                KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                       M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                       M/S
!    *PUD*          U-VELOCITY IN DOWNDRAFTS                     M/S
!    *PVD*          V-VELOCITY IN DOWNDRAFTS                     M/S
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG

!          EXTERNALS
!          ---------
!          *CUADJTQ* TO SPECIFY QS AT HALF LEVELS

!     AUTHOR.
!     -------
!      M.TIEDTKE         E.C.M.W.F.     12/89

!     MODIFICATIONS.
!     --------------
!      P. LOPEZ          E.C.M.W.F.     12/2003
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      P. Lopez      11-Jan-2007 Added LDRAIN switch (1D-Var rain)
!      A. Geer       01-Oct-2008 LDRAIN1D name change to reflect usage
!----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RCPD
USE YOECUMF2 , ONLY : TECUMF2
USE YOEPHLI  , ONLY : TEPHLI

IMPLICIT NONE

TYPE(TECUMF2)     ,INTENT(INOUT) :: YDECUMF2
TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
LOGICAL           ,INTENT(IN)    :: LDRAIN1D
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLWMIN(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLAB(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQENH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQSENH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLU(KLON,KLEV) 
REAL(KIND=JPRB) ::     ZWMAX(KLON)
REAL(KIND=JPRB) ::     ZPH(KLON)
LOGICAL ::  LLFLAG(KLON)

INTEGER(KIND=JPIM) :: ICALL, IK, JK, JL

REAL(KIND=JPRB) :: ZZS, ZRCPD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cuadjtq.intfb.h"
#include "cuadjtqs.intfb.h"

!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUININ2',0,ZHOOK_HANDLE)
ASSOCIATE(NJKT22=>YDECUMF2%NJKT22, &
 & LPHYLIN=>YDEPHLI%LPHYLIN)
! Constant (for optimization)
ZRCPD=1.0_JPRB/RCPD

!*    1.           SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
!*                 ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
!*                 FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
!                  ----------------------------------------------

DO JK=2,KLEV
  DO JL=KIDIA,KFDIA
    PTENH(JL,JK)=(MAX(RCPD*PTEN(JL,JK-1)+PGEO(JL,JK-1),&
     & RCPD*PTEN(JL,JK)+PGEO(JL,JK))-PGEOH(JL,JK))*ZRCPD  
    PQENH(JL,JK)=PQEN(JL,JK-1)
    PQSENH(JL,JK)=PQSEN(JL,JK-1)
    ZPH(JL)=PAPH(JL,JK)
    LLFLAG(JL)=.TRUE.
  ENDDO

!orig   IF(JK.GE.KLEV-1) GO TO 130
  IF(JK >= KLEV-1 .OR. JK < NJKT22) CYCLE
  IK=JK
  IF(LPHYLIN .OR. LDRAIN1D) THEN
    ICALL=0
    CALL CUADJTQS &
     & ( KIDIA,    KFDIA,    KLON,    KLEV,     IK,&
     &   ZPH,      PTENH,    PQSENH,  LLFLAG,   ICALL)  
  ELSE
    ICALL=3
    CALL CUADJTQ &
     & ( YDEPHLI,  KIDIA,    KFDIA,    KLON,    KLEV,      IK,&
     &   ZPH,      PTENH,    PQSENH,   LLFLAG,   ICALL)  
  ENDIF

  DO JL=KIDIA,KFDIA
    PQENH(JL,JK)=MIN(PQEN(JL,JK-1),PQSEN(JL,JK-1))&
     & +(PQSENH(JL,JK)-PQSEN(JL,JK-1))  
    PQENH(JL,JK)=MAX(PQENH(JL,JK),0.0_JPRB)
  ENDDO
!orig  130   continue
ENDDO

DO JL=KIDIA,KFDIA
  PTENH(JL,KLEV)=(RCPD*PTEN(JL,KLEV)+PGEO(JL,KLEV)-PGEOH(JL,KLEV))*ZRCPD
  PQENH(JL,KLEV)=PQEN(JL,KLEV)
  PTENH(JL,1)=PTEN(JL,1)
  PQENH(JL,1)=PQEN(JL,1)
  KLWMIN(JL)=KLEV
  ZWMAX(JL)=0.0_JPRB
ENDDO

DO JK=KLEV-1,2,-1
  DO JL=KIDIA,KFDIA
    ZZS=MAX(RCPD*PTENH(JL,JK)+PGEOH(JL,JK),&
     & RCPD*PTENH(JL,JK+1)+PGEOH(JL,JK+1))  
    PTENH(JL,JK)=(ZZS-PGEOH(JL,JK))*ZRCPD
  ENDDO
ENDDO

DO JK=KLEV,3,-1
!DIR$ IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    IF(PVERVEL(JL,JK) < ZWMAX(JL)) THEN
      ZWMAX(JL)=PVERVEL(JL,JK)
      KLWMIN(JL)=JK
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    2.0          INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
!*                 ---------------------------------------------

DO JK=1,KLEV
  IK=JK-1
  IF(JK == 1) IK=1
  DO JL=KIDIA,KFDIA
    PTU(JL,JK)=PTENH(JL,JK)
    PTD(JL,JK)=PTENH(JL,JK)
    PQU(JL,JK)=PQENH(JL,JK)
    PQD(JL,JK)=PQENH(JL,JK)
    PLU(JL,JK)=0.0_JPRB
    PUU(JL,JK)=PUEN(JL,IK)
    PUD(JL,JK)=PUEN(JL,IK)
    PVU(JL,JK)=PVEN(JL,IK)
    PVD(JL,JK)=PVEN(JL,IK)
    KLAB(JL,JK)=0
  ENDDO
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUININ2',1,ZHOOK_HANDLE)
END SUBROUTINE CUININ2

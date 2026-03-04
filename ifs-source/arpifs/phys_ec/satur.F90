! (C) Copyright 1996- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE SATUR ( KIDIA , KFDIA , KLON  , KTDIA , KLEV, LDPHYLIN, &
 & PAPRSF, PT    , PQSAT , KFLAG)

!***

! **   *SATUR* -  COMPUTES SPECIFIC HUMIDITY AT SATURATION

!       J.F. MAHFOUF       E.C.M.W.F.     15/05/96

!       Modified J. HAGUE          13/01/03 MASS Vector Functions

!       PURPOSE.
!       --------

!       SPECIFIC HUMIDITY AT SATURATION IS USED BY THE
!       DIAGNOSTIC CLOUD SCHEME TO COMPUTE RELATIVE HUMIDITY
!       AND LIQUID WATER CONTENT

!       INTERFACE
!       ---------

!       THIS ROUTINE IS CALLED FROM *CALLPAR*.

!       PARAMETER     DESCRIPTION                                 UNITS
!       ---------     -----------                                 -----
!       INPUT PARAMETERS (INTEGER):

!      *KIDIA*        START POINT
!      *KFDIA*        END POINT
!      *KLON*         NUMBER OF GRID POINTS PER PACKET
!      *KTDIA*        START OF THE VERTICAL LOOP
!      *KLEV*         NUMBER OF LEVELS

!       INPUT PARAMETERS (REAL):

!      *PAPRSF*        PRESSURE ON FULL LEVELS                      PA
!      *PT*            TEMPERATURE AT T-DT                          K

!       INPUT PARAMETERS (INTEGER):

!      *KFLAG*         FLAG TO DETECT CALL FROM

!                      CONVECTION  KFLAG=1
!                      OTHER       KFLAG=2

!       OUTPUT PARAMETER (REAL):

!      *PQSAT*         SATURATION SPECIFIC HUMIDITY                 KG/KG

!-------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RETV     ,RLVTT    ,RLSTT    ,RTT
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                    R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &                    RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
 &                    RTWAT_RTICE_R      ,RTWAT_RTICECU_R

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA
LOGICAL           ,INTENT(IN)    :: LDPHYLIN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSAT(KLON,KLEV)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLAG
INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) :: ZCOR, ZEW, ZFOEEW, ZQMAX, ZQS, ZTARG
REAL(KIND=JPRB) :: ZALFA, ZFOEEWL, ZFOEEWI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!DIR$ VFUNCTION EXPHF

#include "fcttre.func.h"

!----------------------------------------------------------------------

!*    1.           DEFINE CONSTANTS
!                  ----------------

IF (LHOOK) CALL DR_HOOK('SATUR',0,ZHOOK_HANDLE)
ZQMAX=0.5_JPRB

!     *
!----------------------------------------------------------------------

!     *    2.           CALCULATE SATURATION SPECIFIC HUMIDITY
!                       --------------------------------------

IF (LDPHYLIN) THEN
  DO JK=KTDIA,KLEV
    DO JL=KIDIA, KFDIA
      ZTARG = PT(JL,JK)
      ZALFA = FOEALFA(ZTARG)

      ZFOEEWL = R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
      ZFOEEWI = R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
      ZFOEEW = ZALFA*ZFOEEWL+(1.0_JPRB-ZALFA)*ZFOEEWI

      ZQS    = ZFOEEW/PAPRSF(JL,JK)
      IF (ZQS > ZQMAX) THEN
        ZQS=ZQMAX
      ENDIF
      ZCOR = 1.0_JPRB/(1.0_JPRB-RETV*ZQS)
      PQSAT(JL,JK)=ZQS*ZCOR
    ENDDO
  ENDDO
ELSE

  DO JK=KTDIA,KLEV
    DO JL=KIDIA, KFDIA
      IF(KFLAG == 1) THEN
        ZEW  = FOEEWMCU(PT(JL,JK))
      ELSE
        ZEW  = FOEEWM(PT(JL,JK))
      ENDIF
      ZQS  = ZEW/PAPRSF(JL,JK)
      ZQS  = MIN(ZQMAX,ZQS)
      ZCOR = 1.0_JPRB/(1.0_JPRB-RETV*ZQS)
      PQSAT(JL,JK)=ZQS*ZCOR
    ENDDO
  ENDDO

ENDIF

IF (LHOOK) CALL DR_HOOK('SATUR',1,ZHOOK_HANDLE)
END SUBROUTINE SATUR

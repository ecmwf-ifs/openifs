! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE VDFEIS (YDECUMF,KIDIA   , KFDIA   , KLON    , KLEV   , LDCLCOMP,&
                   & PT    , PQ      , PAP     , PGEO   , PLCL,   PEIS)

!     ------------------------------------------------------------------

!          PURPOSE.
!          --------
!          TO PRODUCE ESTIMATE ON PBL INVERSION STRENGTH BASED
!          ON WOOD&BRETHERTON (2006)

!          BASED ON SWITCH LDCLCOMP THE LIFTING CONDENSATION LEVEL (LCL) IS
!          EITHER INPUT OR IF LATTER NOT AVAILABLE SIMPLY ESTIMATED (e.g.
!          FROM NEAR SURFACE TEMPERATURE AND DEWPOINT FOLLOWING EPSY OR BOLTON FORMULA

!          AUTHOR
!          ------
!          Martin Koehler         E.C.M.W.F
!          revised: P. Bechtold, K. Lonitz           Aug. 2014
!        R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : YDCST=>YRCST ! allows use of included functions. REK.
USE YOETHF   , ONLY : YDTHF=>YRTHF ! allows use of included functions. REK.
USE YOECUMF  , ONLY : TECUMF

IMPLICIT NONE

!*         0.1    GLOBAL VARIABLES

TYPE(TECUMF)      ,INTENT(IN)    :: YDECUMF
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
LOGICAL,           INTENT(IN)    :: LDCLCOMP              !use input LCL(.False.)
                                                         !or estimate LCL(.true.)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)        !temperature
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)        !specific humidity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)       !full level pressure
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV)      !geopot full levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLCL(KLON)           !lifting condensation level (m)
                                                         !or boundary layer-height
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEIS(KLON)           !inversion strength (K)

REAL(KIND=JPRB)    :: ZLCL(KLON), ZSTABIL(KLON)
REAL(KIND=JPRB)    :: ZQS, ZTD, ZT850, ZGAMMA850, ZRG, ZRDOCP, ZA
INTEGER(KIND=JPIM) :: JK, JL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

#include "fcttre.func.h"

IF (LHOOK) CALL DR_HOOK('VDFEIS',0,ZHOOK_HANDLE)
ASSOCIATE(NJKT4=>YDECUMF%NJKT4, NJKT6=>YDECUMF%NJKT6, &
 & RG=>YDCST%RG, RCPD=>YDCST%RCPD, RETV=>YDCST%RETV, RLVTT=>YDCST%RLVTT, RLSTT=>YDCST%RLSTT, RTT=>YDCST%RTT, &
 & RV=>YDCST%RV, RD=>YDCST%RD, &
 & R2ES=>YDTHF%R2ES, R3LES=>YDTHF%R3LES, R3IES=>YDTHF%R3IES, R4LES=>YDTHF%R4LES, &
 & R4IES=>YDTHF%R4IES, R5LES=>YDTHF%R5LES, R5IES=>YDTHF%R5IES, R5ALVCP=>YDTHF%R5ALVCP, &
 & R5ALSCP=>YDTHF%R5ALSCP, RALVDCP=>YDTHF%RALVDCP, RALSDCP=>YDTHF%RALSDCP, &
 & RTWAT=>YDTHF%RTWAT, RTICE=>YDTHF%RTICE, RTICECU=>YDTHF%RTICECU, RTWAT_RTICE_R=>YDTHF%RTWAT_RTICE_R, &
 & RTWAT_RTICECU_R=>YDTHF%RTWAT_RTICECU_R)

! optimization
ZRG    = 1.0_JPRB/RG
ZRDOCP = RD/RCPD

DO JL=KIDIA,KFDIA
  PEIS(JL) = -999.0_JPRB
ENDDO

IF (LDCLCOMP) THEN
   !compute dew point temperature and LCL
   DO JL=KIDIA,KFDIA
      ZA = LOG(MAX(1.E-4_JPRB,PQ(JL,KLEV))*PAP(JL,KLEV)/380.04_JPRB)
      ZTD = (R4LES*ZA-R3LES*RTT)/(ZA-R3LES)
      ZLCL(JL) = 125.0_JPRB*(PT(JL,KLEV)-ZTD) ! Epsy equation
   ENDDO
ELSE
   DO JL=KIDIA,KFDIA
      ZLCL(JL) = PLCL(JL)
   ENDDO
ENDIF

!*    Estimated Inversion Strength (EIS) criteria from Wood & Bretherton (2006)

DO JL=KIDIA,KFDIA
    ZSTABIL(JL)    = PT(JL,NJKT6) * ( 1.0E5_JPRB/PAP(JL,NJKT6) ) ** ZRDOCP &
             & - PT(JL,KLEV)  * ( 1.0E5_JPRB/PAP(JL,KLEV) )  ** ZRDOCP  
    ZT850      = 0.5 *( PT(JL,KLEV) + PT(JL,NJKT6) )
    JK         = NJKT4
!     qsat (full level)
    ZQS = FOEEWM(ZT850)/PAP(JL,JK)
    ZQS = MIN(0.5_JPRB,ZQS)
    ZQS = ZQS/(1.0_JPRB-RETV*ZQS)

    ZGAMMA850  = RG/RCPD* (1 - ( 1 + RLVTT   *ZQS / (     RD*ZT850   ) ) &
                           & / ( 1 + RLVTT**2*ZQS / (RCPD*RV*ZT850**2) ) )
    PEIS(JL)   = ZSTABIL(JL) - ZGAMMA850 * ( PGEO(JL,NJKT6)*ZRG - ZLCL(JL) )
      
ENDDO
  
!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VDFEIS',1,ZHOOK_HANDLE)
END SUBROUTINE VDFEIS

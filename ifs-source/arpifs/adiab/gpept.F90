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

SUBROUTINE GPEPT(YDPHY,KPROMA,KSTART,KPROF,KFLEV,PTETA,PTT0,PRPRESF,PTETAE)

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RV       ,RCPD     ,RCPV     ,RETV     ,&
 &                    RCW      ,RCS      ,RLVTT    ,RLSTT    ,RTT      ,&
 &                    RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
 &                    RGAMS    ,RALPD    ,RBETD    ,RGAMD  
USE YOMPHY   , ONLY : TPHY

!**** *GPEPT*    Compute the equivalent potential temperature on model levels

!      PURPOSE.
!      --------
!           Computes the equivalent potential temperature on model levels

!**    INTERFACE.
!      ----------
!           *CALL* GPEPT( ... )

!           EXPLICITE ARGUMENTS.
!           --------------------
!           KPROMA      : horizontal dimension.                     (INPUT)
!           KSTART      : start of work.                            (INPUT)
!           KPROF       : depth of work.                            (INPUT)
!           KFLEV       : number of input pressure levels.          (INPUT)
!           PTETA       : teta                                      (INPUT)
!           PTT0        : model leval temperature                   (INPUT)
!           PRPRESF     : model leval pressure                      (INPUT)
!           PTETAE      : equivalent pot. temperature on model levels (OUTPUT)

!      EXTERNALS.
!      ----------
!          thermodynamic inline functions

!      AUTHOR.
!      -------
!           RYAD EL KHATIB and C. PERIARD  *METEO-FRANCE* 

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 94-04-08
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     -----------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTETA(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTETAE(KPROMA,KFLEV) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZCPS, ZDELTA, ZES, ZESP, ZLH, ZQS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "fcttrm.func.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPEPT',0,ZHOOK_HANDLE)
ASSOCIATE(LNEIGE=>YDPHY%LNEIGE)
!     ------------------------------------------------------------------

!*       1.COMPUTE POTENTIAL EQUIVALENT TEMPERATURE ON MODEL LEVELS
!          --------------------------------------------------------

DO JLEV = 1, KFLEV
  DO JROF = KSTART, KPROF

    IF (LNEIGE) THEN
      ZDELTA = MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT - PTT0(JROF,JLEV) ))
    ELSE
      ZDELTA = 0.0_JPRB
    ENDIF
    ZES = FOEW( PTT0(JROF,JLEV),ZDELTA )
    ZESP = ZES/PRPRESF(JROF,JLEV)
    ZQS = FOQS( ZESP )
    ZLH = FOLH( PTT0(JROF,JLEV), ZDELTA )
    ZCPS = RCPD*(1.0_JPRB-ZQS) + RCPV*ZQS
    PTETAE(JROF,JLEV) = PTETA(JROF,JLEV)*&
     & EXP( ZLH*ZQS/( ZCPS*PTT0(JROF,JLEV) ) )  

  ENDDO
ENDDO

!     -----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPEPT',1,ZHOOK_HANDLE)
END SUBROUTINE GPEPT

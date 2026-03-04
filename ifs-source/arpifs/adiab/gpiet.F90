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

SUBROUTINE GPIET(YDPHY,KPROMA,KSTART,KPROF,KFLEV,PQT0,PTT0,PCP,PRPRESF,PTETAE)

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY :  RD       ,RV       ,RCPD     ,RCPV     ,RETV     ,&
 &                     RCW      ,RCS      ,RLVTT    ,RLSTT    ,RTT      ,RATM     ,&
 &                     RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
 &                     RGAMS    ,RALPD    ,RBETD    ,RGAMD  
USE YOMPHY   , ONLY : TPHY

!**** *GPIET*    Compute the isobaric equivalent temperature on model levels

!      PURPOSE.
!      --------
!           Computes the isobaric equivalent temperature on model levels

!**    INTERFACE.
!      ----------
!           *CALL* GPIET( ... )

!           EXPLICITE ARGUMENTS.
!           --------------------
!           KPROMA      : horizontal dimension.                     (INPUT)
!           KSTART      : start of work.                            (INPUT)
!           KPROF       : depth of work.                            (INPUT)
!           KFLEV       : number of input pressure levels.          (INPUT)
!           PQT0        : model level specific humidity             (INPUT)
!           PTT0        : model level temperature                   (INPUT)
!           PCP         : model level Cp                            (INPUT)
!           PRPRESF     : model level pressure                      (INPUT)
!           PTETAE      : izo. equiv. pot. temp. on model levels    (OUTPUT)

!      EXTERNALS.
!      ----------
!          thermodynamic inline functions

!      AUTHOR.
!      -------
!           FILIP VANA  * LACE *

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 02-03-19
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     -----------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCP(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTETAE(KPROMA,KFLEV) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZDELTA, ZLH
REAL(KIND=JPRB) :: ZUSRATM, ZEXP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "fcttrm.func.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPIET',0,ZHOOK_HANDLE)
ASSOCIATE(LNEIGE=>YDPHY%LNEIGE)
!     ------------------------------------------------------------------

!*       1.COMPUTE ISOBARIC EQUIVALENT POT TEMPERATURE ON MODEL LEVELS
!          -----------------------------------------------------------

ZUSRATM=1.0_JPRB/RATM
ZEXP=-RD/RCPD

DO JLEV = 1, KFLEV
  DO JROF = KSTART, KPROF

    IF (LNEIGE) THEN
      ZDELTA = MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT - PTT0(JROF,JLEV) ))
    ELSE
      ZDELTA = 0.0_JPRB
    ENDIF
    ZLH = FOLH( PTT0(JROF,JLEV), ZDELTA )

    PTETAE(JROF,JLEV) = (PTT0(JROF,JLEV) +&
     & (ZLH*PQT0(JROF,JLEV))/PCP(JROF,JLEV) )*&
     & (PRPRESF(JROF,JLEV)*ZUSRATM)**(ZEXP)  

  ENDDO
ENDDO

!     -----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPIET',1,ZHOOK_HANDLE)
END SUBROUTINE GPIET

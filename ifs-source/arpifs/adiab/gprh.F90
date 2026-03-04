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

SUBROUTINE GPRH(LDWMORH,KPROMA,KSTART,KPROF,KFLEV,&
 & PRHMAX,PRHMIN,PQ,PT,PRESF,PES,PRH,&  
 ! --- OPTIONAL INPUT --------------------------------------------------------
 & LDSONNTAG)

!**** *GPRH* - COMPUTES ES AND RH FROM T AND Q

!     PURPOSE.  COMPUTES SATURATION VAPOUR PRESSURE AND RELATIVE
!     --------  HUMIDITY FROM TEMPERATURE AND SPECIFIC HUMIDITY.

!**   INTERFACE.
!     ----------
!        *CALL* *GPRH(..)*

!        EXPLICIT ARGUMENTS
!        --------------------

!        LDWMORH              - WMO CONSIDER WATER PHASE ONLY          (INPUT)
!        KPROMA               - HORIZONTAL DIMENSIONS.                 (INPUT)
!        KSTART               - START OF WORK                          (INPUT)
!        KPROF                - DEPTH OF WORK                          (INPUT)
!        KFLEV                - NUMBER OF MODEL LEVELS                 (INPUT)
!        PRHMAX               - MAXIMUM RELATIVE HUMIDITY              (INPUT)
!        PRHMIN               - MINIMUM RELATIVE HUMIDITY              (INPUT)
!        PQ(KPROMA,KFLEV)     - SPECIFIC HUMIDITY ON FULL MODEL LEVELS (INPUT)
!        PT(KPROMA,KFLEV)     - TEMPERATURE ON FULL MODEL LEVELS       (INPUT)
!        PRESF(KPROMA,KFLEV)  - MODEL FULL LEVEL PRESSURES             (INPUT)

!        PES(KPROMA,KFLEV)    - SATURATION PRESSURE                    (OUTPUT)
!        PRH(KPROMA,KFLEV)    - RELATIVE HUMIDITY                      (OUTPUT)

!   * OPTIONAL INPUT:
!        LDSONNTAG - Use Sonntag(1994) equation (for comparison with sondes)

!        IMPLICIT ARGUMENTS :   PHYSICAL CONSTANTS FROM COMMONS YOMCST
!        --------------------   AND YOMPHY, FUNCTIONS FROM FCTTRM.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.     NONE.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*
!      ORIGINAL : 88-02-04

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RV       ,RCPV     ,RETV     ,RCW      ,&
 &                    RCS      ,RLVTT    ,RLSTT    ,RTT      ,RALPW    ,&
 &                    RBETW    ,RGAMW    ,RALPS    ,RBETS    ,RGAMS    ,&
 &                    RALPD    ,RBETD    ,RGAMD  
USE YOMPHY   , ONLY : YRPHY
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                    R4IES    ,R5LES    ,R5IES    ,R5ALVCP ,R5ALSCP  ,&
 &                    RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE   ,RTICECU,&
 &                    RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
USE YOMCT0   , ONLY : LECMWF

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
LOGICAL           ,INTENT(IN)    :: LDWMORH 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHMAX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHMIN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PES(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRH(KPROMA,KFLEV) 
LOGICAL,OPTIONAL  ,INTENT(IN)    :: LDSONNTAG

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL, JLEV
LOGICAL :: LLSONNTAG

REAL(KIND=JPRB) :: ZDELTA, ZRH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "fcttre.func.h"
#include "fcttrm.func.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPRH',0,ZHOOK_HANDLE)
ASSOCIATE(LNEIGE=>YRPHY%LNEIGE)
IF (PRESENT(LDSONNTAG)) THEN
  LLSONNTAG=LDSONNTAG
ELSE
  LLSONNTAG=.FALSE.
ENDIF
!     ------------------------------------------------------------------

!*       1.    COMPUTES ES AND RH.
!              -------------------

DO JLEV=1,KFLEV
  IF(LECMWF) THEN
    IF(LLSONNTAG) THEN
      DO JL=KSTART,KPROF
        PES(JL,JLEV)=FOELSON(PT(JL,JLEV))
      ENDDO
    ELSEIF(LDWMORH) THEN
      DO JL=KSTART,KPROF
        PES(JL,JLEV)=(RETV+1.0_JPRB)*FOEEWMO(PT(JL,JLEV))
      ENDDO
    ELSE
      DO JL=KSTART,KPROF
        PES(JL,JLEV)=(RETV+1.0_JPRB)*FOEEWM(PT(JL,JLEV))
      ENDDO
    ENDIF
  ELSE
    IF(LNEIGE) THEN
      DO JL=KSTART,KPROF
        ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-PT(JL,JLEV)))
        PES(JL,JLEV)=FOEW(PT(JL,JLEV),ZDELTA)
      ENDDO
    ELSE
      DO JL=KSTART,KPROF
        ZDELTA=0.0_JPRB
        PES(JL,JLEV)=FOEW(PT(JL,JLEV),ZDELTA)
      ENDDO
    ENDIF
  ENDIF

  DO JL=KSTART,KPROF
    ZRH=(PRESF(JL,JLEV)*PQ(JL,JLEV)*(RETV+1.0_JPRB))&
     & /((1.0_JPRB+RETV*PQ(JL,JLEV))*PES(JL,JLEV))  
    PRH(JL,JLEV)=MAX(PRHMIN,MIN(ZRH,PRHMAX))
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPRH',1,ZHOOK_HANDLE)
END SUBROUTINE GPRH

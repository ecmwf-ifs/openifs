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

SUBROUTINE PPUV(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,KLEVB,&
 & LDBELO,LDBLOW,PRPRES,PRXP,PRXPD,&
 & PUF,PVF,PUPP,PVPP)  

!**** *PPUV* - POST-PROCESS WINDS

!     PURPOSE.
!     --------
!           PERFORMS THE VERTICAL INTERPOLATION OF U AND V TO GIVEN
!       PRESSURE LEVELS

!**   INTERFACE.
!     ----------
!        *CALL* *PPUV(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KLEVB,
!                     LDBELO,LDBLOW,PRPRES,PRXP,PRXPD,
!                     PUF,PVF,PUPP,PVPP)

!        EXPLICIT ARGUMENTS
!        --------------------

!        KPROMA                    - HORIZONTAL DIMENSION.             (INPUT-C)
!        KSTART                    - START OF WORK.                    (INPUT-C)
!        KPROF                     - DEPTH OF WORK.                    (INPUT-C)
!        KFLEV                     - NUMBER OF INPUT PRESSURE LEVELS   (INPUT-C)
!        KLEVP                     - NUMBER OF OUTPUT PRESSURE LEVELS  (INPUT-C)
!        KLOLEV                    - BEGINING FOR THE INTERPOLATION    (INPUT-C)
!        KLEVB(KPROMA,KLEVP,KPPM)  - INPUT LEVEL BELOW PRPRES          (INPUT-C)
!                                    (SEE PPFLEV)
!        KPPM     - Number of interpolation methods in post-processing (INPUT)

!        PRPRES(KPROMA,KLEVP)      - LOG(POST-PROCESSING LEVEL PRES.)  (INPUT-C)
!        LDBELO(KPROMA,KLEVP)      - .TRUE. IF PRESSURE IS UNDER
!                                     LOWEST  MODEL LEVEL              (INPUT-C)
!        LDBLOW(KLEVP)             - .TRUE IF LDBELO(J) IS CONTAINING
!                                    AT LEAST ONE .TRUE.               (INPUT-C)
!        PRXP(KPROMA,0:KFLEV,KPPM) - HALF,FULL AND LN HALF,FULL LEVEL
!                                    PRESSURES (SEE PPINIT)            (INPUT)
!        PRXPD(KPROMA,0:KFLEV,KPPM)- 1./D(P) AND 1./D(LN(P))           (INPUT)
!        PUF(KPROMA,0:KFLEV)       - ZONAL WIND ON FULL INPUT LEVELS   (INPUT)
!        PVF(KPROMA,0:KFLEV)       - MERID. WIND ON FULL INPUT LEVELS  (INPUT)

!        PUPP(KPROMA,KLEVP)        - POST-PROCESSED ZONAL WIND         (OUTPUT)
!        PVPP(KPROMA,KLEVP)        - POST-PROCESSED MERIDIONAL WIND    (OUTPUT)

!        IMPLICIT ARGUMENTS :  NONE.
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.  PPINTP - LINEAR INTERPOLATION
!     ----------  PPITPQ - QUADRATIC INTERPOLATION FOR TOP LEVELS

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-01-26
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K.Yessad (may 2009): merge PPUV_OLD and PPUV
!        A.Geer        24-Jul-2015 Pre-OOPS: 0:NFLEVG deprecated
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPPVI  , ONLY : LRPPUV_CSTEXT,LRPPUV_CALLITPQ

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPPM
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVB(KPROMA,KLEVP,KPPM) 
LOGICAL           ,INTENT(IN)    :: LDBELO(KPROMA,KLEVP) 
LOGICAL           ,INTENT(IN)    :: LDBLOW(KLEVP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRES(KPROMA,KLEVP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXP(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXPD(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUPP(KPROMA,KLEVP) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVPP(KPROMA,KLEVP) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ISLCT, JL, JLEVP
REAL(KIND=JPRB) :: ZUF(KPROMA,0:KFLEV), ZVF(KPROMA,0:KFLEV)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "ppintp.intfb.h"
#include "ppitpq.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PPUV',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    POST-PROCESS WINDS.
!              -------------------

ISLCT=4

!*       1.1   PREPARE FOR EXTRAPOLATION ABOVE TOP MODEL LEVEL.
ZUF(KSTART:KPROF,1:KFLEV)=PUF(KSTART:KPROF,1:KFLEV)
ZVF(KSTART:KPROF,1:KFLEV)=PVF(KSTART:KPROF,1:KFLEV)
IF(LRPPUV_CSTEXT .OR. KFLEV==1) THEN
  DO JL=KSTART,KPROF
    ZUF(JL,0)=PUF(JL,1)
    ZVF(JL,0)=PVF(JL,1)
  ENDDO
ELSE
  DO JL=KSTART,KPROF
    ZUF(JL,0)=PUF(JL,2)+PRXP(JL,2,ISLCT)*PRXPD(JL,1,ISLCT)*(PUF(JL,1)-PUF(JL,2))
    ZVF(JL,0)=PVF(JL,2)+PRXP(JL,2,ISLCT)*PRXPD(JL,1,ISLCT)*(PVF(JL,1)-PVF(JL,2))
  ENDDO
ENDIF

!*       1.2   INTERPOLATE WINDS.

CALL PPINTP(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,KLEVB,ISLCT,&
 & LDBELO,LDBLOW,PRPRES,PRXP,PRXPD,ZUF,PUPP)  
CALL PPINTP(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,KLEVB,ISLCT,&
 & LDBELO,LDBLOW,PRPRES,PRXP,PRXPD,ZVF,PVPP)  

!*       1.3   CONSTANT EXTRAPOLATION BELOW LOWEST MODEL LEVEL.

DO JLEVP=KLOLEV,KLEVP
  IF (LDBLOW(JLEVP)) THEN
    DO JL=KSTART,KPROF
      IF (LDBELO(JL,JLEVP)) THEN
        PUPP(JL,JLEVP)=ZUF(JL,KFLEV)
        PVPP(JL,JLEVP)=ZVF(JL,KFLEV)
      ENDIF
    ENDDO
  ENDIF
ENDDO

!*       1.4   MODIFY INTERPOLATION/EXTRAPOLATION ABOVE SECOND
!*          FULL MODEL LEVEL. FIT QUADRATIC POLYNOMIAL IN LN(P).

IF(LRPPUV_CALLITPQ) THEN
  CALL PPITPQ(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,ISLCT,&
   & PRPRES,PRXP,ZUF,PUPP)
  CALL PPITPQ(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,ISLCT,&
   & PRPRES,PRXP,ZVF,PVPP)
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPUV',1,ZHOOK_HANDLE)
END SUBROUTINE PPUV

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

SUBROUTINE PPT_OLD(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,KLEVB,&
 & PRES,LDBELO,LDBELS,LDBLOW,LDBLES,LDEXTR,&
 & POROG,PRXP,PRXPD,PTSTAR,PT0,PTF,PTPP)  

!**** *PPT* - POST-PROCESS TEMPERATURE.

!     PURPOSE.
!     --------
!           PERFORMS THE VERTICAL INTERPOLATION OF TEMPERATURETO A GIVEN
!       MODEL CO-ORDINATE LEVEL. ALSO EXTRAPOLATES TEMPERATURE
!       BELOW SURFACE AND BOTTOM MODEL LEVEL IF SO REQUESTED.

!**   INTERFACE.
!     ----------
!        EXPLICIT ARGUMENTS
!        --------------------

!        KPROMA                    - HORIZONTAL DIMENSION.             (INPUT-C)
!        KSTART                    - START OF WORK.                    (INPUT-C)
!        KPROF                     - DEPTH OF WORK.                    (INPUT-C)
!        KFLEV                     - NUMBER OF INPUT PRESSURE LEVELS   (INPUT-C)
!        KLEVP                     - NUMBER OF OUTPUT PRESSURE LEVELS  (INPUT-C)
!        KLOLEV                    - BEGINING FOR THE INTERPOLATION    (INPUT-C)
!        KPPM     - Number of interpolation methods in post-processing (INPUT-C)

!        KLEVB(KPROMA,KLEVP,KPPM)  - INPUT LEVEL BELOW PRES          (INPUT-C)
!                                    (SEE PPFLEV)
!        PRES(KPROMA,KFLEV)       - POST-PROCESS   LEVEL PRESSURES.   (INPUT-C)
!        PLNPRT(KPROMA,KLEVP)      - PLNPRT(.,JLEVP)=
!                                     LOG(PRES(.,JLEVP)/PSOL)        (INPUT-C)
!        POROG(KPROMA)             - MODEL OROGRAPHY.                  (INPUT-C)
!        LDBELO(KPROMA,KLEVP)      - .TRUE. IF PRESSURE IS UNDER
!                                     LOWEST (FULL) MODEL LEVEL        (INPUT-C)
!        LDBELS(KPROMA,KLEVP)      - .TRUE. IF PRESSURE IS UNDER
!                                     MODEL SURFACE                    (INPUT-C)
!        LDBLOW                    - .TRUE. IF LDBELO(J) IS CONTAINING
!                                    AT LEAST ONE .TRUE.               (INPUT-C)
!        LDBLES                    - .TRUE. IF LDBELS(J) IS CONTAINING
!                                    AT LEAST ONE .TRUE.               (INPUT-C)
!        LDEXTR                    - .TRUE. IF EXTRAPOLATION REQUESTED (INPUT-C)
!        PRXP(KPROMA,0:KFLEV,KPPM) - HALF,FULL AND LN HALF,FULL LEVEL
!                                    PRESSURES (SEE PPINIT)            (INPUT)
!        PRXPD(KPROMA,0:KFLEV,KPPM)- 1./D(P) AND 1./D(LN(P))           (INPUT)
!        PTSTAR(KPROMA)            - SURFACE TEMPERATURE               (INPUT)
!        PT0(KPROMA)               - STANDARD SURFACE TEMPERATURE      (INPUT)
!        PTF(KPROMA,0:KFLEV)       - TEMPERATURE ON FULL INPUT LEVELS. (INPUT)

!        PTPP(KPROMA,KLEVP)        - POST-PROCESSED TEMPERATURE        (OUTPUT)

!        IMPLICIT ARGUMENTS :  CONSTANTS FROM YOMCST AND YOMSTA.
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.  PPINTP - LINEAR INTERPOLATION
!     ----------

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
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RD
USE YOMSTA   , ONLY : RDTDZ1

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPPM
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVB(KPROMA,KLEVP,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES(KPROMA,KLEVP) 
LOGICAL           ,INTENT(IN)    :: LDBELO(KPROMA,KLEVP) 
LOGICAL           ,INTENT(IN)    :: LDBELS(KPROMA,KLEVP) 
LOGICAL           ,INTENT(IN)    :: LDBLOW(KLEVP) 
LOGICAL           ,INTENT(IN)    :: LDBLES(KLEVP) 
LOGICAL           ,INTENT(IN)    :: LDEXTR 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXP(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXPD(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTAR(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTF(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTPP(KPROMA,KLEVP) 

REAL(KIND=JPRB) :: ZALPHA(KPROMA)

INTEGER(KIND=JPIM) :: ISLCT, JL, JLEVP

REAL(KIND=JPRB) :: ZALFLP2, ZALFLP3, ZALFLPR, ZALPHAC, ZCOEF,&
 & ZCOEFPL, ZDFI, ZFI2000, ZFI2500, ZLNPRT, &
 & ZTPLAT, ZTX, ZUSDFI  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "ppintp.intfb.h"

!     ------------------------------------------------------------------

!*       1.    POST-PROCESS TEMPERATURE.
!              -------------------------

!*       1.1   PREPARE FOR EXTRAPOLATION ABOVE TOP

IF (LHOOK) CALL DR_HOOK('PPT_OLD',0,ZHOOK_HANDLE)
DO JL=KSTART,KPROF
  PTF(JL,0)=PTF(JL,1)
ENDDO

!*       1.2   INTERPOLATE TEMPERATURE.

ISLCT=2
CALL PPINTP(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,&
 & KLEVB,ISLCT,&
 & LDBELO,LDBLOW,PRES,PRXP,PRXPD,&
 & PTF,PTPP)  

!*       1.3   INTERPOLATE BETWEEN LOWEST LEVEL AND SURFACE.

DO JLEVP=KLOLEV,KLEVP
  IF (LDBLOW(JLEVP)) THEN
    DO JL=KSTART,KPROF
      IF (LDBELO(JL,JLEVP).AND..NOT.LDBELS(JL,JLEVP) ) THEN
        ZCOEF=(PRXP(JL,KFLEV,1)-PRES(JL,JLEVP) )/&
         & (PRXP(JL,KFLEV,1)-PRXP(JL,KFLEV,2))  
        PTPP(JL,JLEVP)=PTSTAR(JL)+ZCOEF*(PTF(JL,KFLEV)-PTSTAR(JL))
      ENDIF
    ENDDO
  ENDIF
ENDDO

!*       1.4   EXTRAPOLATE TEMPERATURE BELOW SURFACE.

IF(LDEXTR) THEN

  ZTX=298.0_JPRB
  ZFI2000=2000.0_JPRB*RG
  ZFI2500=2500.0_JPRB*RG
  ZDFI   =ZFI2500-ZFI2000
  ZUSDFI =1.0_JPRB/ZDFI
  ZALPHAC=-RDTDZ1*RD/RG
  DO JL=KSTART,KPROF
    ZTPLAT =MIN(PT0(JL),ZTX)
    ZCOEFPL=MIN(1.0_JPRB,MAX(0.0_JPRB,(POROG(JL)-ZFI2000)*ZUSDFI))
    ZALPHA(JL)=MAX( 0.0_JPRB , ZALPHAC+&
     & (RD*ZCOEFPL*(ZTPLAT-PT0(JL)))/MAX(POROG(JL),ZFI2000) )  
  ENDDO

  DO JLEVP=KLOLEV,KLEVP
    IF (LDBLES(JLEVP)) THEN
      DO JL=KSTART,KPROF
        IF (LDBELS(JL,JLEVP)) THEN
          ZLNPRT=LOG(PRES(JL,JLEVP)/PRXP(JL,KFLEV,1))
          ZALFLPR=ZALPHA(JL)*ZLNPRT
          ZALFLP2=ZALFLPR*ZALFLPR
          ZALFLP3=ZALFLP2*ZALFLPR
          PTPP(JL,JLEVP)=&
           & PTSTAR(JL)*(1.0_JPRB+ZALFLPR+0.5_JPRB*ZALFLP2+1.0_JPRB/6._JPRB*ZALFLP3)  
        ENDIF
      ENDDO
    ENDIF
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPT_OLD',1,ZHOOK_HANDLE)
END SUBROUTINE PPT_OLD

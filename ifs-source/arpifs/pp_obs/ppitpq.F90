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

SUBROUTINE PPITPQ(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,KSL,&
 & PRPRES,PRXP,PFLD,&
 & PZPP)  

!**** *PPITPQ* - POST-PROCESS TOP LEVELS QUADRATICALLY

!     PURPOSE.
!     --------
!           FITS A QUADRATIC POLYNOMIAL TO THE TO THREE LEVELS OF THE
!       MODEL. POST-PROCESSED VALUES ABOVE THE SECOND LEVEL ARE
!       OWER WRITTEN WITH NEW QUADRATICALLY INTERPOLATED VALUES.

!**   INTERFACE.
!     ----------
!        *CALL* *PPITPQ(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KSL
!                       PRPRES,PRXP,PFLD,
!                       PZPP)

!        EXPLICIT ARGUMENTS
!        --------------------

!        KPROMA                    - HORIZONTAL DIMENSION.             (INPUT-C)
!        KSTART                    - START OF WORK.                    (INPUT-C)
!        KPROF                     - DEPTH OF WORK.                    (INPUT-C)
!        KFLEV                     - NUMBER OF INPUT PRESSURE LEVELS   (INPUT-C)
!        KLEVP                     - NUMBER OF OUTPUT PRESSURE LEVELS  (INPUT-C)
!        KLOLEV                    - BEGINING FOR THE INTERPOLATION    (INPUT-C)
!        KPPM     - Number of interpolation methods in post-processing (INPUT-C)

!        KSL                       - 1:  FIELD GIVEN ON HALF LEVELS.
!                                        INTERPOLATE IN P
!                                  - 2:  FIELD GIVEN ON FULL LEVELS.
!                                        INTERPOLATE IN P
!                                  - 3:  FIELD GIVEN ON HALF LEVELS.
!                                        INTERPOLATE IN LN(P)
!                                  - 4:  FIELD GIVEN ON FULL LEVELS.
!                                        INTERPOLATE IN LN(P)          (INPUT-C)
!        PRPRES(KPROMA,KLEVP)      -  LOG(POST-PROCESSING LEVEL PRES.) (INPUT-C)
!        PRXP(KPROMA,0:KFLEV,KPPM) - HALF,FULL AND LN HALF,FULL LEVEL
!                                    PRESSURES (SEE PPINIT)            (INPUT)
!        PFLD(KPROMA,0:KFLEV)      - MODEL LEVEL FIELD VALUES          (INPUT)
!        PZPP(KPROMA,KLEVP)        - POST-PROCESSED GEOPOTENTIAL       (OUTPUT)

!        IMPLICIT ARGUMENTS :  CONSTANTS FROM YOMCST
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.  NONE
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        ERIK ANDERSSON  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 92-04-07
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPPM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRES(KPROMA,KLEVP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXP(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFLD(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZPP(KPROMA,KLEVP) 
INTEGER(KIND=JPIM) :: I1, I2, I3, JL, JLEV

REAL(KIND=JPRB) :: ZDEN1, ZDEN2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    POST-PROCESS TOP LEVELS
!              -----------------------

!*       1.0   SET POINTERS TO TOP LEVELS

IF (LHOOK) CALL DR_HOOK('PPITPQ',0,ZHOOK_HANDLE)
IF(KSL == 1.OR. KSL == 3) THEN
  I1=0
  I2=1
  I3=2
ELSE
  I1=1
  I2=2
  I3=3
ENDIF

!*       1.2   MODIFY INTERPOLATION/EXTRAPOLATION ABOVE SECOND
!*          LEVEL (JLEV=I2). FIT QUADRATIC POLYNOMIAL IN LN(P) TO THE
!*          GEOPOTENTIAL AT THE TOP THREE LEVELS (JLEV=I1,I2,I3) AND
!*          APPLY TO PRESSURES ABOVE THE SECOND LEVEL FROM THE TOP
!*          (JLEV=I2).

DO JLEV=KLOLEV,KLEVP
  DO JL=KSTART,KPROF
    IF(PRPRES(JL,JLEV) < PRXP(JL,I2,KSL)) THEN
      ZDEN1=(PRXP(JL,I2,KSL)-PRXP(JL,I1,KSL))&
       & *(PRXP(JL,I1,KSL)-PRXP(JL,I3,KSL))  
      ZDEN2=(PRXP(JL,I2,KSL)-PRXP(JL,I3,KSL))&
       & *(PRXP(JL,I1,KSL)-PRXP(JL,I3,KSL))  
      PZPP(JL,JLEV)=PFLD(JL,I2)&
       & +(PFLD(JL,I2)-PFLD(JL,I1))&
       & *(PRPRES(JL,JLEV)-PRXP(JL,I2,KSL))&
       & *(PRPRES(JL,JLEV)-PRXP(JL,I3,KSL))&
       & /ZDEN1 &
       & -(PFLD(JL,I2)-PFLD(JL,I3))&
       & *(PRPRES(JL,JLEV)-PRXP(JL,I1,KSL))&
       & *(PRPRES(JL,JLEV)-PRXP(JL,I2,KSL))&
       & /ZDEN2  
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPITPQ',1,ZHOOK_HANDLE)
END SUBROUTINE PPITPQ

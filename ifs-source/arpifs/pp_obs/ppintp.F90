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

SUBROUTINE PPINTP(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,KPPM,&
 & KLEVB,KSLCT,&
 & LDBELO,LDBLOW,PRPRESC,PRXP,PRXPD,&
 & PFLDI,PFLDO)  

!**** *PPINTP* - LINEARLY INTERPOLATE FIELDS TO GIVEN LEVEL

!     PURPOSE.
!     --------
!           PERFORMS THE VERTICAL LINEAR INTERPOLATION OF A FIELD
!           TO A GIVEN LEVEL.
!           IMPLICITLY ALSO EXTRAPOLATES ABOVE THE TOP MODEL LEVEL.

!**   INTERFACE.
!     ----------
!        *CALL* *PPINTP(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KLOLEV,
!                       KLEVB,KSLCT,
!                       LDBELO,LDBLOW,PRPRESC,PRXP,PRXPD,
!                       PFLDI,PFLDO)

!        EXPLICIT ARGUMENTS :
!        --------------------

!        KPROMA                    - HORIZONTAL DIMENSION.             (INPUT-C)
!        KSTART                    - START OF WORK.                    (INPUT-C)
!        KPROF                     - DEPTH OF WORK.                    (INPUT-C)
!        KFLEV                     - NUMBER OF INPUT PRESSURE LEVELS   (INPUT-C)
!        KLEVB(KPROMA,KLEVP,KPPM)  - INPUT LEVEL BELOW PRPRES          (INPUT-C)
!        KLOLEV                    - BEGINING FOR THE INTERPOLATION    (INPUT-C)
!        KPPM     - Number of interpolation methods in post-processing (INPUT-C)

!        KLEVB(KPROMA,KLEVP,3)     - INPUT LEVEL BELOW PRPRES           (INPUT-C)
!                                    (SEE PPFLEV)
!        KSLCT                     - 1:  FIELD GIVEN ON HALF LEVELS.
!                                        INTERPOLATE IN P
!                                  - 2:  FIELD GIVEN ON FULL LEVELS.
!                                        INTERPOLATE IN P
!                                  - 3:  FIELD GIVEN ON HALF LEVELS.
!                                        INTERPOLATE IN LN(P)          (INPUT-C)
!                                  - 4:  FIELD GIVEN ON FULL LEVELS.
!                                        INTERPOLATE IN LN(P)
!        LDBELO(KPROMA,KLEVP)      - .TRUE. IF PRESSURE IS UNDER
!                                     LOWEST  MODEL LEVEL              (INPUT-C)
!        LDBLOW(KLEVP)             - .TRUE. IF LDBELO IS   CONTAINING
!                                    AT LEAST ONE .TRUE.               (INPUT-C)
!        PRPRESC(KPROMA,KLEVP)      - POST-PROCESSING LEVEL
!                                    PRESSURES:
!                                    IF KVINT=1 OR KVINT = 2:
!                                     PRESSURE
!                                    IF KVINT=3:
!                                     LN PRESSURE                      (INPUT-C)
!        PRXP(KPROMA,0:KFLEV,KPPM) - HALF,FULL AND LN HALF,FULL LEVEL
!                                    PRESSURES (SEE PPINIT)            (INPUT)
!        PRXPD(KPROMA,0:KFLEV,KPPM)- 1./D(P) AND 1./D(LN(P))           (INPUT)
!        PFLDI(KPROMA,KLB:KFLEV)   - MODEL LEVEL FIELD VALUES          (INPUT)

!        PFLDO(KPROMA,KLEVP)       - INTERPOLATED FIELD VALUES         (OUTPUT)

!        IMPLICIT ARGUMENTS :  NONE.
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.  NONE.
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
!        R. El Khatib  27-Apr-2007 Optimization
!        T. Wilhelmsson 20-Jul-2012 Exact when pp model level are equal 
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPPM
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVB(KPROMA,KLEVP,KPPM) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSLCT 
LOGICAL           ,INTENT(IN)    :: LDBELO(KPROMA,KLEVP) 
LOGICAL           ,INTENT(IN)    :: LDBLOW(KLEVP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESC(KPROMA,KLEVP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXP(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXPD(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFLDI(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLDO(KPROMA,KLEVP) 
INTEGER(KIND=JPIM) :: IAB, IBL, JL, JLEVP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    LINEAR INTERPOLATION AND UPWARD EXTRAPOLATION.
!              ----------------------------------------------

!*       1.1   INTERPOLATE .
IF (LHOOK) CALL DR_HOOK('PPINTP',0,ZHOOK_HANDLE)

DO JLEVP=KLOLEV,KLEVP
  IF (.NOT.LDBLOW(JLEVP)) THEN
    DO JL=KSTART,KPROF
      IBL=KLEVB(JL,JLEVP,KSLCT)
      IAB=IBL-1
      PFLDO(JL,JLEVP)=PFLDI(JL,IBL)+REAL(PFLDI(JL,IBL)-PFLDI(JL,IAB),JPRD)* &
       & REAL(PRPRESC(JL,JLEVP)-PRXP(JL,IBL,KSLCT),JPRD)*REAL(PRXPD(JL,IAB,KSLCT),JPRD)
    ENDDO
  ELSE
    DO JL=KSTART,KPROF
      IF (.NOT.LDBELO(JL,JLEVP)) THEN
        IBL=KLEVB(JL,JLEVP,KSLCT)
        IAB=IBL-1
        PFLDO(JL,JLEVP)=PFLDI(JL,IBL)+REAL(PFLDI(JL,IBL)-PFLDI(JL,IAB),JPRD)* &
         & REAL(PRPRESC(JL,JLEVP)-PRXP(JL,IBL,KSLCT),JPRD)*REAL(PRXPD(JL,IAB,KSLCT),JPRD)
      ENDIF
    ENDDO
  ENDIF
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPINTP',1,ZHOOK_HANDLE)
END SUBROUTINE PPINTP

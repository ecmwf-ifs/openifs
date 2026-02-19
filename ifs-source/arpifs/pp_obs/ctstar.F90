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

SUBROUTINE CTSTAR(KPROMA,KSTART,KPROF,PTB,PRESBH,PRESBF,POROG,PTSTAR,PT0)

!**** *CTSTAR* - COMPUTES STANDARD SURFACE TEMPERATURE
!                              AND SURFACE TEMPERATURE.

!     PURPOSE.
!     --------
!           COMPUTES THE STANDARD SURFACE TEMPERATURE AND THE SURFACE
!           TEMPERATURE TO BE USED FOR EXTRAPOLATIONS OF TEMPERATURE
!           AND GEOPOTENTIEL.

!**   INTERFACE.
!     ----------
!        *CALL* *CTSTAR(..)*

!        EXPLICIT ARGUMENTS
!        --------------------

!        KPROMA         - HORIZONTAL DIMENSIONS.             (INPUT)
!        KSTART         - START OF WORK                      (INPUT)
!        KPROF          - DEPTH OF WORK                      (INPUT)
!        PTB(KPROMA)    - TEMPERATURE AT NFLEVG-1             (INPUT)
!        PRESBH(KPROMA) - LOWEST MODEL HALF LEVEL PRESSURES  (INPUT)
!        PRESBF(KPROMA) - PRESSURE AT NFLEVG-1                (INPUT)
!        POROG(KPROMA)  - MODEL ORGRAPHY                     (INPUT)

!        PTSTAR(KPROMA) - SURFACE TEMPERATURE                (OUTPUT)
!        PT0(KPROMA)    - STANDARD SURFACE TEMPERATURE       (OUTPUT)

!        IMPLICIT ARGUMENTS :    CONSTANTS FROM YOMSTA,YOMCST.
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.   NONE.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*
!      ORIGINAL : 89-05-02

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RD
USE YOMSTA   , ONLY : RDTDZ1

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTB(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESBH(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESBF(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSTAR(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT0(KPROMA) 
INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) :: ZALPHA, ZDTDZSG
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    COMPUTES SURFACE TEMPERATURE
!*             THEN STANDARD SURFACE TEMPERATURE.

IF (LHOOK) CALL DR_HOOK('CTSTAR',0,ZHOOK_HANDLE)
ZDTDZSG=-RDTDZ1/RG
ZALPHA=ZDTDZSG*RD
DO JL=KSTART,KPROF
  PTSTAR(JL)=PTB(JL)*(1.0_JPRB+ZALPHA*(PRESBH(JL)/PRESBF(JL)-1.0_JPRB))
  PT0(JL)=PTSTAR(JL)+ZDTDZSG*POROG(JL)
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CTSTAR',1,ZHOOK_HANDLE)
END SUBROUTINE CTSTAR

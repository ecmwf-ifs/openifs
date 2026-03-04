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

SUBROUTINE PPINIT(KPROMA,KSTART,KPROF,KFLEV,KPPM,PRPRESH,PRPRESF,PRXP,PRXPD)

!**** *PPINIT* - SET UP HELP-ARRAYS FOR VERTICAL INTERPOLATION

!     PURPOSE.
!     --------
!           SETS UP HELP-ARRAYS FOR PRESSURE LEVEL VERTICAL INTERPOLATION.
!           MUST BE CALLED BEFORE ANY INTEERPOLATION ROUTINES.

!**   INTERFACE.
!     ----------
!        *CALL* *PPINIT(KPROMA,KSTART,KPROF,KFLEV,PRPRESH,PRPRESF,
!                       PRXP,PRXPD)

!        EXPLICIT ARGUMENTS :
!        --------------------

!        KPROMA                  - HORIZONTAL DIMENSIONS.              (INPUT)
!        KSTART                  - START OF WORK                       (INPUT)
!        KPROF                   - DEPTH OF WORK                       (INPUT)
!        KFLEV                   - NUMBER OF INPUT PRESSURE LEVELS     (INPUT)
!        KPPM     - Number of interpolation methods in post-processing (INPUT)

!        PRPRESH(KPROMA,0:KFLEV)  - INPUT HALF LEVEL PRESSURES          (INPUT)
!        PRPRESF(KPROMA,KFLEV)    - INPUT FULL LEVEL PRESSURES          (INPUT)

!        PRXP(KPROMA,0:KFLEV,KPPM)- HALF,FULL AND LN HALF,FULL LEVEL
!                                  PRESSURES                           (OUTPUT)
!        PRXPD(KPROMA,0:KFLEV,KPPM)- 1./D(P) AND 1./D(LN(P))           (OUTPUT)

!        IMPLICIT ARGUMENTS :  NONE.
!        --------------------

!     EXTERNALS.   NONE.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        MATS HAMRUD  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-12-05
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPPM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRXP(KPROMA,0:KFLEV,KPPM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRXPD(KPROMA,0:KFLEV,KPPM) 
INTEGER(KIND=JPIM) :: JFLEV, JL, JSLCT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    COMPUTE HELP VARIABLES RELATED TO INPUT GRID.
!              ---------------------------------------------

!*       1.1   PRESSURE/LN PRESSURE AT HALF AND FULL LEVELS

IF (LHOOK) CALL DR_HOOK('PPINIT',0,ZHOOK_HANDLE)
DO JFLEV=1,KFLEV
  DO JL=KSTART,KPROF
    PRXP(JL,JFLEV,1)=PRPRESH(JL,JFLEV)
    PRXP(JL,JFLEV,2)=PRPRESF(JL,JFLEV)
    PRXP(JL,JFLEV,3)=LOG(PRPRESH(JL,JFLEV))
    PRXP(JL,JFLEV,4)=LOG(PRPRESF(JL,JFLEV))
  ENDDO
ENDDO

DO JL=KSTART,KPROF
  PRXP(JL,0,1)=PRPRESH(JL,0)
  PRXP(JL,0,2)=0.1_JPRB
  PRXP(JL,0,3)=LOG(0.1_JPRB)
  PRXP(JL,0,4)=LOG(0.1_JPRB)
ENDDO

!*       1.2   1/D(PRES) AND 1/D(LN(PRES))

DO JSLCT=1,KPPM
  DO JFLEV=0,KFLEV-1
    DO JL=KSTART,KPROF
      PRXPD(JL,JFLEV,JSLCT)=1.0_JPRB/&
       & (PRXP(JL,JFLEV+1,JSLCT)-PRXP(JL,JFLEV,JSLCT))  
    ENDDO
  ENDDO
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPINIT',1,ZHOOK_HANDLE)
END SUBROUTINE PPINIT

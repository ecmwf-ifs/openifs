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

SUBROUTINE GPPWC(KPROMA,KSTART,KPROF,KFLEV,PWCH,PQ,PRXP)

!**** *GPPWC* - Computes half level PWC.

!     Purpose.
!     --------
!           Computes half level PWC.

!**   Interface.
!     ----------
!        *CALL* *GPPWC(KPROMA,KSTART,KPROF,KFLEV,PWCH,PQ,PRXP)

!        Explicit arguments :
!        --------------------
!        KPROMA                     - HORIZ. DIMENSIONING      (INPUT)
!        KSTART                     - START OF WORK            (INPUT)
!        KPROF                      - DEPTH OF WORK            (INPUT)
!        KFLEV                      - NUMBER OF LEVELS         (INPUT)
!        PWCH (KPROMA,0:KFLEV)      - HALF LEVEL PWC           (OUTPUT)
!        PQ  (KPROMA,KFLEV)         - HUMIDITY                 (INPUT)
!        PRXP(KPROMA,0:KFLEV,NPPM)  - HALF,FULL AND LN HALF LEVEL
!                                     PRESSURES (SEE PPINIT)   (INPUT)

!        Implicit arguments :    None.
!        --------------------

!     Method.
!     -------
!        See documentation.

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Erik Andersson, based on PP-routines by
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 92-06-08
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A.Geer        09-Oct-2015 remove dependency on YRDIMF%NPPM to call from OOPS obsop
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWCH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRXP(KPROMA,0:KFLEV) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL, JLEV

REAL(KIND=JPRB) :: ZDELP, ZRGI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPPWC',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    COMPUTES HALF LEVEL PWC
!              -----------------------

ZRGI=1.0_JPRB/RG
DO JL=KSTART,KPROF
  PWCH(JL,0)=0.0_JPRB
ENDDO
DO JLEV=1,KFLEV
  DO JL=KSTART,KPROF
    ZDELP=PRXP(JL,JLEV)-PRXP(JL,JLEV-1)
    PWCH(JL,JLEV)=PWCH(JL,JLEV-1)+PQ(JL,JLEV)*ZDELP*ZRGI
  ENDDO
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPPWC',1,ZHOOK_HANDLE)
END SUBROUTINE GPPWC

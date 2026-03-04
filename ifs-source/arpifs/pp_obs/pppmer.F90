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

SUBROUTINE PPPMER(KPROMA,KSTART,KPROF,PRPRESS,POROG,PTSTAR,PT0,PMSLPPP,&
             & LDQNH)

!**** *PPPMER* - POST-PROCESS MSL PRESSURE.

!     PURPOSE.
!     --------
!           COMPUTES MSL PRESSURE.

!**   INTERFACE.
!     ----------
!        *CALL* *PPPMER(...)*

!        EXPLICIT ARGUMENTS
!        --------------------

!        KPROMA                    - HORIZONTAL DIMENSION.             (INPUT)
!        KSTART                    - START OF WORK.                    (INPUT)
!        KPROF                     - DEPTH OF WORK.                    (INPUT)
!        PRPRESS(KPROMA)           - SURFACE PRESSURE                  (INPUT)
!        POROG(KPROMA)             - MODEL OROGRAPHY.                  (INPUT)
!        PTSTAR(KPROMA)            - SURFACE TEMPERATURE               (INPUT)
!        PT0(KPROMA)               - STANDARD SURFACE TEMPERATURE      (INPUT)

!        PMSLPPP(KPROMA)           - POST-PROCESSED MSL PRESSURE       (OUTPUT)

!        IMPLICIT ARGUMENTS :  CONSTANTS FROM YOMCST,YOMGEM,YOMSTA.
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
!      MATS HAMRUD AND PHILIPPE COURTIER  *ECMWF*
!      ORIGINAL : 89-01-26

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 17-Jul-2013 FABEC post-processing
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTAR(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMSLPPP(KPROMA)
LOGICAL, OPTIONAL ,INTENT(IN)    :: LDQNH
 
REAL(KIND=JPRB) :: ZTSTAR(KPROMA)
REAL(KIND=JPRB) :: ZALPHA(KPROMA)

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) :: ZDTDZSG, ZOROG, ZT0, ZTX, ZTY, ZX, ZY, ZY2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL :: LLQNH


!     ------------------------------------------------------------------

!*       1.    POST-PROCESS MSL PRESSURE.
!              --------------------------

!*       1.1   COMPUTATION OF MODIFIED ALPHA AND TSTAR.

IF (LHOOK) CALL DR_HOOK('PPPMER',0,ZHOOK_HANDLE)

IF (PRESENT(LDQNH)) THEN
  LLQNH=LDQNH
ELSE
  LLQNH=.FALSE.
ENDIF

ZTX=290.5_JPRB
ZTY=255.0_JPRB
ZDTDZSG=-RDTDZ1/RG

IF (LLQNH) THEN
  ZT0=288.15_JPRB
  DO JL=KSTART,KPROF
    ZTSTAR(JL)= ZT0+(-ZDTDZSG)*POROG(JL)
    ZOROG=SIGN(MAX(1.0_JPRB,ABS(POROG(JL))),POROG(JL))
    ZALPHA(JL)=RD*(ZT0-ZTSTAR(JL))/ZOROG
  ENDDO
ELSE
  DO JL=KSTART,KPROF

    IF(PTSTAR(JL) < ZTY) THEN
      ZTSTAR(JL)=0.5_JPRB*(ZTY+PTSTAR(JL))
    ELSEIF(PTSTAR(JL) < ZTX) THEN
      ZTSTAR(JL)=PTSTAR(JL)
    ELSE
      ZTSTAR(JL)=0.5_JPRB*(ZTX+PTSTAR(JL))
    ENDIF

    ZT0=ZTSTAR(JL)+ZDTDZSG*POROG(JL)
    IF(ZTX > ZTSTAR(JL) .AND. ZT0 > ZTX) THEN
      ZT0=ZTX
    ELSEIF(ZTX <= ZTSTAR(JL) .AND. ZT0 > ZTSTAR(JL)) THEN
      ZT0=ZTSTAR(JL)
    ELSE
      ZT0=PT0(JL)
    ENDIF

    ZOROG=SIGN(MAX(1.0_JPRB,ABS(POROG(JL))),POROG(JL))
    ZALPHA(JL)=RD*(ZT0-ZTSTAR(JL))/ZOROG
  ENDDO
ENDIF

!*       1.2   COMPUTATION OF MSL PRESSURE.

DO JL=KSTART,KPROF
  IF (ABS(POROG(JL)) >= 0.001_JPRB) THEN
    ZX=POROG(JL)/(RD*ZTSTAR(JL))
    ZY=ZALPHA(JL)*ZX
    ZY2=ZY*ZY
    PMSLPPP(JL)=PRPRESS(JL)*EXP(ZX*(1.0_JPRB-0.5_JPRB*ZY+1.0_JPRB/3._JPRB*ZY2))
  ELSE
    PMSLPPP(JL)=PRPRESS(JL)
  ENDIF
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPPMER',1,ZHOOK_HANDLE)
END SUBROUTINE PPPMER

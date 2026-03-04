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

SUBROUTINE PPLTETA(KPROMA,KSTART,KPROF,KFLEV,PRPRESF,PTETA,&
 & PXTETA,PRPRES,PRLNPRES)  

!**** *PPLTETA*    Compute the pressures of the THETA level

!      PURPOSE.
!      --------
!           COMPUTE THE PRESSURES OF THE THETA LEVEL

!**    INTERFACE.
!      ----------
!           *CALL* PPLTETA( ... )

!           EXPLICITE ARGUMENTS.
!           --------------------
!           KPROMA      : horizontal dimension.                     (INPUT)
!           KSTART      : start of work.                            (INPUT)
!           KPROF       : depth of work.                            (INPUT)
!           KFLEV       : number of input pressure levels.          (INPUT)
!           PRPRESF(kproma,kflev) : model full level pressures.     (INPUT)
!           PTETA(kproma,kflev)   : THETA at full levels            (INPUT)
!           PXTETA      : value of the THETA level                  (INPUT)
!           PRPRES(kproma)    : pressures of the TETA level.        (OUTPUT)
!           PRLNPRES(kproma)  : ln-pressure of the TETA level.      (OUTPUT)

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!      EXTERNALS.
!      ----------
!           NONE

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!      AUTHOR.
!      -------
!           C. PERIARD, Klaus Von Der Emde, Ryad El Khatib

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 94-04-08
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!      -----------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTETA(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXTETA 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRPRES(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRLNPRES(KPROMA) 
INTEGER(KIND=JPIM) :: ILEV, JL, JSS

REAL(KIND=JPRB) :: ZQNN, ZSECUR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    COMPUTE PRESSURES
!              -----------------

IF (LHOOK) CALL DR_HOOK('PPLTETA',0,ZHOOK_HANDLE)
ZSECUR=1.E-10_JPRB

DO JL=KSTART,KPROF
!        Find nearest model level above the post-processing level
  ILEV=KFLEV
  DO JSS=KFLEV,1,-1
    IF (PTETA(JL,JSS) < PXTETA) THEN
      ILEV=JSS
    ENDIF
  ENDDO
!        Interpolation
  IF (ILEV == KFLEV) THEN
    PRPRES(JL)  =PRPRESF(JL,KFLEV)
    PRLNPRES(JL)=LOG(PRPRES(JL))
  ELSEIF(ILEV == 1) THEN
    PRPRES(JL)  =PRPRESF(JL,1)
    PRLNPRES(JL)=LOG(PRPRES(JL))
  ELSE
    ZQNN=PTETA(JL,ILEV-1)-PTETA(JL,ILEV)
    IF (ABS(ZQNN) > ZSECUR) THEN
      PRPRES(JL)=(PRPRESF(JL,ILEV-1)-PRPRESF(JL,ILEV))*&
       & (PXTETA-PTETA(JL,ILEV))/ZQNN &
       & + PRPRESF(JL,ILEV)  
    ELSE
      PRPRES(JL)=PRPRESF(JL,ILEV)
    ENDIF
    PRLNPRES(JL)=LOG(PRPRES(JL))
  ENDIF
ENDDO

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPLTETA',1,ZHOOK_HANDLE)
END SUBROUTINE PPLTETA


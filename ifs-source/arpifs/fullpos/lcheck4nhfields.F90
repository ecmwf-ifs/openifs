! (C) Copyright 1989- Meteo-France.

SUBROUTINE LCHECK4NHFIELDS(LDNHDYN,YDTFP,YD3,YD2)

!**** LCHECK4NHFIELDS **  Check for Non-Hydrostatic post-processing

!     Purpose.
!     --------
!           To check if any NH fields is requested in the post-processing 

!**   Interface.
!     ----------
!        *FUNCTION LCHECK4NHFIELDS(LDNHDYN,YD3,YD2)*

!        Explicit arguments :
!        --------------------
!           LDNHDYN : True id NH post-processing
!           YD3 : 3D fields fullpos structure
!           YD2 : 2D fields fullpos structure

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      R. El Khatib *METEO-FRANCE*
!      Original : 06-Jun-2019

!     Modifications.
!     ------------------------------------------------------------------

USE PARKIND1      , ONLY : JPIM,   JPRB
USE YOMHOOK       , ONLY : LHOOK,  DR_HOOK,  JPHOOK
USE YOM4FPOS      , ONLY : TRQ3FP, TRQ2FP
USE YOMAFN        , ONLY : ALL_FULLPOS_TYPES

IMPLICIT NONE

LOGICAL,                 INTENT(INOUT)        :: LDNHDYN
TYPE(ALL_FULLPOS_TYPES), INTENT(IN)           :: YDTFP
TYPE(TRQ3FP),            INTENT(IN), OPTIONAL :: YD3
TYPE(TRQ2FP),            INTENT(IN), OPTIONAL :: YD2

INTEGER(KIND=JPIM), EXTERNAL :: ISRCHEQ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LCHECK4NHFIELDS',0,ZHOOK_HANDLE)

IF(.NOT.LDNHDYN .AND. PRESENT(YD3)) THEN
  LDNHDYN = (YD3%NPPFIELDG > 0) .AND. &
   & ((ISRCHEQ(YD3%NPPFIELDG,YD3%ICOD(:),1,YDTFP%PD%ICOD) <= YD3%NPPFIELDG) .OR. &
   &  (ISRCHEQ(YD3%NPPFIELDG,YD3%ICOD(:),1,YDTFP%VD%ICOD) <= YD3%NPPFIELDG))
ENDIF
IF(.NOT.LDNHDYN .AND. PRESENT(YD2)) THEN
  LDNHDYN = (YD2%NPPFIELDG > 0) .AND. &
   & (ISRCHEQ(YD2%NPPFIELDG,YD2%ICOD(:),1,YDTFP%WWS%ICOD) <= YD2%NPPFIELDG)
ENDIF

IF (LHOOK) CALL DR_HOOK('LCHECK4NHFIELDS',1,ZHOOK_HANDLE)

END SUBROUTINE LCHECK4NHFIELDS

! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPSURF(YDFPOS,KPPEDR,KVCLIX)

        !**** *SUFPSURF*  - CONTRIBUTION OF FULLPOS TO SETUP THE MODEL SURFACE FIELDS

!     PURPOSE.
!     --------
!       Post-processing request may drive certain variables used to setup the model surface fields

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPSURF*

!        EXPLICIT ARGUMENTS
!        --------------------
!         YDFPOS : fullpos object
!         KPPEDR : 1 if EDR should be post-processed (INOUT because it may be initialized in input) ; 0 otherwise
!         KVCLIX : control variable to read auxilary model climatology fields (INOUT because it may be > 0 in input)

!        IMPLICIT ARGUMENTS
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 28-Sep-2021

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE FULLPOS  , ONLY : TFPOS

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TFPOS), INTENT(IN)           :: YDFPOS
INTEGER(KIND=JPIM), INTENT(INOUT) :: KPPEDR
INTEGER(KIND=JPIM), INTENT(INOUT) :: KVCLIX

INTEGER(KIND=JPIM) :: JFLD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUFPSURF',0,ZHOOK_HANDLE)
ASSOCIATE(CFP3DF=>YDFPOS%YNAMFPL%CFP3DF,YDNAMFPSCI=>YDFPOS%YNAMFPSCI)
!     ------------------------------------------------------------------

DO JFLD=1,SIZE(CFP3DF)
  IF(CFP3DF(JFLD) == ' ') THEN
    EXIT
  ELSEIF (CFP3DF(JFLD) == YDFPOS%YAFN%TFP%EDR%CLNAME) THEN
    KPPEDR=1
    EXIT
  ENDIF
ENDDO

IF (YDNAMFPSCI%NFPCLI >= 3) THEN
  IF (YDNAMFPSCI%NFPSWI >= 2) THEN
    KVCLIX=MAX(KVCLIX,4)
  ELSE
    KVCLIX=MAX(KVCLIX,3)
  ENDIF
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUFPSURF',1,ZHOOK_HANDLE)
END SUBROUTINE SUFPSURF

! (C) Copyright 1989- Meteo-France.

SUBROUTINE CHECKINCLIM(YDNAMFPSCI,KPROMA,YDSURF,KEND,POROG,PSD_VX,PSD_VV)

!**** *CHECKINCLIM*  - CHECK CLIMATOLOGY FIELDS AGAINST AN EXTERNAL DATASET

!     PURPOSE.
!     --------
!         To control the consistency of an external dataset on climatology
!         fields against the model fields (including orography or internal
!         climatology fields). 

!**   INTERFACE.
!     ----------
!       *CALL* *CHECKINCLIM*

!        EXPLICIT ARGUMENTS
!        --------------------
!         INPUT:
!          KEND       : last  adress in post-processing buffers to read
!          POROG      : Model orography
!          PSD_VX     : auxilary climatological diagnostic fields
!          PSD_VV     : vegetation diagnostic fields

!        IMPLICIT ARGUMENTS
!        --------------------
!        se #include below.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*
!        ORIGINAL : 19-Sep-2016 from HPOS

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMFPC   , ONLY : TNAMFPSCI
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TNAMFPSCI),  INTENT(IN) :: YDNAMFPSCI
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
TYPE(TSURF)       ,INTENT(IN)    :: YDSURF
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VX(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VV(:,:)

REAL(KIND=JPRB) :: ZDIF(KPROMA)

INTEGER(KIND=JPIM) :: JI, IST

REAL(KIND=JPRB) :: ZDIFMAX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHECKINCLIM',0,ZHOOK_HANDLE)
ASSOCIATE(RFPCORR=>YDNAMFPSCI%RFPCORR, RFPCSAB=>YDNAMFPSCI%RFPCSAB, RFPCD2=>YDNAMFPSCI%RFPCD2, &
 & YSD_VV=>YDSURF%YSD_VV, YSD_VX=>YDSURF%YSD_VX)
!     ------------------------------------------------------------------

!*       1. PRELIMINARY INITIALISATIONS.
!           ----------------------------

IST=1

IF (YSD_VX%YORO%LSET) THEN
  DO JI=IST,KEND
    ZDIF(JI)=ABS(POROG(JI)-PSD_VX(JI,YSD_VX%YORO%MP))
  ENDDO
  ZDIFMAX=MAXVAL(ZDIF(IST:KEND))
  IF (ZDIFMAX > RFPCORR) THEN
    WRITE (NULOUT,'(A,F9.2,A)')&
     & ' MAX. DIFFERENCE BETWEEN SOURCE AND MODEL OROGRAPHIES :  ',&
     & ZDIFMAX,' J/kg '  
    WRITE (NULOUT,'(A,F9.2,A)')' ALLOWED DIFFERENCE : ', RFPCORR,' J/kg '
    WRITE (NULOUT,'(A)')' TOO MUCH DIFFERENCE >> RERUN CONFIGURATION 923',&
     & ' OR INCREASE RFPCORR IN NAMELIST NAMFPC '
    CALL ABOR1(' CHECKINCLIM: SOURCE AND MODEL OROGRAPHIES MISMATCH ')
  ENDIF
ENDIF

IF (YSD_VV%YSAB%LSET.AND.YSD_VX%YSAB%LSET) THEN
  DO JI=IST,KEND
    ZDIF(JI)=ABS(PSD_VV(JI,YSD_VV%YSAB%MP)-PSD_VX(JI,YSD_VX%YSAB%MP))
  ENDDO
  ZDIFMAX=MAXVAL(ZDIF(IST:KEND))
  IF (ZDIFMAX > RFPCSAB) THEN
    WRITE (NULOUT,'(A,F9.6,A)')&
     & ' MAX. DIFFERENCE BETWEEN SOURCE AND MODEL PERCENTAGE OF SAND :  ',&
     & ZDIFMAX,' % '  
    WRITE (NULOUT,'(A,F9.6,A)')' ALLOWED DIFFERENCE : ',RFPCSAB,' % '
    WRITE (NULOUT,'(A)')' TOO MUCH DIFFERENCE >> RERUN CONFIGURATION 923',&
     & ' OR INCREASE RFPCSAB IN NAMELIST NAMFPC ' 
    CALL ABOR1(' CHECKINCLIM: SOURCE AND MODEL SAND MISMATCH ')
  ENDIF
ENDIF

IF (YSD_VV%YD2%LSET.AND.YSD_VX%YXD2%LSET) THEN
  DO JI=IST,KEND
    ZDIF(JI)=ABS(PSD_VV(JI,YSD_VV%YD2%MP)-PSD_VX(JI,YSD_VX%YXD2%MP))
  ENDDO
  ZDIFMAX=MAXVAL(ZDIF(IST:KEND))
  IF (ZDIFMAX > RFPCD2) THEN
    WRITE (NULOUT,'(A,F9.4,A)')&
     & ' MAX. DIFFERENCE BETWEEN SOURCE AND MODEL SOIL DEPTH :  ',&
     & ZDIFMAX,' m ' 
    WRITE (NULOUT,'(A,F9.4,A)')' ALLOWED DIFFERENCE : ',RFPCD2,' m '
    WRITE (NULOUT,'(A)')' TOO MUCH DIFFERENCE >> RERUN CONFIGURATION 923',&
     & ' OR INCREASE RFPCD2 IN NAMELIST NAMFPC ' 
    CALL ABOR1(' CHECKINCLIM: SOURCE AND MODEL SOIL DEPTH MISMATCH ')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHECKINCLIM',1,ZHOOK_HANDLE)
END SUBROUTINE CHECKINCLIM

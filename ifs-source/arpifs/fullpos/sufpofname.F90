! (C) Copyright 1989- Meteo-France.

SUBROUTINE SUFPOFNAME(CDNAME,CDFPFN,KSTEP,CDEXT,KDIGITS,PTSTEP)

!**** *SUFPOFNAME*  - PREPARE TO OPEN A SET OF FILES ARPEGE/ALADIN FOR
!                   POST-PROCESSING PURPOSE

!     PURPOSE.
!     --------
!           To compute the variables needed before opening files for 
!           post-processing purpose

!**   INTERFACE.
!     ----------
!       *CALL* *SUFPOFNAME(CDFPFN)*

!        EXPLICIT ARGUMENTS
!        --------------------
!            CDNAME : partial filename (filename without the time stamp)
!            CDFPFN : complete filename
!            KSTEP  : step stamp. If not present then the output file is stamped 'INIT'
!            CDEXT  : extension after time stamp
!            KDIGITS: number of digits used for the time stamp. 
!                     If negative the KDIGITS is used as abs(KDIGITS) without
!                     the character "+" preceeding it. This facility enables to
!                     make file names like boundary files in some sense.

!        IMPLICIT ARGUMENTS
!        --------------------
!         See modules above

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!        RYAD EL KHATIB *METEO-FRANCE*
!        ORIGINAL : 25-Feb-2016 from ini3wrfp

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*),   INTENT(IN)  :: CDNAME(:) 
CHARACTER(LEN=*),   INTENT(OUT) :: CDFPFN(:) 
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KSTEP
CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: CDEXT
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KDIGITS
REAL(KIND=JPRB),    INTENT(IN), OPTIONAL :: PTSTEP

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JD, IDIGITS, IST
CHARACTER(LEN=16)  :: CLINC ! output filename increment
! LLINC=.T. = compute the increment like a time step ("+XXXX")
! LLINC=.F. = compute the increment like a number ("XXXX")
LOGICAL :: LLINC 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "get_clinc.intfb.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUFPOFNAME',0,ZHOOK_HANDLE)

IF (SIZE(CDFPFN) /= SIZE(CDNAME)) CALL ABOR1('SUFPOFNAME:SIZES OF CDFPFN AND CDNAME MISMATCH')

IF (PRESENT(KDIGITS)) THEN
  LLINC=(KDIGITS > 0)
  IDIGITS=ABS(KDIGITS)
ELSE
  LLINC=.TRUE.
  IDIGITS=3
ENDIF

IF (PRESENT(KSTEP)) THEN
  IF (LLINC) THEN
    IST=6-IDIGITS+1
    IF (PRESENT(PTSTEP)) THEN
      CALL GET_CLINC (CLINC,KSTEP,IDIGITS,PTSTEP)
    ELSE
      CALL GET_CLINC (CLINC,KSTEP,IDIGITS)
    ENDIF
    DO JD=1,SIZE(CDFPFN)
      IF (LEN_TRIM(CDNAME(JD))+1+LEN_TRIM(CLINC) > LEN(CDFPFN(JD))) CALL ABOR1('SUFPOFNAME:LEN(CDFPFN) TOO SHORT')
      CDFPFN(JD)=TRIM(CDNAME(JD))//'+'//TRIM(CLINC)
    ENDDO
  ELSE
    IST=LEN(CLINC)-IDIGITS+1
    WRITE(CLINC,'(I16.16)') KSTEP
    DO JD=1,SIZE(CDFPFN)
      IF (LEN_TRIM(CDNAME(JD))+ IDIGITS > LEN(CDFPFN(JD))) CALL ABOR1('SUFPOFNAME:LEN(CDFPFN) TOO SHORT')
      CDFPFN(JD)=TRIM(CDNAME(JD))//TRIM(CLINC(IST:))
    ENDDO
  ENDIF
ELSE
  DO JD=1,SIZE(CDFPFN)
    IF (LEN_TRIM(CDNAME(JD))+4 > LEN(CDFPFN(JD))) CALL ABOR1('SUFPOFNAME:LEN(CDFPFN) TOO SHORT')
    CDFPFN(JD)=TRIM(CDNAME(JD))//'INIT'
  ENDDO
ENDIF

IF (PRESENT(CDEXT)) THEN
  IF (LEN_TRIM(CDEXT) > 0) THEN
    DO JD=1,SIZE(CDFPFN)
      IF (LEN_TRIM(CDFPFN(JD))+LEN_TRIM(CDEXT) > LEN(CDFPFN(JD))) CALL ABOR1('SUFPOFNAME:LEN(CDFPFN) TOO SHORT')
      CDFPFN(JD)=TRIM(CDFPFN(JD))//TRIM(CDEXT)
    ENDDO
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('SUFPOFNAME',1,ZHOOK_HANDLE)

END SUBROUTINE SUFPOFNAME

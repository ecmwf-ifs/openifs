! (C) Copyright 1989- Meteo-France.

SUBROUTINE FPSPNORMS(KRESOL,PSPEC,KVSET,LDSURF,CDPREF,PLEV,CDNAME,KGRIB,CDCONF)

!**** *FPSPNORMS*  - Fullpos spectral norms

!     Purpose.
!     --------
!        Compute and print norms of spectral post-processed fields

!**   Interface.
!     ----------

!        Explicit arguments :
!        --------------------
!         PSPEC  : distributed data array
!         KVSET  : V-set for each field
!         LDSURF : .TRUE. if surface field (not 3D) for each field
!         CDPREF : prefix / kind of level for each field
!         PLEV   : level for each field
!         CDNAME : FA name for each field
!         KGRIB  : GRIB code name for each field
!         CDCONF : configuration of I/Os

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Ryad El Khatib  *Meteo-France*
!      Original : 19-Nov-2015 from miscellaneous parts

! Modifications
! -------------
! End Modifications
!------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0    ,ONLY : MYPROC

!      -----------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL(:)
REAL(KIND=JPRB),    INTENT(IN) :: PSPEC(:,:)
INTEGER(KIND=JPIM), INTENT(IN) :: KVSET(:)
LOGICAL,            INTENT(IN) :: LDSURF(:)
CHARACTER(LEN=*),   INTENT(IN) :: CDPREF(:)
REAL(KIND=JPRB),    INTENT(IN) :: PLEV(:)
CHARACTER(LEN=*),   INTENT(IN) :: CDNAME(:)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)  :: KGRIB(:)
CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: CDCONF

INTEGER(KIND=JPIM), PARAMETER :: JPLEN=16
CHARACTER(LEN=JPLEN) :: CLEVEL, CLNOMA
INTEGER(KIND=JPIM) :: JFL, INCHIF, ILPRFU, ILSUFU, IOMASTER, ILEV
REAL(KIND=JPRB) :: ZAVEG(SIZE(PLEV))
LOGICAL :: LLETRANS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------

#include "specnorm.h"
#include "especnorm.h"

#include "abor1.intfb.h"
#include "updtrans.intfb.h"

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FPSPNORMS',0,ZHOOK_HANDLE)

IF (SIZE(KRESOL) > 1) CALL ABOR1 ('FPSPNORMS : NOT READY FOR MULTI-SPECTRAL TRANSFORMS')
IF (SIZE(KVSET) /= SIZE(PLEV)) CALL ABOR1('FPSPNORMS: SIZES OF KVSET AND PLEV MISMATCH')
IF (SIZE(LDSURF) /= SIZE(PLEV)) CALL ABOR1('FPSPNORMS: SIZES OF LDSURF AND PLEV MISMATCH')
IF (SIZE(CDPREF) /= SIZE(PLEV)) CALL ABOR1('FPSPNORMS: SIZES OF CDPREF AND PLEV MISMATCH')
IF (SIZE(CDNAME) /= SIZE(PLEV)) CALL ABOR1('FPSPNORMS: SIZES OF CDNAME AND PLEV MISMATCH')
IF (PRESENT(KGRIB)) THEN
  IF (SIZE(KGRIB) /= SIZE(PLEV)) CALL ABOR1('FPSPNORMS: SIZES OF KGRIB AND PLEV MISMATCH')
ENDIF

IOMASTER=1

CALL UPDTRANS(KRESOL(1),LLETRANS)
IF (LLETRANS) THEN
  CALL ESPECNORM(KRESOL=KRESOL(1),PSPEC=PSPEC,KVSET=KVSET,KMASTER=IOMASTER,PNORM=ZAVEG)
ELSE
  CALL SPECNORM(KRESOL=KRESOL(1),PSPEC=PSPEC,KVSET=KVSET,KMASTER=IOMASTER,PNORM=ZAVEG)
ENDIF
IF (MYPROC == IOMASTER)  THEN
!  A security while comparing listings of two experiments :
  WHERE (ZAVEG < EPSILON(1._JPRB)) ZAVEG = 0._JPRB
!  Print title the first time : 
  WRITE (NULOUT,FMT='(/,''FULL-POS SPNORMS'')')
  IF (PRESENT(KGRIB)) WRITE (NULOUT, FMT='(26X,A7,9X,A7)') 'PARAMID','AVERAGE'
  DO JFL= 1, SIZE(ZAVEG)
    ILPRFU=LEN_TRIM(CDPREF(JFL))
    IF (LDSURF(JFL)) THEN
      ILEV=0
      INCHIF=0
    ELSEIF(PRESENT(CDCONF)) THEN
      CALL FA_FIELD_FORMAT(PLEV(JFL),CDCONF,INCHIF,ILEV)
    ELSE
      IF (NINT(PLEV(JFL)) /= 0) THEN
        ILEV=NINT(PLEV(JFL))
      ELSE
!       Most likely PV level => deci-PVU
        ILEV=NINT(PLEV(JFL)*10000000._JPRB)
      ENDIF      
      INCHIF=MAX(INT(LOG10(ABS(PLEV(JFL))))+1,3)
    ENDIF
    IF (INCHIF /= 0) THEN
      WRITE (UNIT=CLEVEL,FMT='(I8.8)') ILEV
      IF (.NOT.PRESENT(KGRIB)) THEN
        ILSUFU=MIN(JPLEN-ILPRFU-INCHIF,LEN_TRIM(CDNAME(JFL)))
        CLNOMA=CDPREF(JFL)(1:ILPRFU)//CLEVEL(9-INCHIF:8)//CDNAME(JFL)(1:ILSUFU)
        WRITE (NULOUT, FMT='(1X,A16,'' : '',E21.15)') CLNOMA,ZAVEG(JFL)
      ELSE
        CLNOMA=CDPREF(JFL)(1:ILPRFU)//CLEVEL(9-INCHIF:8)
        WRITE(NULOUT,FMT='(1X,A8,1X,A16,1X,I6.0,'' : '',E21.15)') ADJUSTL(CLNOMA),CDNAME(JFL),KGRIB(JFL),ZAVEG(JFL)
      ENDIF
    ELSE
      IF (.NOT.PRESENT(KGRIB)) THEN
        ILSUFU=MIN(JPLEN-ILPRFU-INCHIF,LEN_TRIM(CDNAME(JFL)))
        CLNOMA=CDPREF(JFL)(1:ILPRFU)//CDNAME(JFL)(1:ILSUFU)
        WRITE (NULOUT, FMT='(1X,A16,'' : '',E21.15)') CLNOMA,ZAVEG(JFL)
      ELSE
        WRITE(NULOUT,FMT='(1X,A16,1X,I6.0,'' : '',E21.15)') CDNAME(JFL),KGRIB(JFL),ZAVEG(JFL)
      ENDIF
    ENDIF
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('FPSPNORMS',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE FA_FIELD_FORMAT(PLEV,CDCONF,KNCHIF,KLEV)

! Helper to define the format for printing the norms. Inspired from "fanfar".

REAL(KIND=JPRB),    INTENT(IN)  :: PLEV
CHARACTER(LEN=*),   INTENT(IN)  :: CDCONF
INTEGER(KIND=JPIM), INTENT(OUT) :: KNCHIF
INTEGER(KIND=JPIM), INTENT(OUT) :: KLEV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FPSPNORMS:FA_FIELD_FORMAT',0,ZHOOK_HANDLE)

! PLEV   : raw level
! CDNAME : field name
! CDCONF : configuration of I/Os (kind of level)
 
! KNCHIF : number of digits for formatted level
! KLEV   : formatted level as integer

SELECT CASE (CDCONF)
CASE ('V')
!   Potential vorticity levels are written in deci-PVU :
  KLEV=NINT(PLEV*10000000._JPRB)
CASE ('K')
  KLEV=NINT(PLEV)
  IF (KLEV < 0) KLEV=-KLEV
CASE DEFAULT
  KLEV=NINT(PLEV)
END SELECT
IF (KLEV <= 0 ) THEN
  KNCHIF=0
ELSE
  SELECT CASE (CDCONF)
    CASE ('P','H')
      KNCHIF=5
    CASE ('F')
      KNCHIF=4
    CASE DEFAULT
      KNCHIF=3
  END SELECT
ENDIF

IF (LHOOK) CALL DR_HOOK('FPSPNORMS:FA_FIELD_FORMAT',1,ZHOOK_HANDLE)

END SUBROUTINE FA_FIELD_FORMAT

END SUBROUTINE FPSPNORMS

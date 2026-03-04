! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GSTATS_LABEL(KNUM,CTYPE,CDESC)

USE EC_PARKIND  ,ONLY : JPIM

USE YOMGSTATS

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KNUM
CHARACTER(*) CDESC
CHARACTER(*) CTYPE
INTEGER(KIND=JPIM) :: ILEN, ITLEN

IF(KNUM < 0 .OR. KNUM>JPMAXSTAT) CALL ABOR1('GSTATS_LABEL:ILLEGAL KNUM')
ILEN = LEN(CDESC)
ILEN = MIN(ILEN,50)
ITLEN = LEN(CTYPE)
ITLEN = MIN(ILEN,3)
IF(CCDESC(KNUM) == '') THEN
  CCDESC(KNUM) = CDESC(1:ILEN)
  CCTYPE(KNUM) = CTYPE(1:ITLEN)
ELSEIF(CCDESC(KNUM)(1:ILEN) /= CDESC(1:ILEN)) THEN
  WRITE(JPERR,*)'LABEL',KNUM,' USED ',CCDESC(KNUM)
  CALL ABOR1('GSTATS_LABEL:OVERWRITE OF USED LABEL')
ENDIF

END SUBROUTINE GSTATS_LABEL


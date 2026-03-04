! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUEFPBIP(YDCLIMO,KFPCLI,YDAFN,KFIELDS,KCOD,CDCONF,LDBIP,YDRQAUX)
use parkind1 , only:&
 & jpim
USE TYPE_FPOSBUF, ONLY : FPOSBUF
USE YOMFP4L, ONLY : TRQFP
USE YOMAFN, ONLY : TAFN
TYPE (FPOSBUF), INTENT(IN) :: YDCLIMO
INTEGER(KIND=JPIM),INTENT(IN) :: KFPCLI
TYPE (TAFN), INTENT(IN) :: YDAFN
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN) :: KCOD(KFIELDS)
CHARACTER(LEN=1) ,INTENT(IN) :: CDCONF
LOGICAL ,INTENT(OUT) :: LDBIP(KFIELDS)
TYPE(TRQFP), INTENT(IN), OPTIONAL :: YDRQAUX
call abor1("suefpbip.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE SUEFPBIP

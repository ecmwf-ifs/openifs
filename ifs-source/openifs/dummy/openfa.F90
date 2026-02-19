! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE OPENFA (YDGEOMETRY,YDRIP,KFILE, KUNTIN, YDML_LBC, KDEB2, KFIN, KINDEX,&
 & KDATE, CDFILE, CDNOMC, KINITMONTH)
USE YEMLBC_MODEL , ONLY : TELBC_MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
use parkind1 , only:&
 & jpim
USE YOMRIP , ONLY : TRIP
TYPE(GEOMETRY) , INTENT(IN), TARGET :: YDGEOMETRY
TYPE(TRIP) , INTENT(IN) :: YDRIP
INTEGER(KIND=JPIM), INTENT(IN) :: KFILE
INTEGER(KIND=JPIM), INTENT(OUT) :: KUNTIN
TYPE(TELBC_MODEL) , INTENT(IN), OPTIONAL :: YDML_LBC
INTEGER(KIND=JPIM), INTENT(OUT), OPTIONAL :: KDEB2
INTEGER(KIND=JPIM), INTENT(OUT), OPTIONAL :: KFIN
INTEGER(KIND=JPIM), INTENT(OUT), OPTIONAL :: KINDEX
INTEGER(KIND=JPIM), INTENT(OUT), OPTIONAL :: KDATE(11)
CHARACTER(LEN=*) , INTENT(IN), OPTIONAL :: CDFILE
CHARACTER(LEN=*) , INTENT(OUT), OPTIONAL :: CDNOMC
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KINITMONTH
call abor1("openfa.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE OPENFA

! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE WRSCMR(KUNIT,CDNOM,PIN,KLON,KLEN)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KUNIT
REAL(KIND=JPRB) ,INTENT(IN) :: PIN(KLON,KLEN)
CHARACTER(LEN=*) ,INTENT(IN) :: CDNOM
call abor1("wrscmr.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE WRSCMR

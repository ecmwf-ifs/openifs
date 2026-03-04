! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE OPENFAINFO(KFILE,PTSTEP,KUNTIN,CDFILE,CDLEC,KTEST,CDMESS,LDERR)
use parkind1 , only:&
 & jpim,&
 & jprb
INTEGER(KIND=JPIM),INTENT(IN) :: KFILE
REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: PTSTEP
INTEGER(KIND=JPIM),INTENT(OUT), OPTIONAL :: KUNTIN
CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: CDFILE
CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: CDLEC
INTEGER(KIND=JPIM),INTENT(OUT), OPTIONAL :: KTEST
CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: CDMESS
LOGICAL, INTENT(OUT), OPTIONAL :: LDERR
call abor1("openfainfo.F90 should never be called with OpenIFS - EXITING")
END SUBROUTINE OPENFAINFO

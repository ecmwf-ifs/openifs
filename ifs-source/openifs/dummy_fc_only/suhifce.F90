! (C) Copyright 2011- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUHIFCE(PLON0E,PDLONE,PLATE,KLATE,KLONE,&
 & PLONM,PLATM,KNPTS,KNDIM1,KNFLDS,PFER,PFLDM,LDNG,KMASK)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLATE
INTEGER(KIND=JPIM),INTENT(IN) :: KNPTS
INTEGER(KIND=JPIM),INTENT(IN) :: KNDIM1
INTEGER(KIND=JPIM),INTENT(IN) :: KNFLDS
REAL(KIND=JPRB) ,INTENT(IN) :: PLON0E(KLATE)
REAL(KIND=JPRB) ,INTENT(IN) :: PDLONE(KLATE)
REAL(KIND=JPRB) ,INTENT(IN) :: PLATE(KLATE)
INTEGER(KIND=JPIM),INTENT(IN) :: KLONE(KLATE)
REAL(KIND=JPRB) ,INTENT(IN) :: PLONM(KNPTS)
REAL(KIND=JPRB) ,INTENT(IN) :: PLATM(KNPTS)
REAL(KIND=JPRB) ,INTENT(IN) :: PFER(KNDIM1,KNFLDS)
REAL(KIND=JPRB) ,INTENT(OUT) :: PFLDM(KNPTS,KNFLDS)
LOGICAL, OPTIONAL, INTENT(IN) :: LDNG
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KMASK(KNPTS)
call abor1("SUHIFCE should never be called with this build configuration - EXITING")
END SUBROUTINE SUHIFCE

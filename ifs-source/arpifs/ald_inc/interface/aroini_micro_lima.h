! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

INTERFACE
SUBROUTINE AROINI_MICRO_LIMA(KULOUT,KULNAM,PTSTEP,LDWARM,CMICRO,KSPLITR,KSPLITG,CCSEDIM,LDCRIAUTI,&
          PCRIAUTI,PT0CRIAUTI,PCRIAUTC)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
INTEGER(KIND=JPIM), INTENT (IN) :: KULOUT
INTEGER(KIND=JPIM), INTENT (IN) :: KULNAM
REAL(KIND=JPRB), INTENT (IN) :: PTSTEP
LOGICAL, INTENT (IN) :: LDWARM
CHARACTER (LEN=4), INTENT (IN) :: CMICRO
CHARACTER(4), INTENT (IN) :: CCSEDIM
INTEGER(KIND=JPIM), INTENT (OUT) :: KSPLITR
INTEGER(KIND=JPIM), INTENT (OUT) :: KSPLITG
LOGICAL, INTENT (IN) :: LDCRIAUTI
REAL(KIND=JPRB), INTENT (IN) :: PCRIAUTI
REAL(KIND=JPRB), INTENT (IN) :: PT0CRIAUTI
REAL(KIND=JPRB), INTENT (IN) :: PCRIAUTC
END SUBROUTINE AROINI_MICRO_LIMA
END INTERFACE

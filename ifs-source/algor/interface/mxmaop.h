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
SUBROUTINE MXMAOP(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)
USE PARKIND1 ,ONLY : JPIM ,JPRB
REAL(KIND=JPRB) ,INTENT(IN) :: PA(*)
INTEGER(KIND=JPIM),INTENT(IN) :: KA
INTEGER(KIND=JPIM),INTENT(IN) :: KAD
REAL(KIND=JPRB) ,INTENT(IN) :: PB(*)
INTEGER(KIND=JPIM),INTENT(IN) :: KB
INTEGER(KIND=JPIM),INTENT(IN) :: KBD
REAL(KIND=JPRB) ,INTENT(OUT) :: PC(*)
INTEGER(KIND=JPIM),INTENT(IN) :: KC
INTEGER(KIND=JPIM),INTENT(IN) :: KCA
INTEGER(KIND=JPIM),INTENT(IN) :: KAR
INTEGER(KIND=JPIM),INTENT(IN) :: KAC
INTEGER(KIND=JPIM),INTENT(IN) :: KBC
END SUBROUTINE MXMAOP
END INTERFACE

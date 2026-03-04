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
SUBROUTINE SUHER(KLX,KVX,KVXS,PD,PEI,PES,PFI,PFS,PA,PB,PC,PG,PH)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLX
INTEGER(KIND=JPIM),INTENT(IN) :: KVXS
INTEGER(KIND=JPIM),INTENT(IN) :: KVX
REAL(KIND=JPRB) ,INTENT(IN) :: PD(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PEI(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PES(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PFI(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PFS(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PA(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PB(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PC(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PG(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PH(KVXS,KLX)
END SUBROUTINE SUHER
END INTERFACE

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
SUBROUTINE MXTURHD(KLX,KVX,KVXS,KT,LDMT,PA,PB,PY,PX)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLX
INTEGER(KIND=JPIM),INTENT(IN) :: KVX
INTEGER(KIND=JPIM),INTENT(IN) :: KVXS
INTEGER(KIND=JPIM),INTENT(IN) :: KT
LOGICAL ,INTENT(IN) :: LDMT
REAL(KIND=JPRB) ,INTENT(IN) :: PA(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PB(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PY(KVXS,KLX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PX(KVXS,KLX)
END SUBROUTINE MXTURHD
END INTERFACE

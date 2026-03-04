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
SUBROUTINE MXPTMA(KLX,KVX,KVXS,KIX,PA,PBI,PCI,PBS,PCS,PX,PY)
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLX
INTEGER(KIND=JPIM),INTENT(IN) :: KVXS
INTEGER(KIND=JPIM),INTENT(IN) :: KIX
INTEGER(KIND=JPIM),INTENT(IN) :: KVX
REAL(KIND=JPRB) ,INTENT(IN) :: PA(KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PBI(KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PCI(KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PBS(KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PCS(KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PX(KVXS,KLX,KIX)
REAL(KIND=JPRB) ,INTENT(OUT) :: PY(KVXS,KLX,KIX)
END SUBROUTINE MXPTMA
END INTERFACE

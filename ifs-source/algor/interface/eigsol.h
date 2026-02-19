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
SUBROUTINE EIGSOL(KFLEVG,KNFLEVG,PA,PFR,PFI,K,PMO,KWO,PWO,KER)
USE PARKIND1 ,ONLY : JPIM ,JPRB

INTEGER(KIND=JPIM),INTENT(IN)   :: KFLEVG
INTEGER(KIND=JPIM),INTENT(IN)   :: KNFLEVG
REAL(KIND=JPRB),INTENT(IN)      :: PA(*)
REAL(KIND=JPRB),INTENT(OUT)     :: PFR(*)
REAL(KIND=JPRB),INTENT(OUT)     :: PFI(*)
INTEGER(KIND=JPIM),INTENT(IN)   :: K
REAL(KIND=JPRB),INTENT(OUT)     :: PMO(*)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KWO(*)
REAL(KIND=JPRB),INTENT(OUT)     :: PWO(*)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KER

END SUBROUTINE EIGSOL
END INTERFACE

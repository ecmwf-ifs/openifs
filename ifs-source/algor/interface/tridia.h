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
SUBROUTINE TRIDIA(KN,KSYS,KFIRST,KEND,KTYP,PM,PRHS,PSOL)

! -----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! -----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSYS
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIRST
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PM(1+(KTYP-1)*(KSYS-1),KN,-1:1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHS(KSYS,KN)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSOL(KSYS,KN)

! -----------------------------------------------------------------------------

END SUBROUTINE TRIDIA
END INTERFACE

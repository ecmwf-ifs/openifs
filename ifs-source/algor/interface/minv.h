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
SUBROUTINE MINV(PAB,KDIMN,KDBA,PZSCRA,PDET1,PTOL,KDIMM,KMODE)
USE PARKIND1, ONLY : JPIM, JPRB
INTEGER(KIND=JPIM), INTENT(IN)    :: KDIMN
INTEGER(KIND=JPIM), INTENT(IN)    :: KDBA
INTEGER(KIND=JPIM), INTENT(IN)    :: KDIMM
INTEGER(KIND=JPIM), INTENT(IN)    :: KMODE
REAL(KIND=JPRB),    INTENT(IN)    :: PTOL
REAL(KIND=JPRB),    INTENT(OUT)   :: PDET1
REAL(KIND=JPRB),    INTENT(INOUT) :: PAB(KDBA,KDIMN+KDIMM)
REAL(KIND=JPRB),    INTENT(INOUT) :: PZSCRA(2*KDIMN)
END SUBROUTINE MINV
END INTERFACE

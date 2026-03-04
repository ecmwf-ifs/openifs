! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MEMINFO(KOUT,KSTEP)
USE EC_PARKIND, ONLY : JPIM
IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN) :: KOUT, KSTEP
CHARACTER(LEN=32) CLSTEP
INTEGER(KIND=JPIM) :: ICOMM
#include "ec_meminfo.intfb.h"
WRITE(CLSTEP,'(11X,"STEP",I5," :")') KSTEP
ICOMM = -2 ! No headers from EC_MEMINFO by default
IF (KSTEP == 0) ICOMM = -1 ! Print also headers
CALL EC_MEMINFO(KOUT,TRIM(CLSTEP),ICOMM,KBARR=0,KIOTASK=-1,KCALL=-1)
CALL EC_FLUSH(KOUT)
END SUBROUTINE MEMINFO

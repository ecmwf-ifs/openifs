! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM MASTER1S
USE PARKIND1  ,ONLY : JPRD, JPIM
USE YOMLUN1S , ONLY : NULOUT
USE MPL_MODULE
USE OMP_LIB

IMPLICIT NONE

REAL (KIND=JPRD) :: ZTT0,ZTT1
REAL (KIND=JPIM) :: IERR, NPROC, MYPROC


#include "cnt01s.intfb.h"

!  Driver for offline version of surface code

CALL MPL_INIT()

NPROC=MPL_NPROC()
MYPROC=MPL_MYRANK()

ZTT0 = OMP_GET_WTIME()
CALL CNT01S
ZTT1 = OMP_GET_WTIME()

WRITE(NULOUT,'(a22,f12.3)') 'MASTER1s: Time total: ',ZTT1-ZTT0
CALL MPL_END()

END PROGRAM MASTER1S

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

SUBROUTINE PRTIME(YDRIP,CDGREP)

!     Purpose.
!     --------
!       Prints current date and time.

!     Author.
!     -------
!       Y. Tremolet
!       Original : 22-Feb-2005

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
! ----------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB, JPRD
USE YOMHOOK , ONLY: LHOOK, DR_HOOK, JPHOOK
USE YOMRIP0 , ONLY: NINDAT, NSSSSS
USE YOMRIP  , ONLY : TRIP
USE YOMCST  , ONLY: RDAY
USE YOMLUN  , ONLY: NULOUT

! ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP)       ,INTENT(INOUT):: YDRIP
CHARACTER(LEN=*), INTENT(IN) :: CDGREP

! ----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IYY,IMM,IDD,IS,IM,IH
INTEGER(KIND=JPIM) :: IY0,IM0,ID0,INC,IMON(12)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

#include "updcal.intfb.h"
#include "fcttim.func.h"

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PRTIME',0,ZHOOK_HANDLE)
ASSOCIATE(RSTATI=>YDRIP%RSTATI)
! ----------------------------------------------------------------------

IY0=NCCAA(NINDAT)
IM0=NMM(NINDAT)
ID0=NDD(NINDAT)
INC=(NSSSSS+NINT(RSTATI))/NINT(RDAY)
CALL UPDCAL(ID0,IM0,IY0,INC,IDD,IMM,IYY,IMON,-1)

WRITE(NULOUT,'(A,A,I4,2(''-'',I2.2))')CDGREP,' Current date is: ',IYY,IMM,IDD

IS=MOD(NSSSSS+NINT(RSTATI),NINT(RDAY))
IH=IS/3600
IS=IS-3600*IH
IM=IS/60
IS=IS-60*IM

WRITE(NULOUT,'(A,A,3(I2.2,A))')CDGREP,' Current time is: ',&
                             & IH,'h',IM,':',IS,' GMT'

! ----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PRTIME',1,ZHOOK_HANDLE)
END SUBROUTINE PRTIME

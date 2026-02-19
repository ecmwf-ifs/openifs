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

SUBROUTINE FPWRNCF(CDFPNCF,KSTEP,PTSTEP)

!**** *FPWRNCF*  - Write Fullpos status

!        Explicit arguments :
!        --------------------

!        CDFPNCF : status file
!        KSTEP : time step. If absent, 'INIT' will be used

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014
!      R. El Khatib 10-Dec-2015 KSTEP in argument (OOPS)
!      R. El Khatib 17-Oct-2017 pathname (to support multiple objects)


USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0  , ONLY : MYPROC 
USE YOMLUN  , ONLY : NULFPOS, NULOUT

IMPLICIT NONE

CHARACTER(LEN=*),   INTENT(IN) :: CDFPNCF
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KSTEP
REAL (KIND=JPRB),   INTENT(IN), OPTIONAL :: PTSTEP

INTEGER(KIND=JPIM) :: IDIGITS
CHARACTER (LEN=16) :: CLINC

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

#include "get_clinc.intfb.h"


IF (LHOOK) CALL DR_HOOK ('FPWRNCF',0,ZHOOK_HANDLE)


  IF (MYPROC == 1)  THEN
    CALL GSTATS(1941,0)
    IF (PRESENT(KSTEP)) THEN
      IDIGITS=6
      IF (PRESENT(PTSTEP)) THEN
        CALL GET_CLINC (CLINC,KSTEP,IDIGITS,PTSTEP)
      ELSE
        CALL GET_CLINC (CLINC,KSTEP,IDIGITS)
      ENDIF
    ELSE
      CLINC='INIT'
    ENDIF
    OPEN (UNIT=NULFPOS,FILE=CDFPNCF,FORM='FORMATTED')
    REWIND(NULFPOS)
    WRITE(NULFPOS,'(A)') TRIM (CLINC)
    WRITE(NULOUT,'('' POST FULLPOS EVENT : '',A)') TRIM (CLINC)
    CLOSE(UNIT=NULFPOS,STATUS='KEEP')
    CALL GSTATS(1941,1)
  ENDIF

IF (LHOOK) CALL DR_HOOK ('FPWRNCF',1,ZHOOK_HANDLE)

END SUBROUTINE FPWRNCF

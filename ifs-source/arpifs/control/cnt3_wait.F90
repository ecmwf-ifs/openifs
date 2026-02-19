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

SUBROUTINE CNT3_WAIT (KCNT3)

!     Modifications.
!     --------------
!      R. El Khatib  04-Jul-2014 call sualspa1 before sualspa to facilitate the 
!      subsequent allocation of the spectral structure by alloc_spec
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib : 06-Mar-2015 script name for pp server in namelist
!      R. El Khatib  03-Jul-2015 Move setup to cprep3.F90
! ---------------------------------------------------------------

USE PARKIND1      , ONLY : JPRB, JPIM
USE YOMHOOK       , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0        , ONLY : CSCRIPT_PPSERVER
USE YOMMP0        , ONLY : MYPROC
USE MPL_MODULE    , ONLY : MPL_BROADCAST, MPL_BARRIER

IMPLICIT NONE

INTEGER (KIND=JPIM), INTENT (INOUT) :: KCNT3

#include "abor1.intfb.h"

LOGICAL :: LLEXIST
INTEGER (KIND=JPIM) :: IRET

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('CNT3_WAIT',0,ZHOOK_HANDLE)

IF (CSCRIPT_PPSERVER == ' ') THEN

  KCNT3 = 0

ELSE

  CALL MPL_BARRIER (CDSTRING='CNT3_WAIT:')

  IF (MYPROC == 1) THEN

  ! Check for sync script existence

    INQUIRE (FILE=TRIM(CSCRIPT_PPSERVER), EXIST=LLEXIST)

    IF (LLEXIST) THEN

      CALL SPAWN (TRIM(CSCRIPT_PPSERVER), IRET)

      IF (IRET /= 0) THEN
        CALL ABOR1 ('CNT3_WAIT: '//TRIM(CSCRIPT_PPSERVER)//' FAILED')
      ENDIF

      INQUIRE (FILE=TRIM(CSCRIPT_PPSERVER), EXIST=LLEXIST)

      IF (LLEXIST) THEN
        KCNT3 = KCNT3 + 1
      ELSE
        KCNT3 = -1
      ENDIF

    ELSE

      KCNT3 = 0

    ENDIF

  ENDIF

  CALL MPL_BROADCAST (KCNT3, KTAG=0, CDSTRING='CNT3_WAIT:')

ENDIF

IF (LHOOK) CALL DR_HOOK ('CNT3_WAIT',1,ZHOOK_HANDLE)

END SUBROUTINE CNT3_WAIT

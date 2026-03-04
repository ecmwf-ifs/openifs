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

SUBROUTINE SUFAINIT

!**** *SUFAINIT*  - Set FA limits

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014


USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
USE FA_MOD, ONLY : YFA => FA_COM_DEFAULT, NEW_FA, FA_COM_DEFAULT_INIT
USE LFIMOD, ONLY : NEW_LFI_DEFAULT, LFICOM_DEFAULT
USE YOMFA, ONLY : JPXTRO, JPXLAT, JPXNIV, JPNXFA, JPNXCA
USE YOMLUN, ONLY : NULNAM

IMPLICIT NONE

#include "posname.intfb.h"
#include "namfainit.nam.h"

INTEGER (KIND=JPIM) :: ISTAT, IERR

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('SUFAINIT',0,ZHOOK_HANDLE)

JPXTRO=1280
JPXLAT=2560
JPXNIV=200 
JPNXFA=20 
JPNXCA=20 

CALL POSNAME (NULNAM, 'NAMFAINIT', ISTAT)

IF (ISTAT == 0) THEN
  READ (NULNAM, NAMFAINIT)
ENDIF

CALL NEW_LFI_DEFAULT
CALL NEW_FA (YFA, IERR, JPXTRO, JPXLAT, JPXNIV, JPNXFA, JPNXCA)
YFA%LFI => LFICOM_DEFAULT
FA_COM_DEFAULT_INIT = .TRUE.

IF (LHOOK) CALL DR_HOOK ('SUFAINIT',1,ZHOOK_HANDLE)

END SUBROUTINE SUFAINIT


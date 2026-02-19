! (C) Copyright 2005- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE CDRHOOKINIT(KRET)
!-- Makes sure Dr.Hook gets properly initialized from C-main program, too
USE EC_PARKIND  ,ONLY  : JPIM, JPRD
USE YOMHOOK   ,ONLY  : LHOOK, DR_HOOK_INIT
IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(OUT) :: KRET
CALL DR_HOOK_INIT()
IF (LHOOK) THEN
  KRET = 1 ! Dr.Hook is ON
ELSE
  KRET = 0 ! Dr.Hook is OFF
ENDIF
END SUBROUTINE CDRHOOKINIT

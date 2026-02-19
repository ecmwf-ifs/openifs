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

INTEGER (KIND=JPIM) FUNCTION JMKOUT ()

IF ((1 <= KLSIM) .AND. (KLSIM <= 2)) THEN

  JMKOUT=NINT(PLSMOUT(JLEN))

  IF(JMKOUT == 0)THEN
    JMKOUT=1
    IF(PTSOUT(JLEN) <= PTMERGL.AND.KLSIM == 2) JMKOUT=2
  ELSE
    JMKOUT=3
  ENDIF

ELSEIF (KLSIM == 3) THEN

  IF (PRESENT (LDPSL)) THEN
    JMKOUT=4
  ELSEIF (PLSMOUT(JLEN) > 0._JPRB) THEN
    JMKOUT=4
  ELSE
    JMKOUT=5
  ENDIF

ENDIF

END FUNCTION JMKOUT

INTEGER (KIND=JPIM) FUNCTION JMKIN (KOFF, KIN1, KIN2)

INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL :: KOFF
INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL :: KIN1
INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL :: KIN2

#include "abor1.intfb.h"

INTEGER (KIND=JPIM) :: IIND

IF (PRESENT (KOFF)) THEN
  IIND = KL0(JLEN,II1(KOFF)+1) + II2 (KOFF)
ELSEIF (PRESENT (KIN1) .AND. PRESENT (KIN2)) THEN
  IIND = KL0(JLEN,KIN1+1)+KIN2
ELSE
  CALL ABOR1 ('JMKIN: KOFF OR (KIN1+KIN2) REQUIRED')
ENDIF

IF ((1 <= KLSIM) .AND. (KLSIM <= 2)) THEN

  JMKIN=NINT(PLSMIN(IIND))
  
  IF(JMKIN == 0)THEN
    JMKIN=1
    IF(PTSIN(IIND) <= PTMERGL.AND.KLSIM == 2)JMKIN=2
  ELSE
    JMKIN=3
  ENDIF

ELSEIF (KLSIM == 3) THEN

  IF (PRESENT (LDPSL)) THEN
    IF (LDPSL (KOFF)) THEN
      JMKIN=4
    ELSE
! This will affect a 0. weight to this point
      JMKIN=6
    ENDIF
  ELSEIF (PLSMIN(IIND) > 0._JPRB) THEN
    JMKIN=4
  ELSE
! This will affect a 0. weight to this point
    JMKIN=6
  ENDIF

ENDIF

END FUNCTION JMKIN


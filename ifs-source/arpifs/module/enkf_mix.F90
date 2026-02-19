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

MODULE ENKF_MIX

!   Purpose
!   -------
!     Switches for Ensemble Kalman Filter within the IFS

!   Author
!   ------
!     Mats Hamrud

!   Modifications
!   -------------
!     Original    23-Feb-2010
!
USE PARKIND1, ONLY: JPIM, JPRB
USE YOMHOOK , ONLY: LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT   ,NULNAM

IMPLICIT NONE
SAVE

LOGICAL :: LENKF
INTEGER(KIND=JPIM) :: NSIZE_ENSEMBLE
INTEGER(KIND=JPIM) :: MYMEMBER
INTEGER(KIND=JPIM) :: NOFFSET_OBS

CONTAINS
!========================================================================
SUBROUTINE SETUP_ENKF

!   Purpose
!   -------
!    Setup routine for ENKF control

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

NAMELIST/NAMENKF/LENKF,NSIZE_ENSEMBLE,MYMEMBER,NOFFSET_OBS

#include "posnam.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ENKF_MIX:SETUP_ENKF',0,ZHOOK_HANDLE)


LENKF = .FALSE.
NSIZE_ENSEMBLE = 0
MYMEMBER = 0
NOFFSET_OBS = 0

CALL POSNAM(NULNAM,'NAMENKF')
READ(NULNAM,NAMENKF)

WRITE(NULOUT,NAMENKF)

IF (LHOOK) CALL DR_HOOK('ENKF_MIX:SETUP_ENKF',1,ZHOOK_HANDLE)

END SUBROUTINE SETUP_ENKF
!========================================================================
END MODULE ENKF_MIX

! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WVWAMDECOMP

! ----------------------------------------------------------------------

!**** *WVWAMDECOMP* - WAVE MODEL MODEL DECOMPOSITION

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : LOUTMDLDCP
      USE YOWMPP   , ONLY : NPROC    ,MPMAXLENGTH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpdecomp.intfb.h"
#include "outmdldcp.intfb.h"

      INTEGER(KIND=JWIM) :: NPR, MAXLEN

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL :: LLIRANK, LLWVENVI

! ----------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('WVWAMDECOMP',0,ZHOOK_HANDLE)

      NPR=NPROC
      LLIRANK = .FALSE.
      LLWVENVI = .TRUE.
      CALL MPDECOMP(NPR, MAXLEN, LLIRANK, LLWVENVI)
      MPMAXLENGTH=MAXLEN

      ! Write model decomposition to a file
      IF ( LOUTMDLDCP ) CALL OUTMDLDCP

      IF (LHOOK) CALL DR_HOOK('WVWAMDECOMP',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WVWAMDECOMP

! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE STOKESTRN (KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, WSWAVE, WDWAVE, CICOVER, CITHICK, &
& USTOKES, VSTOKES, STRNMS, NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN)

! ----------------------------------------------------------------------

!**** *STOKESTRN* - WRAPPER TO CALL STOKESDRIFT and CIMSSTRN 

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *STOKESTRN (KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, FF_NOW, INTFLDS, WAM2NEMO)

!          *KIJS*    - INDEX OF FIRST GRIDPOINT.
!          *KIJL*    - INDEX OF LAST GRIDPOINT.
!          *FL1*     - SPECTRUM(INPUT).
!          *WAVNUM*  - WAVE NUMBER.
!          *STOKFAC* - STOKES DRIFT FACTOR.
!          *DEPTH*   - WATER DEPTH.
!          *FF_NOW*  - FORCING FIELDS AT CURRENT TIME.
!          *INTFLDS* - INTEGRATED/DERIVED PARAMETERS
!          *WAM2NEMO*- WAVE FIELDS PASSED TO NEMO

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS, INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWCOUP  , ONLY : LWCOU, LWNEMOCOU, LWNEMOCOUSEND, LWNEMOCOUSTK, LWNEMOCOUSTRN
      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "cimsstrn.intfb.h"
#include "stokesdrift.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, STOKFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WSWAVE, WDWAVE, CICOVER, CITHICK
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN


      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('STOKESTRN',0,ZHOOK_HANDLE)

CALL STOKESDRIFT(KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES,VSTOKES)

IF (LWNEMOCOUSTRN) CALL CIMSSTRN(KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRNMS)


IF (LWNEMOCOU .AND.                                    & 
&   ( (LWNEMOCOUSEND .AND. LWCOU) .OR. (.NOT.LWCOU) )  &
&  ) THEN
  IF (LWNEMOCOUSTK) THEN
    NEMOUSTOKES(KIJS:KIJL) = USTOKES(KIJS:KIJL)
    NEMOVSTOKES(KIJS:KIJL) = VSTOKES(KIJS:KIJL)
  ELSE
    NEMOUSTOKES(KIJS:KIJL) = 0.0_JWRO
    NEMOVSTOKES(KIJS:KIJL) = 0.0_JWRO
  ENDIF

  IF (LWNEMOCOUSTRN) NEMOSTRN(KIJS:KIJL) = STRNMS(KIJS:KIJL)
ENDIF


! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('STOKESTRN',1,ZHOOK_HANDLE)

END SUBROUTINE STOKESTRN

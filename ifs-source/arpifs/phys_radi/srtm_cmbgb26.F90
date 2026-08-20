! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
SUBROUTINE SRTM_CMBGB26

!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM , JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOESRTM  , ONLY : NGN
USE YOESRTWN , ONLY : NGC, NGS, RWGT
!USE YOESRTWN , ONLY : NGC, NGS, NGN, RWGT
USE YOESRTA26, ONLY : SFLUXREF, RAYL, &
                    & SFLUXREFC, RAYLC

IMPLICIT NONE

! Local variables
INTEGER(KIND=JPIM) :: IGC, IPR, IPRSM
REAL(KIND=JPRB)    :: ZSUMF1, ZSUMF2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB26',0,ZHOOK_HANDLE)

IPRSM = 0
DO IGC = 1,NGC(11)
  ZSUMF1 = 0.
  ZSUMF2 = 0.
  DO IPR = 1, NGN(NGS(10)+IGC)
    IPRSM = IPRSM + 1
    ZSUMF1 = ZSUMF1 + RAYL(IPRSM)*RWGT(IPRSM+160)
    ZSUMF2 = ZSUMF2 + SFLUXREF(IPRSM)
  ENDDO
  RAYLC(IGC) = ZSUMF1
  SFLUXREFC(IGC) = ZSUMF2
ENDDO

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB26',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_CMBGB26


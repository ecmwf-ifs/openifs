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

MODULE YOMCSGEOM

USE PARKIND1, ONLY : JPIM, JPRB, JPRD

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    * Computational sphere horizontal geometry

!     RCOLON(NGPTOT) cosine of longitude on transformed sphere
!     RSILON(NGPTOT)   sine        "             "         "
!     RINDX (NGPTOT) Longitude index
!     RINDY (NGPTOT) Latitude index
!     RATATH(NGPTOT) RA*TAN(THETA) on real sphere
!     RATATX(NGPTOT) Curvature term for LAM (for u eq.)

TYPE TCSGEOM
  REAL(KIND=JPRD), POINTER :: RCOLON(:) => NULL()
  REAL(KIND=JPRD), POINTER :: RSILON(:) => NULL()
  REAL(KIND=JPRB), POINTER :: RINDX (:) => NULL()
  REAL(KIND=JPRB), POINTER :: RINDY (:) => NULL()
  REAL(KIND=JPRB), POINTER :: RATATH(:) => NULL()
  REAL(KIND=JPRB), POINTER :: RATATX(:) => NULL()
END TYPE TCSGEOM

! define blocked and non-blocked (_NB) structures
! note that the blocked structure YRCSGEOM will be initialised to point into 
! the non-blocked structure YRCSGEOM_NB
!!TYPE(TCSGEOM), POINTER :: YRCSGEOM(:) => NULL()
!!TYPE(TCSGEOM), POINTER :: YRCSGEOM_NB => NULL()

! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

SUBROUTINE DEALLO_TCSGEOM(YDCSGEOM)

USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

TYPE(TCSGEOM), POINTER, INTENT(INOUT) :: YDCSGEOM(:)
INTEGER(KIND=JPIM) :: J
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMCSGEOM:DEALLO_TCSGEOM',0,ZHOOK_HANDLE)

IF (ASSOCIATED(YDCSGEOM)) THEN
  DO J = 1, SIZE(YDCSGEOM)
    NULLIFY(YDCSGEOM(J)%RCOLON)
    NULLIFY(YDCSGEOM(J)%RSILON)
    NULLIFY(YDCSGEOM(J)%RINDX )
    NULLIFY(YDCSGEOM(J)%RINDY )
    NULLIFY(YDCSGEOM(J)%RATATH)
    NULLIFY(YDCSGEOM(J)%RATATX)
  ENDDO
  DEALLOCATE(YDCSGEOM)
ENDIF

IF (LHOOK) CALL DR_HOOK('YOMCSGEOM:DEALLO_TCSGEOM',1,ZHOOK_HANDLE)

END SUBROUTINE DEALLO_TCSGEOM

!     ------------------------------------------------------------------

SUBROUTINE DEALLO_TCSGEOM_NB(YDCSGEOM)

USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
USE DEALLOCATE_IF_ASSOCIATED_MOD, ONLY : DEALLOCATE_IF_ASSOCIATED

TYPE(TCSGEOM), INTENT(INOUT) :: YDCSGEOM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOMCSGEOM:DEALLO_TCSGEOM_NB',0,ZHOOK_HANDLE)

CALL DEALLOCATE_IF_ASSOCIATED(YDCSGEOM%RCOLON)
CALL DEALLOCATE_IF_ASSOCIATED(YDCSGEOM%RSILON)
CALL DEALLOCATE_IF_ASSOCIATED(YDCSGEOM%RINDX )
CALL DEALLOCATE_IF_ASSOCIATED(YDCSGEOM%RINDY )
CALL DEALLOCATE_IF_ASSOCIATED(YDCSGEOM%RATATH)
CALL DEALLOCATE_IF_ASSOCIATED(YDCSGEOM%RATATX)

IF (LHOOK) CALL DR_HOOK('YOMCSGEOM:DEALLO_TCSGEOM_NB',1,ZHOOK_HANDLE)

END SUBROUTINE DEALLO_TCSGEOM_NB

!     ------------------------------------------------------------------

END MODULE YOMCSGEOM

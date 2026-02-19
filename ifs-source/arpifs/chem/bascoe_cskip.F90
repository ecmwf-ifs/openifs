! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_CSKIP( CD , KUNIT )
!--------------------------------------------------------------------
!   ... Skip lines starting with char cd in opened file kunit
!--------------------------------------------------------------------
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARKIND1           , ONLY : JPIM,  JPRB

      IMPLICIT NONE
!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
  CHARACTER(LEN=*), INTENT(IN)   :: CD
  INTEGER(KIND=JPIM), INTENT(IN) :: KUNIT

  REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE
  CHARACTER(LEN=LEN(CD))         :: CL_C

  IF (LHOOK) CALL DR_HOOK('BASCOE_CSKIP',0,ZHOOK_HANDLE )

    DO
       READ(KUNIT,'(a)')CL_C
       IF (CL_C/=CD)  EXIT
    ENDDO
    BACKSPACE(KUNIT)

  IF (LHOOK) CALL DR_HOOK('BASCOE_CSKIP',1,ZHOOK_HANDLE )
END SUBROUTINE BASCOE_CSKIP

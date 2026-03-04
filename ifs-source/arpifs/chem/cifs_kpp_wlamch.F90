! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction


!--------------------------------------------------------------
SUBROUTINE CIFS_KPP_WLAMCH( PROUNDOFF, CD_Char )
!--------------------------------------------------------------
!     returns ZEPSilon machine
!     after LAPACK
!     replace this by the function from the optimized LAPACK implementation:
!          CALL SLAMCH('E') or CALL DLAMCH('E')
!--------------------------------------------------------------
!      USE cifs_kpp_Precision
USE PARKIND1  , ONLY : JPIM,JPRB,JPRD
USE YOMLUN   , ONLY : NULOUT 
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

REAL(KIND=JPRD), INTENT(OUT) :: PROUNDOFF
CHARACTER, INTENT(IN)     :: CD_Char
INTEGER(KIND=JPIM)    :: i
REAL(KIND=JPRD)       ::  ZEPS
REAL(KIND=JPRD)       ::  ZSUMA
REAL(KIND=JPRD), PARAMETER  ::  ZONE=1.0_JPRD, ZHALF=0.5_JPRD
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('CIFS_KPP_WLAMCH',0,ZHOOK_HANDLE )

  ZEPS = ZHALF**(16)
  DO i = 17, 120
    ZEPS = ZEPS*ZHALF
    CALL WLAMCH_ADD(ZONE,ZEPS,ZSUMA)
    IF (ZSUMA <= ZONE) CYCLE
  ENDDO
  IF (ZSUMA <= ZONE) THEN
    ZEPS = ZEPS*2
    i = i-1 
  ELSE
    WRITE(NULOUT,*) 'ERROR IN WLAMCH ZEPS < ',ZEPS
    ! VH nevertheless, use this tiny number!
    PROUNDOFF = ZEPS
    IF (LHOOK) CALL DR_HOOK('CIFS_KPP_WLAMCH',1,ZHOOK_HANDLE )
    RETURN
  ENDIF

PROUNDOFF = ZEPS

IF (LHOOK) CALL DR_HOOK('CIFS_KPP_WLAMCH',1,ZHOOK_HANDLE )

CONTAINS

SUBROUTINE WLAMCH_ADD( PA, PB, PSUMA )
 USE PARKIND1  , ONLY : JPRD
 USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

REAL(KIND=JPRD),INTENT(IN) :: PA, PB
REAL(KIND=JPRD),INTENT(OUT):: PSUMA
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WLAMCH_ADD',0,ZHOOK_HANDLE )

PSUMA = PA + PB

IF (LHOOK) CALL DR_HOOK('WLAMCH_ADD',1,ZHOOK_HANDLE )
END SUBROUTINE WLAMCH_ADD

END SUBROUTINE CIFS_KPP_WLAMCH


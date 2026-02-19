! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE MODIFY_WV_CONTINUUM(NWVCONTINUUM)

! MODIFY_WV_CONTINUUM - Adjust the shortwave continuum coefficients
!
! PURPOSE
! -------
!   The default water vapour continuum model in SRTM is MT_CKD 2.5,
!   but some measurement programmes, notably from the CAVIAR project
!   (Shine et al., J. Mol. Spectrosc., 2016) suggest a much stronger
!   absorption in the near infrared. This routine provides the option
!   to implement an approximate scaling of the shortwave continuum
!   coefficients to match the CAVIAR continuum. Further details on the
!   impact were provided by Hogan et al. (2017, ECMWF Tech. Memo. 816).
!
! INTERFACE
! ---------
!   This routine is called from SUECRAD. If its argument is 0, it does
!   nothing so that the default SRTM continuum is used. If its
!   argument is 1 then it implements the CAVIAR continuum by scaling
!   coefficients within the relevant SRTM modules.
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF
!   Original: 2018-02-21
!
! MODIFICATIONS
! -------------
!
! -----------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
! Load the coefficients for each relevant shortwave band
USE YOESRTA16, ONLY : SELFREF16 => SELFREF, FORREF16 => FORREF
USE YOESRTA17, ONLY : SELFREF17 => SELFREF, FORREF17 => FORREF
USE YOESRTA18, ONLY : SELFREF18 => SELFREF, FORREF18 => FORREF
USE YOESRTA19, ONLY : SELFREF19 => SELFREF, FORREF19 => FORREF
USE YOESRTA20, ONLY : SELFREF20 => SELFREF, FORREF20 => FORREF
USE YOESRTA21, ONLY : SELFREF21 => SELFREF, FORREF21 => FORREF
USE YOESRTA22, ONLY : SELFREF22 => SELFREF, FORREF22 => FORREF
USE YOESRTA23, ONLY : SELFREF23 => SELFREF, FORREF23 => FORREF
USE YOESRTA29, ONLY : SELFREF29 => SELFREF, FORREF29 => FORREF

IMPLICIT NONE

! CAVIAR continuum enhancements
REAL(KIND=JPRB), PARAMETER :: SELF_ENH16(16) = (/ 2.42,  2.42,  2.91,  2.91,  2.52, &
     &  2.52,  2.53,  2.53,  2.51,  2.51,  2.51,  2.51,  2.58,  2.58,  2.58,  2.58 /)
REAL(KIND=JPRB), PARAMETER :: FORE_ENH16(16) = (/ 3.38,  3.38,  3.19,  3.19,  1.21,  &
     &  1.21,  1.09,  1.09,  1.07,  1.07,  1.07,  1.07,  1.12,  1.12,  1.12,  1.12 /)
REAL(KIND=JPRB), PARAMETER :: SELF_ENH17(16) = (/ 2.18,  1.40,  1.09,  1.19,  1.02,  1.00, &
     &  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
REAL(KIND=JPRB), PARAMETER :: FORE_ENH17(16) = (/ 3.17,  3.40,  1.66,  1.00,  1.00,  1.00, &
     &  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
REAL(KIND=JPRB), PARAMETER :: SELF_ENH18(16) = (/ 9.67, 12.36,  9.22,  3.71,  1.12,  1.12, &
     &  0.53,  0.53,  0.49,  0.49,  0.49,  0.49,  0.35,  0.35,  0.35,  0.35 /)
REAL(KIND=JPRB), PARAMETER :: FORE_ENH18(16) = (/ 38.90, 15.37, 16.55, 14.81,  4.91,  4.91, &
     &  2.59,  2.59,  2.21,  2.21,  2.21,  2.21,  1.77,  1.77,  1.77,  1.77 /)
REAL(KIND=JPRB), PARAMETER :: SELF_ENH19(16) = (/ 28.53, 26.12, 19.14, 10.12,  3.69, &
     &  3.69,  1.63,  1.63,  2.52,  2.52,  2.52,  2.52,  2.40,  2.40,  2.40,  2.40 /)
REAL(KIND=JPRB), PARAMETER :: FORE_ENH19(16) = (/ 11.66,  9.78,  9.57,  9.55,  4.96, &
     &  4.96,  2.68,  2.68,  2.61,  2.61,  2.61,  2.61,  2.37,  2.37,  2.37,  2.37 /)
REAL(KIND=JPRB), PARAMETER :: SELF_ENH20(16) = (/ 4.93,  2.76,  1.23,  0.66,  1.41, &
     &  1.11,  1.07,  1.03,  1.03,  1.03,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
REAL(KIND=JPRB), PARAMETER :: FORE_ENH20(16) = (/ 24.16,  9.04,  2.73,  2.17,  1.05, &
     &  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
REAL(KIND=JPRB), PARAMETER :: SELF_ENH21(16) = (/ 9.70,  4.56,  0.99,  1.21,  1.37, &
     &  1.25,  0.94,  0.99,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
REAL(KIND=JPRB), PARAMETER :: FORE_ENH21(16) = (/ 50.84, 19.27,  1.49,  1.16,  0.97, &
     &  1.64,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
REAL(KIND=JPRB), PARAMETER :: SELF_ENH22(16) = (/ 3.37,  3.37,  3.37,  3.37,  3.37, &
     &  3.37,  3.37,  3.37,  1.42,  1.42,  1.42,  1.42,  1.42,  1.42,  1.42,  1.42 /)
REAL(KIND=JPRB), PARAMETER :: FORE_ENH22(16) = (/ 12.31, 12.31, 12.31, 12.31, 12.31, &
     &  12.31, 12.31, 12.31,  3.20,  3.20,  3.20,  3.20,  3.20,  3.20,  3.20,  3.20 /)
REAL(KIND=JPRB), PARAMETER :: SELF_ENH23(16) = (/ 1.00,  1.00,  1.19,  1.19,  1.65, &
     &  1.46,  1.32,  1.07,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
REAL(KIND=JPRB), PARAMETER :: FORE_ENH23(16) = (/ 1.04,  1.04,  1.08,  1.08,  1.12, &
     &  1.10,  1.18,  1.06,  1.01,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
REAL(KIND=JPRB), PARAMETER :: SELF_ENH29(16) = (/ 1.70,  1.00,  1.00,  1.03,  1.19,  &
     &  1.19, 1.43,  1.43,  1.30,  1.30,  1.33,  1.33,  1.28,  1.28,  1.08,  1.23 /)
REAL(KIND=JPRB), PARAMETER :: FORE_ENH29(16) = (/ 107.42,  5.87,  3.26,  2.42,  1.39, &
     &  1.39,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)

INTEGER(KIND=JPIM), INTENT(IN) :: NWVCONTINUUM

INTEGER(KIND=JPIM) :: JG

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MODIFY_WV_CONTINUUM',0,ZHOOK_HANDLE)
! -----------------------------------------------------------------------

IF (NWVCONTINUUM == 1) THEN
  ! Apply CAVIAR continuum enhancements
  DO JG = 1,16
    FORREF16(:,JG)  = FORREF16(:,JG)  * FORE_ENH16(JG)
    SELFREF16(:,JG) = SELFREF16(:,JG) * SELF_ENH16(JG)
  ENDDO
  DO JG = 1,16
    FORREF17(:,JG)  = FORREF17(:,JG)  * FORE_ENH17(JG)
    SELFREF17(:,JG) = SELFREF17(:,JG) * SELF_ENH17(JG)
  ENDDO
  DO JG = 1,16
    FORREF18(:,JG)  = FORREF18(:,JG)  * FORE_ENH18(JG)
    SELFREF18(:,JG) = SELFREF18(:,JG) * SELF_ENH18(JG)
  ENDDO
  DO JG = 1,16
    FORREF19(:,JG)  = FORREF19(:,JG)  * FORE_ENH19(JG)
    SELFREF19(:,JG) = SELFREF19(:,JG) * SELF_ENH19(JG)
  ENDDO
  DO JG = 1,16
    FORREF20(:,JG)  = FORREF20(:,JG)  * FORE_ENH20(JG)
    SELFREF20(:,JG) = SELFREF20(:,JG) * SELF_ENH20(JG)
  ENDDO
  DO JG = 1,16
    FORREF21(:,JG)  = FORREF21(:,JG)  * FORE_ENH21(JG)
    SELFREF21(:,JG) = SELFREF21(:,JG) * SELF_ENH21(JG)
  ENDDO
  DO JG = 1,16
    FORREF22(:,JG)  = FORREF22(:,JG)  * FORE_ENH22(JG)
    SELFREF22(:,JG) = SELFREF22(:,JG) * SELF_ENH22(JG)
  ENDDO
  DO JG = 1,16
    FORREF23(:,JG)  = FORREF23(:,JG)  * FORE_ENH23(JG)
    SELFREF23(:,JG) = SELFREF23(:,JG) * SELF_ENH23(JG)
  ENDDO
  DO JG = 1,16
    FORREF29(:,JG)  = FORREF29(:,JG)  * FORE_ENH29(JG)
    SELFREF29(:,JG) = SELFREF29(:,JG) * SELF_ENH29(JG)
  ENDDO
ENDIF
  
! -----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MODIFY_WV_CONTINUUM',1,ZHOOK_HANDLE)

END SUBROUTINE MODIFY_WV_CONTINUUM

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

MODULE MODULE_RADTC_MIX

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : NUNDEFLD

IMPLICIT NONE

SAVE

! -------------------------------------------------------------------------
! Define structures, and declare variables, for radiation transmission
! coefficients.
! -------------------------------------------------------------------------


! Meaning of roots 'RAB3C' to 'COR', and 'AC':

! RAB3C: Clear-sky parallel/diffuse reflectivity on layers (code 'RAB3C').
! RAB4C: Clear-sky diffuse transmissivity on layers (code 'RAB4C').
! RAB6C: Clear-sky diffuse reflectivity on layers (code 'RAB6C').
! RAB3N: Cloudy parallel/diffuse reflectivity on layers (code 'RAB3N').
! RAB4N: Cloudy diffuse transmissivity on layers (code 'RAB4N').
! RAB6N: Cloudy diffuse reflectivity on layers (code 'RAB6N').
! RAT1C: Clear-sky parallel transmissivity on layers (code 'RAT1C').
! RAT2C: Clear-sky parallel/diffuse transmissivity on layers (code 'RAT2C')
! RAT3C: Clear-sky parallel/diffuse reflectivity on layers (code 'RAT3C').
! RAT4C: Clear-sky diffuse transmissivity on layers (code 'RAT4C').
! RAT5C: Clear-sky diffuse reflectivity on layers (code 'RAT5C').
! RAT1N: Cloudy parallel transmissivity on layers (code 'RAT1N').
! RAT2N: Cloudy parallel/diffuse transmissivity on layers (code 'RAT2N').
! RAT3N: Cloudy parallel/diffuse reflectivity on layers (code 'RAT3N').
! RAT4N: Cloudy diffuse transmissivity on layers (code 'RAT4N').
! RAT5N: Cloudy diffuse reflectivity on layers (code 'RAT5N').
! COR  : Correction for cloudiness for thermal radiation (code 'COR').
! AC   : Transmission coefficients matrix for thermal radiation
!        for clear sky = Curtis matrix (code 'AC').

! Arrays contain equivalent coefficients from top to given layer.

! -------------------------------------------------------------------------

!      1       TYPE DEFINITION

! Attributes for radiation transmission coefficients RAB3C to COR.
TYPE TYPE_RADTC
INTEGER(KIND=JPIM) :: MRAB3C        ! RAB3C
INTEGER(KIND=JPIM) :: MRAB4C        ! RAB4C
INTEGER(KIND=JPIM) :: MRAB6C        ! RAB6C
INTEGER(KIND=JPIM) :: MRAB3N        ! RAB3N
INTEGER(KIND=JPIM) :: MRAB4N        ! RAB4N
INTEGER(KIND=JPIM) :: MRAB6N        ! RAB6N
INTEGER(KIND=JPIM) :: MRAT1C        ! RAT1C
INTEGER(KIND=JPIM) :: MRAT2C        ! RAT2C
INTEGER(KIND=JPIM) :: MRAT3C        ! RAT3C
INTEGER(KIND=JPIM) :: MRAT4C        ! RAT4C
INTEGER(KIND=JPIM) :: MRAT5C        ! RAT5C
INTEGER(KIND=JPIM) :: MRAT1N        ! RAT1N
INTEGER(KIND=JPIM) :: MRAT2N        ! RAT2N
INTEGER(KIND=JPIM) :: MRAT3N        ! RAT3N
INTEGER(KIND=JPIM) :: MRAT4N        ! RAT4N
INTEGER(KIND=JPIM) :: MRAT5N        ! RAT5N
INTEGER(KIND=JPIM) :: MCOR          ! COR
INTEGER(KIND=JPIM) :: NDIM          ! Total number of fields.
CHARACTER(LEN=5)   :: CLNAME(17)    ! Root of transmission coefficient name.
END TYPE TYPE_RADTC

! -------------------------------------------------------------------------

!      2       VARIABLES DECLARATION

! Information about roots of transmission coefficient names.
TYPE(TYPE_RADTC) :: YM_RADTC

! -------------------------------------------------------------------------

!      3       SET-UP

CONTAINS

!=========================================================================

SUBROUTINE SUPTRTC(LDUSE_RADTC,YDM_RADTC)

! Sets-up YM_RADTC
LOGICAL,INTENT(IN)           :: LDUSE_RADTC ! use transmission coeffs if T
TYPE(TYPE_RADTC),INTENT(OUT) :: YDM_RADTC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODULE_RADTC_MIX:SUPTRTC',0,ZHOOK_HANDLE)

IF (LDUSE_RADTC) THEN

  ! transmission coefficients are used (may be stored in grid-point or Fourier space).

  YDM_RADTC%MRAB3C= 1
  YDM_RADTC%MRAB4C= 2
  YDM_RADTC%MRAB6C= 3
  YDM_RADTC%MRAB3N= 4
  YDM_RADTC%MRAB4N= 5
  YDM_RADTC%MRAB6N= 6
  YDM_RADTC%MRAT1C= 7
  YDM_RADTC%MRAT2C= 8
  YDM_RADTC%MRAT3C= 9
  YDM_RADTC%MRAT4C=10
  YDM_RADTC%MRAT5C=11
  YDM_RADTC%MRAT1N=12
  YDM_RADTC%MRAT2N=13
  YDM_RADTC%MRAT3N=14
  YDM_RADTC%MRAT4N=15
  YDM_RADTC%MRAT5N=16
  YDM_RADTC%MCOR  =17
  ! YDM_RADTC%NDIM and NG3SR (YOMRCOEF) should have the same value if LRCOEF=T
  YDM_RADTC%NDIM  =17

  YDM_RADTC%CLNAME( 1)='RAB3C'
  YDM_RADTC%CLNAME( 2)='RAB4C'
  YDM_RADTC%CLNAME( 3)='RAB6C'
  YDM_RADTC%CLNAME( 4)='RAB3N'
  YDM_RADTC%CLNAME( 5)='RAB4N'
  YDM_RADTC%CLNAME( 6)='RAB6N'
  YDM_RADTC%CLNAME( 7)='RAT1C'
  YDM_RADTC%CLNAME( 8)='RAT2C'
  YDM_RADTC%CLNAME( 9)='RAT3C'
  YDM_RADTC%CLNAME(10)='RAT4C'
  YDM_RADTC%CLNAME(11)='RAT5C'
  YDM_RADTC%CLNAME(12)='RAT1N'
  YDM_RADTC%CLNAME(13)='RAT2N'
  YDM_RADTC%CLNAME(14)='RAT3N'
  YDM_RADTC%CLNAME(15)='RAT4N'
  YDM_RADTC%CLNAME(16)='RAT5N'
  YDM_RADTC%CLNAME(17)='COR  '

ELSE

  ! transmission coefficients are not used.

  YDM_RADTC%MRAB3C=NUNDEFLD
  YDM_RADTC%MRAB4C=NUNDEFLD
  YDM_RADTC%MRAB6C=NUNDEFLD
  YDM_RADTC%MRAB3N=NUNDEFLD
  YDM_RADTC%MRAB4N=NUNDEFLD
  YDM_RADTC%MRAB6N=NUNDEFLD
  YDM_RADTC%MRAT1C=NUNDEFLD
  YDM_RADTC%MRAT2C=NUNDEFLD
  YDM_RADTC%MRAT3C=NUNDEFLD
  YDM_RADTC%MRAT4C=NUNDEFLD
  YDM_RADTC%MRAT5C=NUNDEFLD
  YDM_RADTC%MRAT1N=NUNDEFLD
  YDM_RADTC%MRAT2N=NUNDEFLD
  YDM_RADTC%MRAT3N=NUNDEFLD
  YDM_RADTC%MRAT4N=NUNDEFLD
  YDM_RADTC%MRAT5N=NUNDEFLD
  YDM_RADTC%MCOR  =NUNDEFLD

  YDM_RADTC%NDIM  =1

  ! YDM_RADTC%CLNAME is not used in this case.

ENDIF

IF (LHOOK) CALL DR_HOOK('MODULE_RADTC_MIX:SUPTRTC',1,ZHOOK_HANDLE)

END SUBROUTINE SUPTRTC

!=========================================================================

END MODULE MODULE_RADTC_MIX

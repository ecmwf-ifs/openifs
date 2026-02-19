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

SUBROUTINE SUVFE_TESTOPER(CDTAG,CDEX,LDINT_FROM_SURF, &
  & KTYPE,KFLEV_IN,PETA_IN,KTYPE_IN,  &
  & KFLEV_OUT,PETA_OUT, &
  & POPER)

!**** *SUVFE_TESTOPER - basic tests of finite element operator quality

!**   Interface.
!     ----------

!     *CALL* SUVFE_TESTOPER

!     Explicit arguments :
!     --------------------
!      * INPUT:
!        CDTAG                     : name of operator
!        CDEX                      : kind of the test (CST=constant, LIN=linear,
!                                    QUA=guadratic, SIN=periodic function)
!        LDINT_FROM_SURF           : integral operator from surface/top
!        KTYPE                     : TYPE OF OPERATOR:
!                                    KTYPE = -1 INTEGRAL OPERATOR
!                                    KTYPE =  1 FIRTS  ORDER DERIVATIVE OPERATOR
!                                    KTYPE =  2 SECOND ORDER DERIVATIVE OPERATOR
!        KFLEV_IN                  : NUMBER OF INPUT LEVELS OF OPERATOR
!        PETA_IN                   : VALUES OF INPUT ETA LEVELS
!        KTYPE_IN                  : 0 = value , 1=derivative
!        KFLEV_OUT                 : NUMBER OF OUTPUT LEVELS OF OPERATOR
!        PETA_OUT                  : VALUES OF OUTPUT ETA LEVELS
!        POPER                     : VFE OPERATOR (ACCORDING KTYPE)

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ALADIN/LACE documentation on NH dynamics.

!     Author.
!     -------
!        Jozef Vivoda, SHMU/LACE 
!        Original : 2010-09

!     Modifications.
!     --------------
!     F. Vana  14-Jan-2020  Single precision support
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB   ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=15) , INTENT(IN) :: CDTAG
CHARACTER(LEN=3)  , INTENT(IN) :: CDEX
LOGICAL           , INTENT(IN) :: LDINT_FROM_SURF
INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV_IN
REAL   (KIND=JPRB), INTENT(IN) :: PETA_IN(KFLEV_IN)
INTEGER(KIND=JPIM), INTENT(IN) :: KTYPE_IN(KFLEV_IN)
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEV_OUT
REAL   (KIND=JPRB), INTENT(IN) :: PETA_OUT(KFLEV_OUT)
REAL   (KIND=JPRD), INTENT(IN) :: POPER(KFLEV_OUT,KFLEV_IN)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: ZARG, ZBC, ZMAE, ZRMSE

REAL(KIND=JPRD) :: ZFUN(KFLEV_IN)   ! To comply for MATMUL
REAL(KIND=JPRB) :: ZRES(KFLEV_OUT)
REAL(KIND=JPRB) :: ZANA(KFLEV_OUT)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_TESTOPER',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*  1. APPLY THE OPERATOR ON A GIVEN FUNCTION
!   -----------------------------------------

!   ---------------------------------
!*  1.1 Example 1 : constant function
!   ---------------------------------

IF (CDEX=='CST') THEN

  ZBC = 1.0_JPRB
  DO JLEV=1,KFLEV_IN
    IF( KTYPE_IN(JLEV) == 0 )THEN ! value
      ZFUN(JLEV) = 1.0_JPRD
      IF( KTYPE > 0 )THEN
        ZFUN(JLEV) = ZFUN(JLEV) - ZBC
      ENDIF
    ELSE ! first derivative
      ZFUN(JLEV) = 0.0_JPRD
    ENDIF
  ENDDO
  DO JLEV=1,KFLEV_OUT
    IF( KTYPE == -1 )THEN
      IF( .NOT.LDINT_FROM_SURF )THEN
        ZANA(JLEV) = ZBC * PETA_OUT(JLEV)
      ELSE
        ZANA(JLEV) = ZBC * ( 1.0_JPRB - PETA_OUT(JLEV) )
      ENDIF
    ELSEIF( KTYPE == 1 )THEN
      ZANA(JLEV) = 0.0_JPRB
    ELSEIF( KTYPE == 2 )THEN
      ZANA(JLEV) = 0.0_JPRB
    ENDIF
  ENDDO
WRITE(NULOUT,'(A," ZRES :: OPER * 1")') TRIM(CDTAG)

!   -------------------------------
!*  1.2 Example 2 : linear function
!   -------------------------------

ELSEIF (CDEX=='LIN') THEN

  IF( .NOT.LDINT_FROM_SURF )THEN
    ZBC = 2.0_JPRB
  ELSE
    ZBC = 1.0_JPRB
  ENDIF
  DO JLEV=1,KFLEV_IN
    IF( KTYPE_IN(JLEV) == 0 )THEN ! value
      ZFUN(JLEV) = 2.0_JPRD - PETA_IN(JLEV)
      IF( KTYPE > 0 )THEN
        ZFUN(JLEV) = ZFUN(JLEV) - ZBC
      ENDIF
    ELSE ! first derivative
      ZFUN(JLEV) = - 1.0_JPRD
    ENDIF
  ENDDO
  DO JLEV=1,KFLEV_OUT
    IF( KTYPE == -1 )THEN
      ZANA(JLEV) = 2.0_JPRB * PETA_OUT(JLEV) - &
       & 1.0_JPRB/2.0_JPRB*PETA_OUT(JLEV)**2
      IF( LDINT_FROM_SURF )THEN
        ZANA(JLEV) = 3.0_JPRB/2.0_JPRB - ZANA(JLEV)
      ENDIF
    ELSEIF( KTYPE == 1 )THEN
      ZANA(JLEV) = -1.0_JPRB
    ELSEIF( KTYPE == 2 )THEN
      ZANA(JLEV) = 0.0_JPRB
    ENDIF
  ENDDO
  WRITE(NULOUT,'(A," ZRES :: OPER * ( 2 - eta )")') TRIM(CDTAG)

!   ----------------------------------
!*  1.3 Example 3 : quadratic function
!   ----------------------------------

ELSEIF (CDEX=='QUA') THEN

  IF( .NOT.LDINT_FROM_SURF )THEN
    ZBC = 0.0_JPRB
  ELSE
    ZBC = 2.0_JPRB
  ENDIF
  DO JLEV=1,KFLEV_IN
    IF( KTYPE_IN(JLEV) == 0 )THEN ! value
      ZFUN(JLEV) = 2.0_JPRD*PETA_IN(JLEV)**2
      IF( KTYPE > 0 )THEN
        ZFUN(JLEV) = ZFUN(JLEV) - ZBC
      ENDIF
    ELSE ! first derivative
      ZFUN(JLEV) = 4.0_JPRD*PETA_IN(JLEV)
    ENDIF
  ENDDO
  DO JLEV=1,KFLEV_OUT
    IF( KTYPE == -1 )THEN
      ZANA(JLEV) = (2.0_JPRB/3.0_JPRB)*PETA_OUT(JLEV)**3
      IF( LDINT_FROM_SURF )THEN
        ZANA(JLEV) = 2.0_JPRB/3.0_JPRB - ZANA(JLEV)
      ENDIF
    ELSEIF( KTYPE == 1 )THEN
      ZANA(JLEV) = 4.0_JPRB*PETA_OUT(JLEV)
    ELSEIF( KTYPE == 2 )THEN
      ZANA(JLEV) = 4.0_JPRB
    ENDIF
  ENDDO
  WRITE(NULOUT,'(A," ZRES :: OPER * (2 * eta * eta)")') TRIM(CDTAG)

!   ---------------------------------------------------------
!*  1.4 Example 4 : function satisfying BC for x=0,pi
!                   f=0, df/dx=0, d^2f/dx^2=0, int_0^1 f dx=0
!   ---------------------------------------------------------

ELSEIF (CDEX=='SIN') THEN

  IF( .NOT.LDINT_FROM_SURF )THEN
    ZBC = 0.0_JPRB
  ELSE
    ZBC = 0.0_JPRB
  ENDIF
  DO JLEV=1,KFLEV_IN
    ZARG = RPI*PETA_IN(JLEV)
    IF( KTYPE_IN(JLEV) == 0 )THEN ! value
      ZFUN(JLEV) = SIN(ZARG)**3*COS(ZARG)
      IF( KTYPE > 0 )THEN
        ZFUN(JLEV) = ZFUN(JLEV) - ZBC
      ENDIF
    ELSE ! first derivative
      ZFUN(JLEV) = 3.0_JPRD*RPI*SIN(ZARG)**2*COS(ZARG)**2 &
       & - RPI*SIN(ZARG)**4
    ENDIF
  ENDDO
  DO JLEV=1,KFLEV_OUT
    ZARG  = RPI*PETA_OUT(JLEV)
    IF( KTYPE == -1 )THEN
      ZANA(JLEV) = SIN(ZARG)**4/(4.0_JPRB*RPI)
      IF( LDINT_FROM_SURF )THEN
        ZANA(JLEV) = -ZANA(JLEV)
      ENDIF
    ELSEIF( KTYPE == 1 )THEN
      ZANA(JLEV) = 3.0_JPRB*RPI*SIN(ZARG)**2*COS(ZARG)**2 &
       & - RPI*SIN(ZARG)**4
    ELSEIF( KTYPE == 2 )THEN
      ZANA(JLEV) = 6.0_JPRB*RPI*RPI*SIN(ZARG)*COS(ZARG)**3 &
       & - 10.0_JPRB*RPI*RPI*SIN(ZARG)**3*COS(ZARG)
    ENDIF
  ENDDO
  WRITE(NULOUT,'(A," ZRES :: OPER * ( sin(pi*eta)^3*cos(pi*eta) )")') TRIM(CDTAG)

ELSE
  WRITE(NULOUT,'(" SUVFE_TESTOPER : UNDEFINED EXAMPLE ",A)') TRIM(CDEX)
ENDIF

!*  2. LIST THE RESULTS.
!   --------------------

ZRES = MATMUL(POPER,ZFUN)
DO JLEV=1,KFLEV_IN
  WRITE(NULOUT,'(A," TYPE: ",I2," F: ",F15.6," LEV: ",I5)') TRIM(CDTAG), &
   & KTYPE_IN(JLEV), ZFUN(JLEV), JLEV
ENDDO
DO JLEV=1,KFLEV_OUT
  WRITE(NULOUT,'(A," VFE: ",F15.6," ANAL: ",F15.6," LEV: ",I5)') TRIM(CDTAG), &
   & ZRES(JLEV), ZANA(JLEV), JLEV
ENDDO

!*  3. COMPUTE MEAN ABSOLUTE AND RMS ERRORS.
!   ----------------------------------------

ZMAE  = 0.0_JPRB
ZRMSE = 0.0_JPRB

DO JLEV=1,KFLEV_OUT
  ZMAE  = ZMAE  + ABS(ZRES(JLEV) - ZANA(JLEV))
  ZRMSE = ZRMSE + (ZRES(JLEV) - ZANA(JLEV))**2
ENDDO

ZMAE  = ZMAE/REAL(KFLEV_OUT,JPRB)
ZRMSE = SQRT(ZRMSE/REAL(KFLEV_OUT,JPRB))
WRITE(NULOUT,'(A," CALCULATED ERRORS, MAE: ",F15.12," RMSE: ",F15.12)') &
 & TRIM(CDTAG), ZMAE, ZRMSE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_TESTOPER',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE SUVFE_TESTOPER

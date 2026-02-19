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

SUBROUTINE MINV(PAB,KDIMN,KDBA,PZSCRA,PDET1,PTOL,KDIMM,KMODE)

!     Author
!     ------
!       Filip Vana (c) ECMWF, based on previous F77 routine minv.F

!     Modifications.
!     --------------
!     F. Vana  05-Mar-2015  Support for single precision


USE PARKIND1, ONLY : JPRD, JPRB,  JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KDIMN
INTEGER(KIND=JPIM), INTENT(IN)    :: KDBA
INTEGER(KIND=JPIM), INTENT(IN)    :: KDIMM
INTEGER(KIND=JPIM), INTENT(IN)    :: KMODE
REAL(KIND=JPRB),    INTENT(IN)    :: PTOL
REAL(KIND=JPRB),    INTENT(OUT)   :: PDET1
REAL(KIND=JPRB),    INTENT(INOUT) :: PAB(KDBA,KDIMN+KDIMM)
REAL(KIND=JPRB),    INTENT(INOUT) :: PZSCRA(2*KDIMN)  ! not used


!-----------------------------------------------------------------------
REAL(KIND=JPRB)   :: ZRCOND, ZALPHA, ZBETA
REAL(KIND=JPRB)   :: ZSCRATCH(4*KDIMN), ZMAT(KDBA,KDIMN), ZX(KDIMN), ZY(KDBA)
CHARACTER*1        CTRANS
INTEGER(KIND=JPIM) ::  INCX, JL, JC, JSYS
INTEGER(KIND=JPIM) ::  IPVT(KDIMN)
LOGICAL, PARAMETER :: LLDOUBLE = (JPRB == JPRD)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MINV',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

IF (KDBA.NE.KDIMN) CALL ABOR1 ('  ERROR IN MINV -- Matrix MUST be square ')

IF (KDIMM < 0) CALL ABOR1 ('  ERROR IN MINV -- KDIMM MUST BE >= 0  ')

!    Extraction de la matrice ZMAT a factoriser
DO JL = 1,KDBA
  DO JC = 1,KDIMN
    ZMAT(JL,JC) = PAB(JL,JC)
  ENDDO
ENDDO
CALL GECO
IF(ZRCOND <= PTOL) CALL ABOR1(' MINV : MATRIX IS SINGULAR  ')

!    Inversion de ZMAT
CALL GEDI

PDET1 =  1._JPRB 

!    Remplacement de A (ou ZMAT) par son inverse:

IF (KMODE /= 0) THEN
  DO JL = 1,KDBA
    DO JC = 1,KDIMN
      PAB(JL,JC) = ZMAT(JL,JC)
    ENDDO
  ENDDO
ENDIF

!    Resolution des differents systemes lineaires

IF (KDIMM > 0) THEN
  CTRANS = 'N'
  ZALPHA = 1._JPRB
  ZBETA  = 0._JPRB
  INCX = 1
  DO JSYS = 1,KDIMM
    !    Extraction du second membre X
    DO JL = 1,KDBA
      ZX(JL) = PAB(JL,KDIMN+JSYS)
    ENDDO
    ZY = 0._JPRB
    IF (LLDOUBLE) THEN
      CALL DGEMV(CTRANS,KDBA,KDIMN,ZALPHA,ZMAT,KDBA,ZX,INCX,ZBETA,ZY,INCX)
    ELSE
      CALL SGEMV(CTRANS,KDBA,KDIMN,ZALPHA,ZMAT,KDBA,ZX,INCX,ZBETA,ZY,INCX)
    ENDIF
    !    Sauvegarde de la solution
    DO JL = 1,KDBA
      PAB(JL,KDIMN+JSYS) = ZY(JL)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MINV',1,ZHOOK_HANDLE)

CONTAINS
SUBROUTINE GECO

!--- simulate LINPAC routines SGECO/DGECO using LAPACK

INTEGER(KIND=JPIM) ::  J, INFO
INTEGER(KIND=JPIM) ::  IWORK(KDIMN)
REAL(KIND=JPRB) :: ZANORM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MINV:GECO',0,ZHOOK_HANDLE)

ZANORM= 0._JPRB
DO J = 1, KDIMN
  ZANORM = MAX(ZANORM,SUM(ABS(ZMAT(1:KDIMN,J))))
ENDDO

IF (JPRB == JPRD) THEN
  ! Double precision version
  CALL DGETRF (KDIMN,KDIMN,ZMAT,KDBA,IPVT,INFO)
  IF (INFO < 0) CALL ABOR1 ('DGETRF RETURNS NEGATIVE INFO')

  CALL DGECON ('1',KDIMN,ZMAT,KDBA,ZANORM,ZRCOND,ZSCRATCH,IWORK,INFO)
  IF (INFO /= 0) CALL ABOR1 ('DGECON RETURNS NON-ZERO INFO')
ELSE
  ! Single precision version
  CALL SGETRF (KDIMN,KDIMN,ZMAT,KDBA,IPVT,INFO)
  IF (INFO < 0) CALL ABOR1 ('SGETRF RETURNS NEGATIVE INFO')

  CALL SGECON ('1',KDIMN,ZMAT,KDBA,ZANORM,ZRCOND,ZSCRATCH,IWORK,INFO)
  IF (INFO /= 0) CALL ABOR1 ('SGECON RETURNS NON-ZERO INFO')
ENDIF
IF (LHOOK) CALL DR_HOOK('MINV:GECO',1,ZHOOK_HANDLE)
END SUBROUTINE GECO

SUBROUTINE GEDI

!--- simulate LINPAC routines SGEDI/DGEDI using LAPACK

INTEGER(KIND=JPIM) :: INFO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MINV:GEDI',0,ZHOOK_HANDLE)

IF (JPRB == JPRD) THEN
  CALL DGETRI (KDIMN,ZMAT,KDBA,IPVT,ZSCRATCH,SIZE(ZSCRATCH),INFO)
ELSE
  CALL SGETRI (KDIMN,ZMAT,KDBA,IPVT,ZSCRATCH,SIZE(ZSCRATCH),INFO)
ENDIF
IF (INFO /= 0) CALL ABOR1 ('DGETRI RETURNS NON-ZERO INFO')
      
IF (LHOOK) CALL DR_HOOK('MINV:GEDI',1,ZHOOK_HANDLE)
END SUBROUTINE GEDI
      
END SUBROUTINE MINV

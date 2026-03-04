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

SUBROUTINE MINV_8(PAB,KDIMN,KDBA,PZSCRA,PDET1,PTOL,KDIMM,KMODE)

!     Author
!     ------
!       Filip Vana (c) ECMWF, based on previous F77 routine minv.F
!            Double precision only version possibly also being called 
!            from single precision code.

!     Modifications.
!     --------------


USE PARKIND1, ONLY : JPRD, JPRB,  JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

!-----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KDIMN
INTEGER(KIND=JPIM), INTENT(IN)    :: KDBA
INTEGER(KIND=JPIM), INTENT(IN)    :: KDIMM
INTEGER(KIND=JPIM), INTENT(IN)    :: KMODE
REAL(KIND=JPRD),    INTENT(IN)    :: PTOL
REAL(KIND=JPRD),    INTENT(OUT)   :: PDET1
REAL(KIND=JPRD),    INTENT(INOUT) :: PAB(KDBA,KDIMN+KDIMM)
REAL(KIND=JPRD),    INTENT(INOUT) :: PZSCRA(2*KDIMN)  ! not used


!-----------------------------------------------------------------------
REAL(KIND=JPRD)   :: ZRCOND, ZALPHA, ZBETA
REAL(KIND=JPRD)   :: ZSCRATCH(4*KDIMN), ZMAT(KDBA,KDIMN), ZX(KDIMN), ZY(KDBA)
CHARACTER*1        CTRANS
INTEGER(KIND=JPIM) ::  INCX, JL, JC, JSYS
INTEGER(KIND=JPIM) ::  IPVT(KDIMN)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MINV_8',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

IF (KDBA.NE.KDIMN) CALL ABOR1 ('  ERROR IN MINV_8 -- Matrix MUST be square ')

IF (KDIMM < 0) CALL ABOR1 ('  ERROR IN MINV_8 -- KDIMM MUST BE >= 0  ')

!    Extraction de la matrice ZMAT a factoriser
DO JL = 1,KDBA
  DO JC = 1,KDIMN
    ZMAT(JL,JC) = PAB(JL,JC)
  ENDDO
ENDDO
CALL GECO_8
IF(ZRCOND <= PTOL) CALL ABOR1(' MINV_8 : MATRIX IS SINGULAR  ')

!    Inversion de ZMAT
CALL GEDI_8

PDET1 =  1._JPRD 

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
  ZALPHA = 1._JPRD
  ZBETA  = 0._JPRD
  INCX = 1
  DO JSYS = 1,KDIMM
    !    Extraction du second membre X
    DO JL = 1,KDBA
      ZX(JL) = PAB(JL,KDIMN+JSYS)
    ENDDO
    Zy = 0._JPRD
    CALL DGEMV(CTRANS,KDBA,KDIMN,ZALPHA,ZMAT,KDBA,ZX,INCX,ZBETA,ZY,INCX)

    !    Sauvegarde de la solution
    DO JL = 1,KDBA
      PAB(JL,KDIMN+JSYS) = ZY(JL)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MINV_8',1,ZHOOK_HANDLE)

CONTAINS
SUBROUTINE GECO_8

!--- simulate LINPAC routines SGECO/DGECO using LAPACK

INTEGER(KIND=JPIM) ::  J, INFO
INTEGER(KIND=JPIM) ::  IWORK(KDIMN)
REAL(KIND=JPRD) :: ZANORM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MINV_8:GECO',0,ZHOOK_HANDLE)

ZANORM= 0._JPRD
DO J = 1, KDIMN
  ZANORM = MAX(ZANORM,SUM(ABS(ZMAT(1:KDIMN,J))))
ENDDO

CALL DGETRF (KDIMN,KDIMN,ZMAT,KDBA,IPVT,INFO)
IF (INFO < 0) CALL ABOR1 ('DGETRF RETURNS NEGATIVE INFO')

CALL DGECON ('1',KDIMN,ZMAT,KDBA,ZANORM,ZRCOND,ZSCRATCH,IWORK,INFO)
IF (INFO /= 0) CALL ABOR1 ('DGECON RETURNS NON-ZERO INFO')

IF (LHOOK) CALL DR_HOOK('MINV_8:GECO',1,ZHOOK_HANDLE)
END SUBROUTINE GECO_8

SUBROUTINE GEDI_8

!--- simulate LINPAC routines SGEDI/DGEDI using LAPACK

INTEGER(KIND=JPIM) :: INFO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MINV_8:GEDI',0,ZHOOK_HANDLE)

CALL DGETRI (KDIMN,ZMAT,KDBA,IPVT,ZSCRATCH,SIZE(ZSCRATCH),INFO)
IF (INFO /= 0) CALL ABOR1 ('DGETRI RETURNS NON-ZERO INFO')
      
IF (LHOOK) CALL DR_HOOK('MINV_8:GEDI',1,ZHOOK_HANDLE)
END SUBROUTINE GEDI_8
      
END SUBROUTINE MINV_8

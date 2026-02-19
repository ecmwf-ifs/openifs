! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FSPGL_INT_MOD
CONTAINS
SUBROUTINE FSPGL_INT(KM,KMLOC,KF_UV,KF_SCALARS,KF_SCDERS,KF_OUT_LT,&
 & FSPGL_PROC,KFLDPTRUV,KFLDPTRSC)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_TRANS       ,ONLY : FOUBUF_IN, LDIVGP, LVORGP
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_DISTR       ,ONLY : D
USE TPM_FIELDS      ,ONLY : F
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV,KF_SCALARS,KF_SCDERS,KF_OUT_LT
EXTERNAL  FSPGL_PROC
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)

!
! ZFIELD 2nd dimension is extended from 0 to R%NDGL+1, while only 1 to R%NDGL
! is given from the north/south transforms, and only 1 to R%NDGL rows will be
! passed to the east/west transforms.
! the 2 extra rows are used inside the model Fourier space computations
! (outside the transform package - see FSPGLH in Arpege/IFS).
!
REAL(KIND=JPRB) :: ZFIELD(2*KF_OUT_LT,0:R%NDGL+1)


INTEGER(KIND=JPIM) :: ISL, IGLS, JFLD, JGL ,IPROC,  IPROCS
INTEGER(KIND=JPIM) :: IPTRU,IST,J
INTEGER(KIND=JPIM) :: IDGNH,IDGL
INTEGER(KIND=JPIM) :: ISTAN(R%NDGNH),ISTAS(R%NDGNH)
INTEGER(KIND=JPIM) :: IFLDPTRUV(KF_UV),IFLDPTRSC(KF_SCALARS)
!     ------------------------------------------------------------------

IF(PRESENT(KFLDPTRUV)) THEN
  IFLDPTRUV(:) = KFLDPTRUV(1:KF_UV)
  IFLDPTRSC(:) = KFLDPTRSC(1:KF_SCALARS)
ELSE
  DO J=1,KF_UV
    IFLDPTRUV(J) = J
  ENDDO
  DO J=1,KF_SCALARS
    IFLDPTRSC(J) = J
  ENDDO
ENDIF

ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IDGNH = R%NDGNH
IDGL = R%NDGL
DO JGL=ISL,IDGNH
  IPROC = D%NPROCL(JGL)
  ISTAN(JGL) = (D%NSTAGT0B(IPROC) + D%NPNTGTB1(KMLOC,JGL))*2*KF_OUT_LT
  IGLS = IDGL+1-JGL
  IPROCS = D%NPROCL(IGLS)
  ISTAS(JGL) = (D%NSTAGT0B(IPROCS) + D%NPNTGTB1(KMLOC,IGLS))*2*KF_OUT_LT
ENDDO

DO JGL=ISL,IDGNH
  IGLS = IDGL+1-JGL
  DO JFLD=1,2*KF_OUT_LT
    ZFIELD(JFLD,JGL)  = FOUBUF_IN(ISTAN(JGL)+JFLD)
    ZFIELD(JFLD,IGLS) = FOUBUF_IN(ISTAS(JGL)+JFLD)
  ENDDO
ENDDO

IST = 1
IF(LVORGP) THEN
  IST = IST+2*KF_UV
ENDIF
IF(LDIVGP) THEN
  IST = IST+2*KF_UV
ENDIF
IPTRU = IST




CALL FSPGL_PROC(KM,ISL,IDGL,KF_OUT_LT,F%R1MU2,ZFIELD,&
 &   IPTRU,KF_UV,KF_SCALARS,&
 &   IFLDPTRUV)

DO JGL=ISL,IDGNH
  IGLS = IDGL+1-JGL
!OCL      NOVREC
  DO JFLD=1,2*KF_OUT_LT
    FOUBUF_IN(ISTAN(JGL)+JFLD) = ZFIELD(JFLD,JGL)
    FOUBUF_IN(ISTAS(JGL)+JFLD) = ZFIELD(JFLD,IGLS)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FSPGL_INT
END MODULE FSPGL_INT_MOD

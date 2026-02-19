! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIR_CTLAD_MOD
CONTAINS
SUBROUTINE LTDIR_CTLAD(KF_FS,KF_UV,KF_SCALARS, &
 & PSPVOR,PSPDIV,PSPSCALAR, &
 & PSPSC3A,PSPSC3B,PSPSC2, &
 & KFLDPTRUV,KFLDPTRSC)

!**** *LTDIR_CTLAD* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL LTDIR_CTLAD(...)

!     Explicit arguments :
!     --------------------
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : LALLOPERM
!USE TPM_DIM
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D

USE LTDIRAD_MOD     ,ONLY : LTDIRAD
USE TRMTOL_MOD      ,ONLY : TRMTOL
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)   :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)   :: KFLDPTRSC(:)

INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILED2

!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

CALL GSTATS(105,0)
IBLEN = D%NLENGT0B*2*KF_FS
IF (ALLOCATED(FOUBUF_IN)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF_IN)) THEN
    DEALLOCATE(FOUBUF_IN)
    ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
ENDIF
IF (ALLOCATED(FOUBUF)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
    DEALLOCATE(FOUBUF)
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF(MAX(1,IBLEN)))
ENDIF

! Direct Legendre transform

ILED2 = 2*KF_FS
CALL GSTATS(1646,0)
IF(KF_FS > 0) THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL LTDIRAD(IM,JM,KF_FS,KF_UV,KF_SCALARS,ILED2, &
     & PSPVOR,PSPDIV,PSPSCALAR,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & KFLDPTRUV,KFLDPTRSC)
  ENDDO
!$OMP END PARALLEL DO
ENDIF
CALL GSTATS(1646,1)

CALL GSTATS(105,1)

CALL GSTATS(181,0)
CALL TRMTOL(FOUBUF,FOUBUF_IN,2*KF_FS)
CALL GSTATS(181,1)
IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF)
!     ------------------------------------------------------------------

END SUBROUTINE LTDIR_CTLAD
END MODULE LTDIR_CTLAD_MOD

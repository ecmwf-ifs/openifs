! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SPNORM_CTL_MOD
CONTAINS
SUBROUTINE SPNORM_CTL(PSPEC,KFLD,KFLD_G,KVSET,KMASTER,PMET,PNORM)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D, MYPROC, MYSETV

USE SPNORMD_MOD     ,ONLY : SPNORMD
USE SPNORMC_MOD     ,ONLY : SPNORMC
!

IMPLICIT NONE

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KMASTER
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PMET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PNORM(:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFLD,KFLD_G
INTEGER(KIND=JPIM) :: IVSET(KFLD_G)
REAL(KIND=JPRB)    :: ZMET(0:R%NSMAX)
REAL(KIND=JPRB)    :: ZSM(KFLD,D%NUMP)
REAL(KIND=JPRB)    :: ZGM(KFLD_G,0:R%NSMAX)

!     ------------------------------------------------------------------

IF(PRESENT(KVSET)) THEN
  IVSET(:) = KVSET(:)
ELSE
  IVSET(:) = MYSETV
ENDIF

IF(PRESENT(PMET)) THEN
  ZMET(:) = PMET(:)
ELSE
  ZMET(:) = 1.0_JPRB
ENDIF

CALL SPNORMD(PSPEC,KFLD,ZMET,ZSM)

CALL SPNORMC(ZSM,KFLD_G,IVSET,KMASTER,R%NSMAX,ZGM)

IF(MYPROC == KMASTER) THEN
  PNORM(1:KFLD_G) = SUM(ZGM,DIM=2)
  PNORM(1:KFLD_G) = SQRT(PNORM(1:KFLD_G))
ENDIF
!     ------------------------------------------------------------------

END SUBROUTINE SPNORM_CTL
END MODULE SPNORM_CTL_MOD

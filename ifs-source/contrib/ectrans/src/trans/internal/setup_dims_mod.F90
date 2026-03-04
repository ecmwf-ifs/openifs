! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SETUP_DIMS_MOD
CONTAINS
SUBROUTINE SETUP_DIMS

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_FLT         ,ONLY : S
!

IMPLICIT NONE

INTEGER(KIND=JPIM) :: JM,JN,ISPOLEG

!     ------------------------------------------------------------------

ISPOLEG = 0
DO JM=0,R%NSMAX
  DO JN=JM,R%NTMAX+1
    ISPOLEG = ISPOLEG+1
  ENDDO
ENDDO
R%NSPOLEG = ISPOLEG

R%NSPEC_G = (R%NSMAX+1)*(R%NSMAX+2)/2
R%NSPEC2_G = R%NSPEC_G*2

R%NDGNH = (R%NDGL+1)/2

R%NLEI1 = R%NSMAX+4+MOD(R%NSMAX+4+1,2)
R%NLEI3 = R%NDGNH+MOD(R%NDGNH+2,2)
IF (S%LSOUTHPNM) R%NLEI3=2*R%NLEI3 

R%NLED3 = R%NTMAX+2+MOD(R%NTMAX+3,2)
R%NLED4 = R%NTMAX+3+MOD(R%NTMAX+4,2)

!     ------------------------------------------------------------------

END SUBROUTINE SETUP_DIMS
END MODULE SETUP_DIMS_MOD

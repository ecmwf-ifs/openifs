! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRE_SULEG_MOD
CONTAINS
SUBROUTINE PRE_SULEG
USE PARKIND1  ,ONLY : JPRD, JPIM
USE PARKIND2  ,ONLY : JPRH
USE TPM_GEN
USE TPM_DIM
USE TPM_CONSTANTS
USE TPM_DISTR
USE TPM_FIELDS

INTEGER(KIND=JPIM) :: INM, IM, ICOUNT,JMLOC,JN
LOGICAL :: LLP1,LLP2


LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1

ICOUNT = 0
DO JMLOC=1,D%NUMP
  IM = D%MYMS(JMLOC)
  DO JN=IM,R%NTMAX+2
    ICOUNT = ICOUNT+1
  ENDDO
ENDDO

ALLOCATE(F%REPSNM(ICOUNT))
IF (LLP2) WRITE(NOUT,9) 'F%REPSNM  ',SIZE(F%REPSNM ),SHAPE(F%REPSNM )
ALLOCATE(F%RN(-1:R%NTMAX+3))
IF (LLP2) WRITE(NOUT,9) 'F%RN      ',SIZE(F%RN     ),SHAPE(F%RN     ) 
ALLOCATE(F%RLAPIN(-1:R%NSMAX+2))
IF (LLP2) WRITE(NOUT,9) 'F%RLAPIN  ',SIZE(F%RLAPIN ),SHAPE(F%RLAPIN ) 
ALLOCATE(F%NLTN(-1:R%NTMAX+3))
IF (LLP2) WRITE(NOUT,9) 'F%NLTN    ',SIZE(F%NLTN ),SHAPE(F%NLTN ) 

ICOUNT = 0
DO JMLOC=1,D%NUMP
  IM = D%MYMS(JMLOC)
  DO JN=IM,R%NTMAX+2
    ICOUNT = ICOUNT+1
    F%REPSNM(ICOUNT) = SQRT(REAL(JN*JN-IM*IM,JPRD)/&
     &REAL(4*JN*JN-1,JPRD))
  ENDDO
ENDDO

DO JN=-1,R%NTMAX+3
  F%RN(JN) = REAL(JN,JPRD)
  F%NLTN(JN) = R%NTMAX+2-JN
ENDDO
F%RLAPIN(:)  = 0.0_JPRD
F%RLAPIN(0)  = 0.0_JPRD
F%RLAPIN(-1) = 0.0_JPRD
DO JN=1,R%NSMAX+2
  F%RLAPIN(JN)=-(REAL(RA,JPRD)*REAL(RA,JPRD)/REAL(JN*(JN+1),JPRD))
ENDDO

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END SUBROUTINE PRE_SULEG
END MODULE PRE_SULEG_MOD

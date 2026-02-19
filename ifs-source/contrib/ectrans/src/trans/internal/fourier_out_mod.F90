! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_OUT_MOD
CONTAINS
SUBROUTINE FOURIER_OUT(PREEL,KFIELDS,KGL)

!**** *FOURIER_OUT* - Copy fourier data from local array to buffer

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUT(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_GEOMETRY    ,ONLY : G
!

IMPLICIT NONE

REAL(KIND=JPRB),   INTENT(IN) :: PREEL(:,:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL

INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,IR,II,ISTA

!     ------------------------------------------------------------------

IGLG = D%NPTRLS(MYSETW)+KGL-1
DO JM=0,G%NMEN(IGLG)
  IPROC = D%NPROCM(JM)
  IR    = 2*JM+1+D%NSTAGTF(KGL)
  II    = 2*JM+2+D%NSTAGTF(KGL)
  ISTA  = (D%NSTAGT1B(D%MSTABF(IPROC))+D%NPNTGTB0(JM,KGL))*2*KFIELDS
  DO JF=1,KFIELDS
    FOUBUF_IN(ISTA+2*JF-1) = PREEL(JF,IR)
    FOUBUF_IN(ISTA+2*JF  ) = PREEL(JF,II)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_OUT
END MODULE FOURIER_OUT_MOD


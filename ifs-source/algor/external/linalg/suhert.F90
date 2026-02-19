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

SUBROUTINE SUHERT(KLX,KVX,KVXS,PD,PEI,PES,PA,PB,PG)

!**** *SUHERT*   - Decomposition of a set of tridiagonal matrices.

!     Purpose.    Decomposition of a set of tridiagonal matrices into
!     --------    products of two bidiagonal triangular matrices.

!    Example, for KLX=5, matrix number l:

!    I PD(l,1)  PES(l,1)   0         0        0     I
!    I PEI(l,1) PD(l,2)  PES(l,2)    0        0     I
!    I   0      PEI(l,2) PD(l,3)  PES(l,3)    0     I
!    I   0        0      PEI(l,3) PD(l,4)  PES(l,4) I
!    I   0        0        0      PEI(l,4) PD(l,5)  I

!                      EQUAL

!    I PA(l,1)    0       0       0       0    I
!    I PB(l,1) PA(l,2)    0       0       0    I
!    I    0    PB(l,2) PA(l,3)    0       0    I
!    I    0       0    PB(l,3) PA(l,4)    0    I
!    I    0       0       0    PB(l,4) PA(l,5) I

!                      TIMES

!    I   1   PG(l,1)    0       0       0    I
!    I   0      1    PG(l,2)    0       0    I
!    I   0      0       1    PG(l,3)    0    I
!    I   0      0       0       1    PG(l,4) I
!    I   0      0       0       0       1    I

!**   Interface.
!     ----------
!        *CALL* *SUHERT(KLX,KVX,KVXS,PD,PEI,PES,PA,PB,PG)

!        Explicit arguments :
!        --------------------
!         KLX:            - Dimension of the system.                (input)
!         KVX             - Number of matrices to be factorised.    (input)
!         KVXS            - Surdimension associated to KVX.         (input)
!         PD,PEI,PES      - non-zero diagonals of the initial       (input)
!                           matrix (see figure above).
!                           All minor determinants of this
!                           matrix must be non zero.
!         PA,PB:          - non-zero diagonals of the lower         (output)
!                           triangular matrix (see figure above).
!         PG:             - non-zero diagonals of the upper         (output)
!                           triangular matrix (see figure above).

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD: DECEMBER 2003.

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!      ----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVXS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEI(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PES(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PA(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PG(KVXS,KLX) 

!      ----------------------------------------------------------------

INTEGER(KIND=JPIM) :: IL, JL, JV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHERT',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

!*       1.    COMPUTATION OF PA,PB,PG.
!              ------------------------

IF (KLX >= 2) THEN

  DO JV=1,KVX
    PB(JV,1)=PEI(JV,1)
    PA(JV,1)= PD(JV,1)
    PG(JV,1)=PES(JV,1)/PA(JV,1)
  ENDDO

  IF (KLX >= 3) THEN
    DO JL=2,KLX-1
      DO JV=1,KVX
        PB(JV,JL)=PEI(JV,JL)
        PA(JV,JL)= PD(JV,JL)-PB(JV,JL-1)*PG(JV,JL-1)
        PG(JV,JL)=PES(JV,JL)/PA(JV,JL)
      ENDDO
    ENDDO
  ENDIF

  IL=KLX
  DO JV=1,KVX
    PA(JV,IL)=PD(JV,IL)-PB(JV,IL-1)*PG(JV,IL-1)
  ENDDO

ELSEIF (KLX == 1) THEN

  DO JV=1,KVX
    PA(JV,1)=PD(JV,1)
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUHERT',1,ZHOOK_HANDLE)
END SUBROUTINE SUHERT


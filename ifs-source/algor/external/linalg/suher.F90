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

SUBROUTINE SUHER(KLX,KVX,KVXS,PD,PEI,PES,PFI,PFS,PA,PB,PC,PG,PH)

!**** *SUHER*   - Decomposition of a set of pentadiagonal matrixes.

!     Purpose.    Decomposition of a set of pentadiagonal matrixes into
!     --------    products of two tridiagonal triangular matrixes.

!    Example, for KLX=5, matrix number l:

!    I PD(l,1)  PES(l,1) PFS(l,1)    0        0     I
!    I PEI(l,1) PD(l,2)  PES(l,2)    0        0     I
!    I PFI(l,1) PEI(l,2) PD(l,3)  PES(l,3) PFS(l,3) I
!    I   0      PFI(l,2) PEI(l,3) PD(l,4)  PES(l,4) I
!    I   0        0      PFI(l,3) PEI(l,4) PD(l,5)  I

!                      EQUAL

!    I PA(l,1)    0       0       0       0    I
!    I PB(l,1) PA(l,2)    0       0       0    I
!    I PC(l,1) PB(l,2) PA(l,3)    0       0    I
!    I    0    PC(l,2) PB(l,3) PA(l,4)    0    I
!    I    0       0    PC(l,3) PB(l,4) PA(l,5) I

!                      TIMES

!    I   1   PG(l,1) PH(l,1)    0       0    I
!    I   0      1    PG(l,2) PH(l,2)    0    I
!    I   0      0       1    PG(l,3) PH(l,3) I
!    I   0      0       0       1    PG(l,4) I
!    I   0      0       0       0       1    I

!**   Interface.
!     ----------
!        *CALL* *SUHER(KLX,KVX,KVXS,PD,PEI,PES,PFI,PFS,PA,PB,PC,PG,PH)

!        Explicit arguments :
!        --------------------
!         KLX:            - Dimension of the system.                (input)
!         KVX             - Number of matrixes to be factorised.    (input)
!         KVXS            - Surdimension associated to KVX.         (input)
!         PD,PEI,PFI      - non-zero diagonals of the initial       (input)
!           ,PES,PFS:       matrix (see figure above).
!                           All minor determinants of this
!                           matrix must be non zero.
!         PA,PB,PC:       - non-zero diagonals of the lower         (output)
!                           triangular matrix (see figure above).
!         PG,PH:          - non-zero diagonals of the upper         (output)
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
!        K. YESSAD: OCTOBER 1993.

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVXS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEI(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PES(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFI(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFS(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PA(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PC(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PG(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PH(KVXS,KLX) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IL, JL, JV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHER',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF PA,PB,PC,PG,PH.
!              ------------------------------

IF (KLX >= 4) THEN

  DO JV=1,KVX
    PC(JV,1)=PFI(JV,1)
    PB(JV,1)=PEI(JV,1)
    PA(JV,1)= PD(JV,1)
    PG(JV,1)=PES(JV,1)/PA(JV,1)
    PH(JV,1)=PFS(JV,1)/PA(JV,1)
  ENDDO

  DO JV=1,KVX
    PC(JV,2)=PFI(JV,2)
    PB(JV,2)=PEI(JV,2)-PC(JV,1)*PG(JV,1)
    PA(JV,2)= PD(JV,2)-PB(JV,1)*PG(JV,1)
    PG(JV,2)=(PES(JV,2)-PB(JV,1)*PH(JV,1))/PA(JV,2)
    PH(JV,2)=PFS(JV,2)/PA(JV,2)
  ENDDO

  IF (KLX >= 5) THEN
    DO JL=3,KLX-2
      DO JV=1,KVX
        PC(JV,JL)=PFI(JV,JL)
        PB(JV,JL)=PEI(JV,JL)-PC(JV,JL-1)*PG(JV,JL-1)
        PA(JV,JL)= PD(JV,JL)-PC(JV,JL-2)*PH(JV,JL-2)-PB(JV,JL-1)*PG(JV,JL-1)
        PG(JV,JL)=(PES(JV,JL)-PB(JV,JL-1)*PH(JV,JL-1))/PA(JV,JL)
        PH(JV,JL)=PFS(JV,JL)/PA(JV,JL)
      ENDDO
    ENDDO
  ENDIF

  IL=KLX-1
  DO JV=1,KVX
    PB(JV,IL)=PEI(JV,IL)-PC(JV,IL-1)*PG(JV,IL-1)
    PA(JV,IL)= PD(JV,IL)-PC(JV,IL-2)*PH(JV,IL-2)-PB(JV,IL-1)*PG(JV,IL-1)
    PG(JV,IL)=(PES(JV,IL)-PB(JV,IL-1)*PH(JV,IL-1))/PA(JV,IL)
  ENDDO

  IL=KLX
  DO JV=1,KVX
    PA(JV,IL)=PD(JV,IL)-PC(JV,IL-2)*PH(JV,IL-2)-PB(JV,IL-1)*PG(JV,IL-1)
  ENDDO

ELSEIF (KLX == 3) THEN

  DO JV=1,KVX
    PC(JV,1)=PFI(JV,1)
    PB(JV,1)=PEI(JV,1)
    PA(JV,1)= PD(JV,1)
    PG(JV,1)=PES(JV,1)/PA(JV,1)
    PH(JV,1)=PFS(JV,1)/PA(JV,1)
    PB(JV,2)=PEI(JV,2)-PC(JV,1)*PG(JV,1)
    PA(JV,2)= PD(JV,2)-PB(JV,1)*PG(JV,1)
    PG(JV,2)=(PES(JV,2)-PB(JV,1)*PH(JV,1))/PA(JV,2)
    PA(JV,3)=PD(JV,3)-PC(JV,1)*PH(JV,1)-PB(JV,2)*PG(JV,2)
  ENDDO

ELSEIF (KLX == 2) THEN

  DO JV=1,KVX
    PB(JV,1)=PEI(JV,1)
    PA(JV,1)= PD(JV,1)
    PG(JV,1)=PES(JV,1)/PA(JV,1)
    PA(JV,2)=PD(JV,2)-PB(JV,1)*PG(JV,1)
  ENDDO

ELSEIF (KLX == 1) THEN

  DO JV=1,KVX
    PA(JV,1)=PD(JV,1)
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUHER',1,ZHOOK_HANDLE)
END SUBROUTINE SUHER


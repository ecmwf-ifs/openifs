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

SUBROUTINE SUHES(KLX,KVX,KVXS,PD,PE,PF,PA,PB,PC)

!**** *SUHES*   - Decomposition of a set of symmetric pentadiagonal matrixes.

!     Purpose.    Decomposition of a set of symmetric pentadiagonal matrixes
!     --------    into products of two tridiagonal triangular matrixes.

!    Example, for KLX=5, matrix number l:

!    I PD(l,1) PE(l,1) PF(l,1)    0       0    I
!    I PE(l,1) PD(l,2) PE(l,2)    0       0    I
!    I PF(l,1) PE(l,2) PD(l,3) PE(l,3) PF(l,3) I
!    I    0    PF(l,2) PE(l,3) PD(l,4) PE(l,4) I
!    I    0       0    PF(l,3) PE(l,4) PD(l,5) I

!                      EQUAL

!    I PA(l,1)    0       0       0       0    I
!    I PB(l,1) PA(l,2)    0       0       0    I
!    I PC(l,1) PB(l,2) PA(l,3)    0       0    I
!    I    0    PC(l,2) PB(l,3) PA(l,4)    0    I
!    I    0       0    PC(l,3) PB(l,4) PA(l,5) I

!                      TIMES

!    I   1   PB(l,1)/PA(l,1) PC(l,1)/PA(l,1)        0               0        I
!    I   0          1        PB(l,2)/PA(l,2) PC(l,2)/PA(l,2)        0        I
!    I   0          0               1        PB(l,3)/PA(l,3) PC(l,3)/PA(l,3) I
!    I   0          0               0               1        PB(l,4)/PA(l,4) I
!    I   0          0               0               0               1        I

!**   Interface.
!     ----------
!        *CALL* *SUHES(KLX,KVX,KVXS,PD,PE,PF,PA,PB,PC)

!        Explicit arguments :
!        --------------------
!         KLX:       - Dimension of the system.                     (input)
!         KVX        - Number of matrixes to be factorised.         (input)
!         KVXS       - Surdimension associated to KVX.              (input)
!         PD,PE,PF:  - non-zero diagonals of the initial matrixes   (input)
!                      (see figure above). All minor determinants
!                      of these matrixes must be non zero.
!         PA,PB,PC:  - non-zero diagonals of the lower triangular   (output)
!                      matrixes (see figure above).

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PE(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PF(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PA(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PC(KVXS,KLX) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IL, JL, JV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUHES',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF PA,PB,PC.
!              ------------------------

IF (KLX >= 4) THEN

  DO JV=1,KVX
    PC(JV,1)=PF(JV,1)
    PB(JV,1)=PE(JV,1)
    PA(JV,1)=PD(JV,1)
  ENDDO

  DO JV=1,KVX
    PC(JV,2)=PF(JV,2)
    PB(JV,2)=PE(JV,2)-PC(JV,1)*PB(JV,1)/PA(JV,1)
    PA(JV,2)=PD(JV,2)-PB(JV,1)*PB(JV,1)/PA(JV,1)
  ENDDO

  IF (KLX >= 5) THEN
    DO JL=3,KLX-2
      DO JV=1,KVX
        PC(JV,JL)=PF(JV,JL)
        PB(JV,JL)=PE(JV,JL)-PC(JV,JL-1)*PB(JV,JL-1)/PA(JV,JL-1)
        PA(JV,JL)=PD(JV,JL)-PC(JV,JL-2)*PC(JV,JL-2)/PA(JV,JL-2)&
         & -PB(JV,JL-1)*PB(JV,JL-1)/PA(JV,JL-1)  
      ENDDO
    ENDDO
  ENDIF

  IL=KLX-1
  DO JV=1,KVX
    PB(JV,IL)=PE(JV,IL)-PC(JV,IL-1)*PB(JV,IL-1)/PA(JV,IL-1)
    PA(JV,IL)=PD(JV,IL)-PC(JV,IL-2)*PC(JV,IL-2)/PA(JV,IL-2)&
     & -PB(JV,IL-1)*PB(JV,IL-1)/PA(JV,IL-1)  
  ENDDO

  IL=KLX
  DO JV=1,KVX
    PA(JV,IL)=PD(JV,IL)-PC(JV,IL-2)*PC(JV,IL-2)/PA(JV,IL-2)&
     & -PB(JV,IL-1)*PB(JV,IL-1)/PA(JV,IL-1)  
  ENDDO

ELSEIF (KLX == 3) THEN

  DO JV=1,KVX
    PC(JV,1)=PF(JV,1)
    PB(JV,1)=PE(JV,1)
    PA(JV,1)=PD(JV,1)
    PB(JV,2)=PE(JV,2)-PC(JV,1)*PB(JV,1)/PA(JV,1)
    PA(JV,2)=PD(JV,2)-PB(JV,1)*PB(JV,1)/PA(JV,1)
    PA(JV,3)=PD(JV,3)-PC(JV,1)*PC(JV,1)/PA(JV,1)-PB(JV,2)*PB(JV,2)/PA(JV,2)
  ENDDO

ELSEIF (KLX == 2) THEN

  DO JV=1,KVX
    PB(JV,1)=PE(JV,1)
    PA(JV,1)=PD(JV,1)
    PA(JV,2)=PD(JV,2)-PB(JV,1)*PB(JV,1)/PA(JV,1)
  ENDDO

ELSEIF (KLX == 1) THEN

  DO JV=1,KVX
    PA(JV,1)=PD(JV,1)
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUHES',1,ZHOOK_HANDLE)
END SUBROUTINE SUHES


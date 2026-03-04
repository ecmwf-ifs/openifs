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

SUBROUTINE MXPTMA(KLX,KVX,KVXS,KIX,PA,PBI,PCI,PBS,PCS,PX,PY)

!**** *MXPTMA*   - Multiplication of a pentadiagonal matrix by a matrix.

!     Purpose.    Multiplication of a pentadiagonal matrix by a matrix.
!     --------

!    Example, for KLX=5,KVX=2,KIX=1.

!    I PA (1) PBS(1) PCS(1)   0      0    I   I PX(1,1) PX(2,1) I
!    I PBI(1) PA (2) PBS(2) PCS(2)   0    I   I PX(1,2) PX(2,2) I
!    I PCI(1) PBI(2) PA (3) PBS(3) PCS(3) I * I PX(1,3) PX(2,3) I
!    I   0    PCI(2) PBI(3) PA (4) PBS(4) I   I PX(1,4) PX(2,4) I
!    I   0      0    PCI(3) PBI(4) PA (5) I   I PX(1,5) PX(2,5) I

!      I PY(1,1) PY(2,1) I
!      I PY(1,2) PY(2,2) I
!    = I PY(1,3) PY(2,3) I
!      I PY(1,4) PY(2,4) I
!      I PY(1,5) PY(2,5) I

!**   Interface.
!     ----------
!        *CALL* *MXPTMA(KLX,KVX,KVXS,KIX,PA,PBI,PCI,PBS,PCS,PX,PY)

!        Explicit arguments :
!        --------------------
!         KLX:         - Dimension of the matrix.                   (input)
!         KVX,KIX:     - Number of vectors to be multiplied is      (input)
!                        KVX*KIX.
!         KVXS:        - Surdimension corresponding to KVX.         (input)
!         PA:          - Diagonal of the matrix.                    (input)
!         PBI,PCI:     - Lower diagonals of the matrix.             (input)
!         PBS,PCS:     - Upper diagonals of the matrix.             (input)
!         PX:          - Initial vector:                            (input)
!         PY:          - Final vector:                              (output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None.

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
INTEGER(KIND=JPIM),INTENT(IN)    :: KIX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBI(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBS(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCS(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX(KVXS,KLX,KIX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PY(KVXS,KLX,KIX) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JI, JL, JV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MXPTMA',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF PY.
!              ------------------

IF (KLX >= 4) THEN

  DO JI=1,KIX
    DO JV=1,KVX
      PY(JV,1,JI) = PA (1)*PX(JV,1,JI)+PBS(1)*PX(JV,2,JI)+PCS(1)*PX(JV,3,JI)
      PY(JV,2,JI) = PBI(1)*PX(JV,1,JI)&
       & +PA (2)*PX(JV,2,JI)&
       & +PBS(2)*PX(JV,3,JI)&
       & +PCS(2)*PX(JV,4,JI)  
    ENDDO
  ENDDO

  DO JI=1,KIX
    DO JV=1,KVX
      DO JL=3,KLX-2
        PY(JV,JL,JI) = PCI(JL-2)*PX(JV,JL-2,JI)&
         & +PBI(JL-1)*PX(JV,JL-1,JI)&
         & +PA (JL  )*PX(JV,JL  ,JI)&
         & +PBS(JL  )*PX(JV,JL+1,JI)&
         & +PCS(JL  )*PX(JV,JL+2,JI)  
      ENDDO
    ENDDO
  ENDDO

  DO JI=1,KIX
    DO JV=1,KVX
      PY(JV,KLX-1,JI) = PCI(KLX-3)*PX(JV,KLX-3,JI)&
       & +PBI(KLX-2)*PX(JV,KLX-2,JI)&
       & +PA (KLX-1)*PX(JV,KLX-1,JI)&
       & +PBS(KLX-1)*PX(JV,KLX  ,JI)  
      PY(JV,KLX,JI) = PCI(KLX-2)*PX(JV,KLX-2,JI)&
       & +PBI(KLX-1)*PX(JV,KLX-1,JI)&
       & +PA (KLX  )*PX(JV,KLX  ,JI)  
    ENDDO
  ENDDO

ELSEIF (KLX == 3) THEN

  DO JI=1,KIX
    DO JV=1,KVX
      PY(JV,1,JI) = PA (1)*PX(JV,1,JI)+PBS(1)*PX(JV,2,JI)+PCS(1)*PX(JV,3,JI)
      PY(JV,2,JI) = PBI(1)*PX(JV,1,JI)+PA (2)*PX(JV,2,JI)+PBS(2)*PX(JV,3,JI)
      PY(JV,3,JI) = PCI(1)*PX(JV,1,JI)+PBI(2)*PX(JV,2,JI)+PA (3)*PX(JV,3,JI)
    ENDDO
  ENDDO

ELSEIF (KLX == 2) THEN

  DO JI=1,KIX
    DO JV=1,KVX
      PY(JV,1,JI) = PA (1)*PX(JV,1,JI)+PBS(1)*PX(JV,2,JI)
      PY(JV,2,JI) = PBI(1)*PX(JV,1,JI)+PA (2)*PX(JV,2,JI)
    ENDDO
  ENDDO

ELSEIF (KLX == 1) THEN

  DO JI=1,KIX
    DO JV=1,KVX
      PY(JV,1,JI) = PA (1)*PX(JV,1,JI)
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MXPTMA',1,ZHOOK_HANDLE)
END SUBROUTINE MXPTMA


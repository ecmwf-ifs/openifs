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

SUBROUTINE MXTURE(KLX,KVX,KVXS,KIX,KT,LDMT,PA,PB,PC,PY,PX)

!**** *MXTURE*   - Resolution of a set of triangular tridiagonal systems.

!     Purpose.    Resolution of a set of triangular tridiagonal systems.
!     --------

!    Example, for KLX=4:

!    PY is known, PX is unknown. We want to find PX, solution of the
!     following system, for l=1 to KVX, i=1 to KIX.

!    * KT = -2:

!    I PA(l,1)    0       0       0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I PB(l,1) PA(l,2)    0       0    I   I PX(l,2,i) I   I PY(l,2,i) I
!    I PC(l,1) PB(l,2) PA(l,3)    0    I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0    PC(l,2) PB(l,3) PA(l,4) I   I PX(l,4,i) I   I PY(l,4,i) I

!    * KT = -3:

!    I    1       0       0       0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I PB(l,1)    1       0       0    I   I PX(l,2,i) I   I PY(l,2,i) I
!    I PC(l,1) PB(l,2)    1       0    I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0    PC(l,2) PB(l,3)    1    I   I PX(l,4,i) I   I PY(l,4,i) I

!    Dummy array PA is not used in this case.

!    * KT = 1:

!    I    1    ZB(l,1) ZC(l,1)    0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I    0       1    ZB(l,2) ZC(l,2) I   I PX(l,2,i) I   I PY(l,2,i) I
!    I    0       0       1    ZB(l,3) I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0       0       0       1    I   I PX(l,4,i) I   I PY(l,4,i) I

!    where ZB(l,j)=PB(l,j)/PA(l,j); ZC(l,j)=PC(l,j)/PA(l,j) .

!    * KT = 2:

!    I PA(l,1) PB(l,1) PC(l,1)    0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I    0    PA(l,2) PB(l,2) PC(l,2) I   I PX(l,2,i) I   I PY(l,2,i) I
!    I    0       0    PA(l,3) PB(l,3) I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0       0       0    PA(l,4) I   I PX(l,4,i) I   I PY(l,4,i) I

!    * KT = 3:

!    I    1    PB(l,1) PC(l,1)    0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I    0       1    PB(l,2) PC(l,2) I   I PX(l,2,i) I   I PY(l,2,i) I
!    I    0       0       1    PB(l,3) I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0       0       0       1    I   I PX(l,4,i) I   I PY(l,4,i) I

!    Dummy array PA is not used in this case.

!**   Interface.
!     ----------
!        *CALL* *MXTURE(KLX,KVX,KVXS,KIX,KT,LDMT,PA,PB,PC,PY,PX)

!        Explicit arguments :
!        --------------------
!         KLX:        - Dimension of the system.                    (input)
!         KVX:        - Number of variables (second                 (input)
!                       dimension of PX and PY).
!         KVXS:       - Surdimension corresponding to KVX.          (input)
!         KIX:        - Number of variables multiplied by the same
!                       matrix.                                     (input)
!         KT:         - Type of matrix (see figures above).         (input)
!         LDMT:       - .T.: copy of PX on PY at the end of subr.   (input)
!                       .F.: no copy.
!         PA,PB,PC:   - non-zero diagonals of the triangular        (input)
!                       matrixes (see figures above).
!                       Caution! All PA coefficients must be non 0.
!         PY:         - known vector.                               (input)
!         PX:         - unknown vector.                             (output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None.
!        Called by SPC,SPCAD,SPC2,SPC2AD,MXTURS.

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
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVXS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KT 
LOGICAL           ,INTENT(IN)    :: LDMT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KVX,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(KVX,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC(KVX,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PY(KVXS,KLX,KIX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PX(KVXS,KLX,KIX) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IIX, IT, JI, JL, JV

REAL(KIND=JPRB) :: ZBB, ZCC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MXTURE',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

IIX=KIX-MOD(KIX,2)
IT=MIN(3,MAX(-3,KT))
IF (KT == 0.OR.KT == -1) IT=-2

!      ----------------------------------------------------------------

!*       2.    COMPUTATION OF PX.
!              ------------------

IF (IT == -3) THEN

  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,1,JI)=PY(JV,1,JI)
      PX(JV,1,JI+1)=PY(JV,1,JI+1)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,1,JI)=PY(JV,1,JI)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        ZBB=PB(JV,1)
        PX(JV,2,JI)=PY(JV,2,JI)-ZBB*PX(JV,1,JI)
        PX(JV,2,JI+1)=PY(JV,2,JI+1)-ZBB*PX(JV,1,JI+1)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        ZBB=PB(JV,1)
        PX(JV,2,JI)=PY(JV,2,JI)-ZBB*PX(JV,1,JI)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=3,KLX
        DO JV=1,KVX
          ZBB=PB(JV,JL-1)
          ZCC=PC(JV,JL-2)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL-1,JI)-ZCC*PX(JV,JL-2,JI)
          PX(JV,JL,JI+1)=PY(JV,JL,JI+1)&
           & -ZBB*PX(JV,JL-1,JI+1)-ZCC*PX(JV,JL-2,JI+1)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=3,KLX
        DO JV=1,KVX
          ZBB=PB(JV,JL-1)
          ZCC=PC(JV,JL-2)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL-1,JI)-ZCC*PX(JV,JL-2,JI)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT == -2) THEN

  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,1,JI)=PY(JV,1,JI)/PA(JV,1)
      PX(JV,1,JI+1)=PY(JV,1,JI+1)/PA(JV,1)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,1,JI)=PY(JV,1,JI)/PA(JV,1)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        PX(JV,2,JI)=(PY(JV,2,JI)-PB(JV,1)*PX(JV,1,JI))/PA(JV,2)
        PX(JV,2,JI+1)=(PY(JV,2,JI+1)-PB(JV,1)*PX(JV,1,JI+1))/PA(JV,2)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        PX(JV,2,JI)=(PY(JV,2,JI)-PB(JV,1)*PX(JV,1,JI))/PA(JV,2)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=3,KLX
        DO JV=1,KVX
          PX(JV,JL,JI)=(PY(JV,JL,JI)-PC(JV,JL-2)*PX(JV,JL-2,JI)&
           & -PB(JV,JL-1)*PX(JV,JL-1,JI))/PA(JV,JL)  
          PX(JV,JL,JI+1)=(PY(JV,JL,JI+1)-PC(JV,JL-2)*PX(JV,JL-2,&
           & JI+1)&
           & -PB(JV,JL-1)*PX(JV,JL-1,JI+1))/PA(JV,JL)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=3,KLX
        DO JV=1,KVX
          PX(JV,JL,JI)=(PY(JV,JL,JI)-PC(JV,JL-2)*PX(JV,JL-2,JI)&
           & -PB(JV,JL-1)*PX(JV,JL-1,JI))/PA(JV,JL)  
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT == 1) THEN

  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)
      PX(JV,KLX,JI+1)=PY(JV,KLX,JI+1)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        ZBB=PB(JV,KLX-1)/PA(JV,KLX-1)
        PX(JV,KLX-1,JI)=PY(JV,KLX-1,JI)-ZBB*PX(JV,KLX,JI)
        PX(JV,KLX-1,JI+1)=PY(JV,KLX-1,JI+1)-ZBB*PX(JV,KLX,JI+1)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        ZBB=PB(JV,KLX-1)/PA(JV,KLX-1)
        PX(JV,KLX-1,JI)=PY(JV,KLX-1,JI)-ZBB*PX(JV,KLX,JI)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          ZBB=PB(JV,JL)/PA(JV,JL)
          ZCC=PC(JV,JL)/PA(JV,JL)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL+1,JI)-ZCC*PX(JV,JL+2,JI)
          PX(JV,JL,JI+1)=PY(JV,JL,JI+1)&
           & -ZBB*PX(JV,JL+1,JI+1)-ZCC*PX(JV,JL+2,JI+1)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          ZBB=PB(JV,JL)/PA(JV,JL)
          ZCC=PC(JV,JL)/PA(JV,JL)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL+1,JI)-ZCC*PX(JV,JL+2,JI)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT == 2) THEN

  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)/PA(JV,KLX)
      PX(JV,KLX,JI+1)=PY(JV,KLX,JI+1)/PA(JV,KLX)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)/PA(JV,KLX)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        PX(JV,KLX-1,JI)=&
         & (PY(JV,KLX-1,JI)-PB(JV,KLX-1)*PX(JV,KLX,JI))/PA(JV,KLX-1)  
        PX(JV,KLX-1,JI+1)=&
         & (PY(JV,KLX-1,JI+1)-PB(JV,KLX-1)*PX(JV,KLX,JI+1))/PA(JV,&
         & KLX-1)  
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        PX(JV,KLX-1,JI)=&
         & (PY(JV,KLX-1,JI)-PB(JV,KLX-1)*PX(JV,KLX,JI))/PA(JV,KLX-1)  
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          PX(JV,JL,JI)=(PY(JV,JL,JI)-PB(JV,JL)*PX(JV,JL+1,JI)&
           & -PC(JV,JL)*PX(JV,JL+2,JI))/PA(JV,JL)  
          PX(JV,JL,JI+1)=(PY(JV,JL,JI+1)-PB(JV,JL)*PX(JV,JL+1,JI+&
           & 1)&
           & -PC(JV,JL)*PX(JV,JL+2,JI+1))/PA(JV,JL)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          PX(JV,JL,JI)=(PY(JV,JL,JI)-PB(JV,JL)*PX(JV,JL+1,JI)&
           & -PC(JV,JL)*PX(JV,JL+2,JI))/PA(JV,JL)  
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT == 3) THEN

  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)
      PX(JV,KLX,JI+1)=PY(JV,KLX,JI+1)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        ZBB=PB(JV,KLX-1)
        PX(JV,KLX-1,JI)=PY(JV,KLX-1,JI)-ZBB*PX(JV,KLX,JI)
        PX(JV,KLX-1,JI+1)=PY(JV,KLX-1,JI+1)-ZBB*PX(JV,KLX,JI+1)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        ZBB=PB(JV,KLX-1)
        PX(JV,KLX-1,JI)=PY(JV,KLX-1,JI)-ZBB*PX(JV,KLX,JI)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          ZBB=PB(JV,JL)
          ZCC=PC(JV,JL)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL+1,JI)-ZCC*PX(JV,JL+2,JI)
          PX(JV,JL,JI+1)=PY(JV,JL,JI+1)&
           & -ZBB*PX(JV,JL+1,JI+1)-ZCC*PX(JV,JL+2,JI+1)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          ZBB=PB(JV,JL)
          ZCC=PC(JV,JL)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL+1,JI)-ZCC*PX(JV,JL+2,JI)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ENDIF
!      ----------------------------------------------------------------

!*       3.    FINAL MEMORY TRANSFER.
!              ----------------------

IF (LDMT) THEN
  DO JI=1,KIX
    DO JL=1,KLX
      DO JV=1,KVX
        PY(JV,JL,JI)=PX(JV,JL,JI)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MXTURE',1,ZHOOK_HANDLE)
END SUBROUTINE MXTURE


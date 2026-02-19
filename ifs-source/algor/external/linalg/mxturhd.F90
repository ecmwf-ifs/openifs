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

SUBROUTINE MXTURHD(KLX,KVX,KVXS,KT,LDMT,PA,PB,PY,PX)

!**** *MXTURHD*   - Resolution of a set of triangular bidiagonal systems.

!     Purpose.    Resolution of a set of triangular bidiagonal systems.
!     --------

!    Example, for KLX=4:

!    PY is known, PX is unknown. We want to find PX, solution of the
!     following system, for l=1 to KVX.

!    * KT = -2:

!    I PA(l,1)    0       0       0    I   I PX(l,1) I   I PY(l,1) I
!    I PB(l,1) PA(l,2)    0       0    I   I PX(l,2) I   I PY(l,2) I
!    I    0    PB(l,2) PA(l,3)    0    I * I PX(l,3) I = I PY(l,3) I
!    I    0       0    PB(l,3) PA(l,4) I   I PX(l,4) I   I PY(l,4) I

!    * KT = -3:

!    I    1       0       0       0    I   I PX(l,1) I   I PY(l,1) I
!    I PB(l,1)    1       0       0    I   I PX(l,2) I   I PY(l,2) I
!    I    0    PB(l,2)    1       0    I * I PX(l,3) I = I PY(l,3) I
!    I    0       0    PB(l,3)    1    I   I PX(l,4) I   I PY(l,4) I

!    Dummy array PA is not used in this case.

!    * KT = 1:

!    I    1    ZB(l,1)    0       0    I   I PX(l,1) I   I PY(l,1) I
!    I    0       1    ZB(l,2)    0    I   I PX(l,2) I   I PY(l,2) I
!    I    0       0       1    ZB(l,3) I * I PX(l,3) I = I PY(l,3) I
!    I    0       0       0       1    I   I PX(l,4) I   I PY(l,4) I

!    where ZB(l,j)=PB(l,j)/PA(l,j)

!    * KT = 2:

!    I PA(l,1) PB(l,1)    0       0    I   I PX(l,1) I   I PY(l,1) I
!    I    0    PA(l,2) PB(l,2)    0    I   I PX(l,2) I   I PY(l,2) I
!    I    0       0    PA(l,3) PB(l,3) I * I PX(l,3) I = I PY(l,3) I
!    I    0       0       0    PA(l,4) I   I PX(l,4) I   I PY(l,4) I

!    * KT = 3:

!    I    1    PB(l,1)    0       0    I   I PX(l,1) I   I PY(l,1) I
!    I    0       1    PB(l,2)    0    I   I PX(l,2) I   I PY(l,2) I
!    I    0       0       1    PB(l,3) I * I PX(l,3) I = I PY(l,3) I
!    I    0       0       0       1    I   I PX(l,4) I   I PY(l,4) I

!    Dummy array PA is not used in this case.

!**   Interface.
!     ----------
!        *CALL* *MXTURHD(KLX,KVX,KVXS,KT,LDMT,PA,PB,PY,PX)

!        Explicit arguments :
!        --------------------
!         KLX:        - Dimension of the system.                    (input)
!         KVX:        - Number of variables (based on  NFLEVL )     (input)
!         KVXS:       - Surdimension based on NFLEVG for PA,PB      (input)
!         KT:         - Type of matrix (see figures above).         (input)
!                       Only values -2 , 1 and 3 are considered
!         LDMT:       - .T.: copy of PX on PY at the end of subr.   (input)
!                       .F.: no copy.
!         PA,PB:      - non-zero diagonals of the triangular        (input)
!                       matrixes (see figures above).
!                       Caution! All PA coefficients must be non 0.
!                       Initialisation in SUHDU (based on NFLEVG)
!         PY:         - known vector.                               (input/output)
!                       If LDMT=.F. PY is an input array only,
!                       otherwise it is updated at the end.
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

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        (MXTURE) K. YESSAD: OCTOBER 1993.
!        MXTURHD : Jean Latour (CRI) September 1996

!     Modifications.
!     --------------
!     L.Isaksen 96-12-12 : More MXTURE options added, cleaned
!     M.Hamrud      01-Oct-2003 CY28 Cleaning
!     K.Yessad  03-12-02 : Bidiagonal matrices instead of tridiagonal ones.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!      ----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVXS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KT 
LOGICAL           ,INTENT(IN)    :: LDMT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KVX,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(KVX,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PY(KVXS,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PX(KVXS,KLX) 

!      ----------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL, JV

REAL(KIND=JPRB) :: ZBB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

#include "abor1.intfb.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MXTURHD',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

!*       1.    COMPUTATION OF PX.
!              ------------------

IF (KT == -3) THEN

  PX(:,1)=PY(:,1)

  IF (KLX >= 2) THEN
    DO JL=2,KLX
      DO JV=1,KVX
        ZBB=PB(JV,JL-1)
        PX(JV,JL)=PY(JV,JL)-ZBB*PX(JV,JL-1)
      ENDDO
    ENDDO
  ENDIF

ELSEIF (KT == -2) THEN

  DO JV=1,KVX
    PX(JV,1)=PY(JV,1)/PA(JV,1)
  ENDDO

  IF (KLX >= 2) THEN
    DO JL=2,KLX
      DO JV=1,KVX
        PX(JV,JL)=(PY(JV,JL)-PB(JV,JL-1)*PX(JV,JL-1))/PA(JV,JL)  
      ENDDO
    ENDDO
  ENDIF

ELSEIF (KT == 1) THEN

  PX(:,KLX)=PY(:,KLX)

  IF (KLX >= 2) THEN
    DO JL=KLX-1,1,-1
      DO JV=1,KVX
        ZBB=PB(JV,JL)/PA(JV,JL)
        PX(JV,JL)=PY(JV,JL)-ZBB*PX(JV,JL+1)
      ENDDO
    ENDDO
  ENDIF

ELSEIF (KT == 2) THEN

  DO JV=1,KVX
    PX(JV,KLX)=PY(JV,KLX)/PA(JV,KLX)
  ENDDO

  IF (KLX >= 2) THEN
    DO JL=KLX-1,1,-1
      DO JV=1,KVX
        PX(JV,JL)=(PY(JV,JL)-PB(JV,JL)*PX(JV,JL+1))/PA(JV,JL)  
      ENDDO
    ENDDO
  ENDIF

ELSEIF (KT == 3) THEN

  PX(:,KLX)=PY(:,KLX)

  IF (KLX >= 2) THEN
    DO JL=KLX-1,1,-1
      DO JV=1,KVX
        ZBB=PB(JV,JL)
        PX(JV,JL)=PY(JV,JL)-ZBB*PX(JV,JL+1)
      ENDDO
    ENDDO
  ENDIF

ELSE
  CALL ABOR1('KT value in MXTURHD NOT IMPLEMENTED')
ENDIF
!      ----------------------------------------------------------------

!*       2.    FINAL MEMORY TRANSFER.
!              ----------------------

IF (LDMT) THEN
  PY(:,:)=PX(:,:)
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MXTURHD',1,ZHOOK_HANDLE)
END SUBROUTINE MXTURHD


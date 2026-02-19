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

SUBROUTINE MXTURS(KLX,KVX,KVXS,KIX,PA,PB,PC,PY,PX)

!**** *MXTURS*   - Resolution of a set of pentadiagonal symmetric systems.

!     Purpose.    Resolution of a set of pentadiagonal symmetric systems,
!     --------    once known the LU decomposition of these symmetric matrixes.

!    Example, for KLX=4:

!    PY is known, PX is unknown. We want to find PX, solution of the
!    following system, for l=1 to KVX, i=1 to KIX.
!                            S(l)*PX(l,i)=PY(l,i)
!    where S(l) are symmetric pentadiagonal matrixes, the LU decomposition
!    of which yields:

!           I PA(l,1)    0       0       0    I
!    L(l) = I PB(l,1) PA(l,2)    0       0    I
!           I PC(l,1) PB(l,2) PA(l,3)    0    I
!           I    0    PC(l,2) PB(l,3) PA(l,4) I

!           I    1    ZB(l,1) ZC(l,1)    0    I
!    U(l) = I    0       1    ZB(l,2) ZC(l,2) I
!           I    0       0       1    ZB(l,3) I
!           I    0       0       0       1    I

!    where ZB(l,j)=PB(l,j)/PA(l,j); ZC(l,j)=PC(l,j)/PA(l,j) .

!    S(l)=L(l)*U(l)

!    We call routine MXTURE to inverse L, then U.

!**   Interface.
!     ----------
!        *CALL* *MXTURS(KLX,KVX,KVXS,KIX,PA,PB,PC,PY,PX)

!        Explicit arguments :
!        --------------------
!         KLX:        - Dimension of the system.                    (input)
!         KVX:        - Number of variables (second                 (input)
!                       dimension of PX and PY).
!         KVXS:       - Surdimension corresponding to KVX.          (input)
!         KIX:        - Number of variables multiplied by the same
!                       matrix.                                     (input)
!         PA,PB,PC:   - non-zero diagonals of the triangular        (input)
!                       matrixes L (see figures above).
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
!        Calls MXTURE.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD: MARCH 1994.

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KVX,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(KVX,KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC(KVX,KLX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PY(KVXS,KLX,KIX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PX(KVXS,KLX,KIX) 

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxture.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MXTURS',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    INVERSION OF THE TWO TRIANGULAR TRIDIAGONAL MATRIXES.
!              -----------------------------------------------------

CALL MXTURE(KLX,KVX,KVXS,KIX,-2,.TRUE. ,PA,PB,PC,PY,PX)
CALL MXTURE(KLX,KVX,KVXS,KIX, 1,.FALSE.,PA,PB,PC,PY,PX)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MXTURS',1,ZHOOK_HANDLE)
END SUBROUTINE MXTURS


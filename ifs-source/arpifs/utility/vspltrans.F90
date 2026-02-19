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

SUBROUTINE VSPLTRANS(YDVSPLIP,KPROMA,KST,KPROF,KFLEV,KFLSA,KFLEN,PX)

!**** *VSPLTRANS* Transforms a field to B-spline space in the vertical

!     Purpose.
!     --------
!        * Transforms a field to B-spline space in the vertical for
!          easy vertical interpolation with cubic splines

!**   Interface.
!     ----------
!        *CALL* *VSPLTRANS(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA  - horizontal dimension.
!          KST     - first element of work.
!          KPROF   - depth of work.
!          KFLEV   - number of layers.
!          KFLSA   - lower number of layer for arrays with extra-layers.
!          KFLEN   - upper number of layer for arrays with extra-layers.

!        INPUT/OUTPUT:

!          PX      - on input: field to be transformed
!                    on output: result of transform

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           TRIDIA

!           Called by LATTEX.

!     Author.
!     -------
!        A. UNTCH 

!     Modifications.
!     --------------
!        Original : AUGUST 2001.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K.Yessad (Nov 2008): replace TRIDIAVSPL by TRIDIA + cleanings
!     ------------------------------------------------------------------

USE YOMVSPLIP , ONLY : TVSPLIP
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVSPLIP)     ,INTENT(IN)    :: YDVSPLIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLSA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PX(KPROMA,KFLSA:KFLEN) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSOL(KST:KPROF,KFLEV)
REAL(KIND=JPRB) :: ZRHS(KST:KPROF,KFLEV)
REAL(KIND=JPRB) :: ZDERI_T(KST:KPROF)
REAL(KIND=JPRB) :: ZDERI_B(KST:KPROF)

INTEGER(KIND=JPIM) :: JLEV, IDIM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "tridia.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('VSPLTRANS',0,ZHOOK_HANDLE)
ASSOCIATE(RVSPTRI=>YDVSPLIP%RVSPTRI, RVSPC=>YDVSPLIP%RVSPC)
!     ------------------------------------------------------------------

IDIM=KPROF-KST+1

! estimate first derivatives of PX at top and bottom
!!ZDERI_T(KST:KPROF)=(PX(KST:KPROF,2)-PX(KST:KPROF,1))*RVSPC(9)
ZDERI_T(KST:KPROF)=0.0_JPRB
!!ZDERI_B(KST:KPROF)=(PX(KST:KPROF,KFLEV)-PX(KST:KPROF,KFLEV-1))*RVSPC(10)
ZDERI_B(KST:KPROF)=0.0_JPRB

ZRHS(KST:KPROF,1)=PX(KST:KPROF,1)+RVSPC(4)*ZDERI_T(KST:KPROF)
DO JLEV=2,KFLEV-1
  ZRHS(KST:KPROF,JLEV)=PX(KST:KPROF,JLEV)
ENDDO
ZRHS(KST:KPROF,KFLEV)=PX(KST:KPROF,KFLEV)+RVSPC(8)*ZDERI_B(KST:KPROF)

CALL TRIDIA(KFLEV,IDIM,KST,KPROF,1,RVSPTRI,ZRHS,ZSOL)

DO JLEV=1,KFLEV
  PX(KST:KPROF,JLEV)=ZSOL(KST:KPROF,JLEV)
ENDDO
PX(KST:KPROF,      0)=RVSPC(1)*ZDERI_T(KST:KPROF)&
 & +RVSPC(2)*PX(KST:KPROF,1)+RVSPC(3)*PX(KST:KPROF,2)  
PX(KST:KPROF,KFLEV+1)=RVSPC(5)*ZDERI_B(KST:KPROF)&
 & +RVSPC(6)*PX(KST:KPROF,KFLEV-1)+RVSPC(7)*PX(KST:KPROF,KFLEV)  

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VSPLTRANS',1,ZHOOK_HANDLE)
END SUBROUTINE VSPLTRANS

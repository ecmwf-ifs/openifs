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

MODULE LARCHE_HLP
IMPLICIT NONE
!**** *LARCHE_HLP -  semi-LAgrangian scheme:
!                    interface with routines of library MASS;
!                    these routines being called for example in routine LARCHE.

!     Purpose.
!     --------
!      Routines JFH_VMOD,JFH_VIF2,VEXP_,VASIN_

!**   Interface.
!     ----------

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        ???

!     Modifications.
!     --------------
!        Original : SEPTEMBER 2003.
!        Modified : 22 Sep 2003 K. YESSAD: rewrite in DOCTOR norm.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

CONTAINS 

! ======================================================================

#ifdef RS6K
@PROCESS OPT(3),STRICT
#endif

SUBROUTINE JFH_VMOD(PF,PC,KN)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!     USE XLF_FP_UTIL
!     ------------------------------------------------------------------
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KN 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PF(KN) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC 
INTEGER(KIND=JPIM) :: JI
REAL(KIND=JPRB) :: ZT,ZX

REAL(KIND=JPRB) :: PPRND
PARAMETER(PPRND=2.0_JPRB**52+2.0_JPRB**51)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
!     MODE=SET_ROUND_MODE(FP_RND_RN)
IF (LHOOK) CALL DR_HOOK('JFH_VMOD',0,ZHOOK_HANDLE)
DO JI=1,KN
  ZT=PF(JI)*(1.0_JPRB/PC)
  ZX=ZT-SIGN(0.5_JPRB,ZT)
  ZT=ZT-((PPRND+ZX)-PPRND)
  PF(JI)=PC*ZT
ENDDO
IF (LHOOK) CALL DR_HOOK('JFH_VMOD',1,ZHOOK_HANDLE)

!     MODE=SET_ROUND_MODE(MODE)
!     ------------------------------------------------------------------
END SUBROUTINE JFH_VMOD
! ======================================================================

#ifdef RS6K
@PROCESS OPT(3),STRICT
#endif

SUBROUTINE JFH_VMOD_JPRD(PF,PC,KN)
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!     USE XLF_FP_UTIL
!     ------------------------------------------------------------------
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KN 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PF(KN) 
REAL(KIND=JPRD)   ,INTENT(IN)    :: PC 
INTEGER(KIND=JPIM) :: JI
REAL(KIND=JPRD) :: ZT,ZX

REAL(KIND=JPRD) :: PPRND
PARAMETER(PPRND=2.0_JPRD**52+2.0_JPRD**51)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
!     MODE=SET_ROUND_MODE(FP_RND_RN)
IF (LHOOK) CALL DR_HOOK('JFH_VMOD',0,ZHOOK_HANDLE)
DO JI=1,KN
  ZT=PF(JI)*(1.0_JPRD/PC)
  ZX=ZT-SIGN(0.5_JPRD,ZT)
  ZT=ZT-((PPRND+ZX)-PPRND)
  PF(JI)=PC*ZT
ENDDO
IF (LHOOK) CALL DR_HOOK('JFH_VMOD',1,ZHOOK_HANDLE)

!     MODE=SET_ROUND_MODE(MODE)
!     ------------------------------------------------------------------
END SUBROUTINE JFH_VMOD_JPRD

! ======================================================================

#ifdef RS6K
@PROCESS OPT(3),NOSTRICT
#endif
SUBROUTINE JFH_VIF2(PA,KN)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
!     ------------------------------------------------------------------
IMPLICIT NONE
! * DUMMY:
INTEGER(KIND=JPIM),INTENT(IN) :: KN
REAL(KIND=JPRB),INTENT(INOUT) :: PA(KN)
INTEGER(KIND=JPIM) :: J
!     ------------------------------------------------------------------
DO J=1,KN
  PA(J)=MIN(1.0_JPRB,MAX(PA(J),-1.0_JPRB))
ENDDO
!     ------------------------------------------------------------------
END SUBROUTINE JFH_VIF2

! ======================================================================

SUBROUTINE JFH_VIF2_JPRD(PA,KN)
USE PARKIND1  ,ONLY : JPIM     ,JPRD
!     ------------------------------------------------------------------
IMPLICIT NONE
! * DUMMY:
INTEGER(KIND=JPIM),INTENT(IN) :: KN
REAL(KIND=JPRD),INTENT(INOUT) :: PA(KN)
INTEGER(KIND=JPIM) :: J
!     ------------------------------------------------------------------
DO J=1,KN
  PA(J)=MIN(1.0_JPRD,MAX(PA(J),-1.0_JPRD))
ENDDO
!     ------------------------------------------------------------------
END SUBROUTINE JFH_VIF2_JPRD

! ======================================================================
END MODULE LARCHE_HLP

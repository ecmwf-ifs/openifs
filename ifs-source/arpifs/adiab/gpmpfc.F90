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

SUBROUTINE GPMPFC(YDGMV,YDML_GCONF,YDDYN,YDDYNA,KPROMA,KFLEV,KST,KEN,KFLAG,PGM,PGMV,PGMVS,PGFL,&
                & P0U, P0V, P0DIV, P0TL, P0TM, P9U, P9V, P0UL, P0VL, P0VOR, P0SPDL, P0SPDM, &
                & P0SVDL, P0SVDM, P0NHXL, P0NHXM, P9DIV, P9TL, P9TM, P9SPDL, P9SPDM, P9SVDL,&
                & P9SVDM, P0SPL, P0SPM, P9SPL, P9SPM)

!**** *GPMPFC* - Apply map factor to convert 
!                reduced variables -> geographical variables if kflag=0
!                or
!                geographical variables -> reduced variables if kflag=1

!     Purpose.
!     --------
!           Multiply or divide by map factor.

!**   Interface.
!     ----------
!        *CALL* *GPMPFC(...)

!        Explicit arguments :
!        --------------------

!      INPUT:
!      ------
!       KPROMA    - horizontal dimensioning
!       KFLEV     - number of layers
!       KST       - start of work
!       KEN       - depth of work
!       KFLAG     - 0 -> multiply, 1-> divide
!       PGM       - map factor

!      INPUT/OUTPUT:
!      -------------
!       PGMV      - GMV variables
!       PGMVS     - GMVS variables
!       PGFL      - GFL variables

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*
!      Original : 1994-01-18

! Modifications
! -------------
!   Modified 2002-07-02 by C. Fischer  : rename NHS variables T0/T9
!   Modified 2002-11-13 by K. YESSAD   : some cleanings + improve vectorization
!   Modified 2003-08    by M. HAMRUD   : GFL
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   08-Jun-2004 J. Masek   NH cleaning (LVSLWBC)
!   01-Jul-2004 K. Yessad  Make clearer the tests for PC scheme.
!   09-Feb-2006 M. Deque : Dasux compilance
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Jun 2011): new dataflow, GMVS too.
!   K. Yessad (Nov 2012): simplify testings.
!   K. Yessad (July 2014): Move some variables.
!   H. Petithomme (Dec 2020): optimisation
! End Modifications
!     ------------------------------------------------------------------
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMGMV                 , ONLY : TGMV
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYNA                , ONLY : TDYNA
USE YOMDYN                 , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TGMV) , INTENT(INOUT) :: YDGMV
TYPE(TDYN)  ,INTENT(IN)    :: YDDYN
TYPE(TDYNA) ,INTENT(IN)    :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLAG 
REAL(KIND=JPRB),TARGET,INTENT(IN) :: PGM(KPROMA)
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(INOUT) :: PGMV(KPROMA,KFLEV,YDGMV%NDIMGMV)
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(INOUT) :: PGMVS(KPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(KPROMA,KFLEV,YDML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(INOUT),DIMENSION(KPROMA,KFLEV) :: P0U, P0V, P0DIV, P0TL, P0TM, P9U, P9V, P0UL, P0VL, P0VOR, P0SPDL, P0SPDM
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(INOUT),DIMENSION(KPROMA,KFLEV) :: P0SVDL, P0SVDM, P0NHXL, P0NHXM, P9DIV, P9TL, P9TM, P9SPDL, P9SPDM, P9SVDL, P9SVDM
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(INOUT),DIMENSION(KPROMA) :: P0SPL, P0SPM, P9SPL, P9SPM
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV,JGFL,JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZGM(:)
REAL(KIND=JPRB),TARGET :: ZGM0(KPROMA)
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0U, Z0V, Z0DIV, Z0TL, Z0TM, Z9U, Z9V, Z0UL, Z0VL, Z0VOR, Z0SPDL, Z0SPDM
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0SVDL, Z0SVDM, Z0NHXL, Z0NHXM, Z9DIV, Z9TL, Z9TM, Z9SPDL, Z9SPDM, Z9SVDL, Z9SVDM
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:) :: Z0SPL, Z0SPM, Z9SPL, Z9SPM

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPMPFC',0,ZHOOK_HANDLE)
ASSOCIATE(NDIM=>YDML_GCONF%YGFL%NDIM, NUMFLDS=>YDML_GCONF%YGFL%NUMFLDS, YCOMP=>YDML_GCONF%YGFL%YCOMP, &
 & LUVDER=>YDML_GCONF%YRDIMF%LUVDER, LVOR=>YDML_GCONF%YRDIMF%LVOR, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT9=>YDGMV%YT9)
!     ------------------------------------------------------------------

!*       1. APPLY MAP FACTOR.
!           -----------------

Z0U    => NULL ()
Z0V    => NULL ()
Z0DIV  => NULL ()
Z0TL   => NULL ()
Z0TM   => NULL ()
Z9U    => NULL ()
Z9V    => NULL ()
Z0UL   => NULL ()
Z0VL   => NULL ()
Z0VOR  => NULL ()
Z0SPDL => NULL ()
Z0SPDM => NULL ()
Z0SVDL => NULL ()
Z0SVDM => NULL ()
Z0NHXL => NULL ()
Z0NHXM => NULL ()
Z9DIV  => NULL ()
Z9TL   => NULL ()
Z9TM   => NULL ()
Z9SPDL => NULL ()
Z9SPDM => NULL ()
Z9SVDL => NULL ()
Z9SVDM => NULL ()
Z0SPL  => NULL ()
Z0SPM  => NULL ()
Z9SPL  => NULL ()
Z9SPM  => NULL ()

IF (PRESENT (PGMV)) THEN
  IF (YT0%MU    > 0) Z0U    => PGMV(:,:,YT0%MU)
  IF (YT0%MV    > 0) Z0V    => PGMV(:,:,YT0%MV)
  IF (YT0%MDIV  > 0) Z0DIV  => PGMV(:,:,YT0%MDIV)
  IF (YT0%MTL   > 0) Z0TL   => PGMV(:,:,YT0%MTL)
  IF (YT0%MTM   > 0) Z0TM   => PGMV(:,:,YT0%MTM)
  IF (YT9%MU    > 0) Z9U    => PGMV(:,:,YT9%MU)
  IF (YT9%MV    > 0) Z9V    => PGMV(:,:,YT9%MV)
  IF (YT0%MUL   > 0) Z0UL   => PGMV(:,:,YT0%MUL)
  IF (YT0%MVL   > 0) Z0VL   => PGMV(:,:,YT0%MVL)
  IF (YT0%MVOR  > 0) Z0VOR  => PGMV(:,:,YT0%MVOR)
  IF (YT0%MSPDL > 0) Z0SPDL => PGMV(:,:,YT0%MSPDL)
  IF (YT0%MSPDM > 0) Z0SPDM => PGMV(:,:,YT0%MSPDM)
  IF (YT0%MSVDL > 0) Z0SVDL => PGMV(:,:,YT0%MSVDL)
  IF (YT0%MSVDM > 0) Z0SVDM => PGMV(:,:,YT0%MSVDM)
  IF (YT0%MNHXL > 0) Z0NHXL => PGMV(:,:,YT0%MNHXL)
  IF (YT0%MNHXM > 0) Z0NHXM => PGMV(:,:,YT0%MNHXM)
  IF (YT9%MDIV  > 0) Z9DIV  => PGMV(:,:,YT9%MDIV)
  IF (YT9%MTL   > 0) Z9TL   => PGMV(:,:,YT9%MTL)
  IF (YT9%MTM   > 0) Z9TM   => PGMV(:,:,YT9%MTM)
  IF (YT9%MSPDL > 0) Z9SPDL => PGMV(:,:,YT9%MSPDL)
  IF (YT9%MSPDM > 0) Z9SPDM => PGMV(:,:,YT9%MSPDM)
  IF (YT9%MSVDL > 0) Z9SVDL => PGMV(:,:,YT9%MSVDL)
  IF (YT9%MSVDM > 0) Z9SVDM => PGMV(:,:,YT9%MSVDM)
ENDIF

IF (PRESENT (PGMVS)) THEN
 IF (YT0%MSPL > 0)  Z0SPL  => PGMVS(:,YT0%MSPL)
 IF (YT0%MSPM > 0)  Z0SPM  => PGMVS(:,YT0%MSPM)
 IF (YT9%MSPL > 0)  Z9SPL  => PGMVS(:,YT9%MSPL)
 IF (YT9%MSPM > 0)  Z9SPM  => PGMVS(:,YT9%MSPM)
ENDIF

IF (PRESENT (P0U   )) Z0U    => P0U   
IF (PRESENT (P0V   )) Z0V    => P0V   
IF (PRESENT (P0DIV )) Z0DIV  => P0DIV 
IF (PRESENT (P0TL  )) Z0TL   => P0TL  
IF (PRESENT (P0TM  )) Z0TM   => P0TM  
IF (PRESENT (P9U   )) Z9U    => P9U   
IF (PRESENT (P9V   )) Z9V    => P9V   
IF (PRESENT (P0UL  )) Z0UL   => P0UL  
IF (PRESENT (P0VL  )) Z0VL   => P0VL  
IF (PRESENT (P0VOR )) Z0VOR  => P0VOR 
IF (PRESENT (P0SPDL)) Z0SPDL => P0SPDL
IF (PRESENT (P0SPDM)) Z0SPDM => P0SPDM
IF (PRESENT (P0SVDL)) Z0SVDL => P0SVDL
IF (PRESENT (P0SVDM)) Z0SVDM => P0SVDM
IF (PRESENT (P0NHXL)) Z0NHXL => P0NHXL
IF (PRESENT (P0NHXM)) Z0NHXM => P0NHXM
IF (PRESENT (P9DIV )) Z9DIV  => P9DIV 
IF (PRESENT (P9TL  )) Z9TL   => P9TL  
IF (PRESENT (P9TM  )) Z9TM   => P9TM  
IF (PRESENT (P9SPDL)) Z9SPDL => P9SPDL
IF (PRESENT (P9SPDM)) Z9SPDM => P9SPDM
IF (PRESENT (P9SVDL)) Z9SVDL => P9SVDL
IF (PRESENT (P9SVDM)) Z9SVDM => P9SVDM
IF (PRESENT (P0SPL )) Z0SPL  => P0SPL 
IF (PRESENT (P0SPM )) Z0SPM  => P0SPM 
IF (PRESENT (P9SPL )) Z9SPL  => P9SPL 
IF (PRESENT (P9SPM )) Z9SPM  => P9SPM 

IF(KFLAG == 0) THEN
  ZGM => PGM(:)
ELSEIF(KFLAG == 1) THEN
  ZGM0(KST:KEN) = 1.0_JPRB/PGM(KST:KEN)
  ZGM => ZGM0(:)
ELSE
  CALL ABOR1('GPMPFC: ILLEGAL KFLAG')
ENDIF

! * GMV (unconditional t0 and t9 merged):

DO JLEV=1,KFLEV
  DO JROF=KST,KEN
    ! mapping of gmv t0
    Z0U(JROF,JLEV)=Z0U(JROF,JLEV)*ZGM(JROF)
    Z0V(JROF,JLEV)=Z0V(JROF,JLEV)*ZGM(JROF)
    Z0DIV(JROF,JLEV)=Z0DIV(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
    Z0TL(JROF,JLEV)=Z0TL(JROF,JLEV)*ZGM(JROF)
    Z0TM(JROF,JLEV)=Z0TM(JROF,JLEV)*ZGM(JROF)

    ! mapping of gmv t9
    Z9U(JROF,JLEV)=Z9U(JROF,JLEV)*ZGM(JROF)
    Z9V(JROF,JLEV)=Z9V(JROF,JLEV)*ZGM(JROF)
  ENDDO

  IF(LUVDER) THEN
    DO JROF=KST,KEN
      Z0UL(JROF,JLEV)=Z0UL(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
      Z0VL(JROF,JLEV)=Z0VL(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
    ENDDO
  ENDIF

  IF(LVOR) THEN
    DO JROF=KST,KEN
      Z0VOR(JROF,JLEV)=Z0VOR(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
    ENDDO
  ENDIF

  IF(YDDYNA%LNHDYN) THEN
    DO JROF=KST,KEN
      Z0SPDL(JROF,JLEV)=Z0SPDL(JROF,JLEV)*ZGM(JROF)
      Z0SPDM(JROF,JLEV)=Z0SPDM(JROF,JLEV)*ZGM(JROF)
      Z0SVDL(JROF,JLEV)=Z0SVDL(JROF,JLEV)*ZGM(JROF)
      Z0SVDM(JROF,JLEV)=Z0SVDM(JROF,JLEV)*ZGM(JROF)
    ENDDO

    IF (YDDYNA%LNHXDER) THEN
      DO JROF=KST,KEN
        Z0NHXL(JROF,JLEV)=Z0NHXL(JROF,JLEV)*ZGM(JROF)
        Z0NHXM(JROF,JLEV)=Z0NHXM(JROF,JLEV)*ZGM(JROF)
      ENDDO
    ENDIF
  ENDIF
ENDDO

IF(.NOT.YDDYNA%LTWOTL) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KEN
      Z9DIV(JROF,JLEV)=Z9DIV(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
      Z9TL(JROF,JLEV)=Z9TL(JROF,JLEV)*ZGM(JROF)
      Z9TM(JROF,JLEV)=Z9TM(JROF,JLEV)*ZGM(JROF)
    ENDDO

    IF(YDDYNA%LNHDYN) THEN
      DO JROF=KST,KEN
        Z9SPDL(JROF,JLEV)=Z9SPDL(JROF,JLEV)*ZGM(JROF)
        Z9SPDM(JROF,JLEV)=Z9SPDM(JROF,JLEV)*ZGM(JROF)
        Z9SVDL(JROF,JLEV)=Z9SVDL(JROF,JLEV)*ZGM(JROF)
        Z9SVDM(JROF,JLEV)=Z9SVDM(JROF,JLEV)*ZGM(JROF)
      ENDDO
    ENDIF
  ENDDO
ENDIF

! * GMVS:

DO JROF=KST,KEN
  Z0SPL(JROF)=Z0SPL(JROF)*ZGM(JROF)
  Z0SPM(JROF)=Z0SPM(JROF)*ZGM(JROF)
ENDDO

IF (.NOT.YDDYNA%LTWOTL.OR.NCURRENT_ITER > 0.AND.YDDYNA%LPC_FULL) THEN
  DO JROF=KST,KEN
    Z9SPL(JROF)=Z9SPL(JROF)*ZGM(JROF)
    Z9SPM(JROF)=Z9SPM(JROF)*ZGM(JROF)
  ENDDO
ENDIF

! * GFL:

DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%LCDERS) THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KEN
        PGFL(JROF,JLEV,YCOMP(JGFL)%MPL)=PGFL(JROF,JLEV,YCOMP(JGFL)%MPL)*ZGM(JROF)  
        PGFL(JROF,JLEV,YCOMP(JGFL)%MPM)=PGFL(JROF,JLEV,YCOMP(JGFL)%MPM)*ZGM(JROF)  
      ENDDO
    ENDDO
  ENDIF
ENDDO
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GPMPFC',1,ZHOOK_HANDLE)
END SUBROUTINE GPMPFC


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

SUBROUTINE GPMPFC(YDGMV,YDML_GCONF,YDDYN,YDDYNA,KPROMA,KFLEV,KST,KEN,KFLAG,PGM,PGMV,PGMVS,PGFL)

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(KPROMA,KFLEV,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(KPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(KPROMA,KFLEV,YDML_GCONF%YGFL%NDIM) 
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZGM(KPROMA)
LOGICAL :: LL9
INTEGER(KIND=JPIM) :: JLEV,JGFL,JROF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

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

IF (YDDYNA%LTWOTL) THEN
  LL9=(NCURRENT_ITER > 0 .AND. YDDYNA%LPC_FULL)
ELSE
  LL9=.TRUE.
ENDIF

IF(KFLAG == 0) THEN
  ZGM(KST:KEN) = PGM(KST:KEN)
ELSEIF(KFLAG == 1) THEN
  ZGM(KST:KEN) = 1.0_JPRB/PGM(KST:KEN)
ELSE
  CALL ABOR1('GPMPFC: ILLEGAL KFLAG')
ENDIF

! * GMV:

DO JLEV=1,KFLEV
  DO JROF=KST,KEN
    PGMV(JROF,JLEV,YT0%MU)=PGMV(JROF,JLEV,YT0%MU)*ZGM(JROF)
    PGMV(JROF,JLEV,YT0%MV)=PGMV(JROF,JLEV,YT0%MV)*ZGM(JROF)
    PGMV(JROF,JLEV,YT0%MDIV)=PGMV(JROF,JLEV,YT0%MDIV)*ZGM(JROF)*ZGM(JROF)
    PGMV(JROF,JLEV,YT0%MTL)=PGMV(JROF,JLEV,YT0%MTL)*ZGM(JROF)
    PGMV(JROF,JLEV,YT0%MTM)=PGMV(JROF,JLEV,YT0%MTM)*ZGM(JROF)
    IF(LUVDER) THEN
      PGMV(JROF,JLEV,YT0%MUL)=PGMV(JROF,JLEV,YT0%MUL)*ZGM(JROF)*ZGM(JROF)
      PGMV(JROF,JLEV,YT0%MVL)=PGMV(JROF,JLEV,YT0%MVL)*ZGM(JROF)*ZGM(JROF)
    ENDIF
    IF(LVOR) THEN
      PGMV(JROF,JLEV,YT0%MVOR)=PGMV(JROF,JLEV,YT0%MVOR)*ZGM(JROF)*ZGM(JROF)
    ENDIF
  ENDDO
ENDDO

IF(YDDYNA%LNHDYN) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KEN
      PGMV(JROF,JLEV,YT0%MSPDL)=PGMV(JROF,JLEV,YT0%MSPDL)*ZGM(JROF)
      PGMV(JROF,JLEV,YT0%MSPDM)=PGMV(JROF,JLEV,YT0%MSPDM)*ZGM(JROF)
      PGMV(JROF,JLEV,YT0%MSVDL)=PGMV(JROF,JLEV,YT0%MSVDL)*ZGM(JROF)
      PGMV(JROF,JLEV,YT0%MSVDM)=PGMV(JROF,JLEV,YT0%MSVDM)*ZGM(JROF)
      IF (YDDYNA%LNHXDER) THEN
        PGMV(JROF,JLEV,YT0%MNHXL)=PGMV(JROF,JLEV,YT0%MNHXL)*ZGM(JROF)
        PGMV(JROF,JLEV,YT0%MNHXM)=PGMV(JROF,JLEV,YT0%MNHXM)*ZGM(JROF)
      ENDIF
    ENDDO
  ENDDO
ENDIF

DO JLEV=1,KFLEV
  DO JROF=KST,KEN
    PGMV(JROF,JLEV,YT9%MU)=PGMV(JROF,JLEV,YT9%MU)*ZGM(JROF)
    PGMV(JROF,JLEV,YT9%MV)=PGMV(JROF,JLEV,YT9%MV)*ZGM(JROF)
    IF(.NOT.YDDYNA%LTWOTL) THEN
      PGMV(JROF,JLEV,YT9%MDIV)=PGMV(JROF,JLEV,YT9%MDIV)*ZGM(JROF)*ZGM(JROF)
      PGMV(JROF,JLEV,YT9%MTL)=PGMV(JROF,JLEV,YT9%MTL)*ZGM(JROF)
      PGMV(JROF,JLEV,YT9%MTM)=PGMV(JROF,JLEV,YT9%MTM)*ZGM(JROF)
    ENDIF
  ENDDO
ENDDO

IF(YDDYNA%LNHDYN .AND. .NOT.YDDYNA%LTWOTL) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KEN
      PGMV(JROF,JLEV,YT9%MSPDL)=PGMV(JROF,JLEV,YT9%MSPDL)*ZGM(JROF)
      PGMV(JROF,JLEV,YT9%MSPDM)=PGMV(JROF,JLEV,YT9%MSPDM)*ZGM(JROF)
      PGMV(JROF,JLEV,YT9%MSVDL)=PGMV(JROF,JLEV,YT9%MSVDL)*ZGM(JROF)
      PGMV(JROF,JLEV,YT9%MSVDM)=PGMV(JROF,JLEV,YT9%MSVDM)*ZGM(JROF)
    ENDDO
  ENDDO
ENDIF

! * GMVS:

DO JROF=KST,KEN
  PGMVS(JROF,YT0%MSPL)=PGMVS(JROF,YT0%MSPL)*ZGM(JROF)
  PGMVS(JROF,YT0%MSPM)=PGMVS(JROF,YT0%MSPM)*ZGM(JROF)
ENDDO

IF (LL9) THEN
  DO JROF=KST,KEN
    PGMVS(JROF,YT9%MSPL)=PGMVS(JROF,YT9%MSPL)*ZGM(JROF)
    PGMVS(JROF,YT9%MSPM)=PGMVS(JROF,YT9%MSPM)*ZGM(JROF)
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

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPMPFC',1,ZHOOK_HANDLE)
END SUBROUTINE GPMPFC

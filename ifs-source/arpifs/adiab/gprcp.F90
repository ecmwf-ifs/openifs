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

SUBROUTINE GPRCP(KPROMA,KSTART,KPROF,KFLEV,PQ,PQI,PQL,PQR,PQS,PQG,&
 & PCP,PR,PKAP,PGFL,KGFLTYP,LDTHERMACT)

!**** *GPRCP* - Computes Cp, R and R/Cp from Q

!     Purpose.
!     --------
!           COMPUTES CP AND R  AND R/CP FROM Q

!**   Interface.
!     ----------
!        *CALL* *GPRCP(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA               - dimensioning.
!          KSTART               - start of work.
!          KPROF                - depth of work.
!          KFLEV                - number of layers.
!          PQ(KPROMA,KFLEV)     - specific humidity.
!          PQI(KPROMA,KFLEV)    - ice.
!          PQL(KPROMA,KFLEV)    - liquid water.
!          PQR(KPROMA,KFLEV)    - rain.
!          PQS(KPROMA,KFLEV)    - snow.
!          PQG(KPROMA,KFLEV)    - graupel.

!        OUTPUT:
!          PCP(KPROMA,KFLEV)    - CP
!          PR(KPROMA,KFLEV)     - R
!          PKAP(KPROMA,KFLEV)   - KAPPA

!        Implicit arguments :  Physical constants from YOMCST
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
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-02-04

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified by Y.Seity  04-02-13 (Rain, Snow and Graupel)
!      M.Hamrud  15-Jan-2006  Revised GPRCP
!      K. Yessad (Jan 2011): more compact rewriting.
!      R. El Khatib 28-Aug-2014 Optimizations :
!       - compute R or CP only if required
!       - loop collapsing whenever possible, through pure array syntax
!      A. Geer      01-Oct-2015    For independence of observation operator in OOPS, 
!                                  allow calls without YGFL initialised. Removal
!                                  of all YGFL references will have to wait.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RD, RV, RCPD, RCPV, RCW, RCS
USE YOM_YGFL , ONLY : YGFL

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQ(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQI(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQL(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQR(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQS(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQG(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PGFL(KPROMA,KFLEV,YGFL%NDIM) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(OUT)   :: PCP(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(OUT)   :: PR(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(OUT)   :: PKAP(KPROMA,KFLEV) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KGFLTYP
LOGICAL,OPTIONAL,INTENT(IN) :: LDTHERMACT   ! To allow calls independent of YGFL 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZCP(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZR(KPROMA,KFLEV)
REAL(KIND=JPRB), ALLOCATABLE :: ZGFL_R(:), ZGFL_CP(:), ZGFL_X(:)
REAL(KIND=JPRB) ::  ZX

INTEGER(KIND=JPIM) :: JGFL, K, J
INTEGER(KIND=JPIM), ALLOCATABLE :: IACT(:),IPT(:)
INTEGER(KIND=JPIM) :: INUMACT,IGFLTYP
LOGICAL :: LLGFL,LLQ,LLQL,LLQI,LLQR,LLQS,LLQG, LLR, LLCP
LOGICAL :: LLQ_THERMACT,LLQL_THERMACT,LLQI_THERMACT,LLQR_THERMACT,LLQS_THERMACT,LLQG_THERMACT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#ifdef __INTEL_COMPILER
INTEGER, PARAMETER :: JPREFETCH = 3
#else
INTEGER, PARAMETER :: JPREFETCH = 0
#endif

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPRCP',0,ZHOOK_HANDLE)

ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & YG=>YGFL%YG, YI=>YGFL%YI, YL=>YGFL%YL, YQ=>YGFL%YQ, YR=>YGFL%YR, &
 & YS=>YGFL%YS)

!     ------------------------------------------------------------------

!*       1.    COMPUTES R AND CP AND KAPPA.
!              ----------------------------

LLGFL = PRESENT(PGFL)
IF(LLGFL) THEN
  ALLOCATE(ZGFL_R(NUMFLDS))
  ALLOCATE(ZGFL_CP(NUMFLDS))
  ALLOCATE(ZGFL_X(NUMFLDS))
  ALLOCATE(IACT(NUMFLDS))
  ALLOCATE(IPT(NUMFLDS))
ENDIF
IF(.NOT. LLGFL) THEN
  LLQ  = PRESENT(PQ)
  LLQL = PRESENT(PQL)
  LLQI = PRESENT(PQI)
  LLQR = PRESENT(PQR)
  LLQS = PRESENT(PQS)
  LLQG = PRESENT(PQG)
  IF(PRESENT(LDTHERMACT)) THEN
    LLQ_THERMACT=LDTHERMACT
    LLQL_THERMACT=LDTHERMACT
    LLQI_THERMACT=LDTHERMACT
    LLQR_THERMACT=LDTHERMACT
    LLQS_THERMACT=LDTHERMACT
    LLQG_THERMACT=LDTHERMACT
  ELSE
    LLQ_THERMACT=YQ%LTHERMACT
    LLQL_THERMACT=YL%LTHERMACT
    LLQI_THERMACT=YI%LTHERMACT
    LLQR_THERMACT=YR%LTHERMACT
    LLQS_THERMACT=YS%LTHERMACT
    LLQG_THERMACT=YG%LTHERMACT
  ENDIF
ENDIF

LLR=PRESENT(PR).OR.PRESENT(PKAP)
LLCP=PRESENT(PCP).OR.PRESENT(PKAP)

! * compute IPT:
IF(LLGFL) THEN
  IGFLTYP = 0
  IF(PRESENT(KGFLTYP)) IGFLTYP=KGFLTYP
  INUMACT = 0
  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LTHERMACT) THEN
      INUMACT = INUMACT+1
      IACT(INUMACT) = JGFL
      IF(IGFLTYP == 0) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP
      ELSEIF(IGFLTYP == 1) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP1
      ELSEIF(IGFLTYP == 5) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP5
      ELSEIF(IGFLTYP == 9) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP9_PH
      ELSEIF(IGFLTYP == 101) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP_SL1
      ELSE
        CALL ABOR1('GPRCP:UNKNOWN GFL TYPE')
      ENDIF
    ENDIF
#ifdef __INTEL_COMPILER
    IF (INUMACT < JPREFETCH) THEN
       ! Prefetch first few layers of PGFL to L2
       DO J=1,KFLEV
          DO K=KSTART,KPROF,8
             CALL MM_PREFETCH(PGFL(K,J,IPT(INUMACT)), 0)
          ENDDO
       ENDDO
    ENDIF
#endif
  ENDDO
ENDIF

! * compute ZR,ZCP:
IF(LLGFL) THEN
  IF(LLR) THEN
     IF(INUMACT == 0) THEN
        ZR(KSTART:KPROF,:) = RD
     ELSE
        !     Does not vectorize:
        DO JGFL=1,INUMACT
           ZGFL_X(JGFL) = YCOMP(IACT(JGFL))%R
        ENDDO
        !     Vectorizes:
        ZGFL_R(1:INUMACT) = ZGFL_X(1:INUMACT)-RD
        DO J=1,KFLEV
           ZR(KSTART:KPROF,J) = ZGFL_R(1)*PGFL(KSTART:KPROF,J,IPT(1))
           DO JGFL=2,INUMACT
              ! CALL ACCUMULATE(ZGFL_R(JGFL),PGFL(1,1,IPT(JGFL)),ZR)
              DO K=KSTART,KPROF
                 ZR(K,J) = ZR(K,J)+ZGFL_R(JGFL)*PGFL(K,J,IPT(JGFL))
              ENDDO
           ENDDO
           ZR(KSTART:KPROF,J) = RD+ZR(KSTART:KPROF,J)
        ENDDO
     ENDIF
  ENDIF
  IF(LLCP) THEN
     IF(INUMACT == 0) THEN
        ZCP(KSTART:KPROF,:) = RCPD
     ELSE
        !     Does not vectorize:
        DO JGFL=1,INUMACT
           ZGFL_X(JGFL) = YCOMP(IACT(JGFL))%RCP
        ENDDO
        !     Vectorizes:
        ZGFL_CP(1:INUMACT) = ZGFL_X(1:INUMACT)-RCPD
        DO J=1,KFLEV
           ZCP(KSTART:KPROF,J) = ZGFL_CP(1)*PGFL(KSTART:KPROF,J,IPT(1))
           DO JGFL=2,INUMACT
              !        CALL ACCUMULATE(ZGFL_CP(JGFL),PGFL(1,1,IPT(JGFL)),ZCP)
              DO K=KSTART,KPROF
                 ZCP(K,J) = ZCP(K,J)+ZGFL_CP(JGFL)*PGFL(K,J,IPT(JGFL))
              ENDDO
           ENDDO
           ZCP(KSTART:KPROF,J) = RCPD+ZCP(KSTART:KPROF,J)
        ENDDO
     ENDIF
  ENDIF
ELSE
  IF(LLR) THEN
    ZR(KSTART:KPROF,:)  = RD
  ENDIF 
  IF(LLCP) THEN
    ZCP(KSTART:KPROF,:) = RCPD
  ENDIF 
  IF(LLQ .AND. LLQ_THERMACT) THEN
    IF(LLR) THEN
      ZX=RV-RD
      CALL ACCUMULATE(ZX,PQ,ZR)
    ENDIF 
    IF(LLCP) THEN
      ZX=RCPV-RCPD
      CALL ACCUMULATE(ZX,PQ,ZCP)
    ENDIF 
  ENDIF
  IF(LLQL .AND. LLQL_THERMACT) THEN
    IF(LLR) THEN
      ZX=-RD
      CALL ACCUMULATE(ZX,PQL,ZR)
    ENDIF 
    IF(LLCP) THEN
      ZX=RCW-RCPD
      CALL ACCUMULATE(ZX,PQL,ZCP)
    ENDIF 
  ENDIF
  IF(LLQI .AND. LLQI_THERMACT) THEN
    IF(LLR) THEN
      ZX=-RD
      CALL ACCUMULATE(ZX,PQI,ZR)
    ENDIF 
    IF(LLCP) THEN
      ZX=RCS-RCPD
      CALL ACCUMULATE(ZX,PQI,ZCP)
    ENDIF 
  ENDIF
  IF(LLQR .AND. LLQR_THERMACT) THEN
    IF(LLR) THEN
      ZX=-RD
      CALL ACCUMULATE(ZX,PQR,ZR)
    ENDIF 
    IF(LLCP) THEN
      ZX=RCW-RCPD
      CALL ACCUMULATE(ZX,PQR,ZCP)
    ENDIF 
  ENDIF
  IF(LLQS .AND. LLQS_THERMACT) THEN
    IF(LLR) THEN
      ZX=-RD
      CALL ACCUMULATE(ZX,PQS,ZR)
    ENDIF 
    IF(LLCP) THEN
      ZX=RCS-RCPD
      CALL ACCUMULATE(ZX,PQS,ZCP)
    ENDIF 
  ENDIF
  IF(LLQG .AND. LLQG_THERMACT) THEN
    IF(LLR) THEN
      ZX=-RD
      CALL ACCUMULATE(ZX,PQG,ZR)
    ENDIF 
    IF(LLCP) THEN
      ZX=RCS-RCPD
      CALL ACCUMULATE(ZX,PQG,ZCP)
    ENDIF 
  ENDIF
ENDIF

! * fill PR,PCP,PKAP:
IF (PRESENT(PR)) PR(KSTART:KPROF,:)=ZR(KSTART:KPROF,:)
IF (PRESENT(PCP)) PCP(KSTART:KPROF,:)=ZCP(KSTART:KPROF,:)
IF (PRESENT(PKAP)) PKAP(KSTART:KPROF,:) = ZR(KSTART:KPROF,:)/ZCP(KSTART:KPROF,:)

IF(LLGFL) THEN
 DEALLOCATE(ZGFL_R,ZGFL_CP,ZGFL_X,IACT,IPT)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPRCP',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE ACCUMULATE(PX,PGFL,POUT)
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
REAL(KIND=JPRB), INTENT(IN) :: PX
REAL(KIND=JPRB), INTENT(IN) :: PGFL(KPROMA,KFLEV)
REAL(KIND=JPRB), INTENT(INOUT) :: POUT(KPROMA,KFLEV)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GPRCP:ACCUMULATE',0,ZHOOK_HANDLE)


IF (KSTART==1 .AND. KPROF==KPROMA) THEN
! Implicit loop collapsing
  POUT(:,:) = POUT(:,:)+PX*PGFL(:,:)
ELSE
! Nested loops : avoid uninitialized PGFL on the edge
  POUT(KSTART:KPROF,:) = POUT(KSTART:KPROF,:)+PX*PGFL(KSTART:KPROF,:)
ENDIF
IF (LHOOK) CALL DR_HOOK('GPRCP:ACCUMULATE',1,ZHOOK_HANDLE)
END SUBROUTINE ACCUMULATE

END SUBROUTINE GPRCP

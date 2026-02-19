! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE QSUPERSATCLIP(YDECLDP,KIDIA,KFDIA,KLON,KLEV,PT,PA,PAP,PQ)
!**** *QSUPERSATCLIP * - Remove unphysical supersaturation

!     PURPOSE.
!     --------
!       Perturbed initial conditions may be supersaturated unphysically.
!       This supersaturation is removed at the initial time step only.

!**   INTERFACE.
!     ----------
!        *CALL* *QSUPERSATCLIP(...)*

!     INPUT ARGUMENTS.

!     KIDIA   : start of horizontal loop
!     KFDIA   : start of horizontal loop
!     KLON    : horizontal dimension
!     KLEV    : end of vertical loop and vertical dimension
!     PT      : temperature
!     PA      : cloud fraction
!     PAP     : pressure on full levels
!     PQ      : spec. humidity

!     OUTPUT ARGUMENTS.

!     PQ      : spec. humidity

!     EXTERNALS.  NONE
!     ---------  

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!     19-09-2008 M. Leutbecher

!     MODIFICATIONS.
!     --------------
!     DD-MM-YYYY 

!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RETV     ,RTT      ,RLSTT    ,RLVTT
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                    R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &                    RALVDCP  ,RALSDCP  ,RTWAT    ,&
 &                    RTICE    ,RTICECU  ,&
 &                    RTWAT_RTICE_R      ,RTWAT_RTICECU_R  ,&
 &                    RKOOP1   ,RKOOP2
USE YOECLDP  , ONLY : TECLDP

IMPLICIT NONE

TYPE(TECLDP)      ,INTENT(INOUT) :: YDECLDP
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQ(KLON,KLEV)
INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) :: ZQP1, ZQS, ZTP1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcttre.func.h"
#include "fccld.func.h"

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('QSUPERSATCLIP',0,ZHOOK_HANDLE)
ASSOCIATE(NSSOPT=>YDECLDP%NSSOPT)
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

!     Supersaturation check

    ZTP1=PT(JL,JK)
    ZQP1=PQ(JL,JK)

    ZQS=FOEEWM(ZTP1)/PAP(JL,JK)
    ZQS=MIN(0.5_JPRB,ZQS)
    ZQS=ZQS/(1.0_JPRB-RETV*ZQS)

    !----------------------------
    ! supersaturation adjustments
    !----------------------------
    IF (ZTP1<RTICE.AND.NSSOPT>0) &
    & ZQS=ZQS*(PA(JL,JK)+FOKOOP(ZTP1)*(1.0_JPRB-PA(JL,JK)))

    IF (ZQP1 > ZQS) THEN
      PQ(JL,JK)= ZQS 
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('QSUPERSATCLIP',1,ZHOOK_HANDLE)
END SUBROUTINE QSUPERSATCLIP

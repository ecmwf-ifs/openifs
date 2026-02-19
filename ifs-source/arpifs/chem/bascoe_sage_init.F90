! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_SAGE_INIT
!**   DESCRIPTION 
!     ----------
!
!   Part of BASCOE / TM5 routines for IFS chemistry: 
!     AUTHOR.
!     -------
!        Coded in C-IFS by VINCENT HUIJNEN    *KNMI*
!        Original code from BASCOE_CTM v4s09, simonc@oma.be, June 2008
!
!-----------------------------------------------------------------------
USE BASCOE_MODULE      , ONLY :  IALT_SAGE,ILAT_SAGE,Y2AREA_LAT,LAT_SAGE, AREA_LAT, &
        & E_SAGE,ZTROP_SAGE, ALT_SAGE
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK


IMPLICIT NONE
!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM)             :: JK,JL
REAL(KIND=JPRB)                :: ZLATWORK(ILAT_SAGE)
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_INIT',0,ZHOOK_HANDLE )



DO JK=1,IALT_SAGE
  DO JL=1,ILAT_SAGE
      !!  S(Z,ZT,E0)=((0.8*(Z-ZT)/20.0)+0.7)*E0*1.0E-6
      AREA_LAT(JL,JK)=((0.8_JPRB*(ALT_SAGE(JK)-ZTROP_SAGE(JL))/20.0_JPRB)+0.7_JPRB)* &
                     & E_SAGE(JL,JK)*1.0E-6_JPRB
  ENDDO
  CALL SPLINE(LAT_SAGE,AREA_LAT(1,JK),ILAT_SAGE, &
     &        ZLATWORK,Y2AREA_LAT(1,JK))
ENDDO


IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_INIT',1,ZHOOK_HANDLE )

CONTAINS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE SPLINE(PX,PY,KN,PWORK,PY2)
!     Routine to calculate 2.nd derivatives of tabulated function
!     Y(i)=Y(Xi), to be used for cubic spline calculation.

USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK


IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KN
REAL(KIND=JPRB), INTENT(IN)    :: PX(KN),PY(KN)
REAL(KIND=JPRB), INTENT(OUT)   :: PWORK(KN),PY2(KN)
!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
!VH INTEGER(KIND=JPIM)             :: KLO,KHI
INTEGER(KIND=JPIM)             :: JI
REAL(KIND=JPRB)                :: ZSIG,ZP,ZQN,ZUN,ZYP1,ZYPN
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_INIT:SPLINE',0,ZHOOK_HANDLE )
!
      ZYP1=(PY(2)-PY(1))/(PX(2)-PX(1))
      ZYPN=(PY(KN)-PY(KN-1))/(PX(KN)-PX(KN-1))
      IF(ZYP1 > 99.0E+30_JPRB) THEN
          PY2(1)=0.0
          PWORK(1)=0.0
      ELSE
          PY2(1)=-0.5E0
          PWORK(1)=(3.0E0/(PX(2)-PX(1)))*((PY(2)-PY(1))/(PX(2)-PX(1))-ZYP1)
      ENDIF
      DO JI=2,KN-1
          ZSIG=(PX(JI)-PX(JI-1))/(PX(JI+1)-PX(JI-1))
          ZP=ZSIG*PY2(JI-1)+2.0E0
          PY2(JI)=(ZSIG-1.0E0)/ZP
          PWORK(JI)=(6.0E0*((PY(JI+1)-PY(JI))/(PX(JI+1)-PX(JI))-(PY(JI)-PY(JI-1)) &
     &             /(PX(JI)-PX(JI-1)))/(PX(JI+1)-PX(JI-1))-ZSIG*PWORK(JI-1))/ZP
      ENDDO
      IF(ZYPN > 99.0E+30_JPRB) THEN
          ZQN=0.0
          ZUN=0.0
      ELSE
          ZQN=0.5E0_JPRB
          ZUN=(3.0E0_JPRB/(PX(KN)-PX(KN-1)))*(ZYPN-(PY(KN)-PY(KN-1))/(PX(KN)-PX(KN-1)))
      ENDIF
      PY2(KN)=(ZUN-ZQN*PWORK(KN-1))/(ZQN*PY2(KN-1)+1.0E0_JPRB)
      DO JI=KN-1,1,-1
          PY2(JI)=PY2(JI)*PY2(JI+1)+PWORK(JI)
      ENDDO
!

IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_INIT:SPLINE',1,ZHOOK_HANDLE )
END SUBROUTINE SPLINE


END SUBROUTINE BASCOE_SAGE_INIT

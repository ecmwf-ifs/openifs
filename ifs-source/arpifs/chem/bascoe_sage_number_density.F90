! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_SAGE_NUMBER_DENSITY(PAIR,PLAT,PMEDIAN,PSTDEV, PND)
!
!     Subroutine used to calculate aerosol number densities consistent
!     with SAGE II averaged extinction measurements, assuming constant
!     lognormal size distribution median radius and geometric standard
!     deviation.
!     Data source: Hitchman et al.,  JGR 99, 20689, 1994.
!
!    Routine is integrated with SAGE_AREA from BASCOE...
!
!     Input:
!               PAIR:    (Pa)        Air pressure
!               LATITUDE (deg)       Latitude [-90,+90]
!               MEDIAN   (m)         Lognormal median radius
!               GSTDEV               Lognormal geometric standar deviation
!     Intermediate:
!               AREA:    (m**2/m**3) Sulfate aerosol surface area density
!     Output:
!               ND       (m**-3)     Aerosol number density
!
USE PARKIND1  ,    ONLY : JPIM,   JPRB
USE YOMHOOK   ,    ONLY : LHOOK,  DR_HOOK, JPHOOK
USE BASCOE_MODULE, ONLY : ILAT_SAGE,IALT_SAGE,LAT_SAGE,ALT_SAGE, &
            & AREA_LAT,Y2AREA_LAT

IMPLICIT NONE
REAL(KIND=JPRB), INTENT(IN)  :: PAIR,PLAT,PMEDIAN,PSTDEV
REAL(KIND=JPRB), INTENT(OUT) :: PND

! Local
INTEGER(KIND=JPIM) :: JA
REAL(KIND=JPRB), PARAMETER :: ZH0=6.5, ZP0=1000.0E2_JPRB
REAL(KIND=JPRB)    :: ZFACTOR,ZAREA,ZHEIGHT
REAL(KIND=JPRB)    :: Z2AREA_ALT(IALT_SAGE),ZALTWORK(IALT_SAGE)
REAL(KIND=JPRB)    :: ZAREA_ALT(IALT_SAGE)
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE


IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_NUMBER_DENSITY',0,ZHOOK_HANDLE )
!----------------------------------------------------------------------------

ZFACTOR=4.0*4.0*ATAN(1.0)*(PMEDIAN**2)*EXP(2.0*(LOG(PSTDEV))**2)


      ZHEIGHT = ZH0 * LOG(ZP0/PAIR)
      ZHEIGHT = MAX(ALT_SAGE(1),MIN(ZHEIGHT,ALT_SAGE(IALT_SAGE)))

      DO JA=1,IALT_SAGE
          CALL SPLINT(LAT_SAGE,AREA_LAT(1,JA),Y2AREA_LAT(1,JA),ILAT_SAGE, &
     &                PLAT,ZAREA_ALT(JA))
      ENDDO
      CALL SPLINE(ALT_SAGE,ZAREA_ALT,IALT_SAGE,ZALTWORK,Z2AREA_ALT)
      CALL SPLINT(ALT_SAGE,ZAREA_ALT,Z2AREA_ALT,IALT_SAGE,ZHEIGHT,ZAREA)


PND=ZAREA/ZFACTOR

IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_NUMBER_DENSITY',1,ZHOOK_HANDLE )

CONTAINS

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
IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_NUMBER_DENSITY:SPLINE',0,ZHOOK_HANDLE )
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

IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_NUMBER_DENSITY:SPLINE',1,ZHOOK_HANDLE )
END SUBROUTINE SPLINE


SUBROUTINE SPLINT(PXA,PYA,PY2A,KN,PX,PY)
!******************************************************************************
!     Cubic spline calculation
!
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN):: KN
REAL(KIND=JPRB), INTENT(IN)   :: PXA(KN),PYA(KN),PY2A(KN), PX
REAL(KIND=JPRB), INTENT(OUT)  :: PY
!* Local
REAL(KIND=JPRB)               :: ZH,ZA,ZB
INTEGER(KIND=JPIM)            :: JKLO,JKHI,JK
REAL(KIND=JPHOOK)               :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_NUMBER_DENSITY:SPLINT',0,ZHOOK_HANDLE )

      JKLO=1_JPIM
      JKHI=KN
      DO WHILE(JKHI-JKLO > 1_JPIM) 
          JK=(JKHI+JKLO)/2_JPIM
          IF(PXA(JK) > PX) THEN
              JKHI=JK
          ELSE
              JKLO=JK
          ENDIF
          !GOTO 1
      ENDDO
      ZH=PXA(JKHI)-PXA(JKLO)
      ZA=(PXA(JKHI)-PX)/ZH
      ZB=(PX-PXA(JKLO))/ZH
      PY=ZA*PYA(JKLO)+ZB*PYA(JKHI)+ &
     &        ((ZA**3-ZA)*PY2A(JKLO)+(ZB**3-ZB)*PY2A(JKHI))*(ZH**2)/6.0E0_JPRB
!
IF (LHOOK) CALL DR_HOOK('BASCOE_SAGE_NUMBER_DENSITY:SPLINT',1,ZHOOK_HANDLE )
END SUBROUTINE SPLINT

END SUBROUTINE BASCOE_SAGE_NUMBER_DENSITY


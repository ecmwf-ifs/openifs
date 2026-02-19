! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DUCTDIA ( KIDIA, KFDIA, KLON, KLEV, &
                 &   PT, PQ, PAP, PGEOM1, &
                 &   PREFRAC, PDNDZ_MIN, PDNDZ_AVG, &
                 &   PZDCBOT, PZTLBOT, PZTLTOP)   
     
!**** COMPUTES REFRACTIVITY PROFILE FROM T, Q AND PS INFORMATION
!     AND MINIMUM REFRACTIVITY GRADIENT IN THE LOWER TROPOSPHERE.


!     P. LOPEZ   ECMWF    OCTOBER 2008

!     DUCTDIA IS CALLED FROM CALLPAR

!-----------------------------------------------------------------------------------

! - 3D INPUT
!     PT         : TEMPERATURE                            K
!     PQ         : SPECIFIC HUMIDITY                      KG/KG
!     PAP        : PRESSURE ON FULL LEVELS                Pa
!     PGEOM1     : GEOPOTENTIAL ON FULL LEVELS            M^2/S^2

! - 2D OUTPUT
!     PDNDZ_MIN  : MINIMUM REFRACTIVITY GRADIENT 
!                  INSIDE TRAPPING LAYER                  M^-1
!     PDNDZ_AVG  : AVERAGE REFRACTIVITY GRADIENT 
!                  INSIDE TRAPPING LAYER                  M^-1
!     PZDCBOT    : DUCT BASE HEIGHT                       M
!     PZTLBOT    : TRAPPING LAYER BASE HEIGHT             M
!     PZTLTOP    : TRAPPING LAYER TOP HEIGHT              M

! - 3D OUTPUT
!     PREFRAC    : ATMOSPHERIC REFRACTIVITY ON MODEL LEVELS

!-----------------------------------------------------------------------------------
!     METHOD:

!     ATMOSPHERIC REFRACTIVITY, N, IS COMPUTED AS A FUNCTION OF PRESSURE (Pa), 
!     TEMPERATURE (K) AND WATER VAPOUR PARTIAL PRESSURE, E (Pa):

!                       N = 0.776 P/T + 3730. E/T^2 

!     MODIFIED REFRACTIVITY, M, IS DEFINED AS:

!                         M = N + 0.157 Z

!     WHERE Z IS THE ALTITUDE (M), TO ACCOUNT FOR THE EARTH'S CURVATURE.

!     A TRAPPING LAYER IS ASSUMED TO BE PRESENT WHENEVER THE VERTICAL GRADIENT 
!     OF ATMOSPHERIC REFRACTIVITITY (DN/DZ) IS LOWER THAN -0.157 M^-1 
!     (STARTING FROM THE LOWEST MODEL LEVEL AND UP).

!     TRAPPING LAYER BASE HEIGHT (PZTLBOT) IS DEFINED BY THE LOWEST MODEL LEVEL 
!     ABOVE WHICH DN/DZ <= -0.157 M^-1.
!     NOTE: PZTLBOT IS RESET TO ZERO (=SURFACE) WHEN TRAPPING LAYER BASE IS EQUAL TO
!     THE MODEL LEVEL JUST ABOVE THE SURFACE.

!     ONCE A TRAPPING LAYER IS IDENTIFIED, ITS TOP HEIGHT (PZTLTOP) IS DEFINED 
!     AS THE BOTTOM OF THE FIRST MODEL LAYER INSIDE WHICH DN/DZ > -0.157 M^-1.

!     THE DETECTION OF TRAPPING LAYERS IS RESTRICTED TO THOSE WHOSE BASE LIES 
!     WITHIN THE LOWEST 2.5 KM OF THE TROPOSPHERE (TO ACCOUNT FOR VERTICAL 
!     RESOLUTION DEGRADATION ABOVE). 
!     NOTE THAT TRAPPING LAYER TOP HEIGHT CAN BE ABOVE 2.5 KM.

!     THE "DUCT" IS DEFINED AS THE LAYER THAT EXTENDS BETWEEN THE TOP OF THE 
!     TRAPPING LAYER (PZTLTOP) AND THE MODEL LEVEL (HEIGHT PZDCBOT) DEFINED BY:

!       IF ZM_SURF > ZM_MIN, PZDCBOT = 0,
!       ELSE PZDCBOT IS EQUAL TO THE LOWEST MODEL LEVEL AT WHICH ZM = ZM_MIN. 

!     WHERE ZM_MIN IS THE VALUE OF M AT THE TOP OF THE TRAPPING LAYER AND 
!     ZM_SURF ITS VALUE AT THE LOWEST MODEL LEVEL.

!     REFERENCE: TURTON ET AL. (METEOROL. MAG., 1988).

!-----------------------------------------------------------------------------------
!     MODIFICATIONS:

!-----------------------------------------------------------------------------------
!**** MODULES

!* KIND     
USE PARKIND1     , ONLY: JPIM, JPRB
USE YOMHOOK      , ONLY: LHOOK, DR_HOOK, JPHOOK

!* IFS
USE YOMCST       , ONLY: RG, RETV
             
IMPLICIT NONE                         

!**** VARIABLES

!* INPUT
INTEGER(KIND=JPIM), INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM), INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM), INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM), INTENT(IN)    :: KFDIA 
REAL (KIND=JPRB)  , INTENT(IN)    :: PT(KLON,KLEV) 
REAL (KIND=JPRB)  , INTENT(IN)    :: PQ(KLON,KLEV) 
REAL (KIND=JPRB)  , INTENT(IN)    :: PAP(KLON,KLEV) 
REAL (KIND=JPRB)  , INTENT(IN)    :: PGEOM1(KLON,KLEV) 
!* OUTPUT
REAL (KIND=JPRB)  , INTENT(OUT)   :: PREFRAC(KLON,KLEV)
REAL (KIND=JPRB)  , INTENT(OUT)   :: PDNDZ_MIN(KLON)
REAL (KIND=JPRB)  , INTENT(OUT)   :: PDNDZ_AVG(KLON)
REAL (KIND=JPRB)  , INTENT(OUT)   :: PZDCBOT(KLON)
REAL (KIND=JPRB)  , INTENT(OUT)   :: PZTLBOT(KLON)
REAL (KIND=JPRB)  , INTENT(OUT)   :: PZTLTOP(KLON)

!* LOCAL
REAL (KIND=JPRB) :: ZZMAX, ZE, ZN, ZN_OLD, ZMISSH, ZMISSN
REAL (KIND=JPRB) :: ZDNDZ, ZDNDZ_MIN, ZRG, ZM_MIN, ZM_SURF
REAL (KIND=JPRB) :: ZM(KLEV), ZZ(KLEV)

INTEGER (KIND=JPIM) :: JL, JK, ITLBOT, ITLTOP

LOGICAL :: LLDUCT

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK ('DUCTDIA',0,ZHOOK_HANDLE)

! Constants
ZRG=1.0_JPRB/RG
ZZMAX=2500._JPRB
!LLL ZMISSH=-9999.0_JPRB
!LLL ZMISSN=-9999.0_JPRB
ZMISSH=-1.0_JPRB
ZMISSN=-1.0_JPRB

! Initializations
PREFRAC(:,:)=0.0_JPRB
PDNDZ_MIN(:)=ZMISSN
PDNDZ_AVG(:)=ZMISSN
PZTLBOT  (:)=ZMISSH
PZTLTOP  (:)=ZMISSH
PZDCBOT  (:)=ZMISSH

DO JL = KIDIA, KFDIA

! Initializations
  ZDNDZ_MIN = 1.E+12_JPRB
  ZM_MIN = 1.E+12_JPRB
  LLDUCT=.FALSE.
  ZM(:)=-9999._JPRB

! Model level heights
  DO JK=1,KLEV
    ZZ(JK) = PGEOM1(JL,JK)*ZRG
  ENDDO

! Water vapour partial pressure at lowest model level
  ZE = PQ(JL,KLEV) * PAP(JL,KLEV) *(RETV + 1.0_JPRB) / (1.0_JPRB + RETV * PQ(JL,KLEV))
! Refractivity at 2m level (see e.g. Steiner and Smith, 2002)
  ZN_OLD = 0.776_JPRB * PAP(JL,KLEV) / PT(JL,KLEV) + 3730._JPRB * ZE / PT(JL,KLEV)**2
  PREFRAC(JL,KLEV) = ZN_OLD
! Modified refractivity
  ZM(KLEV) = ZN_OLD + 0.157_JPRB * ZZ(KLEV)

  DO JK = KLEV-1, 1, -1
! Water vapour partial pressure
    ZE = PQ(JL,JK) * PAP(JL,JK) *(RETV + 1.0_JPRB) / (1.0_JPRB + RETV * PQ(JL,JK))
! Refractivity (see e.g. Steiner and Smith, 2002)
    ZN = 0.776_JPRB * PAP(JL,JK) / PT(JL,JK) + 3730._JPRB * ZE / PT(JL,JK)**2
! Store refractivity profile
    PREFRAC(JL,JK) = ZN
! Refractivity gradient
    ZDNDZ = (ZN - ZN_OLD)/(ZZ(JK) - ZZ(JK+1))
! Modified refractivity
    ZM(JK) = ZN + 0.157_JPRB * ZZ(JK)
! Bottom height of trapping layer (only if below ZZMAX)
    IF (ZZ(JK) <= ZZMAX .AND. ZDNDZ <= -0.157_JPRB .AND. .NOT.LLDUCT) THEN
      PZTLBOT(JL) = ZZ(JK+1)
      ITLBOT = JK+1
! Reset trapping layer base height to zero if lowest model level
      IF (ITLBOT == KLEV) PZTLBOT(JL) = 0.0_JPRB
      LLDUCT=.TRUE.
    ENDIF
! Top height of trapping layer
    IF (ZDNDZ > -0.157_JPRB .AND. LLDUCT .AND. PZTLTOP(JL) == ZMISSH) THEN
      PZTLTOP(JL) = ZZ(JK+1)
      ITLTOP = JK+1
    ENDIF
! Minimum refractivity gradient and minimum modified refractivity inside trapping layer
    IF (LLDUCT .AND. PZTLBOT(JL) >= 0.0_JPRB .AND. PZTLTOP(JL) == ZMISSH) THEN
      ZDNDZ_MIN = MIN(ZDNDZ_MIN, ZDNDZ)
      ZM_MIN = MIN(ZM_MIN, ZM(JK))
    ENDIF
! Store current value of refractivity for use at level above
    ZN_OLD = ZN
  ENDDO

  ZM_SURF = ZM(KLEV)

! Bottom height of duct
  IF (LLDUCT) THEN
    PZDCBOT(JL)=0.0_JPRB
    IF (ZM_SURF < ZM_MIN) THEN
      DO JK=KLEV-1,ITLBOT,-1
        IF (ZM(JK) > ZM_MIN .AND. ZM(JK+1) <= ZM_MIN .AND. PZDCBOT(JL) == 0.0_JPRB) THEN
          PZDCBOT(JL) = ZZ(JK+1) + (ZZ(JK)-ZZ(JK+1))*(ZM_MIN-ZM(JK+1))/(ZM(JK)-ZM(JK+1))
          PZDCBOT(JL) = MAX(0.0_JPRB,PZDCBOT(JL))
        ENDIF
      ENDDO
    ENDIF
! Store minimum gradient below trapping layer top (-dN/dz is stored; always positive)
    PDNDZ_MIN(JL) = -ZDNDZ_MIN
! Mean refractivity gradient inside trapping layer (-dN/dz is stored; always positive)
    PDNDZ_AVG(JL) = -(PREFRAC(JL,ITLTOP)-PREFRAC(JL,ITLBOT))/(PZTLTOP(JL)-PZTLBOT(JL))
  ENDIF

ENDDO 

IF (LHOOK) CALL DR_HOOK ('DUCTDIA',1,ZHOOK_HANDLE)
END SUBROUTINE DUCTDIA

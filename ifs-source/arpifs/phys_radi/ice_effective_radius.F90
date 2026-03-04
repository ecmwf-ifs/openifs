! (C) Copyright 2016- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
SUBROUTINE ICE_EFFECTIVE_RADIUS &
     & (YDERAD,YDSPP_CONFIG,KIDIA, KFDIA, KLON, KLEV, &
     &  PPRESSURE, PTEMPERATURE, PCLOUD_FRAC, PQ_ICE, PQ_SNOW, PGEMU, &
     &  PRE_UM, PPERT)

! ICE_EFFECTIVE_RADIUS
!
! PURPOSE
! -------
!   Calculate effective radius of ice clouds
!
! AUTHOR
! ------
!   Robin Hogan, ECMWF (using code extracted from radlswr.F90)
!   Original: 2016-02-24
!
! MODIFICATIONS
! -------------
!
!   M Leutbecher  (Oct 2020) SPP abstraction
!
! -------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOERAD   , ONLY : TERAD
USE YOMLUN   , ONLY : NULERR
USE YOMCST   , ONLY : RD, RTT
USE SPP_MOD  , ONLY : TSPP_CONFIG
USE SPP_GEN_MOD, ONLY:SPP_PERT

! -------------------------------------------------------------------

IMPLICIT NONE

! INPUT ARGUMENTS

! *** Array dimensions and ranges
TYPE(TERAD)       ,INTENT(INOUT):: YDERAD
TYPE(TSPP_CONFIG) ,INTENT(IN) :: YDSPP_CONFIG
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA    ! Start column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA    ! End column to process
INTEGER(KIND=JPIM),INTENT(IN) :: KLON     ! Number of columns
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV     ! Number of levels

! *** Variables on model levels
REAL(KIND=JPRB),   INTENT(IN) :: PPRESSURE(KLON,KLEV)    ! (Pa)
REAL(KIND=JPRB),   INTENT(IN) :: PTEMPERATURE(KLON,KLEV) ! (K)
REAL(KIND=JPRB),   INTENT(IN) :: PCLOUD_FRAC(KLON,KLEV)  ! (kg/kg)
REAL(KIND=JPRB),   INTENT(IN) :: PQ_ICE(KLON,KLEV)       ! (kg/kg)
REAL(KIND=JPRB),   INTENT(IN) :: PQ_SNOW(KLON,KLEV)      ! (kg/kg)

! *** Single level variable
REAL(KIND=JPRB),   INTENT(IN) :: PGEMU(KLON) ! Sine of latitude

! OUTPUT ARGUMENT
! Effective radius
REAL(KIND=JPRB),  INTENT(OUT) :: PRE_UM(KLON,KLEV) ! (microns)

! OPTIONAL INPUT ARGUMENT
! SPP perturbation pattern
REAL(KIND=JPRB),  INTENT(IN), OPTIONAL :: PPERT(KLON,YDSPP_CONFIG%SM%NRFTOTAL_RADGRID)

! LOCAL VARIABLES

REAL(KIND=JPRB) :: ZIWC_INCLOUD_GM3 ! In-cloud ice+snow water content in g m-3
REAL(KIND=JPRB) :: ZAIR_DENSITY_GM3 ! Air density in g m-3

REAL(KIND=JPRB) :: ZTEMPERATURE_C   ! Temperature, degrees Celcius
REAL(KIND=JPRB) :: ZAIWC, ZBIWC     ! Factors in empirical relationship
REAL(KIND=JPRB) :: ZDEFAULT_RE_UM   ! Default effective radius in microns 
REAL(KIND=JPRB) :: ZDIAMETER_UM     ! Effective diameter in microns

! Min effective diameter in microns; may vary with latitude
REAL(KIND=JPRB) :: ZMIN_DIAMETER_UM(KLON)

INTEGER(KIND=JPIM) :: JL, JK

! SPP variables
LOGICAL            :: LLUSE_SPP
INTEGER(KIND=JPIM) :: IPN
TYPE(SPP_PERT)     :: PN1

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -------------------------------------------------------------------

#include "abor1.intfb.h"

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ICE_EFFECTIVE_RADIUS',0,ZHOOK_HANDLE)

! -------------------------------------------------------------------

! Do we apply SPP perturbations?
IF (YDSPP_CONFIG%LSPP) THEN
  IPN = YDSPP_CONFIG%PPTR%ZRADEFFI
  LLUSE_SPP = IPN > 0
  IF (LLUSE_SPP) THEN
    PN1 = YDSPP_CONFIG%SM%PN(IPN)
  ENDIF
ELSE
  LLUSE_SPP = .FALSE.
ENDIF

SELECT CASE(YDERAD%NRADIP)
CASE(0)
  ! Ice effective radius fixed at 40 microns
  PRE_UM(KIDIA:KFDIA,:) = 40.0_JPRB  

CASE(1,2)
  ! Ice effective radius from Liou and Ou (1994)
  DO JK = 1,KLEV
    DO JL = KIDIA,KFDIA
      ! Convert Kelvin to Celcius, preventing positive numbers
      ZTEMPERATURE_C = MIN(PTEMPERATURE(JL,JK) - RTT, -0.1_JPRB)
      ! Liou and Ou's empirical formula
      PRE_UM(JL,JK) = 326.3_JPRB + ZTEMPERATURE_C * (12.42_JPRB&
           &  + ZTEMPERATURE_C * (0.197_JPRB + ZTEMPERATURE_C * 0.0012_JPRB))
      IF (YDERAD%NRADIP == 1) THEN
        ! Original Liou and Ou (1994) bounds of 40-130 microns
        PRE_UM(JL,JK) = MAX(PRE_UM(JL,JK), 40.0_JPRB)
        PRE_UM(JL,JK) = MIN(PRE_UM(JL,JK),130.0_JPRB)
      ELSE
        ! Formulation following Jakob, Klein modifications to ice
        ! content
        PRE_UM(JL,JK) = MAX(PRE_UM(JL,JK), 30.0_JPRB)
        PRE_UM(JL,JK) = MIN(PRE_UM(JL,JK), 60.0_JPRB)
      ENDIF
    ENDDO
  ENDDO

CASE(3)
  ! Ice effective radius = f(T,IWC) from Sun and Rikus (1999), revised
  ! by Sun (2001)

  ! Default effective radius is computed from an effective diameter of
  ! 80 microns; note that multiplying by re2de actually converts from
  ! effective diameter to effective radius.
  ZDEFAULT_RE_UM = 80.0_JPRB * YDERAD%RRE2DE

  ! Minimum effective diameter may vary with latitude
  IF (YDERAD%NMINICE == 0) THEN
    ! Constant effective diameter
    ZMIN_DIAMETER_UM(KIDIA:KFDIA) = YDERAD%RMINICE
  ELSE
    ! Ice effective radius varies with latitude, smaller at poles
    DO JL = KIDIA,KFDIA
      ZMIN_DIAMETER_UM(JL) = 20.0_JPRB + (YDERAD%RMINICE - 20.0_JPRB)&
           &                          * COS(ASIN(PGEMU(JL)))
    ENDDO
  ENDIF

  DO JK = 1,KLEV
    DO JL = KIDIA,KFDIA
      IF (PCLOUD_FRAC(JL,JK) > 0.001_JPRB&
           &  .AND. (PQ_ICE(JL,JK)+PQ_SNOW(JL,JK)) > 0.0_JPRB) THEN
        ZAIR_DENSITY_GM3 = 1000.0_JPRB * PPRESSURE(JL,JK) / (RD*PTEMPERATURE(JL,JK))
        ZIWC_INCLOUD_GM3 = ZAIR_DENSITY_GM3 * (PQ_ICE(JL,JK) + PQ_SNOW(JL,JK))&
             &           / PCLOUD_FRAC(JL,JK)
        ZTEMPERATURE_C = PTEMPERATURE(JL,JK) - RTT
        ! Sun, 2001 (corrected from Sun & Rikus, 1999)
        ZAIWC = 45.8966_JPRB * ZIWC_INCLOUD_GM3**0.2214_JPRB
        ZBIWC = 0.7957_JPRB  * ZIWC_INCLOUD_GM3**0.2535_JPRB
        ZDIAMETER_UM = (1.2351_JPRB + 0.0105_JPRB * ZTEMPERATURE_C)&
             & * (ZAIWC + ZBIWC*(PTEMPERATURE(JL,JK) - 83.15_JPRB))

        ! SPP perturbations
        IF (LLUSE_SPP) THEN
           ZDIAMETER_UM = ZDIAMETER_UM * EXP(PN1%MU(1) &
                &  +PN1%XMAG(1)*PPERT(JL,PN1%MP_RADGRID))
        ENDIF

        ZDIAMETER_UM = MIN ( MAX( ZDIAMETER_UM, ZMIN_DIAMETER_UM(JL)), 155.0_JPRB)
        PRE_UM(JL,JK) = ZDIAMETER_UM * YDERAD%RRE2DE
      ELSE
        PRE_UM(JL,JK) = ZDEFAULT_RE_UM
      ENDIF
    ENDDO
  ENDDO
  
CASE DEFAULT
  WRITE(NULERR,'(A,I0,A)') 'ICE EFFECTIVE RADIUS OPTION NRADLP=',YDERAD%NRADIP,' NOT AVAILABLE'
  CALL ABOR1('ERROR IN ICE_EFFECTIVE_RADIUS')

END SELECT

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ICE_EFFECTIVE_RADIUS',1,ZHOOK_HANDLE)
  
END SUBROUTINE ICE_EFFECTIVE_RADIUS

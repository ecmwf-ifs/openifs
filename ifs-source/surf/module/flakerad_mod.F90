MODULE FLAKERAD_MOD
CONTAINS
SUBROUTINE FLAKERAD&
  & (KIDIA        , KFDIA           , LDLAKEPOINT     ,  &
  &  PDEPTH_W     , POPTICPAR_WATER , POPTICPAR_ICE   ,  &
  &  YDFLAKE      , &
  &  PI_ATM_FLK   , PH_ICE_P_FLK    , PH_ML_P_FLK     ,  &

  &  PI_ICE_FLK   , PI_BOT_FLK      , PI_W_FLK        ,  &
  &  PI_H_FLK     , PI_INTM_0_H_FLK , PI_INTM_H_D_FLK    )

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the radiation fluxes 
!  at the ice-water, air-water, 
!  mixed layer-thermocline and water column-bottom sediment interfaces,
!  the mean radiation flux over the mixed layer,
!  and the mean radiation flux over the thermocline.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! <Modifications> 
! !VERSION!  !DATE!     <Your name>
! 1.01T    27-Feb-2008  V. M. Stepanenko
! The code relevant to snow is dropped
! <End modifications>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_FLAKE, ONLY : TFLAKE, ROPTICPAR_MEDIUM

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

! Input (procedure arguments)

INTEGER (KIND = JPIM), INTENT(IN) :: KIDIA
INTEGER (KIND = JPIM), INTENT(IN) :: KFDIA

LOGICAL, INTENT(IN):: LDLAKEPOINT (:)

REAL (KIND = JPRD), INTENT(IN) :: PDEPTH_W     (:)  ! The lake depth [m]

REAL (KIND = JPRD), INTENT(IN) :: PI_ATM_FLK   (:)
REAL (KIND = JPRD), INTENT(IN) :: PH_ICE_P_FLK (:)
REAL (KIND = JPRD), INTENT(IN) :: PH_ML_P_FLK  (:)                   

TYPE (ROPTICPAR_MEDIUM), INTENT(IN) :: POPTICPAR_WATER ! Optical characteristics of water
TYPE (ROPTICPAR_MEDIUM), INTENT(IN) :: POPTICPAR_ICE   ! Optical characteristics of ice

TYPE(TFLAKE), INTENT(IN) :: YDFLAKE

REAL (KIND = JPRD), INTENT(OUT) :: PI_ICE_FLK      (:)
REAL (KIND = JPRD), INTENT(OUT) :: PI_BOT_FLK      (:)
REAL (KIND = JPRD), INTENT(OUT) :: PI_W_FLK        (:)
REAL (KIND = JPRD), INTENT(OUT) :: PI_H_FLK        (:)
REAL (KIND = JPRD), INTENT(OUT) :: PI_INTM_0_H_FLK (:)
REAL (KIND = JPRD), INTENT(OUT) :: PI_INTM_H_D_FLK (:)

!  Local variables of type REAL
REAL (KIND = JPRD):: ZDEPTH_W_MAX    ! Maximum lake depth allowed
REAL (KIND = JPRD):: ZDEPTH_W_MIN    ! Minimum lake depth allowed
REAL (KIND = JPRD):: ZDEPTH_W
!  Local variables of type INTEGER
INTEGER (KIND = JPIM) :: & ! Help variable(s)
  J,JPOINT                 ! DO loop index

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FLAKERAD_MOD:FLAKERAD',0,ZHOOK_HANDLE)
ASSOCIATE(RH_ICE_MIN_FLK=>YDFLAKE%RH_ICE_MIN_FLK, &
 & RH_ML_MIN_FLK=>YDFLAKE%RH_ML_MIN_FLK)

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

ZDEPTH_W_MAX = 50.0_JPRD              ! Maximum lake depth simulated by FLAKE
ZDEPTH_W_MIN =  2.0_JPRD              ! Minimum lake depth simulated by FLAKE

DO JPOINT = KIDIA, KFDIA
  LAKEPOINT: IF (LDLAKEPOINT(JPOINT)) THEN ! The calculations are performed only for lake
    ZDEPTH_W=MIN(ZDEPTH_W_MAX,MAX(ZDEPTH_W_MIN,PDEPTH_W(JPOINT)))
    IF(PH_ICE_P_FLK(JPOINT).GE.RH_ICE_MIN_FLK) THEN            ! Ice exists
      PI_ICE_FLK(JPOINT) = PI_ATM_FLK(JPOINT)
      PI_BOT_FLK(JPOINT) = 0._JPRD
      DO J=1, POPTICPAR_ICE%NBAND_OPTIC
        PI_BOT_FLK(JPOINT) = PI_BOT_FLK(JPOINT) +      & 
       & POPTICPAR_ICE%RFRAC_OPTIC(J)*                 &
	   & EXP(-POPTICPAR_ICE%REXTINCOEF_OPTIC(J)*PH_ICE_P_FLK(JPOINT)) 
      END DO 
      PI_W_FLK(JPOINT) = PI_ICE_FLK(JPOINT)*PI_BOT_FLK(JPOINT)
    ELSE                                                      ! No ice cover
      PI_ICE_FLK(JPOINT) = PI_ATM_FLK(JPOINT)
      PI_W_FLK  (JPOINT) = PI_ATM_FLK(JPOINT)
    END IF 

    IF(PH_ML_P_FLK(JPOINT).GE.RH_ML_MIN_FLK) THEN             ! Radiation flux at the bottom of the mixed layer
      PI_BOT_FLK(JPOINT) = 0._JPRD
      DO J=1, POPTICPAR_WATER%NBAND_OPTIC
        PI_BOT_FLK(JPOINT) = PI_BOT_FLK(JPOINT) +            & 
       & POPTICPAR_WATER%RFRAC_OPTIC(J)*                     &
	   & EXP(-POPTICPAR_WATER%REXTINCOEF_OPTIC(J)*PH_ML_P_FLK(JPOINT)) 
      END DO 
      PI_H_FLK(JPOINT) = PI_W_FLK(JPOINT)*PI_BOT_FLK(JPOINT)
    ELSE                                                      ! Mixed-layer depth is less then a minimum value
      PI_H_FLK(JPOINT) = PI_W_FLK(JPOINT)
    END IF

    PI_BOT_FLK(JPOINT) = 0._JPRD                              ! Radiation flux at the lake bottom
    DO J=1, POPTICPAR_WATER%NBAND_OPTIC
      PI_BOT_FLK(JPOINT) = PI_BOT_FLK(JPOINT) +              & 
     & POPTICPAR_WATER%RFRAC_OPTIC(J)*                       &
     & EXP(-POPTICPAR_WATER%REXTINCOEF_OPTIC(J)*ZDEPTH_W) 
    END DO 
    PI_BOT_FLK(JPOINT) = PI_W_FLK(JPOINT)*PI_BOT_FLK(JPOINT)

    IF(PH_ML_P_FLK(JPOINT).GE.RH_ML_MIN_FLK) THEN           ! Integral-mean radiation flux over the mixed layer
      PI_INTM_0_H_FLK(JPOINT) = 0._JPRD
      DO J=1, POPTICPAR_WATER%NBAND_OPTIC
        PI_INTM_0_H_FLK(JPOINT) = PI_INTM_0_H_FLK(JPOINT) +                   &
       & POPTICPAR_WATER%RFRAC_OPTIC(J)/POPTICPAR_WATER%REXTINCOEF_OPTIC(J)*  &
       & (1._JPRD - EXP(-POPTICPAR_WATER%REXTINCOEF_OPTIC(J)*PH_ML_P_FLK(JPOINT)) )
      END DO 
      PI_INTM_0_H_FLK(JPOINT) = PI_W_FLK(JPOINT)*PI_INTM_0_H_FLK(JPOINT)/     &
     & PH_ML_P_FLK(JPOINT)
    ELSE
      PI_INTM_0_H_FLK(JPOINT) = PI_H_FLK(JPOINT)
    END IF
    
    IF(PH_ML_P_FLK(JPOINT).LE.(ZDEPTH_W-RH_ML_MIN_FLK)) THEN   ! Integral-mean radiation flux over the thermocline
      PI_INTM_H_D_FLK(JPOINT) = 0._JPRD 
      DO J=1, POPTICPAR_WATER%NBAND_OPTIC
        PI_INTM_H_D_FLK(JPOINT) = PI_INTM_H_D_FLK(JPOINT) +                   &
       & POPTICPAR_WATER%RFRAC_OPTIC(J)/POPTICPAR_WATER%REXTINCOEF_OPTIC(J)*  &
       & ( EXP(-POPTICPAR_WATER%REXTINCOEF_OPTIC(J)*PH_ML_P_FLK(JPOINT))      &
       & - EXP(-POPTICPAR_WATER%REXTINCOEF_OPTIC(J)*ZDEPTH_W) )
      END DO 
      PI_INTM_H_D_FLK(JPOINT) = PI_W_FLK(JPOINT)*PI_INTM_H_D_FLK(JPOINT)/     &
     & (ZDEPTH_W-PH_ML_P_FLK(JPOINT))
    ELSE
      PI_INTM_H_D_FLK(JPOINT) = PI_H_FLK(JPOINT)
    END IF
  ENDIF LAKEPOINT
ENDDO

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FLAKERAD_MOD:FLAKERAD',1,ZHOOK_HANDLE)

END SUBROUTINE FLAKERAD
END MODULE FLAKERAD_MOD

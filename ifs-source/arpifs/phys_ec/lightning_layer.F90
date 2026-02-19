! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LIGHTNING_LAYER(&
 ! Input quantities
 &  YDMODEL, KDIM, PAUX, LLKEYS, STATE, PDIAG, FLUX,&
 ! Output quantity
 &  VNOEMI)

!**** *LIGHTNING_LAYER* - Layer routine calling lightning scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variables for general auxiliary quantities
! LLKEYS   : Derived variable with keys
! state    : Derived variable for model state
! tendency : D. V. for model tendencies (entering convection) from processes before 
! PDIAG    : Derived variable for diagnostic quantities
! FLUX     : Derived variable for fluxes

!    ==== Output tendencies from convection ====
! PNOEMI   :  3D NO Emissions in            kg/m*2/s


!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-11-22  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     P. Lopez 08/10/2015 : Added separate routine for lightning NOx computations.
!     P. Lopez 25/03/2021 : Extended usage to the case of explicit convection (i.e. at very high resolutions).

!-----------------------------------------------------------------------

USE TYPE_MODEL,ONLY : MODEL
USE PARKIND1  ,ONLY : JPIM ,   JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOECLDP   ,ONLY : NCLDQL, NCLDQI, NCLDQR, NCLDQS
USE YOMCST    ,ONLY : RD
USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
                    & AUX_DIAG_TYPE, KEYS_LOCAL_TYPE, FLUX_TYPE
USE VARIABLE_MODULE, ONLY: VARIABLE_3D

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE (MODEL)                   , INTENT (INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)    :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)    :: PAUX
TYPE (KEYS_LOCAL_TYPE)         , INTENT (IN)    :: LLKEYS
TYPE (STATE_TYPE)              , INTENT (IN)    :: STATE
TYPE (AUX_DIAG_TYPE)           , INTENT (INOUT) :: PDIAG
TYPE (FLUX_TYPE)               , INTENT (IN)    :: FLUX
TYPE (VARIABLE_3D)             , INTENT (INOUT) :: VNOEMI

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZRHO
REAL(KIND=JPRB) :: ZLIGH_EMI(KDIM%KLON), ZWMFU(KDIM%KLON)
REAL(KIND=JPRB) :: ZCTOPH(KDIM%KLON), ZPRECMX(KDIM%KLON), ZICETOT(KDIM%KLON), ZCDEPTH(KDIM%KLON)
REAL(KIND=JPRB) :: ZLU(KDIM%KLON,KDIM%KLEV), ZFPLCL(KDIM%KLON,0:KDIM%KLEV), ZFPLCN(KDIM%KLON,0:KDIM%KLEV)
REAL(KIND=JPRB) :: ZQPFROZ(KDIM%KLON,KDIM%KLEV)

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

LOGICAL :: LLLINOX(KDIM%KLON)

!-----------------------------------------------------------------------

#include "culight.intfb.h"
#include "culinox.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LIGHTNING_LAYER',0,ZHOOK_HANDLE)

ASSOCIATE(LCHEM_LIGHT=>YDMODEL%YRML_CHEM%YRCHEM%LCHEM_LIGHT, &
        & LMFPEN=>YDMODEL%YRML_PHY_EC%YRECUMF%LMFPEN)

!     ------------------------------------------------------------------

! Decide which precipitation fluxes and cloud condensate amounts 
! to use, depending on whether convection is treated explicitly or implicitly.
IF (LMFPEN) THEN
  DO JK = 0, KDIM%KLEV
    DO JL = KDIM%KIDIA, KDIM%KFDIA
      ZFPLCL(JL,JK) = FLUX%PFPLCL(JL,JK)
      ZFPLCN(JL,JK) = FLUX%PFPLCN(JL,JK)
    ENDDO
  ENDDO
  DO JK = 1, KDIM%KLEV
    DO JL = KDIM%KIDIA, KDIM%KFDIA
      ZLU(JL,JK) = PDIAG%ZLU(JL,JK)
      ZQPFROZ(JL,JK) = 0.0_JPRB     ! Not used when LMFPEN=T.
    ENDDO
  ENDDO
ELSE
  ! Case when deep convection is treated explicitly.
  ! Use resolved precipitation fluxes and hydrometeor contents as input to lightning parameterization.
  DO JL = KDIM%KIDIA, KDIM%KFDIA
    ZFPLCL(JL,0) = 0.0_JPRB
    ZFPLCN(JL,0) = 0.0_JPRB
  ENDDO
  DO JK = 1, KDIM%KLEV
    DO JL = KDIM%KIDIA, KDIM%KFDIA
      ZRHO = PAUX%PRSF1(JL,JK) / (RD * STATE%T(JL,JK))
      ! Convert rain and snow contents into fluxes (assuming fixed fall speeds).
      ZFPLCL(JL,JK) = STATE%CLD(JL,JK,NCLDQR) * ZRHO * 3.0_JPRB
      ZFPLCN(JL,JK) = STATE%CLD(JL,JK,NCLDQS) * ZRHO * 0.5_JPRB
    ENDDO
  ENDDO
  DO JK = 1, KDIM%KLEV
    DO JL = KDIM%KIDIA, KDIM%KFDIA
      ZLU(JL,JK) = STATE%CLD(JL,JK,NCLDQL) + STATE%CLD(JL,JK,NCLDQI)
      ZQPFROZ(JL,JK) = STATE%CLD(JL,JK,NCLDQS)
    ENDDO
  ENDDO
ENDIF

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CULIGHT

CALL CULIGHT ( YDMODEL%YRML_PHY_EC%YREPHY, YDMODEL%YRML_GCONF%YGFL, YDMODEL%YRML_PHY_EC%YRECUMF, &
  &    KDIM%KIDIA,      KDIM%KFDIA,     KDIM%KLON,    KDIM%KLEV, &
  &    PAUX%PGAW,       PAUX%PGELAT,&
  ! Inputs
  &    PAUX%PRSF1,      PAUX%PRS1,      PAUX%PAPHI,   PAUX%PAPHIF,  LLKEYS%LLLAND,&
  &    STATE%T,         ZLU,            PDIAG%PMFU,   PDIAG%PCAPE,&
  &    ZFPLCL,          ZFPLCN,         ZQPFROZ, &
  &    LLKEYS%LLCUM_LIG, PDIAG%ICBOT_LIG, PDIAG%ICTOP_LIG,&
  ! Outputs
  &    LLLINOX,   PDIAG%PLIGH_TOT, PDIAG%PLIGH_CTG, ZCTOPH,&
  &    ZPRECMX,   ZICETOT,   ZCDEPTH,  ZWMFU,  PDIAG%PCHARGE_LIG)


!*         2.     CALL CULINOX (LIGHTNING NOx COMPUTATIONS)

IF (LCHEM_LIGHT) THEN 

  CALL CULINOX ( YDMODEL%YRML_CHEM%YRCHEM, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, &
    &    PAUX%PGELAT,     PAUX%PRS1,    PAUX%PAPHI,   PAUX%PAPHIF, &
    &    LLKEYS%LLLAND,   LLLINOX, &
    &    STATE%T,         PDIAG%ICTOP,&
    &    PDIAG%PLIGH_TOT, PDIAG%PLIGH_CTG, &
    ! Outputs
    &    VNOEMI%P, ZLIGH_EMI)

ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LIGHTNING_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE LIGHTNING_LAYER

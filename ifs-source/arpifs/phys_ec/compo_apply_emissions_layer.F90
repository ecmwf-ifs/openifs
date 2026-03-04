! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE COMPO_APPLY_EMISSIONS_LAYER(YDSURF, &
  & YDMODEL,KDIM,PAUX,PGFL,PSURF,GEMSL)

!-----------------------------------------------------------------------------
! COMPO_APPLY_EMISSIONS_LAYER
!
!   This is the "layer" routing calling the code for applying prescribed
!   emissions to the atmospheric composition tracers.
!
! AUTHOR: Zak Kipling
!-----------------------------------------------------------------------------

USE SURFACE_FIELDS_MIX ,ONLY : TSURF
USE TYPE_MODEL         ,ONLY : MODEL

USE PARKIND1           ,ONLY : JPIM, JPRB
USE YOMHOOK            ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMPHYDER          ,ONLY : DIMENSION_TYPE, AUX_TYPE, &
 &                             SURF_AND_MORE_TYPE, GEMS_LOCAL_TYPE
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE (TSURF)                   , INTENT(IN)    :: YDSURF
TYPE (MODEL)                   , INTENT(IN)    :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT(IN)    :: KDIM
TYPE (AUX_TYPE)                , INTENT(IN)    :: PAUX
REAL (KIND=JPRB)               , INTENT(IN)    :: PGFL(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (GEMS_LOCAL_TYPE)         , INTENT(INOUT) :: GEMSL
!-----------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "compo_apply_emissions.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('COMPO_APPLY_EMISSIONS_LAYER',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
  CALL COMPO_APPLY_EMISSIONS(YDMODEL, &
   &                         KDIM%KIDIA, KDIM%KFDIA, KDIM%KLEV, KDIM%KLON, GEMSL%ITRAC, &
   &                         PAUX%PGELAT, PAUX%PGELAM, PAUX%PDELP, &
   &                         GEMSL%ZCEN, GEMSL%ZTENC, GEMSL%ZCFLX, &
   & PEMIS2D=PSURF%PSD_VF(:,YDSURF%YSD_VF%YEMIS2D(1)%MP:YDSURF%YSD_VF%YEMIS2D(YDMODEL%YRML_GCONF%YGFL%NEMIS2D)%MP), &
   & PEMIS2DAUX=PSURF%PSD_VF(:,YDSURF%YSD_VF%YEMIS2DAUX(1)%MP:YDSURF%YSD_VF%YEMIS2DAUX(YDMODEL%YRML_GCONF%YGFL%NEMIS2DAUX)%MP), &
   & PEMIS3D=PGFL(:,:,YDMODEL%YRML_GCONF%YGFL%YEMIS3D(1)%MP:YDMODEL%YRML_GCONF%YGFL%YEMIS3D(YDMODEL%YRML_GCONF%YGFL%NEMIS3D)%MP), &
   &                         KAERO=GEMSL%IAERO, KCHEM=GEMSL%ICHEM, &
   &                         PEXTRA=PSURF%PSD_XA, &
   &                         PAPHIF=PAUX%PAPHIF, &
   &                         PLSM=PSURF%PSD_VF(:,YDSURF%YSD_VF%YLSM%MP), &
   &                         PBLH=PSURF%PSD_VD(:,YDSURF%YSD_VD%YBLH%MP))

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COMPO_APPLY_EMISSIONS_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE COMPO_APPLY_EMISSIONS_LAYER

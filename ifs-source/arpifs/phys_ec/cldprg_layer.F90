! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CLDPRG_LAYER(YDSURF, &
  & YDMODEL,KDIM,PAUX,STATE,PSURF,PRAD)

!**** *CLDPRG_LAYER* - Layer routine calling cloud computation for radiation when prognostics cloud scheme is on

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variables for general auxiliary quantities
! state    : Derived variable for model state

!     ==== Input/output ====
! PSURF    : Derived variables for general surface quantities
! PRAD     : Derived variables for variables used in radiation scheme


!        IMPLICIT ARGUMENTS :   NONE
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
!      Original : 2012-12-04  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     F. Vana   11-Sep-2020 : Pressure replaced by provisional values at t+dt

!-----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_RAD_TYPE, AUX_TYPE, &
   &                  SURF_AND_MORE_TYPE
USE YOECLDP  , ONLY : NCLDQI,NCLDQL

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)               , INTENT(INOUT) :: YDSURF
TYPE(MODEL)               , INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)     , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)           , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)         , INTENT (IN)   :: STATE
TYPE (SURF_AND_MORE_TYPE) , INTENT(INOUT) :: PSURF
TYPE (AUX_RAD_TYPE)       , INTENT(INOUT) :: PRAD
!-----------------------------------------------------------------------
LOGICAL ::  LLCLDCOVER

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cldpp.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLDPRG_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VD=>YDSURF%YSD_VD)

!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CLDPP

LLCLDCOVER=.TRUE.

CALL  CLDPP(YDMODEL%YRML_PHY_RAD%YRADIATION,YDMODEL%YRML_PHY_RAD%YRERAD, &
 &          YDMODEL%YRML_PHY_RAD%YRERDI,YDMODEL%YRML_PHY_EC%YRECLD,YDMODEL%YRML_GCONF%YRRIP, &
 &          KDIM%KIDIA,KDIM%KFDIA,KDIM%KLON,KDIM%KLEV,PAUX%PRS1,PAUX%PRSF1, &
 &          STATE%A,STATE%CLD(:,:,NCLDQL),STATE%CLD(:,:,NCLDQI),STATE%T,PAUX%PGELAT,PAUX%PGELAM,LLCLDCOVER, &
 &          PSURF%PSD_VD(:,YSD_VD%YCCC%MP),PSURF%PSD_VD(:,YSD_VD%YHCC%MP),PSURF%PSD_VD(:,YSD_VD%YLCC%MP), &
 &          PSURF%PSD_VD(:,YSD_VD%YMCC%MP),PSURF%PSD_VD(:,YSD_VD%YTCC%MP),PRAD%PNEB,PRAD%PQICE,PRAD%PQLI)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLDPRG_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CLDPRG_LAYER

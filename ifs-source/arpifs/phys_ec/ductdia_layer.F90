! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DUCTDIA_LAYER(YDSURF, &
 ! Input quantities
  & KDIM, PAUX, STATE,  &
 ! Input/Output quantities
  & PSURF)

!**** *DUCTDIA_LAYER* - Layer routine calling ducting diagnostics

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM         : Derived variable for dimensions
! PAUX         : Derived variables for general auxiliary quantities
! state        : Derived variable for model state

!     ==== Input/output ====
! PSURF        : D.V. for surface and other mostly 2D fields...


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
!      Original : 2012-11-29  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   & SURF_AND_MORE_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
!-----------------------------------------------------------------------
REAL(KIND=JPRB) :: ZREFRAC(KDIM%KLON,KDIM%KLEV)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "ductdia.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DUCTDIA_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VD=>YDSURF%YSD_VD)

!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL DUCTDIA

CALL DUCTDIA &
 &( KDIM%KIDIA   , KDIM%KFDIA , KDIM%KLON   , KDIM%KLEV   , &
 &  STATE%T      , STATE%Q    , PAUX%PAPRSF , PAUX%PGEOM1 , &
 &  ZREFRAC , PSURF%PSD_VD(:,YSD_VD%YDNDZN%MP), PSURF%PSD_VD(:,YSD_VD%YDNDZA%MP), PSURF%PSD_VD(:,YSD_VD%YDCTB%MP), &
 &  PSURF%PSD_VD(:,YSD_VD%YTPLB%MP), PSURF%PSD_VD(:,YSD_VD%YTPLT%MP))

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DUCTDIA_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE DUCTDIA_LAYER

! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CLDDIA_LAYER(YDSURF, &
  & YDEPHLI,YDECLD,KDIM,PAUX,STATE,AUXL,PDIAG, PSURF,PRAD)

!**** *CLDDIA_LAYER* - Layer routine calling cloud computation for radiation when prognostics cloud scheme is on

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
! AUXL     : D. V for local auxiliary variables
! PDIAG    : Derived variable for diagnostics quantities

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

!-----------------------------------------------------------------------

USE YOECLD             , ONLY : TECLD
USE YOEPHLI            , ONLY : TEPHLI
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_RAD_TYPE, AUX_TYPE, &
   & SURF_AND_MORE_TYPE, AUX_DIAG_LOCAL_TYPE, AUX_DIAG_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TECLD) ,INTENT(INOUT) :: YDECLD
TYPE(TEPHLI),INTENT(INOUT) :: YDEPHLI
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT (IN)   :: AUXL
TYPE (AUX_DIAG_TYPE)           , INTENT (IN)   :: PDIAG
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (AUX_RAD_TYPE)            , INTENT(INOUT) :: PRAD
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "ccloud.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLDDIA_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VD=>YDSURF%YSD_VD, YSD_VN=>YDSURF%YSD_VN)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CLOUD

CALL  CCLOUD &
  & ( YDEPHLI, YDECLD,  KDIM%KIDIA , KDIM%KFDIA , KDIM%KLON   , KDIM%KLEV,&
  & AUXL%IBASC , AUXL%ITOPC,&
  & PAUX%PRSF1 , PSURF%PSD_VN(:,YSD_VN%YACPR%MP),&
  & STATE%Q    , PDIAG%PQSAT , STATE%T     , PAUX%PVERVEL,&
  & PSURF%PSD_VD(:,YSD_VD%YCCC%MP)  , PSURF%PSD_VD(:,YSD_VD%YHCC%MP), &
  & PSURF%PSD_VD(:,YSD_VD%YLCC%MP)  , PSURF%PSD_VD(:,YSD_VD%YMCC%MP),&
  & PRAD%PNEB  , PSURF%PSD_VD(:,YSD_VD%YTCC%MP)  , PRAD%PQICE  , PRAD%PQLI   ) 

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLDDIA_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CLDDIA_LAYER

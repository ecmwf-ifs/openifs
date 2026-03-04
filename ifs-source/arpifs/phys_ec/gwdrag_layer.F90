! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GWDRAG_LAYER(YDSURF, &
 ! Input quantities
  & YDEPHLI,YDEGWD,KDIM, PAUX, STATE, &
 ! Input/Output quantities
  & PSURF, AUXL)

!**** *GWDRAG_LAYER* - Layer routine calling GWD 

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
! PSURF        : Derived variables for general surface quantities
! AUXL         : Derived variables for local auxiliary quantities

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
!      Original : 2012-11-26  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE YOEGWD             , ONLY : TEGWD
USE YOEPHLI            , ONLY : TEPHLI
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   & AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TEGWD) ,INTENT(INOUT) :: YDEGWD
TYPE(TEPHLI),INTENT(INOUT) :: YDEPHLI
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "gwdrag.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GWDRAG_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VF=>YDSURF%YSD_VF)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL GWDRAG

CALL GWDRAG &
 & ( YDEPHLI,YDEGWD, KDIM%KIDIA   , KDIM%KFDIA  , KDIM%KLON , KDIM%KLEV,&
 & PAUX%PAPRS   , PAUX%PAPRSF , PAUX%PGEOM1,&
 & STATE%T      , STATE%U     , STATE%V,&
 & PSURF%PHSTD  , PSURF%PSD_VF(:,YSD_VF%YVRLAN%MP), PSURF%PSD_VF(:,YSD_VF%YVRLDI%MP), &
 & PSURF%PSD_VF(:,YSD_VF%YSIG%MP),&
 ! TENDENCY COEFFICIENTS (OUTPUT)
 & AUXL%ZSOTEU, AUXL%ZSOTEV, AUXL%ZSOBETA)  

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GWDRAG_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE GWDRAG_LAYER

! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GWDRAG_LAYER(YDMODEL,YDSURF, &
 ! Input quantities
  & KDIM, PAUX, STATE, TENDENCY, &
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
! tendency : Derived variable for updated model tendencies

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

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   & AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)                    , INTENT(INOUT) :: YDSURF
TYPE(MODEL)                    , INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE), TARGET        , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZTSTATE(KDIM%KLON,KDIM%KLEV), &
 & ZUSTATE(KDIM%KLON,KDIM%KLEV), ZVSTATE(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB), POINTER :: ZPRES(:,:), ZPRESF(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "gwdrag.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GWDRAG_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VF=>YDSURF%YSD_VF,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,&
  & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,&
  & YDEPHLI=>YDMODEL%YRML_PHY_SLIN%YREPHLI,YDEGWD=>YDMODEL%YRML_PHY_EC%YREGWD)
ASSOCIATE(KIDIA=>KDIM%KIDIA,KFDIA=>KDIM%KFDIA,KLEV=>KDIM%KLEV, KLON=>KDIM%KLON)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL GWDRAG

! Update state and pressure to remain consistent with VDIFF
! Tendencies must not be used
ZTSTATE(KIDIA:KFDIA,1:KLEV)=STATE%T(KIDIA:KFDIA,1:KLEV)
ZUSTATE(KIDIA:KFDIA,1:KLEV)=STATE%U(KIDIA:KFDIA,1:KLEV)
ZVSTATE(KIDIA:KFDIA,1:KLEV)=STATE%V(KIDIA:KFDIA,1:KLEV)
ZPRES =>PAUX%PAPRS
ZPRESF=>PAUX%PAPRSF

CALL GWDRAG &
 & ( YDEPHLI,YDEGWD, KIDIA, KFDIA, KLON, KLEV,&
 & ZPRES , ZPRESF , PAUX%PGEOM1,&
 & ZTSTATE, ZUSTATE, ZVSTATE,&
 & PSURF%PHSTD  , PSURF%PSD_VF(:,YSD_VF%YVRLAN%MP), PSURF%PSD_VF(:,YSD_VF%YVRLDI%MP), &
 & PSURF%PSD_VF(:,YSD_VF%YSIG%MP),&
 ! TENDENCY COEFFICIENTS (OUTPUT)
 & AUXL%ZSOTEU, AUXL%ZSOTEV, AUXL%ZSOBETA)  

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GWDRAG_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE GWDRAG_LAYER

! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADIATION_LAYER(YDDIMV,YDSURF,YDMODEL,&
  & KDIM,STATE,PDIAG,PRAD,PAUX,AUXL,PSURF,SURFL)

!**** *RADIATION_LAYER* - Layer routine calling RADIATION scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! state    : Derived variable for model state
! PDIAG    : Derived variable for diagnostic quantities

!     ==== Input/output ====
! PRAD     : Derived variables for variables used in radiation scheme
! PAUX     : Derived variables for general auxiliary quantities
! AUXL     : Derived variables for local auxiliary quantities
! PSURF    : Derived variables for general surface quantities
! SURFL    : Derived variables for local surface quantities


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
!      Original : 2012-11-22  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!      June 2014  R. Hogan   Pass surface LW down flux (PSRLWD) to RADINA

!-----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE YOMDIMV            , ONLY : TDIMV
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER          , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_RAD_TYPE, AUX_TYPE,&
   &                            AUX_DIAG_TYPE, AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, SURF_AND_MORE_LOCAL_TYPE
USE YOECLDP            , ONLY : NCLDQR,NCLDQS

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV)                    , INTENT(IN)    :: YDDIMV
TYPE(TSURF)                    , INTENT(INOUT) :: YDSURF
TYPE(MODEL)                    , INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (AUX_DIAG_TYPE)           , INTENT (IN)   :: PDIAG
TYPE (AUX_RAD_TYPE)            , INTENT(INOUT) :: PRAD
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (SURF_AND_MORE_LOCAL_TYPE), INTENT(INOUT) :: SURFL
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "radina.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RADIATION_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VF=>YDSURF%YSD_VF, YSP_RR=>YDSURF%YSP_RR)


!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL RADIATION


CALL RADINA&
  & (YDDIMV, YDMODEL, KDIM%KIDIA  , KDIM%KFDIA , KDIM%KLON  , KDIM%KLEV   ,&
  & SURFL%ZALBD  , SURFL%ZALBP  , PAUX%PAPRS  , PAUX%PAPRSF ,&
  & AUXL%ZCCNL  , AUXL%ZCCNO  , PRAD%PNEB   ,&
  & PAUX%PGELAM , PAUX%PCLON  , PAUX%PSLON  ,&
  & PAUX%PDELP  , SURFL%ZEMIR  , SURFL%ZEMIW  ,&
  & PAUX%PGEMU  , PAUX%PMU0   , STATE%Q   , PDIAG%PQSAT  ,&
  & PRAD%PQICE  , PRAD%PQLI   , STATE%CLD(:,:,NCLDQR), STATE%CLD(:,:,NCLDQS),&
  & PSURF%PSD_VF(:,YSD_VF%YLSM%MP), STATE%T, PSURF%PSP_RR(:,YSP_RR%YT%MP9),&
  & PRAD%PEMTD  , PRAD%PTRSW,&
  & AUXL%ZEMIT  , AUXL%ZTHT   , AUXL%ZCTRSO , AUXL%ZCEMTR , AUXL%ZTRSOD, PRAD%PSRLWD(:,1),&
  & AUXL%ZSUDU  , PRAD%PISUND , PRAD%PDSRP&
  & )  

! ecRad radiation scheme has capability to produce longwave
! downwelling fluxes in the same spectral intervals as emissivity
! values; RADINA does not so we pad the rest with zeros.
IF (UBOUND(PRAD%PSRLWD,2) > 1) THEN
  PRAD%PSRLWD(KDIM%KIDIA:KDIM%KFDIA,2:) = 0.0_JPRB
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RADIATION_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE RADIATION_LAYER

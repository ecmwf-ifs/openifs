! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AERDIAG_LAYER(YDSURF,YDECLDP,YDERAD,YDML_GCONF,YDPHY2,KDIM, GEMSL, PAUX, STATE, PDIAG, AUXL, &
 ! Input/Output quantities
 & PSURF)

!**** *AERDIAG_LAYER* - Routine called when radiation scheme is by-passed

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! GEMSL        : Derived variables for local GEMS quantities
! PAUX     : Derived variables for general auxiliary quantities
! state    : Derived variable for model state
! PDIAG    : Derived variable for diagnostics quantities
! AUXL     : Derived variables for local auxiliary quantities

!     ==== Input/Output ====
! PSURF    : Derived variables for general surface quantities

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

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOECLDP                , ONLY : TECLDP, NCLDQI, NCLDQL
USE YOERAD                 , ONLY : TERAD
USE YOMPHY2                , ONLY : TPHY2
USE SURFACE_FIELDS_MIX     , ONLY : TSURF
USE PARKIND1               , ONLY : JPRB
USE YOMHOOK                , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMPHYDER              , ONLY : DIMENSION_TYPE, AUX_DIAG_LOCAL_TYPE, &
   &                                SURF_AND_MORE_TYPE, AUX_DIAG_TYPE, STATE_TYPE, &
   &                                AUX_TYPE, GEMS_LOCAL_TYPE
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)                    , INTENT(INOUT) :: YDSURF
TYPE(TECLDP)                   , INTENT(INOUT) :: YDECLDP
TYPE(TERAD)                    , INTENT(INOUT) :: YDERAD
TYPE(MODEL_GENERAL_CONF_TYPE)  , INTENT(INOUT) :: YDML_GCONF
TYPE(TPHY2)                    , INTENT(INOUT) :: YDPHY2
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (GEMS_LOCAL_TYPE)         , INTENT (IN)   :: GEMSL
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (AUX_DIAG_TYPE)           , INTENT (IN)   :: PDIAG
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "aer_diag1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('AERDIAG_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YSD_VF=>YDSURF%YSD_VF, NACTAERO=>YDML_GCONF%YGFL%NACTAERO)

!     ------------------------------------------------------------------

! Extra diag removed from callpar
AUXL%IEXT3D = 0

!*         1.     Unroll derived structures and call aer_diag1

CALL AER_DIAG1 &
  &( YDERAD,YDECLDP,YDML_GCONF,YDPHY2, &
  & KDIM%KIDIA , KDIM%KFDIA, KDIM%KLON  , KDIM%KTDIA, KDIM%KLEV, KDIM%KLEVX, KDIM%KFLDX, AUXL%IEXT3D, &
  &  GEMSL%ZAEROP(:,:,1:NACTAERO), PAUX%PAPRS, PAUX%PAPRSF, STATE%A  , STATE%CLD(:,:,NCLDQI)   , STATE%CLD(:,:,NCLDQL)  , &
  &  PSURF%PSD_VF(:,YSD_VF%YLSM%MP)  , STATE%Q   , PDIAG%PQSAT , STATE%T  , AUXL%ZWND , &
  &  PSURF%PSD_XA)


WHERE (STATE%A(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV) >  0.001_JPRB)
  PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,AUXL%IEXT3D+11)=&
    &  STATE%CLD(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,NCLDQI)/STATE%A(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
  PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,AUXL%IEXT3D+12)=&
    &  STATE%CLD(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,NCLDQL)/STATE%A(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
  PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,AUXL%IEXT3D+13)= STATE%A(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV)
ELSEWHERE
  PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,AUXL%IEXT3D+11)= 0._JPRB
  PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,AUXL%IEXT3D+12)= 0._JPRB
  PSURF%PSD_XA(KDIM%KIDIA:KDIM%KFDIA,1:KDIM%KLEV,AUXL%IEXT3D+13)= 0._JPRB
ENDWHERE

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AERDIAG_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE AERDIAG_LAYER

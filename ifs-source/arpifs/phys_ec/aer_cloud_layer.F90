! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_CLOUD_LAYER(&
 ! Input quantities
 &  YDMODEL,KDIM, PAUX, STATE, PDIAG, GEMSL,&
 ! Input/Output quantities
 &  AUXL)

!**** *AER_CLOUD_LAYER* - Layer routine calling cloud aerosol interactions

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
! PDIAG        : Derived variable for diagnostic quantities
! GEMSL        : Derived variables for local GEMS quantities

!     ==== Input/output ====
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
!      Original : 2012-12-05  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER, ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
  &  AUX_DIAG_TYPE, AUX_DIAG_LOCAL_TYPE, GEMS_LOCAL_TYPE
USE YOECLDP  , ONLY : NCLDQI, NCLDQL

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)                     ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (AUX_DIAG_TYPE)           , INTENT (IN)   :: PDIAG
TYPE (GEMS_LOCAL_TYPE)         , INTENT (IN)   :: GEMSL
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL

!-----------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "aer_clcld.intfb.h"
#include "aer_cld.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('AER_CLOUD_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & NACTAERO=>YGFL%NACTAERO)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL APPROPRIATE SCHEME

IF (NACTAERO == 0) THEN
 
  ! the state used for this has to be consistent with qsat
  CALL AER_CLCLD( &
    !      input
    & YDMODEL, KDIM%KIDIA,    KDIM%KFDIA,    KDIM%KLON,    KDIM%KLEV,&
    & TSPHY,&
    & state%T,  state%q,   PDIAG%PQSAT, &
    & PAUX%PAPRS,    PAUX%PAPRSF,&
    & PAUX%PGELAM,   PAUX%PGEMU,    PAUX%PCLON,   PAUX%PSLON,& 
    & state%cld(:,:,NCLDQL), state%cld(:,:,NCLDQI),  state%a, &
    !      output
    & AUXL%ZLCRIT_AER,AUXL%ZICRIT_AER,&
    & AUXL%ZRE_LIQ,  AUXL%ZRE_ICE,&
    & AUXL%ZCCN,     AUXL%ZNICE )
    !      diagnostics
    !       & PSURF%PSD_XA,   KDIM%KFLDX)

ELSE  ! NACTAERO /= 0

  CALL AER_CLD( &
    !      input
    & YDMODEL%YRML_PHY_EC%YRECLDP,   KDIM%KIDIA,    KDIM%KFDIA,    KDIM%KLON,    KDIM%KLEV,&
    & NACTAERO, &
    & GEMSL%ZAEROP(:,:,1:NACTAERO),   state%T,&
    & PAUX%PAPRSF,&
    & state%cld(:,:,NCLDQL), state%cld(:,:,NCLDQI),  state%a, &
    !      output
    & AUXL%ZLCRIT_AER,AUXL%ZICRIT_AER,&
    & AUXL%ZRE_LIQ,  AUXL%ZRE_ICE,&
    & AUXL%ZCCN,     AUXL%ZNICE &
    & )
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_CLOUD_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE AER_CLOUD_LAYER

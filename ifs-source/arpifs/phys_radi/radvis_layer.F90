! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE RADVIS_LAYER(YDECLDP,YDERAD,YGFL,YDPHY2,KDIM,PAUX,STATE,GEMSL,PDIAG)

!**** *RADVIS_LAYER* - Layer routine calling visibility diagnostics

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variable for general auxiliary quantities
! state    : Derived variable for model state

!     ==== Input/output ====
! GEMSL    : Derived variable for local GEMS quantities
! PDIAG    : Derived variable for diagn. quantities


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
!      Original : 11-Feb-2012  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE YOECLDP  , ONLY : TECLDP, NCLDQR, NCLDQS, NCLDQI, NCLDQL
USE YOERAD   , ONLY : TERAD
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMPHYDER, ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
  &                   GEMS_LOCAL_TYPE, AUX_DIAG_TYPE
USE YOMPHY2  , ONLY : TPHY2
USE YOM_YGFL , ONLY : TYPE_GFLD

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TECLDP)                   , INTENT(INOUT) :: YDECLDP
TYPE(TERAD)                    , INTENT(INOUT) :: YDERAD
TYPE(TPHY2)                    , INTENT(INOUT) :: YDPHY2
TYPE(TYPE_GFLD)                , INTENT(INOUT) :: YGFL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (GEMS_LOCAL_TYPE)         , INTENT(INOUT) :: GEMSL
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB) :: ZVISIPR(KDIM%KLON), ZVISPAE(KDIM%KLON)
REAL(KIND=JPRB) :: ZVISCAE(KDIM%KLON), ZVISCLD(KDIM%KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "radvis.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RADVIS_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & NACTAERO=>YGFL%NACTAERO)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND EVENTUALLY CALL RADVIS

IF (NACTAERO == 0) THEN
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PDIAG%PVISIH(JL)=GEMSL%ZVISICL(JL)
  ENDDO
ELSE
  GEMSL%ZVISICL=0.0_JPRB
  ZVISIPR=0.0_JPRB
  ZVISCAE=0.0_JPRB
  ZVISPAE=0.0_JPRB
  ZVISCLD=0.0_JPRB

  CALL RADVIS ( YDERAD, KDIM%KIDIA  , KDIM%KFDIA  , KDIM%KLON   , KDIM%KLEV , &
    &  PAUX%PRSF1  , STATE%A    , STATE%CLD(:,:,NCLDQL), STATE%CLD(:,:,NCLDQI), &
    &  STATE%CLD(:,:,NCLDQR), STATE%CLD(:,:,NCLDQS)   , STATE%T    , TSPHY, &
    &  GEMSL%ZCLAERS, GEMSL%ZPRAERS, &
    &  GEMSL%ZVISICL, ZVISIPR, ZVISCAE, ZVISPAE, ZVISCLD )

  DO JL=KDIM%KIDIA,KDIM%KFDIA
    PDIAG%PVISIH(JL)=ZVISIPR(JL)
  ENDDO
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RADVIS_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE RADVIS_LAYER

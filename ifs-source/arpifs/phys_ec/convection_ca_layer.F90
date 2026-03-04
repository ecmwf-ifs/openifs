! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CONVECTION_CA_LAYER(&
 ! Input quantities
 &  YDECUCONVCA,YDECUMF,YDPHY2,KDIM, PAUX, STATE,&
 ! Input/Output quantities
 &  TENDENCY, PPERT, PERTL)

!**** *CONVECTION_CA_LAYER* - Layer routine calling CA scheme in convection

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
! tendency : D. V. for model tendencies (entering convection) 
! PPERT    : Derived variables for local perturbartions
! PERTL    : Derived variables for local perturbartions


!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!       Note: Output directly applied to the tendencies which are then IN/OUT here

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-11-23  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM ,   JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
  &  PERTURB_TYPE, PERTURB_LOCAL_TYPE
USE YOMPHY2  , ONLY : TPHY2
USE YOECUMF  , ONLY : TECUMF
USE YOE_CUCONVCA, ONLY : TECUCONVCA


!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TECUCONVCA)                ,INTENT(INOUT) :: YDECUCONVCA
TYPE(TECUMF)                    ,INTENT(INOUT) :: YDECUMF
TYPE(TPHY2)                     ,INTENT(INOUT) :: YDPHY2
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY
TYPE (PERTURB_TYPE)            , INTENT(INOUT) :: PPERT
TYPE (PERTURB_LOCAL_TYPE)      , INTENT(INOUT) :: PERTL

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB) :: ZFAC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "ca_profpert.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CONVECTION_CA_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(NLIVES=>YDECUCONVCA%NLIVES, &
 & LMFPROFP=>YDECUMF%LMFPROFP, &
 & TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CUCALLN


IF (LMFPROFP) THEN
  CALL CA_PROFPERT(YDECUCONVCA,KDIM%KLON,KDIM%KLEV,KDIM%KIDIA,KDIM%KFDIA,TSPHY,STATE%Q,&
    &     TENDENCY%T,TENDENCY%Q,PAUX%PAPRSF,PPERT%PCUCONVCA,PPERT%PNLCONVCA)
ENDIF

! scaling of CA values before used in convection to modulate base mass flux
ZFAC=3.0_JPRB/REAL(NLIVES)
DO JL=KDIM%KIDIA,KDIM%KFDIA
  PERTL%ZCUCONVCA(JL)=PPERT%PCUCONVCA(JL)*ZFAC
ENDDO



!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CONVECTION_CA_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECTION_CA_LAYER

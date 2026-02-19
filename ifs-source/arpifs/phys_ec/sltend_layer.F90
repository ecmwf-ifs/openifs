! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SLTEND_LAYER(&
 ! Input quantities
 &  YDMODEL,KDIM, PAUX, STATE_T0, &
 &  TENDENCY_DYN, TENDENCY_VDF, TENDENCY_SATADJ, TENDENCY_CML, &
 &  PSLPHY9,&
 ! Input/Output quantities
 &  PSAVTEND, PGFLSLP, PSURF,&
 ! Output tendencies
 &  TENDENCY_LOC)

!**** *SLTEND_LAYER* - Layer routine calling SLTEND

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM         : Derived variable for dimensions
! PAUX         : Derived variable for general auxiliary quantities
! STATE_T0     : Derived variable for model state at t0
! TENDENCY_DYN : D. V. for model tendencies from explicit dynamics
! TENDENCY_VDF : D. V. for model tendencies from turbulence scheme
! TENDENCY_SATADJ : D. V. for model tendencies from saturation adjustment
! TENDENCY_CML : D. V. for cumulated tendencies from all processes
! PSLPHY9      : D. V. for 0.5* phys tendencies from (t) at departure point

!     ==== Input/output ====
! PSAVTEND  : ARRAY OF GMV TENDENCIES + AUX. SUPERSATURATION TO BE SAVED FOR NEXT TIME STEP
! PGFLSLP   : ARRAY OF GFL TENDENCIES + AUX. SUPERSATURATION TO BE SAVED FOR NEXT TIME STEP

!    ==== Output tendencies from convection ====
! TENDENCY_LOC :  Output process tendencies


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
!     R. Forbes 2020-11-15 Added input/output for cloud condensate in sltend call
!     R. Forbes 2021-07-01 Added input of saturation adjustment tendencies (fast process) 
!
!-----------------------------------------------------------------------

USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1  , ONLY : JPRB
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMPHYDER , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
  &                    MODEL_STATE_TYPE, SURF_AND_MORE_TYPE

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)                    , INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE_T0
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_DYN
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_VDF
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_SATADJ
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_CML
TYPE (MODEL_STATE_TYPE)        , INTENT (IN)   :: PSLPHY9
REAL(KIND=JPRB)   ,INTENT(OUT) :: PSAVTEND(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_PHY_G%YRSLPHY%NVTEND) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PGFLSLP(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIMSLP)
TYPE (SURF_AND_MORE_TYPE)      , INTENT (INOUT) :: PSURF
TYPE (STATE_TYPE)              , INTENT (INOUT) :: TENDENCY_LOC

!-----------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "sltend.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SLTEND_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY)
ASSOCIATE(NDIMSLP=>YGFL%NDIMSLP, &
 & NVTEND=>YDSLPHY%NVTEND)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL SLTEND

CALL SLTEND(YDMODEL, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, NVTEND, &
  & STATE_T0%T, STATE_T0%Q, STATE_T0%A, STATE_T0%CLD, PAUX%PRSF1, PAUX%PRS1, &
  & TENDENCY_DYN%U, TENDENCY_DYN%V, TENDENCY_DYN%T, TENDENCY_DYN%Q, TENDENCY_DYN%O3, TENDENCY_DYN%A, TENDENCY_DYN%CLD, &
  & TENDENCY_CML%U, TENDENCY_CML%V, TENDENCY_CML%T, TENDENCY_CML%Q, TENDENCY_CML%O3, TENDENCY_CML%A, TENDENCY_CML%CLD, &
  & TENDENCY_LOC%U, TENDENCY_LOC%V, TENDENCY_LOC%T, TENDENCY_LOC%Q, TENDENCY_LOC%O3, TENDENCY_LOC%A, TENDENCY_LOC%CLD, &
  & PSLPHY9%U, PSLPHY9%V, PSLPHY9%T, PSLPHY9%GFL, &
  & TENDENCY_VDF%U, TENDENCY_VDF%V, TENDENCY_VDF%T, TENDENCY_VDF%Q, TENDENCY_VDF%A, &
  & TENDENCY_SATADJ%T, TENDENCY_SATADJ%Q, TENDENCY_SATADJ%A, TENDENCY_SATADJ%CLD, & 
  & PSAVTEND, PGFLSLP, PSURF%PSD_XA, KDIM%KFLDX ) 


!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SLTEND_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE SLTEND_LAYER

! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE ECRADFR(YDDIM,YDML_PHY_RAD,YDDYNA,YDRIP)

!**** *ECRADFR*   - MODIFY FREQUENCY OF FULL-RADIATION COMPUTATIONS

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *ECRADFR* FROM *CNT4*
!              -------        ----

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMONS YOMDIM, YOMCT3, YOERAD

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        ORIGINAL : 92-11-25

!     MODIFICATIONS.
!     --------------
!        010129 JJMorcrette clean-up LERAD1H, NLNGR1H
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        N. Semane+P.Bechtold   replace 3600s by RHOUR for small planet
!        K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_RADIATION_MOD , ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE YOMDIM   , ONLY : TDIM
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCT3   , ONLY : NSTEP
USE YOMRIP   , ONLY : TRIP
USE YOMDYNA  , ONLY : TDYNA
USE YOMCST   , ONLY : RHOUR

!      ----------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM) , INTENT(IN) :: YDDIM
TYPE(MODEL_PHYSICS_RADIATION_TYPE),INTENT(INOUT):: YDML_PHY_RAD
TYPE(TDYNA), INTENT(IN) :: YDDYNA
TYPE(TRIP)  ,INTENT(INOUT):: YDRIP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

#include "updtier.intfb.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ECRADFR',0,ZHOOK_HANDLE)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & LERAD1H=>YDML_PHY_RAD%YRERAD%LERAD1H, NLNGR1H=>YDML_PHY_RAD%YRERAD%NLNGR1H, NRADFR=>YDML_PHY_RAD%YRERAD%NRADFR, &
 & NRADNFR=>YDML_PHY_RAD%YRERAD%NRADNFR, NRADSFR=>YDML_PHY_RAD%YRERAD%NRADSFR, &
 & NRADUV=>YDML_PHY_RAD%YREUVRAD%NRADUV, &
 & TSTEP=>YDRIP%TSTEP)
!      ----------------------------------------------------------------

IF (NSMAX >= 63 .AND. LERAD1H) THEN
  IF (NSTEP*TSTEP < NLNGR1H *RHOUR) THEN
    NRADFR=NRADSFR
    NRADUV=NRADSFR
  ELSE
    NRADFR=NRADNFR
    NRADUV=NRADNFR
  ENDIF
ENDIF

CALL UPDTIER(YDML_PHY_RAD,YDDYNA,YDRIP,NSTEP,NRADFR,TSTEP)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ECRADFR',1,ZHOOK_HANDLE)
END SUBROUTINE ECRADFR

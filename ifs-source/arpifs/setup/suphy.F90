! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE SUPHY(YDGEOMETRY,YDMODEL,KULOUT)

!**** *SUPHY*   - Calls physic specific set-up routines

!     Purpose.
!     --------
!           Calls set-up routines specific to the different physics
!           packages that can be used in the IFS/ARPEGE model

!**   Interface.
!     ----------
!        *CALL* *SUPHY(...)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY, YOEPHY

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        or
!        Documentation ARPEGE (depending on which physics will be used)

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      J.-F. Geleyn for the ARPEGE rewriting.
!      Original : 87-10-15

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      B.Sass        01-June-2006 (call setup for HIRLAM physics)
!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LR2D


IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)    :: YDGEOMETRY
TYPE(MODEL)    ,INTENT(INOUT)    :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "suphec.intfb.h"
#include "suphmf.intfb.h"
#include "sumts.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUPHY',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    Call initialization of specific physics' commons.
!              -------------------------------------------------

!*       1.1   Meteo-France Physics
!              --------------------

CALL SUPHMF(YDGEOMETRY,YDMODEL,KULOUT)

!*       1.2   ECMWF Physics
!              -------------

IF (.NOT.LR2D) THEN
  CALL SUPHEC(YDGEOMETRY,YDMODEL,KULOUT)
ENDIF

!     ------------------------------------------------------------------

!*       3.    Initialize "model to satellite" RTTOV parameters.
!              ------------------------------------------------

CALL SUMTS(YDMODEL%YRML_GCONF)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPHY',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHY

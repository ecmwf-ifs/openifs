! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE STEPO1C(YDGEOMETRY,YDMODEL,YDFIELDS)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE FIELDS_MOD , ONLY : FIELDS
USE YOMCT0   , ONLY : NFRPOS
USE YOMCT3   , ONLY : NSTEP
USE YOMLOG1C , ONLY : OUTFORM

#ifdef DOC

!**** *STEPO1C*  - Controls integration job at lowest level

!     Purpose.
!     --------
!        CONTROLS THE INTEGRATION

!**   Interface.
!     ----------
!        *CALL* *STEPO1C

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!                 WRTP1C -  Write out prognostic variables
!                 CPG1C  -  Grid point computations

!        Called by CNT41C

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the single column model

!     Author.
!     -------
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original      94-01-06
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),    INTENT(INOUT) :: YDMODEL
TYPE(FIELDS)   ,INTENT(INOUT) :: YDFIELDS


!     ------------------------------------------------------------------
#include "wrtp1c_nc.intfb.h"
#include "cpg1c.intfb.h"
#include "wrtp1c.intfb.h"
!     ------------------------------------------------------------------


!*       1.    WRITE OUT PROGNOSTIC VARIABLES.
!              -------------------------------

IF (MOD(NSTEP,NFRPOS) == 0) THEN
  IF (OUTFORM .EQ. 'netcdf') THEN
    CALL WRTP1C_NC(YDGEOMETRY,YDMODEL,YDFIELDS%YRSURF)
  ELSE
    CALL WRTP1C(YDGEOMETRY,YDMODEL,YDFIELDS%YRSURF)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       2.    GRID POINT COMPUTATIONS.
!              ------------------------

CALL CPG1C(YDGEOMETRY,YDMODEL,YDFIELDS)

!     ------------------------------------------------------------------

END SUBROUTINE STEPO1C

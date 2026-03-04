! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CNT1C

#ifdef DOC

!**** *CNT1C*  - Routine which controls the SCM.

!     Purpose.
!     --------
!          Controls the SCM.

!***  Interface.
!     ----------
!        *CALL* *CNT1C

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        SU0YOM1C    - main initialization
!        SU2YOM      - initialize level 2 commons (few)
!        SU3YOM      - initialize level 3 commons (few)
!        SUINIF1C_NC - read initial condition
!        CNT41C      - integration

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the one column model

!     Author.
!     -------
!        Martin Koehler  *ECMWF*

!     Modifications.
!     --------------
!        Origin:        01-02-22 combination of cnt01c, cnt21c, cnt31c 
!        Martin Koehler 00-09-10 NetCDF adoptation
!        M. Ko"hler     6-6-2006 Single Column Model integration within IFS 
!     ------------------------------------------------------------------
#endif

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE MODEL_MOD, ONLY : MODEL, MODEL_SET
USE FIELDS_MOD,ONLY : FIELDS
!USE MTRAJ_MOD, ONLY : MTRAJ, MTRAJ_SET

!     ------------------------------------------------------------------
IMPLICIT NONE

TYPE(GEOMETRY)     :: YRGEOMETRY
TYPE(MODEL),TARGET :: YRMODEL
TYPE(FIELDS)       :: YRFIELDS
!TYPE(MTRAJ)    :: YRMTRAJ

!     ------------------------------------------------------------------
#include "ifs_init.intfb.h"
#include "su0yom1c.intfb.h"
#include "su2yom.intfb.h"
#include "su3yom.intfb.h"
#include "suinif1c_nc.intfb.h"
#include "cnt41c.intfb.h"
#include "cntend.intfb.h"
!     ------------------------------------------------------------------

!  Copy from CNT0:
!*      0.   Initialize level 0 commons.
!            ---------------------------

CALL MODEL_SET(YRMODEL)


!        1.    Initialize commons.
!              -------------------
CALL IFS_INIT()

CALL SU0YOM1C(YRGEOMETRY,YRMODEL,YRFIELDS)

YRFIELDS%STATE_MODEL => YRMODEL

!        2.    INITIALIZE YOMCT2.
!              ------------------

CALL SU2YOM(YRMODEL%YRML_GCONF%YRRIP)

!        3.    INITIALIZE LEVEL 3 COMMONS.
!              ---------------------------

CALL SU3YOM

!        4.    READ INPUT FILE.
!              ----------------

CALL SUINIF1C_NC(YRGEOMETRY,YRMODEL,YRFIELDS%YRSURF)

!        5.    INTEGRATION JOB.
!              ----------------

CALL CNT41C(YRGEOMETRY,YRMODEL,YRFIELDS)

!        6.    Close datafiles.
!              ----------------

CALL CNTEND


END SUBROUTINE CNT1C

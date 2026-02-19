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

SUBROUTINE CSTA(YDGEOMETRY,YDFIELDS,YDMODEL,KINITMONTH)

!**** *CSTA*  - START OF THE MODEL

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *CSTA

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.    SUINIF - initialize model fields
!     ----------    NOR    - change orography

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.Yessad  : June 2006   Possibility to read purely gp ARPEGE files
!      M.Drusch      17-Jan-2007 introduce nconf 302
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      K. Yessad (July 2014): Move some variables.
!      O. Marsden  Aug 2016 : Removed use of SPA3
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE FIELDS_MOD         , ONLY : FIELDS
USE YOMGFL             , ONLY : TGFL
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0             , ONLY : NPRINTLEV
USE YOMCT0             , ONLY : NCONF, LARPEGEF_RDGP_INIT
USE YOMLUN             , ONLY : NULOUT
IMPLICIT NONE

!      -----------------------------------------------------------

TYPE(GEOMETRY),INTENT(INOUT) :: YDGEOMETRY
TYPE(FIELDS)  ,INTENT(INOUT) :: YDFIELDS
TYPE(MODEL)   ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KINITMONTH

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------

#include "suinif.intfb.h"

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CSTA',0,ZHOOK_HANDLE)

!      -----------------------------------------------------------

IF (NCONF/100==0.OR.NCONF/100==2.OR.NCONF==302) THEN
  CALL SUINIF(YDGEOMETRY,YDFIELDS%YRGFL,YDFIELDS%YRSURF,YDFIELDS%YRSPEC,YDFIELDS%YRCFU,YDFIELDS%YRXFU, &
 &            YDMODEL,0,LDRDGRIDSP=LARPEGEF_RDGP_INIT,YDMCUF=YDFIELDS%YMCUF,KINITMONTH=KINITMONTH)
  IF (NPRINTLEV >=2) THEN
    WRITE(NULOUT,*) ' CSTA: SUINIF(0) called '
    CALL FLUSH(NULOUT)
  ENDIF
ENDIF

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CSTA',1,ZHOOK_HANDLE)
END SUBROUTINE CSTA

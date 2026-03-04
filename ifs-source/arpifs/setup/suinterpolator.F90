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

SUBROUTINE SUINTERPOLATOR(YDGEOMETRY,YDDYNA,YDSLINT)

!**** *SUVERT*  - Routine to initialize vertical interpolator

!     Purpose.
!     --------
!           Initialize vertical interpolator

!**   Interface.
!     ----------

!     *CALL* SUINTERPOLATOR
!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        see the modules used above.

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Tomas Wilhelmsson from SUVERT.
!      Original : 2013-08-22

!     Modifications.
!     --------------
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMSLINT , ONLY : TSLINT
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LRPLANE
USE YOMMP0   , ONLY : NPRINTLEV, LOUTPUT
USE YOMDYNA  , ONLY : TDYNA

!-------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDYNA),    INTENT(IN) :: YDDYNA
TYPE(TSLINT),   INTENT(INOUT) :: YDSLINT
LOGICAL :: LLVERBOSE,LLEQMER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "suvsleta.intfb.h"
#include "suvsplip.intfb.h"
#include "suhslmer.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUINTERPOLATOR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NDGENH=>YDDIM%NDGENH, NDGSAH=>YDDIM%NDGSAH, NDGSUR=>YDDIM%NDGSUR, &
 & NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

!*            1. SET UP VERTICAL INTERPOLATOR.
!             --------------------------------

IF (NFLEVG > 1) THEN
  ! * externalisable part of interpolator: computes RVSPTRI, RVSPC, RFVV.
  CALL SUVSPLIP(YDGEOMETRY,YDSLINT%YRVSPLIP)

  ! * externalisable part of interpolator: attributes of YRVSLETA.
  LLVERBOSE=LOUTPUT.AND.NPRINTLEV >= 1
  CALL SUVSLETA(YDGEOMETRY,YDDYNA,YDSLINT%YRVSLETA,LLVERBOSE)
ENDIF

! * externalisable part of interpolator: intermediate quantities for meridian cubic interpolations.
IF (LRPLANE .OR. (NDGSUR>=2 .AND. .NOT.LRPLANE)) THEN
  LLEQMER=LRPLANE
  CALL SUHSLMER(YDGEOMETRY,NDGSAH,NDGENH,LLEQMER,YDDYNA%SLHDEPSH,YDSLINT%YRHSLMER)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUINTERPOLATOR',1,ZHOOK_HANDLE)
END SUBROUTINE SUINTERPOLATOR

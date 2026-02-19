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

SUBROUTINE SUVERT2(YDGEOMETRY)

!**** *SUVERT2*  - Routine to initialize vertical coordinate
!                  Simplified version of SUVERT for 2D models and more generally
!                  for all configurations using NIOLEVG=NFLEVG=1.

!     Purpose.
!     --------
!           Initialize the hybrid-cordinate system of the model in 2D models.

!**   Interface.
!     ----------

!     *CALL* SUVERT2

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
!      K. Yessad (July 2012) after SUVERT

!     Modifications.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
! End Modifications
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVV1   , ONLY : DVALH, DVBH

!-------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVERT2',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
  & NIOLEVG=>YDGEOMETRY%YRDIMV%NIOLEVG,YDDIMV=>YDGEOMETRY%YRDIMV)!     ------------------------------------------------------------------

IF (NIOLEVG /= 1 .OR. NFLEVG /= 1) THEN
  CALL ABOR1('SUVERT2 may be called only if NIOLEVG=1 and NFLEVG=1.')
ENDIF

! * fill DVALH and DVBH.
DVALH(0)=0.0_JPRB
DVALH(1)=0.0_JPRB
DVBH(0)=0.0_JPRB
DVBH(1)=1.0_JPRB

! * fill YRVETA.
YDVETA%VETAF(0)=0._JPRB
YDVETA%VETAF(1)=0.5_JPRB
YDVETA%VETAF(2)=1.0_JPRB
YDVETA%VETAH(0)=0._JPRB
YDVETA%VETAH(1)=1.0_JPRB

! * fill YRVAB (attributes VC, VDELA, VDELB, VAF, VBF are not used in this case).
YDVAB%VAH(0)=0._JPRB
YDVAB%VAH(1)=0._JPRB
YDVAB%VALH(0)=0._JPRB
YDVAB%VALH(1)=0._JPRB
YDVAB%VBH(0)=0._JPRB
YDVAB%VBH(1)=1._JPRB

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERT2',1,ZHOOK_HANDLE)
END SUBROUTINE SUVERT2

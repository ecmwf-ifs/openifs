! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUGEOMETRY1C(YDGEOMETRY,YDML_DYN)

!**** *SUGEOMETRY1C*  - SETUP MODEL GEOMETRY

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        *CALL* *SUGEOMETRY1C*

!        EXPLICIT ARGUMENTS
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS
!        --------------------
!        NONE

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      Filip Vana  fter SUGEOMETRY
!       updated version 29-Jan-2018

!     MODIFICATIONS.
!     --------------

!     ------------------------------------------------------------------

USE TYPE_GEOMETRY,      ONLY : GEOMETRY
USE MODEL_DYNAMICS_MOD, ONLY : MODEL_DYNAMICS_TYPE
USE PARKIND1,           ONLY : JPRB
USE YOMHOOK,            ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN,             ONLY : NULOUT
USE YOMDYNCORE,         ONLY : RPLRADI
USE YOMCVER,            ONLY : SUCVER_GEOM,PRT_CVER_GEOM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), TARGET,    INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE), INTENT(INOUT) :: YDML_DYN

CHARACTER(LEN=35)  ::  CLINE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "sudim.intfb.h"
#include "suvv1.intfb.h"
#include "sugem_naml.intfb.h"
#include "suvert.intfb.h"
#include "susta.intfb.h"
#include "sugem1c_nc.intfb.h"
#include "suinterpolator.intfb.h"


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUGEOMETRY1C',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,  &
 & YDCSGLEG=>YDGEOMETRY%YRCSGLEG,NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT,NGPBLKS=>YDGEOMETRY%YRDIM%NGPBLKS)

!     ------------------------------------------------------------------

CLINE='----------------------------------'

YDGEOMETRY%LNONHYD_GEOM=>YDGEOMETRY%YRVERT_GEOM%LNONHYD_GEOM
YDGEOMETRY%YRCVER => YDGEOMETRY%YRVERT_GEOM%YRCVER

!*    Initialize geometry dimensions
WRITE(NULOUT,*) '------ Set up geometry dimensions  ------',CLINE
CALL SUDIM(YDGEOMETRY,KSUPERSEDE=0)

CALL SUGEM_NAML(YDGEOMETRY)

WRITE(NULOUT,*) '---- Set up cver -------',CLINE
CALL SUCVER_GEOM(YDGEOMETRY%YRCVER,YDGEOMETRY%LNONHYD_GEOM)
CALL PRT_CVER_GEOM(YDGEOMETRY%YRCVER)

!*    Allocate and initialize vertical coordinate
WRITE(NULOUT,*) '---- Set up vertical coordinate -------',CLINE
CALL SUVV1(YDGEOMETRY)

!*    Initialize geometry 
CALL SUGEM1C_NC(YDGEOMETRY)

!*    Initialize vertical geometry
WRITE(NULOUT,*) '---- Set up vertical geometry -------',CLINE
CALL SUVERT(YDGEOMETRY)

!*    Initialize standard atmosphere
WRITE(NULOUT,*) '---- Set up standard atmosphere -----',CLINE
CALL SUSTA(YDGEOMETRY)

!!*    Initialize geometry (second part)
!!     in particuliar, allocate and fill YRGSGEOM and YRCSGEOM
!WRITE(NULOUT,*) '---- Set up geometry, part 2 for global model --------',CLINE
!CALL SUGEM2(YDGEOMETRY)

!!*    Initialize geometry (third part)
!!     in particuliar, reallocate and fill NLOEN, NMEN, YDCSGLEG%RLATI
!WRITE(NULOUT,*) '---- Set up geometry, part 3 --------',CLINE
!CALL SUGEM3(YDGEOMETRY)

!*    Initialize vertical interpolator
WRITE(NULOUT,*) '---- Set up vertical interpolator -------',CLINE
CALL SUINTERPOLATOR(YDGEOMETRY,YDML_DYN%YRDYNA, YDML_DYN%YRSLINT)
IF (YDML_DYN%YRSLINT%YRVSLETA%NRLEVX > 200000) THEN
  WRITE(NULOUT,*) ' The minimum distance of model levels is not enough to ensure accurate '
  WRITE(NULOUT,*) ' performance for SL advection. Use the LREGETA=.t. or change vertical discretization.'
  ! push it artificially to NRLEVX=200000
  CALL ABOR1('SUGEOMETRY1C: Assumed ZDVETAH is too small.')
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGEOMETRY1C',1,ZHOOK_HANDLE)
END SUBROUTINE SUGEOMETRY1C

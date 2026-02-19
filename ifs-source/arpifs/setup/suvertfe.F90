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

SUBROUTINE SUVERTFE(YDGEOMETRY)

!**** *SUVERTFE*  - Setup VERTical Finite Element scheme

!     Purpose.
!     --------
!           Call of initialisation routines for the different 
!           versions of the vertical finite element scheme

!**   Interface.
!     ----------

!     *CALL* SUVERTFE

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!        see below

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mariano Hortal ECMWF
!      Original : 2000-05

!     Modifications.
!     --------------
!      K.Yessad and J.Vivoda: 28-08-2007 Set-up VFE for NH model.
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      J. Vivoda (Oct 2013): new options for VFE-NH
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPRB ,JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRINTLEV, LOUTPUT
USE YOMLUN   , ONLY : NULERR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "sunh_vertfe1d.intfb.h"
#include "sunh_vertfe1dd.intfb.h"
#include "sunh_vertfe3d.intfb.h"
#include "sunh_vertfe3dbc.intfb.h"
#include "sunh_vertfe3dd.intfb.h"
#include "suvertfe1.intfb.h"
#include "suvertfe3.intfb.h"
#include "suvertfe3d.intfb.h"
#include "sunh_vertfespline.intfb.h"
#include "sunh_vertfespline_half.intfb.h"
#include "sunh_vertfespline_inv.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUVERTFE',0,ZHOOK_HANDLE)
ASSOCIATE(YDVFE=>YDGEOMETRY%YRVFE,YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

IF(YDCVER%LVERTFE) THEN
  IF (YDCVER%LVFE_ECMWF) THEN

    IF(YDCVER%NVFE_TYPE == 1) THEN

      ! setup linear finite element scheme
      CALL SUVERTFE1(YDGEOMETRY)
      IF (YDGEOMETRY%LNONHYD_GEOM) THEN
        CALL SUNH_VERTFE1D(YDGEOMETRY)
        ! ----- SUNH_VERTFE1DBC not yet coded ----------------------
        ! CALL SUNH_VERTFE1DBC
        ! Provisionally RDERB=RDERI with additional zeroed columns.
        YDVFE%RDERB(1:NFLEVG,1)=0.0_JPRB
        YDVFE%RDERB(1:NFLEVG,2:NFLEVG+1)=YDVFE%RDERI(1:NFLEVG,1:NFLEVG)
        YDVFE%RDERB(1:NFLEVG,NFLEVG+2)=0.0_JPRB
        ! ----------------------------------------------------------
        CALL SUNH_VERTFE1DD(YDGEOMETRY)
      ENDIF

    ELSEIF (YDCVER%NVFE_TYPE == 3) THEN

      ! setup cubic finite element scheme
      CALL SUVERTFE3(YDGEOMETRY)
      IF (YDGEOMETRY%LNONHYD_GEOM) THEN
        CALL SUNH_VERTFE3D(YDGEOMETRY)
        CALL SUNH_VERTFE3DBC(YDGEOMETRY)
        CALL SUNH_VERTFE3DD(YDGEOMETRY)
      ELSE
        CALL SUVERTFE3D(YDGEOMETRY)
      ENDIF
    ELSE
      WRITE(NULERR,*) ' SUVERTFE: INVALID VALUE OF NVFE_TYPE : ',YDCVER%NVFE_TYPE
      WRITE(NULERR,*) ' FOR LVFE_ECMWF.'
    ENDIF

  ELSE

    ! setup B-spline finite element scheme
    CALL SUNH_VERTFESPLINE(YDGEOMETRY)
    IF (YDGEOMETRY%LNONHYD_GEOM) THEN
      IF ((YDCVER%LVFE_LAPL.AND.YDCVER%LVFE_LAPL_HALF).OR.YDCVER%LVFE_DELNHPRE) THEN
        CALL SUNH_VERTFESPLINE_HALF(YDGEOMETRY)
      ENDIF
      IF (YDCVER%LVFE_GW.OR.YDCVER%LVFE_GW_HALF) THEN
        CALL SUNH_VERTFESPLINE_INV(YDGEOMETRY)
      ENDIF
    ENDIF

  ENDIF

ELSE
  CALL ABOR1(' IN SUVERTFE LVERTFE SHOULD BE TRUE')
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERTFE',1,ZHOOK_HANDLE)
END SUBROUTINE SUVERTFE

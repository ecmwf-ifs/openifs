! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUGEM1C_NC(YDGEOMETRY)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI      ,ROMEGA
USE YOMLUN   , ONLY : NINIGG

#ifdef DOC

!**** *SUGEM1C_NC*  - Initialize geometry parameters

!     Purpose.
!     --------
!           Initialize geometry parameters of the SCM

!**   Interface.
!     ----------

!     *CALL* SUGEM1C_NC

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the SCM

!     Author.
!     -------
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original      94-01-25
!        M. Ko"hler  2000-12-04  modification to NetCDF
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
REAL(KIND=JPRM)    :: ZLAT, ZLON
INTEGER(KIND=JPIM) :: INCID, VARID, ISTATUS

#include "netcdf.inc"

!     ------------------------------------------------------------------
#include "handle_err_nc.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUGEM1C_NC',0,ZHOOK_HANDLE)

!*         1.  OPEN INPUT FILE.
!              ----------------

ISTATUS = NF_OPEN ('scm_in.nc', NF_NOWRITE, INCID)
CALL HANDLE_ERR_NC(ISTATUS)
WRITE(*,*) 'NETCDF-FILE scm_in.nc OPENED ON UNIT ',INCID


!*         2.  READ LATITUDE AND LONGITUDE.
!              ----------------------------

ISTATUS = NF_INQ_VARID     (INCID, 'lat', VARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_GET_VARA_REAL (INCID, VARID, 1, 1, ZLAT)
CALL HANDLE_ERR_NC(ISTATUS)

ISTATUS = NF_INQ_VARID     (INCID, 'lon', VARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_GET_VARA_REAL (INCID, VARID, 1, 1, ZLON)
CALL HANDLE_ERR_NC(ISTATUS)


!*        3.  CONVERSION FROM DEGREES TO RADIANS.
!             -----------------------------------
! Create YRSGEOM_NB
ALLOCATE(YDGEOMETRY%YRGSGEOM_NB%GELAT(1))
ALLOCATE(YDGEOMETRY%YRGSGEOM_NB%GELAM(1))
ALLOCATE(YDGEOMETRY%YRGSGEOM_NB%GEMU(1))
ALLOCATE(YDGEOMETRY%YRGSGEOM_NB%GSQM2(1))
ALLOCATE(YDGEOMETRY%YRGSGEOM_NB%GECLO(1))
ALLOCATE(YDGEOMETRY%YRGSGEOM_NB%GESLO(1))
ALLOCATE(YDGEOMETRY%YRGSGEOM_NB%RCORI(1))
ALLOCATE(YDGEOMETRY%YRGSGEOM_NB%NGPLAT(1))

! Allocate YRGSGEOM:
ALLOCATE(YDGEOMETRY%YRGSGEOM(1))
YDGEOMETRY%YRGSGEOM(1)%GELAT=> YDGEOMETRY%YRGSGEOM_NB%GELAT
YDGEOMETRY%YRGSGEOM(1)%GELAM=> YDGEOMETRY%YRGSGEOM_NB%GELAM
YDGEOMETRY%YRGSGEOM(1)%GEMU => YDGEOMETRY%YRGSGEOM_NB%GEMU
YDGEOMETRY%YRGSGEOM(1)%GSQM2=> YDGEOMETRY%YRGSGEOM_NB%GSQM2
YDGEOMETRY%YRGSGEOM(1)%GECLO=> YDGEOMETRY%YRGSGEOM_NB%GECLO
YDGEOMETRY%YRGSGEOM(1)%GESLO=> YDGEOMETRY%YRGSGEOM_NB%GESLO
YDGEOMETRY%YRGSGEOM(1)%RCORI=> YDGEOMETRY%YRGSGEOM_NB%RCORI
YDGEOMETRY%YRGSGEOM(1)%NGPLAT=> YDGEOMETRY%YRGSGEOM_NB%NGPLAT


YDGEOMETRY%YRGSGEOM_NB%NGPLAT(1)=1

YDGEOMETRY%YRGSGEOM_NB%GELAT(1)=ZLAT*RPI/180._JPRB
YDGEOMETRY%YRGSGEOM_NB%GELAM(1)=ZLON*RPI/180._JPRB

IF (YDGEOMETRY%YRGSGEOM_NB%GELAM(1) < 0.0_JPRB) THEN
  YDGEOMETRY%YRGSGEOM_NB%GELAM(1)=YDGEOMETRY%YRGSGEOM_NB%GELAM(1)+2.0_JPRB*RPI
ENDIF


!*        4.  COMPUTATION OF GEOMETRY PARAMETERS.
!             -----------------------------------

YDGEOMETRY%YRGSGEOM_NB%GEMU=SIN(YDGEOMETRY%YRGSGEOM_NB%GELAT)
YDGEOMETRY%YRGSGEOM_NB%GSQM2=COS(YDGEOMETRY%YRGSGEOM_NB%GELAT)

YDGEOMETRY%YRGSGEOM_NB%GESLO=SIN(YDGEOMETRY%YRGSGEOM_NB%GELAM)
YDGEOMETRY%YRGSGEOM_NB%GECLO=COS(YDGEOMETRY%YRGSGEOM_NB%GELAM)

YDGEOMETRY%YRGSGEOM_NB%RCORI=2.0_JPRB*ROMEGA*YDGEOMETRY%YRGSGEOM_NB%GEMU


!*        5.  CLOSE INPUT FILE.
!             -----------------

ISTATUS = NF_CLOSE (INCID)
CALL HANDLE_ERR_NC(ISTATUS)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUGEM1C_NC',1,ZHOOK_HANDLE)
END SUBROUTINE SUGEM1C_NC

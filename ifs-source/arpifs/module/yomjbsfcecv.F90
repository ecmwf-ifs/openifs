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

MODULE YOMJBSFCECV

USE PARKIND1, ONLY: JPIM, JPRB

IMPLICIT NONE

SAVE


TYPE SKTECV_CONFIG
  !   ------------------------------------------------------------------
  !   JP_NSKTECV_MAX
  !   CP_KNOWN
  !   ------------------------------------------------------------------
  !
  INTEGER(KIND=JPIM)            :: JP_NSKTECV_MAX = 3
  CHARACTER(LEN=2),DIMENSION(3) :: CP_KNOWN = (/'ir', 'mw', 'hi'/)
END TYPE SKTECV_CONFIG

TYPE SKTECV_DATA
  !   ------------------------------------------------------------------
  !   CSENSOR
  !   Z_ECV_DSKT_IR
  !   ------------------------------------------------------------------
  !
  CHARACTER(LEN=2)  , ALLOCATABLE :: CSENSOR(:)
  REAL(KIND=JPRB)   , ALLOCATABLE :: Z_ECV_DSKT_IR(:,:)
END TYPE SKTECV_DATA

TYPE SSHECV_CONFIG
  !   ------------------------------------------------------------------
  !   ------------------------------------------------------------------
  !
END TYPE SSHECV_CONFIG


TYPE SSHECV_DATA
  !   ------------------------------------------------------------------
  !   LFIXEDZTD     True if zero TL/AD for ZTD
  !   LFIXEDSSH     True if zero TL/AD for SSH
  !   LTESTSSH      True for testing when the ocean is not coupled
  !   PSQRTMTCOR    Time correlation matrix
  !   ------------------------------------------------------------------
  !
  LOGICAL                      :: LFIXEDZTD
  LOGICAL                      :: LFIXEDSSH
  LOGICAL                      :: LTESTSSH
END TYPE SSHECV_DATA

TYPE TSLECV_CONFIG
  !   ------------------------------------------------------------------
  !   JP_NSKTECV_MAX
  !   CP_KNOWN
  !   ------------------------------------------------------------------
  !
  INTEGER(KIND=JPIM)            :: JP_NTSLECV_MAX = 4
  CHARACTER(LEN=4),DIMENSION(4) :: CP_KNOWN = (/'skin', 'snow', 'soil', 'ice '/)
  REAL(KIND=JPRB)               :: PCILIM = 0.01_JPRB
  REAL(KIND=JPRB)               :: PSDLIM = 0.01_JPRB
END TYPE TSLECV_CONFIG

TYPE TSLECV_DATA
  !   ------------------------------------------------------------------
  !   CSFCNAME
  !   ------------------------------------------------------------------
  !
  CHARACTER(LEN=4)  , ALLOCATABLE :: CTSLNAME(:)
END TYPE TSLECV_DATA

TYPE SFCECVBAL
  !   ------------------------------------------------------------------
  !   CNAME         Name of the blance type
  !   CFILE         Name of the file with the stored matrix
  !   LHASBAL       True is balance active
  !   PSQRTMTCOR    Correlation matrix
  !   ------------------------------------------------------------------
  !
  CHARACTER(LEN=20)            :: CNAME
  CHARACTER(LEN=20)            :: CFILE
  LOGICAL                      :: LHASBAL
  REAL(KIND=JPRB), ALLOCATABLE :: PSQRTM(:,:,:)
END TYPE SFCECVBAL

TYPE TJBSFCECV
  !   ------------------------------------------------------------------
  !   NSFCECV       Number of field
  !   NDIM2_SFCECV  Size of the second dimension
  !   LSFCRESCALE   True to rescale std using info from wavelet file
  !   PSSFCSTEP     Time step between fields
  !   LBALTIME      True to have time correlation
  !   LBALTEMP      True to have T balance (SKT only)
  !
  !   ------------------------------------------------------------------
  !
  INTEGER(KIND=JPIM)            :: NSFCECV
  INTEGER(KIND=JPIM)            :: NDIM2_SFCECV = 0
  INTEGER(KIND=JPIM)            :: NSFCBAL
  LOGICAL                       :: LSFCRESCALE
  REAL(KIND=JPRB)               :: PSSFCSTEP
  TYPE(SFCECVBAL), ALLOCATABLE  :: YRBAL(:)
END TYPE TJBSFCECV

TYPE(TJBSFCECV), TARGET         :: YRSKTECV
TYPE(TJBSFCECV), TARGET         :: YRSSHECV
TYPE(TJBSFCECV), TARGET         :: YRTSLECV
!
TYPE(SKTECV_CONFIG)             :: YRSKTECVCFG
TYPE(SSHECV_CONFIG)             :: YRSSHECVCFG
TYPE(TSLECV_CONFIG)             :: YRTSLECVCFG
!
TYPE(SKTECV_DATA)               :: YRSKTECVDATA
TYPE(SSHECV_DATA)               :: YRSSHECVDATA
TYPE(TSLECV_DATA)               :: YRTSLECVDATA

#include "abor1.intfb.h"
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

SUBROUTINE SET_SKT_ECV_INC(YDGEOMETRY,YDSURF,YDFIELD_ECV,KMODSEC)

USE GEOMETRY_MOD          , ONLY : GEOMETRY
USE YOMHOOK               , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMJBECV              , ONLY : LSKTECV
USE YOMVAR                , ONLY : LECV
USE YOMLUN                , ONLY : NULOUT
USE FIELD_CONTAINER_GP_MOD, ONLY : FIELD_CONTAINER_GP
USE SURFACE_FIELDS_MIX    , ONLY : TSURF
USE ECV_DEFINITIONS       , ONLY : VID

IMPLICIT NONE

TYPE(GEOMETRY)          ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)             ,INTENT(IN)    :: YDSURF
TYPE(FIELD_CONTAINER_GP),INTENT(INOUT) :: YDFIELD_ECV
INTEGER(KIND=JPIM)      ,INTENT(IN)    :: KMODSEC
  call abor1("oifs/fc-only - SET_SKT_ECV_INC should never be called")


END SUBROUTINE SET_SKT_ECV_INC


SUBROUTINE GET_SFC_ECV(YDGEOMETRY,YDFIELD_ECV,KID,KFIELDS,KBLOCK,KNUMB,KMODSEC,PSFC_ECV,LDFOUND)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMJBECV     , ONLY : YRECVDATA, NDIAECV
USE YOMVAR       , ONLY : LECV
USE YOMLUN       , ONLY : NULOUT
USE ECV_DEFINITIONS       , ONLY : ECVFIELD_ACCESS,VID
USE FIELD_CONTAINER_GP_MOD, ONLY : FIELD_CONTAINER_GP
USE FIELD_CONTAINER_OPER_MOD, ONLY : FIELD_CONTAINER_GPNORM

IMPLICIT NONE

TYPE(GEOMETRY)            ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_CONTAINER_GP)  ,INTENT(INOUT) :: YDFIELD_ECV
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KID
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KBLOCK
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KNUMB
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KMODSEC
REAL(KIND=JPRB)           ,INTENT(INOUT) :: PSFC_ECV(YDGEOMETRY%YRDIM%NPROMA,KFIELDS)
LOGICAL                   ,INTENT(OUT)   :: LDFOUND
  call abor1("oifs/fc-only -  GET_SFC_ECV should never be called")


END SUBROUTINE GET_SFC_ECV

!-----------------------------------------------------------------------
END MODULE YOMJBSFCECV

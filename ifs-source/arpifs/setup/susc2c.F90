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

SUBROUTINE SUSC2C(YDGEOMETRY,YDEPHY,YDML_GCONF,YDPHY,YDDYNA,LDTENC,YDVARIABLES,YDGFL,YDGMV,YDSURF)

!**** *SUSC2C*  Allocate fields in gridpoint space

!     Purpose.
!     -------
!        Extracted allocations from SUSC2B which belong to the OOPS
!        fieldset object.

!     Interface.
!     ---------
!        *CALL*  *SUSC2C*

!        Explicit arguments : None
!        ------------------

!        Implicit arguments : see above "USE MODULE"
!        ------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        ALLO_SURF, SETUP_GMV

!     Author.
!     ------
!        Tomas Wilhelmsson  * ECMWF *
!        Original  : 11-08-30

!     Modifications.
!     -------------
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE YOMPHY       , ONLY : TPHY
USE YOMDYNA      , ONLY : TDYNA
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE VARIABLES_MOD, ONLY : VARIABLES, HAS_MODEL_FIELDS
USE SURFACE_FIELDS_MIX , ONLY : TSURF, ALLO_SURF
USE YOMGFL   , ONLY : TGFL
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LALLOPR
USE GMV_SUBS_MOD , ONLY : SETUP_GMV

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),  INTENT(INOUT) :: YDGEOMETRY
TYPE(TEPHY)     ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TPHY)      ,INTENT(INOUT) :: YDPHY
TYPE(TDYNA)     ,INTENT(INOUT) :: YDDYNA
LOGICAL         ,INTENT(IN)    :: LDTENC
TYPE(VARIABLES), INTENT(IN)    :: YDVARIABLES
TYPE(TGFL) ,     INTENT(INOUT) :: YDGFL
TYPE(TGMV) ,     INTENT(INOUT) :: YDGMV
TYPE(TSURF),     INTENT(INOUT) :: YDSURF

LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!    -------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUSC2C',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(NDIM=>YGFL%NDIM, NDIMPC=>YGFL%NDIMPC, NDIMPT=>YGFL%NDIMPT, &
 & NDIMSLP=>YGFL%NDIMSLP, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG)
!    -------------------------------------------------------------------

CALL GSTATS(1939,0)

LLP = NPRINTLEV >= 1.OR. LALLOPR

!*        1.  BUFFER ALLOCATIONS.
!         -----------------------

!         1.1   SETUP GMV STRUCTURE

CALL SETUP_GMV(YDGEOMETRY,YDGMV,YDEPHY,YDML_GCONF%YRDIMF,YDPHY,YDDYNA,LDTENC)

!*        1.2   SURFACE FIELDS

IF (.NOT. YDVARIABLES%LINEAR) THEN
  CALL ALLO_SURF(YDGEOMETRY%YRDIM,YDSURF)
ENDIF

!*        1.3   ALLOCATE GFL ARRAYS

! NDIM should depend on YDVARIABLES
ALLOCATE(YDGFL%GFL(NPROMA,NFLEVG,NDIM,NGPBLKS))
IF (LLP) WRITE(NULOUT,9) 'GFL      ',SIZE(YDGFL%GFL),SHAPE(YDGFL%GFL)

!IF (HAS_MODEL_FIELDS(YDVARIABLES)) THEN
  ALLOCATE(YDGFL%GFLSLP(NPROMA,NFLEVG,NDIMSLP,NGPBLKS))
  IF (LLP) WRITE(NULOUT,9) 'GFLSLP   ',SIZE(YDGFL%GFLSLP),SHAPE(YDGFL%GFLSLP)

  WRITE(NULOUT,*)'NDIMPT=',NDIMPT
  ALLOCATE(YDGFL%GFLPT(NPROMA,NFLEVG,NDIMPT,NGPBLKS))
  IF (LLP) WRITE(NULOUT,9) 'GFLPT    ',SIZE(YDGFL%GFLPT),SHAPE(YDGFL%GFLPT)

  WRITE(NULOUT,*)'NDIMPC=',NDIMPC
  ALLOCATE(YDGFL%GFLPC(NPROMA,NFLEVG,NDIMPC,NGPBLKS))
  IF (LLP) WRITE(NULOUT,9) 'GFLPC    ',SIZE(YDGFL%GFLPC),SHAPE(YDGFL%GFLPC)
!ENDIF
!    -------------------------------------------------------------------

9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

CALL GSTATS(1939,1)
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSC2C',1,ZHOOK_HANDLE)
END SUBROUTINE SUSC2C

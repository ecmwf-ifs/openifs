! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SUGFL1C(YDGEOMETRY,YDMODEL)

USE GEOMETRY_MOD, ONLY : GEOMETRY
USE TYPE_MODEL  , ONLY : MODEL
USE PARKIND1    , ONLY : JPIM, JPRB
USE YOMHOOK     , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL    , ONLY : JPGFL
USE YOMGPD1C    , ONLY : GFLSLP
USE YOMLUN      , ONLY : NULOUT

!**** *SUGFL1C*  - Initialize definition of unified_treatment grid_point fields

!     Purpose.
!     --------
!           Initialize definition of unified_treatment fields (GFL)
!           SCM version adopted from sugfl.F90 and gfl_subs_mod.F90.  This
!           sets variables not used in SCM callpar and allows for
!           unchanged version of callpar.

!**   Interface.
!     ----------
!        *CALL* *SUGFL1C

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        MODULE YOM_YGFL

!     Method.
!     -------
!        Structure definitons (e.g. YA) are in yom_ygfl.F90, which 
!        are based on types defined in type_gfls.F90. 

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the single column model

!     Author.
!     -------
!        Martin Ko"hler  *ECMWF*

!     Modifications.
!     --------------
!        Original    2006-03-13
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!        F. Vana     22-Jun-2014  further update including LSLPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL),TARGET, INTENT(INOUT) :: YDMODEL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: IGFLPTR, JL

!-------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUGFL1C',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM,YDMP=>YDGEOMETRY%YRMP, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL,YDVAB=>YDGEOMETRY%YRVAB, &
 & LSLPHY=>YDMODEL%YRML_PHY_EC%YREPHY%LSLPHY, &
 & LSLAG=>YDMODEL%YRML_DYN%YRDYNA%LSLAG)
!-------------------------------------------------------------------------

YGFL%YQ   => YGFL%YCOMP(JPGFL)
YGFL%YL   => YGFL%YCOMP(JPGFL-1)
YGFL%YI   => YGFL%YCOMP(JPGFL-2)
YGFL%YS   => YGFL%YCOMP(JPGFL-3)
YGFL%YR   => YGFL%YCOMP(JPGFL-4)
YGFL%YG   => YGFL%YCOMP(JPGFL-5)
YGFL%YA   => YGFL%YCOMP(JPGFL-7)
YGFL%YO3  => YGFL%YCOMP(JPGFL-8)
YGFL%YNOGW=> YGFL%YCOMP(JPGFL-10:JPGFL-9)  ! attention non-orog. grav waves not set up
YGFL%YTKE => YGFL%YCOMP(JPGFL-11)

YGFL%YQ%MP=1
YGFL%YL%MP=2
YGFL%YI%MP=3
YGFL%YR%MP=4
YGFL%YS%MP=5
YGFL%YA%MP=6
YGFL%YNOGW(1)%MP=7
YGFL%YNOGW(2)%MP=8
YGFL%YO3%MP=9
YGFL%YTKE%MP=10

YGFL%YQ%MP9_PH=1
YGFL%YL%MP9_PH=2
YGFL%YI%MP9_PH=3
YGFL%YR%MP9_PH=4
YGFL%YS%MP9_PH=5
YGFL%YA%MP9_PH=6
YGFL%YNOGW(1)%MP9_PH=7
YGFL%YNOGW(2)%MP9_PH=8
YGFL%YO3%MP9_PH=9
YGFL%YTKE%MP9_PH=10

YGFL%NDIM= 10

YGFL%YQ%MP1= 1
YGFL%YL%MP1= 2
YGFL%YI%MP1= 3
YGFL%YA%MP1= 4
YGFL%YR%MP1= 5
YGFL%YS%MP1= 6
YGFL%YO3%MP1=7
YGFL%YTKE%MP1=8
YGFL%NDIM1    = 8

YGFL%NUMFLDS= YGFL%NDIM
YGFL%NDIMSLP= YGFL%NDIM   ! bit waste of space but who cares in 1D model

!YGFL%NDIM     = 5                ! q,l,i,a,o3

YGFL%YQ%LACTIVE    = .TRUE.
YGFL%YA%LACTIVE    = .TRUE.
YGFL%YL%LACTIVE    = .TRUE.
YGFL%YI%LACTIVE    = .TRUE.
YGFL%YR%LACTIVE    = .TRUE.
YGFL%YS%LACTIVE    = .TRUE.
YGFL%YG%LACTIVE    = .FALSE.
YGFL%YO3%LACTIVE   = .FALSE.          ! namelist consistency?
YGFL%YTKE%LACTIVE  = .FALSE.

YGFL%NAERO         = 0                ! number of aerosol types
IGFLPTR=YGFL%NUMFLDS+1
YGFL%YAERO=>YGFL%YCOMP(IGFLPTR:IGFLPTR+YGFL%NAERO-1)
YGFL%YAERO(:)%LPHY = .FALSE.
YGFL%YAERO(:)%LADV = .FALSE.
YGFL%YAERO(:)%MP1  = -HUGE(JPGFL)

YGFL%NAERO_WVL_DIAG = 1
YGFL%NAERO_WVL_DIAG_TYPES = 1

YGFL%NGFL_EXT      = 0                ! number of advected extra variables
YGFL%YEXT=>YGFL%YCOMP(IGFLPTR:IGFLPTR+YGFL%NGFL_EXT-1)
YGFL%YEXT(:)%LPHY  = .FALSE.
YGFL%YEXT(:)%LADV  = .FALSE.
YGFL%YEXT(:)%MP9   = -HUGE(JPGFL)
YGFL%YEXT(:)%MP1   = -HUGE(JPGFL)

IF (LSLPHY.AND.LSLAG) THEN
  ALLOCATE(GFLSLP(YDDIM%NPROMA,YDDIMV%NFLEVG,YGFL%NDIMSLP))
  WRITE(NULOUT,*) 'GFLSLP  ',SIZE(GFLSLP  ),SHAPE(GFLSLP  )
  GFLSLP(:,:,:) = 0.0_JPRB
  YGFL%YQ%LPHY=.TRUE.
  YGFL%YA%LPHY=.TRUE.
  YGFL%YO3%LPHY=.FALSE.
  YGFL%YQ%MPSLP=YGFL%YQ%MP9_PH
  YGFL%YA%MPSLP=YGFL%YA%MP9_PH
ELSE
  ! Minimal allocation 
  ALLOCATE(GFLSLP(1,1,1))
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGFL1C',1,ZHOOK_HANDLE)

!-------------------------------------------------------------------------

END SUBROUTINE SUGFL1C

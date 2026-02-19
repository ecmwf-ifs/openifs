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

SUBROUTINE GPINOZST(YDGEOMETRY,YDOZO,YDDPHY,KSTART,KPROF,KSTGLO,PKOZO)

!**** *GPINOZST* - Ozone physico-chemical properties
!                  and subgrid surface temperature.

!     Purpose. 
!     -------- 

!**   Interface.
!     ----------
!        *CALL* *GPINOZST(...)

!        Explicit arguments :
!        --------------------    

!        INPUT:
!         KSTART       : start of work
!         KPROF        : depth of work
!         KSTGLO       : global offset

!        OUTPUT:
!         PKOZO        : fields for photochemistery of ozone

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None. 
!     ----------
!      Called by CPG_GP.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        ????

! Modifications.
! --------------
!   Original : ????
!   K. Yessad (Sep 2008): add missing comments, cleanings.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
! End Modifications
!------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDPHY      , ONLY : TDPHY
USE YOMMP0       , ONLY : MY_REGION_NS, MY_REGION_EW
USE YOMOZO       , ONLY : TOZO

! ----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDPHY)       ,INTENT(IN)    :: YDDPHY
TYPE(TOZO)        ,INTENT(IN)    :: YDOZO
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKOZO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG*YDDPHY%NVCLIS) 

! ----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IGRILAT(YDGEOMETRY%YRGEM%NGPTOT)

INTEGER(KIND=JPIM) :: ICH, ILAT, IROF, IROFST, ISTART, ISTOP, JFLEV,&
 & JGL, JLON, JROF, JVCLI  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPINOZST',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NTOZ2D=>YDDPHY%NTOZ2D, NTOZ3D=>YDDPHY%NTOZ3D, NVCLIS=>YDDPHY%NVCLIS, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & NFRSTLAT=>YDMP%NFRSTLAT, NLSTLAT=>YDMP%NLSTLAT, NONL=>YDMP%NONL, &
 & NPTRFLOFF=>YDMP%NPTRFLOFF, &
 & TOZ2DL=>YDOZO%TOZ2DL, TOZ3DBL=>YDOZO%TOZ3DBL)
! ----------------------------------------------------------------------------

IROF=0
DO JGL=1,NLSTLAT(MY_REGION_NS)-NFRSTLAT(MY_REGION_NS)+1
  DO JLON=1,NONL(NPTRFLOFF+JGL,MY_REGION_EW)
    IROF=IROF+1
    IGRILAT(IROF)=JGL
  ENDDO
ENDDO
IROFST=KSTGLO
ISTART=KSTART
ISTOP =KPROF

IF (NTOZ3D == 1) THEN
  DO JFLEV=1,NFLEVG
    DO JVCLI=1,NVCLIS
      ICH=(JVCLI-1)*NFLEVG+JFLEV
      IROF=IROFST
      DO JROF=ISTART,ISTOP
        PKOZO(JROF,ICH)=TOZ3DBL(IROF,ICH)
        IROF=IROF+1
      ENDDO
    ENDDO
  ENDDO
ENDIF

IF (NTOZ2D == 1) THEN
  DO JFLEV=1,NFLEVG
    DO JVCLI=1,NVCLIS
      ICH=(JVCLI-1)*NFLEVG+JFLEV
      IROF=IROFST
      DO JROF=ISTART,ISTOP
        ILAT=IGRILAT(IROF)
        PKOZO(JROF,ICH)=TOZ2DL(ILAT,ICH)
        IROF=IROF+1
      ENDDO
    ENDDO
  ENDDO
ENDIF

! ----------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPINOZST',1,ZHOOK_HANDLE)
END SUBROUTINE GPINOZST

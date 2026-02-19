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

SUBROUTINE STATE2SPEC(YDGEOMETRY,YDGFL,YDML_GCONF,YDSP_IN,YDSPEC)

!    Purpose.
!    --------
!      Copy hard coded model state into spectral field.

!    Author.
!    -------
!      Y.Tremolet

!    Modifications.
!    --------------
!      Original:  03-Aug-2005
!      O. Marsden August 2016 Removed use of SPA3, replaced by explicit YDSP_IN argument
!      Y. Michel Sept 2018 : LAM case
!-----------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMGFL                 , ONLY : TGFL
USE PARKIND1               , ONLY : JPRB, JPIM
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SPECTRAL_FIELDS_MOD    , ONLY : ASSIGNMENT(=), SPECTRAL_FIELD
USE YOMCT0      , ONLY : LELAM

IMPLICIT NONE
TYPE(GEOMETRY)               , INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)                   , INTENT(INOUT) :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE), INTENT(INOUT) :: YDML_GCONF
TYPE(SPECTRAL_FIELD)         , INTENT(IN)    :: YDSP_IN
TYPE(SPECTRAL_FIELD)         , INTENT(INOUT) :: YDSPEC
INTEGER(KIND=JPIM) :: JF, JGFL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
#include "dir_trans.h"
#include "edir_trans.h"
#include "abor1.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('STATE2SPEC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDDIMF=>YDML_GCONF%YRDIMF, &
 & YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(NS3D=>YDDIMF%NS3D, NFD2D=>YDDIMF%NFD2D, &
 & NPROMA=>YDDIM%NPROMA, NSMAX=>YDDIM%NSMAX, NRESOL=>YDDIM%NRESOL, &
 & NPSP=>YDMP%NPSP, NBSETLEV=>YDMP%NBSETLEV, &
 & GFL=>YDGFL%GFL)
!-----------------------------------------------------------------------
IF (YDSPEC%NSMAX/=NSMAX) CALL ABOR1('STATE2SPEC: Error resolution')

YDSPEC%SP3D(:,:,1:NS3D)=YDSP_IN%SP3D(:,:,:)
IF (NPSP==1) YDSPEC%SP2D(:,1:NFD2D)=YDSP_IN%SP2D(:,1:NFD2D)

DO JF=1,YDSPEC%NS3D
  DO JGFL=1,YGFL%NUMFLDS
     IF(YDSPEC%NGRIB3(JF)==YGFL%YCOMP(JGFL)%IGRBCODE .AND. YGFL%YCOMP(JGFL)%LGP) THEN
        IF (LELAM) THEN
           CALL EDIR_TRANS(PGP=GFL(:,:,YGFL%YCOMP(JGFL)%MP,:),PSPSCALAR=YDSPEC%SP3D(:,:,JF),&
                & KVSETSC=NBSETLEV,KRESOL=NRESOL,KPROMA=NPROMA)
        ELSE
           CALL DIR_TRANS(PGP=GFL(:,:,YGFL%YCOMP(JGFL)%MP,:),PSPSCALAR=YDSPEC%SP3D(:,:,JF),&
                & KVSETSC=NBSETLEV,KRESOL=NRESOL,KPROMA=NPROMA)
        ENDIF
    ENDIF
  ENDDO
ENDDO
!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('STATE2SPEC',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------------
END SUBROUTINE STATE2SPEC

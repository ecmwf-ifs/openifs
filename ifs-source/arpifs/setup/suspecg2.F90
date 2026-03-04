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

SUBROUTINE SUSPECG2(YDGEOMETRY,YDGFL,YDML_GCONF,PSPA3,PSPA2)

!**** *SUSPECG2*  - Initialize the spectral fields for N3DINI=2
!                   no file need to be read (NSUPERSEDE=0 could be used)

!     Purpose.
!     --------
!           Initialize the spectral fields
!           For analytical initial conditions

!**   Interface.
!     ----------
!        *CALL* *SUSPECG2(.....)*

!        Explicit arguments :
!        --------------------
!               PSPA3, PSPA2 - spectral fields 
!               + possibly spectral fields in GFL

!        Implicit arguments :
!        --------------------
!        The spectral fields of the model.
!        The boundary condition fields of the model.
!         
!     Method.
!     -------
!        See documentation

!     Externals.  
!     ----------  

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Nils Wedi *ECMWF*

!     Modifications.
!     --------------
!        Original : 07-07-20 Nils Wedi
!        K. Yessad (Jan 2010): remove useless variables.
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        S. Malardel : large cleaning (routine now called only if N3DINI=2)
!                      + main DCMIP12 test cases
!     ------------------------------------------------------------------
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMGFL                 , ONLY : TGFL
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN                 , ONLY : NULOUT, NULERR
USE YOM_GRIB_CODES         , ONLY : NGRBZ, NGRBLNSP, NGRBVO, NGRBD, NGRBT, NGRB118, NGRB119
USE YOMDYNCORE             , ONLY : CTEST_FAMILY, NTESTCASE

IMPLICIT NONE

TYPE(GEOMETRY),         INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL),             INTENT(INOUT) :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN) :: YDML_GCONF
REAL(KIND=JPRB),TARGET, INTENT(OUT)   :: PSPA3(:,:,:) 
REAL(KIND=JPRB),        INTENT(OUT)   :: PSPA2(:,:)

!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZTEMP(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZDIV(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZVOR(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZLNSP(YDGEOMETRY%YRDIM%NSPEC2),ZOROG(YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZPRDPD(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)

INTEGER(KIND=JPIM) :: IGRIB2D(YDML_GCONF%YRDIMF%NS2D),IGRIB3D(YDML_GCONF%YRDIMF%NS3D+YDML_GCONF%YGFL%NUMFLDS)
INTEGER(KIND=JPIM) :: I2D,I3D,JGFL,JK

!     ------------------------------------------------------------------
#include "sudcmip12_spec.intfb.h"
#include "sudcmip16_spec.intfb.h"
#include "sumisc_spec.intfb.h"
#include "abor1.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUSPECG2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIMV=>YDGEOMETRY%YRDIMV, YDMP=>YDGEOMETRY%YRMP, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDVAB=>YDGEOMETRY%YRVAB, &
 & YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(NFTHER=>YDDIMF%NFTHER, NGRBSP3=>YDDIMF%NGRBSP3, NGRBSP2=>YDDIMF%NGRBSP2, &
 & NFLEVL=>YDDIMV%NFLEVL, &
 & NPSP=>YDMP%NPSP)
!     ------------------------------------------------------------------

CALL GSTATS(18,0)

WRITE(UNIT=NULOUT,FMT='('' INITIALISING ARTIFICIAL SPECTRAL DATA FOR: '')')
WRITE(UNIT=NULOUT,FMT='('' CTEST_FAMILY = '',A)') CTEST_FAMILY
WRITE(UNIT=NULOUT,FMT='('' NTESTCASE    = '',I3)') NTESTCASE 

!==================================================================
! Initial setup
ZTEMP = 0.0_JPRB
ZDIV = 0.0_JPRB
ZVOR = 0.0_JPRB
ZLNSP = 0.0_JPRB
ZOROG = 0.0_JPRB
ZPRDPD = 0.0_JPRB

!==================================================================
! Initialize according to CTEST_FAMILY and NTESTCASE

IF      (CTEST_FAMILY == "DCMIP12     ") THEN
  CALL SUDCMIP12_SPEC(YDGEOMETRY,NTESTCASE,ZTEMP,ZDIV,ZVOR,ZLNSP,ZOROG)

ELSEIF  (CTEST_FAMILY == "DCMIP16     ") THEN
  CALL SUDCMIP16_SPEC(YDGEOMETRY,NTESTCASE,ZTEMP,ZDIV,ZVOR,ZLNSP,ZOROG)

ELSEIF  (CTEST_FAMILY == "MISC        ") THEN
  CALL SUMISC_SPEC(YDGEOMETRY,NTESTCASE,ZTEMP,ZDIV,ZVOR,ZLNSP,ZOROG)

ELSE
  WRITE(NULERR,'('' INVALID SETUP for CTEST_FAMILY IN SUSPECG2'')')
  CALL ABOR1(' INVALID SETUP CTEST_FAMILY IN SUSPECG2 ')
ENDIF

!==================================================================
! initialize 2D spectral fields from values computed above
!==================================================================
! 2d variables  
I2D = 1
IGRIB2D(I2D) = NGRBLNSP
I2D = I2D+1
IGRIB2D(I2D) = NGRBZ

! DBG : SMALL PLANET FIX
IF (NPSP == 1) THEN
 DO JGFL=1,I2D
   IF( IGRIB2D(JGFL) == NGRBLNSP ) THEN
     PSPA2(:,JGFL) = ZLNSP(:)
   ENDIF    
   IF( IGRIB2D(JGFL) == NGRBZ ) THEN
     PSPA2(:,JGFL) = ZOROG(:)
   ENDIF
 ENDDO
ENDIF

!==================================================================
! initialize 3D spectral fields from values computed above
!==================================================================
! GMV 
I3D = 2+NFTHER
IGRIB3D(1:I3D) = NGRBSP3(1:I3D)

DO JGFL=1,I3D
  DO JK=1,NFLEVL
    IF ( IGRIB3D(JGFL) == NGRBVO) THEN 
      PSPA3(JK,:,JGFL) = ZVOR(JK,:)
    ENDIF
    IF ( IGRIB3D(JGFL) == NGRBD) THEN
      PSPA3(JK,:,JGFL) = ZDIV(JK,:)
    ENDIF
    IF ( IGRIB3D(JGFL) == NGRBT) THEN
      PSPA3(JK,:,JGFL) = ZTEMP(JK,:)
    ENDIF
    IF ( IGRIB3D(JGFL) == NGRB118) THEN
      PSPA3(JK,:,JGFL) = ZPRDPD(JK,:)
    ENDIF
    IF ( IGRIB3D(JGFL) == NGRB119) THEN
      PSPA3(JK,:,JGFL) = 0._JPRB
    ENDIF
  ENDDO
ENDDO


CALL GSTATS(18,1)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSPECG2',1,ZHOOK_HANDLE)
END SUBROUTINE SUSPECG2



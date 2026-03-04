! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GEMS_INIT(YDSURF, YDML_CHEM,YDEPHY,YGFL,YDPHY2, YDRIP, KIDIA, KFDIA, KLEV, KLON, KTRAC, &
    & PGELAM, PGELAT, PSD_VF, KAERO, KCHEM, PLRCH4, &
    & PGFL, PTENGFL, PCEN, PTENC, PTENC_SKF, PCFLX, PSCAV )

!**   INTERFACE.
!     ----------
!          *GEMS_INIT* IS CALLED FROM *CALLPAR/TL/AD*.

! INPUTS:
!  -------

! INPUTS/OUTPUTS:

!-----------------------------------------------------------------------

!     Externals.
!     ---------


!     Author
!    --------
!         2008-10-15, R. Engelen 

!     Modifications :
!    ----------------
!         2009-03-07, P. de Rosnay offline jacobians in SEKF surface analysis
!         2009-10-07, P. Bechtold:  Initialisation of scavenging coefficients
!                                   for wash out/ wet deposition    
!         2011-10-26, L. Jones: Fluxes no longer passed in individually
!         2011-11-11  J. Flemming, chemistry fluxes and dry deposition and diagnostics added
!         K. Yessad (July 2014): Move some variables.
!         2014-04-01  J. Flemming, setting values of chemistry fluxes moved to chem_initflux.F90
!         2018-05-17  J. Flemming, diurnal cycle for fire emissions 

!-----------------------------------------------------------------------

USE MODEL_CHEM_MOD     , ONLY : MODEL_CHEM_TYPE
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL           , ONLY : TYPE_GFLD
USE YOM_GRIB_CODES     , ONLY : NGRBGHG
USE YOMPHY2            , ONLY : TPHY2
USE YOEPHY             , ONLY : TEPHY
USE YOMRIP             , ONLY : TRIP

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)          ,INTENT(INOUT) :: YDSURF
TYPE(TEPHY)          ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_CHEM_TYPE),INTENT(INOUT) :: YDML_CHEM
TYPE(TPHY2)          ,INTENT(INOUT) :: YDPHY2
TYPE(TYPE_GFLD)      ,INTENT(INOUT) :: YGFL
TYPE(TRIP), INTENT(IN):: YDRIP
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KLON 

REAL(KIND=JPRB)  , INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB)  , INTENT(IN)    :: PGELAT(KLON)

! General array of 2D flux fields
REAL(KIND=JPRB)      ,INTENT(IN), TARGET :: PSD_VF(KLON,YDSURF%YSD_VFD%NDIM)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PGFL(KLON,KLEV,YGFL%NDIM), PTENGFL(KLON,KLEV,YGFL%NDIM1)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PLRCH4(KLON,KLEV)
! In the declaration below the INTENT attribute has been removed to comply 
! strict f95 standards. The attribute INTENT(INOUT) should be put back after all 
! compiler (especially NEC) support this extension. - R. El Khatib 04-Jun-2009
REAL(KIND=JPRB), POINTER, INTENT(INOUT) :: PCEN(:,:,:),PTENC(:,:,:),PCFLX(:,:), &
       & PSCAV(:), PTENC_SKF(:,:,:)
INTEGER(KIND=JPIM), INTENT(INOUT) :: KTRAC
INTEGER(KIND=JPIM), INTENT(OUT) :: KAERO(YGFL%NAERO), KCHEM(YGFL%NCHEM)
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: ITRC, ICO2, ICH4
INTEGER(KIND=JPIM) :: JEXT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GEMS_INIT',0,ZHOOK_HANDLE)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & LCHEM_DIA=>YDML_CHEM%YRCOMPO%LCHEM_DIA, LCHEM_DDFLX=>YDML_CHEM%YRCOMPO%LCHEM_DDFLX, &
 & LCHEM_CONVSCAV=>YDML_CHEM%YRCHEM%LCHEM_CONVSCAV, &
 & YAERO=>YGFL%YAERO, YGHG=>YGFL%YGHG, &
 & NCHEM=>YGFL%NCHEM, YLRCH4=>YGFL%YLRCH4, &
 & NGHG=>YGFL%NGHG, NAERO=>YGFL%NAERO, &
 & NACTAERO=>YGFL%NACTAERO, NCHEM_DV=>YGFL%NCHEM_DV, YCHEM=>YGFL%YCHEM, &
 & YSD_VFD=>YDSURF%YSD_VFD, YSD_VF=>YDSURF%YSD_VF, &
 & NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1)
KTRAC=0
ICO2=0
ICH4=0
KAERO(:)=0
KCHEM(:)=0

! Allocate general variables for tracer transport:
! IMPORTANT: Tracer order is : CO2 - -other tracers - react Gases - Aerosol - extra GFL

DO JEXT=1,NGHG
  KTRAC=KTRAC+1
  IF(YGHG(JEXT)%IGRBCODE == NGRBGHG(1)) ICO2=KTRAC
  IF(YGHG(JEXT)%IGRBCODE == NGRBGHG(2)) ICH4=KTRAC
ENDDO
DO JEXT=1,NAERO
 KTRAC=KTRAC+1
 KAERO(JEXT)=KTRAC
ENDDO
DO JEXT=1,NCHEM
 KTRAC=KTRAC+1
 KCHEM(JEXT)=KTRAC
ENDDO

ALLOCATE( PCEN(KLON,KLEV,KTRAC) )
ALLOCATE( PTENC(KLON,KLEV,KTRAC) )
ALLOCATE( PTENC_SKF(KLON,KLEV,KTRAC) )
ALLOCATE( PCFLX(KLON,KTRAC) )
ALLOCATE( PSCAV(KTRAC) )
PCFLX(KIDIA:KFDIA,1:KTRAC)=0.0_JPRB
PCEN (KIDIA:KFDIA,1:KLEV,1:KTRAC)=0.0_JPRB
PTENC(KIDIA:KFDIA,1:KLEV,1:KTRAC)=0.0_JPRB
PSCAV(1:KTRAC)=0.0_JPRB


! to cp this in GEMSL array would not be needed as data are in gfl array anyway 
PLRCH4(KIDIA:KFDIA,1:KLEV) = 0.0_JPRB 
IF(NGHG > 0 .AND. YLRCH4%LGP) THEN
  PLRCH4(KIDIA:KFDIA,1:KLEV) = PGFL(KIDIA:KFDIA,1:KLEV,YLRCH4%MP9_PH)
ENDIF

ITRC=0
DO JEXT=1,NGHG
    ITRC=ITRC+1
    PCEN(KIDIA:KFDIA,1:KLEV,ITRC) =PGFL(KIDIA:KFDIA,1:KLEV,YGHG(JEXT)%MP9_PH)
    PTENC(KIDIA:KFDIA,1:KLEV,ITRC)=PTENGFL(KIDIA:KFDIA,1:KLEV,YGHG(JEXT)%MP1)
ENDDO
DO JEXT=1,NAERO
  ITRC=ITRC+1
  PCEN(KIDIA:KFDIA,1:KLEV,ITRC) =PGFL(KIDIA:KFDIA,1:KLEV,YAERO(JEXT)%MP9_PH)
  PTENC(KIDIA:KFDIA,1:KLEV,ITRC)=PTENGFL(KIDIA:KFDIA,1:KLEV,YAERO(JEXT)%MP1)
ENDDO
DO JEXT=1,NCHEM
  ITRC=ITRC+1
  PCEN(KIDIA:KFDIA,1:KLEV,ITRC) =PGFL(KIDIA:KFDIA,1:KLEV,YCHEM(JEXT)%MP9_PH)
  IF (YCHEM(JEXT)%LADV ) THEN
    PTENC(KIDIA:KFDIA,1:KLEV,ITRC)=PTENGFL(KIDIA:KFDIA,1:KLEV,YCHEM(JEXT)%MP1)
  ELSE
    PTENC(KIDIA:KFDIA,1:KLEV,ITRC)=0.0_JPRB
  ENDIF
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GEMS_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE GEMS_INIT

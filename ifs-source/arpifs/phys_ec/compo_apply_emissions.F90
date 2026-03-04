! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE COMPO_APPLY_EMISSIONS(YDMODEL, &
    & KIDIA, KFDIA, KLEV, KLON, KTRAC, &
    & PGELAT, PGELAM, PDELP, &
    & PCEN, PTENC, PCFLX, &
    & PEMIS2D, PEMIS3D, PEMIS2DAUX, &
    & KAERO, KCHEM, &
    & PEXTRA, &
    & PAPHIF, PBLH, &
    & PLSM)

!-----------------------------------------------------------------------------
! COMPO_APPLY_EMISSIONS
!
!   This routine iterates over the emission fields for the composition species
!   (aerosols, chemistry, GHG) and applies them either to the surface flux or
!   tendency arrays as appropriate.
!
! AUTHOR: Zak Kipling
!-----------------------------------------------------------------------------

USE TYPE_MODEL,ONLY : MODEL

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG
USE YOMCHEM  , ONLY : IEXTR_EM
USE YOMCOMPO , ONLY : TCOMPO_EMIS
USE YOM_YGFL , ONLY : TYPE_GFL_COMP

IMPLICIT NONE

! Standard model objects
TYPE(MODEL), INTENT(IN) :: YDMODEL

! Standard array dimensions
INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA, KFDIA, KLEV, KLON, KTRAC

! General geophysical/meteorological inputs
REAL(KIND=JPRB), INTENT(IN) :: PGELAT(KLON), PGELAM(KLON)
REAL(KIND=JPRB), INTENT(IN) :: PDELP(KLON,KLEV)

! Tracer state, flux and tendency arrays
REAL(KIND=JPRB), INTENT(IN) :: PCEN(KLON,KLEV,KTRAC)
REAL(KIND=JPRB), INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB), INTENT(INOUT) :: PCFLX(KLON,KTRAC) 

! Emission field inputs
REAL(KIND=JPRB), INTENT(IN) :: PEMIS2D(KLON, YDMODEL%YRML_GCONF%YGFL%NEMIS2D)
REAL(KIND=JPRB), INTENT(IN) :: PEMIS3D(KLON, KLEV, YDMODEL%YRML_GCONF%YGFL%NEMIS3D)
REAL(KIND=JPRB), INTENT(IN), TARGET :: PEMIS2DAUX(KLON, YDMODEL%YRML_GCONF%YGFL%NEMIS2DAUX)

! Optional diagnostic output arrays
REAL(KIND=JPRB),INTENT(INOUT) :: PEXTRA(:,:,:)

! Tracer index arrays
INTEGER(KIND=JPIM), INTENT(IN) :: KCHEM(YDMODEL%YRML_GCONF%YGFL%NCHEM)
INTEGER(KIND=JPIM), INTENT(IN) :: KAERO(YDMODEL%YRML_GCONF%YGFL%NAERO)

! Optional inputs only needed for vertical profiles on 2D emissions
REAL(KIND=JPRB), INTENT(IN) :: PAPHIF(KLON,KLEV), PBLH(KLON)

! Legacy auxiliary inputs
REAL(KIND=JPRB), INTENT(IN), OPTIONAL         :: PLSM(KLON)


INTEGER(KIND=JPIM) :: JK, JSPECIES, JEMIS
INTEGER(KIND=JPIM) :: IFOUND

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZSFCFR_LAND, ZSFCFR_SEA, ZTANDEC
REAL(KIND=JPRB), POINTER :: ZINJF(:), ZVOLCALTI(:), ZPARAM(:)
REAL(KIND=JPRB) :: ZPROFILE(KLON, KLEV+1), ZPROFILE_LAND(KLON, KLEV+1), ZPROFILE_SEA(KLON, KLEV+1)
REAL(KIND=JPRB) :: ZFLUX(KLON)
REAL(KIND=JPRB) , DIMENSION (KLON) :: ZDIURN, ZDIURNTMP

LOGICAL :: LLFIREPRESENT(KLON), LLHAVEFIREPRESENT

!
!-----------------------------------------------------------------------
#include "chem_inext.intfb.h"
#include "compo_diurnal.intfb.h"
#include "compo_injection_profile.intfb.h"
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('COMPO_APPLY_EMISSIONS',0,ZHOOK_HANDLE)
ASSOCIATE(YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDERIP=>YDMODEL%YRML_PHY_RAD%YRERIP)
ASSOCIATE(NAERO=>YGFL%NAERO, NACTAERO=>YGFL%NACTAERO, NGHG=>YGFL%NGHG, NCHEM=>YGFL%NCHEM, &
 & YCHEM=>YGFL%YCHEM, YAERO=>YGFL%YAERO, YGHG=>YGFL%YGHG, &
 & NEMIS2D_DESC=>YDCOMPO%NEMIS2D_DESC, YEMIS2D_DESC=>YDCOMPO%YEMIS2D_DESC, &
 & NEMIS3D_DESC=>YDCOMPO%NEMIS3D_DESC, YEMIS3D_DESC=>YDCOMPO%YEMIS3D_DESC, &
 & LCHEM_DIA=>YDCOMPO%LCHEM_DIA, &
 & RSIDECM=>YDERIP%RSIDECM,RCODECM=>YDERIP%RCODECM, &
 & TSPHY=>YDPHY2%TSPHY)

LLHAVEFIREPRESENT = .FALSE.
NULLIFY(ZINJF)
! Look at OM and CO2 emissions to distinguish (non-negligible)  fires
DO JEMIS=1,NEMIS2D_DESC
  IF (IAND(YEMIS2D_DESC(JEMIS)%LEGACY_CHEM_OVERRIDE,16) /= 0) THEN
    IF (YEMIS2D_DESC(JEMIS)%VERTICAL_PARAM_INDEX > 0) THEN
      ZINJF => PEMIS2DAUX(:,YEMIS2D_DESC(JEMIS)%VERTICAL_PARAM_INDEX)
    ENDIF
    LLFIREPRESENT(KIDIA:KFDIA) = ( PEMIS2D(KIDIA:KFDIA,YEMIS2D_DESC(JEMIS)%PARAM_INDEX) > 0.0_JPRB )
    LLHAVEFIREPRESENT = .TRUE.
    EXIT
  ENDIF
ENDDO

NULLIFY(ZVOLCALTI)
DO JEMIS=1,NEMIS2D_DESC
  IF (IAND(YEMIS2D_DESC(JEMIS)%LEGACY_CHEM_OVERRIDE,32) /= 0) THEN
    IF (YEMIS2D_DESC(JEMIS)%VERTICAL_PARAM_INDEX <= 0) THEN
      CALL ABOR1('COMPO_APPLY_EMISSIONS: VERTICAL_PARAM_INDEX NOT SET FINDING LEGACY VOLCALTI')
    ENDIF
    ZVOLCALTI => PEMIS2DAUX(:,YEMIS2D_DESC(JEMIS)%VERTICAL_PARAM_INDEX)
    EXIT
  ENDIF
ENDDO

ZTANDEC = RSIDECM / MAX(RCODECM, 1.0E-12_JPRB)

! New-style flexible 2D emissions processing for chemistry species
DO JEMIS=1,NEMIS2D_DESC
  IF (.NOT. YEMIS2D_DESC(JEMIS)%NON_SIMPLE_TRACER) THEN
    IFOUND = 0

    CALL APPLY_2D_EMISSION( NCHEM, YCHEM, YEMIS2D_DESC(JEMIS), LDDIAGFLUX=LCHEM_DIA, LDDIAGTEND=LCHEM_DIA, &
                          & KDIAGSHIFT=0, KFOUND=IFOUND, KINDEX=KCHEM )

    CALL APPLY_2D_EMISSION( NACTAERO, YAERO, YEMIS2D_DESC(JEMIS), LDDIAGFLUX=.FALSE., LDDIAGTEND=LCHEM_DIA, &
                          & KDIAGSHIFT=NCHEM, KFOUND=IFOUND, KINDEX=KAERO )
    ! No index for GHG because they lead tracer array
    CALL APPLY_2D_EMISSION( NGHG, YGHG, YEMIS2D_DESC(JEMIS), LDDIAGFLUX=.FALSE., LDDIAGTEND=LCHEM_DIA, &
                          & KDIAGSHIFT=NCHEM+NACTAERO, KFOUND=IFOUND )

    IF (IFOUND == 0) THEN
      CALL ABOR1('COMPO_APPLY_EMISSIONS: (2D) No tracer found for species '//TRIM(YEMIS2D_DESC(JEMIS)%SPECIES))
    ELSEIF (IFOUND > 1) THEN
      CALL ABOR1('COMPO_APPLY_EMISSIONS: (2D) Multiple tracers found for species '//TRIM(YEMIS2D_DESC(JEMIS)%SPECIES))
    ENDIF
  ENDIF
ENDDO

! New-style flexible 3D emissions processing for chemistry species
DO JEMIS=1,NEMIS3D_DESC
  IF (.NOT. YEMIS3D_DESC(JEMIS)%NON_SIMPLE_TRACER) THEN
    IFOUND = 0

    CALL APPLY_3D_EMISSION( NCHEM, YCHEM, YEMIS3D_DESC(JEMIS), LDDIAGTEND=LCHEM_DIA, &
                          & KFOUND=IFOUND, KINDEX=KCHEM )

    CALL APPLY_3D_EMISSION( NACTAERO, YAERO, YEMIS3D_DESC(JEMIS), LDDIAGTEND=LCHEM_DIA, &
                          & KFOUND=IFOUND, KINDEX=KAERO )

    ! No index for GHG because they lead tracer array
    CALL APPLY_3D_EMISSION( NGHG, YGHG, YEMIS3D_DESC(JEMIS), &
                          & LDDIAGTEND=LCHEM_DIA, KFOUND=IFOUND )

    IF (IFOUND == 0) THEN
      CALL ABOR1('COMPO_APPLY_EMISSIONS: (3D) No tracer found for species '//TRIM(YEMIS3D_DESC(JEMIS)%SPECIES))
    ELSEIF (IFOUND > 1) THEN
      CALL ABOR1('COMPO_APPLY_EMISSIONS: (3D) Multiple tracers found for species '//TRIM(YEMIS3D_DESC(JEMIS)%SPECIES))
    ENDIF
  ENDIF
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('COMPO_APPLY_EMISSIONS',1,ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE APPLY_2D_EMISSION(KSPECIES, YDSPECIES, YDEMIS_DESC, LDDIAGFLUX, LDDIAGTEND, KDIAGSHIFT, KFOUND, KINDEX)
  
    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN)              :: KSPECIES            ! = NCHEM, NAERO or NGHG
    TYPE(TYPE_GFL_COMP), INTENT(IN)             :: YDSPECIES(KSPECIES) ! = YCHEM, YAERO or YGHG
    TYPE(TCOMPO_EMIS), INTENT(IN)               :: YDEMIS_DESC         ! descriptor for this emission
    LOGICAL, INTENT(IN)                         :: LDDIAGFLUX          ! do diags here for sfc flux (for CHEM, not AERO/GHG)
    LOGICAL, INTENT(IN)                         :: LDDIAGTEND          ! do diags here for tendencies (in 2D, for all)
    INTEGER(KIND=JPIM), INTENT(IN)              :: KDIAGSHIFT          ! index shift for diags in PEXTRA (reordered cf tracers!)
    INTEGER(KIND=JPIM), INTENT(INOUT)           :: KFOUND              ! incremented when emission applied to a field
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL    :: KINDEX(KSPECIES)    ! = KCHEM, KAERO or nothing for GHG

    INTEGER(KIND=JPIM) :: ITRAC     ! index of current tracer in all-tracers arrays (looked up via KINDEX)
    INTEGER(KIND=JPIM) :: IDIAGSLOT ! index of current tracer in PEXTRA diagnostics (after adding KDIAGSHIFT)

    REAL(KIND=JPRB) :: ZSFCFLUX(KLON) ! actual surface flux after elevated emissions applied

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('COMPO_APPLY_EMISSIONS:APPLY_2D_EMISSION',0,ZHOOK_HANDLE)
    ASSOCIATE(YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
            & YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2)
      ASSOCIATE(TSPHY=>YDPHY2%TSPHY)
        DO JSPECIES=1,KSPECIES
          IF (YDSPECIES(JSPECIES)%CNAME == YDEMIS_DESC%SPECIES) THEN
            IF (PRESENT(KINDEX)) THEN
              ITRAC = KINDEX(JSPECIES)
            ELSE
              ITRAC = JSPECIES
            ENDIF

            IDIAGSLOT = KDIAGSHIFT + JSPECIES

            CALL COMPO_DIURNAL(YDRIP, KIDIA, KFDIA, KLON, YDEMIS_DESC%DIURNAL_CYCLE_TYPE, &
             &                 PGELAM, PGELAT, ZDIURN, &
             &                 PAMPLITUDE=YDEMIS_DESC%DIURNAL_AMPLITUDE, &
             &                 PHOURPEAK=YDEMIS_DESC%DIURNAL_PEAK_HOUR, &
             &                 PBASELINE=YDEMIS_DESC%DIURNAL_BASELINE, &
             &                 PTANDEC=ZTANDEC)
      
            ! Legacy support: override diurnal cycle where fire present
            IF (IAND(YDEMIS_DESC%LEGACY_CHEM_OVERRIDE,1) /= 0) THEN
              IF (.NOT. LLHAVEFIREPRESENT) THEN
                CALL ABOR1('COMPO_APPLY_EMISSIONS: FIRE PRESENCE FLAG REQUIRED FOR LEGACY_CHEM_OVERRIDE=1')
              ENDIF
              CALL COMPO_DIURNAL( YDRIP, KIDIA, KFDIA, KLON, 'GFAS', PGELAM, PGELAT, ZDIURNTMP, &
                                & PHOURPEAK=13.5_JPRB, PBASELINE=0.2_JPRB )
              WHERE(LLFIREPRESENT(KIDIA:KFDIA)) ZDIURN(KIDIA:KFDIA) = ZDIURNTMP(KIDIA:KFDIA)
            ENDIF
      
            ! This is the surface / total column flux after scaling and diurnal cycle application
            ZFLUX(KIDIA:KFDIA) = - PEMIS2D(KIDIA:KFDIA,YDEMIS_DESC%PARAM_INDEX) * ZDIURN(KIDIA:KFDIA) * YDEMIS_DESC%SCALING
      
            IF ( YDEMIS_DESC%VERTICAL_PROFILE_TYPE == 'HeightMap' .OR. &
               & YDEMIS_DESC%VERTICAL_PROFILE_TYPE == 'AltitudeMap' .OR. &
               & YDEMIS_DESC%VERTICAL_PROFILE_TYPE == 'GFAS' ) THEN
              IF (YDEMIS_DESC%VERTICAL_PARAM_INDEX <= 0) THEN
                CALL ABOR1('COMPO_APPLY_EMISSIONS: VERTICAL_PARAM_INDEX NOT SET')
              ENDIF
            ENDIF

            NULLIFY(ZPARAM)
            IF (YDEMIS_DESC%VERTICAL_PARAM_INDEX > 0) ZPARAM => PEMIS2DAUX(:,YDEMIS_DESC%VERTICAL_PARAM_INDEX)

            CALL COMPO_INJECTION_PROFILE(KIDIA, KFDIA, KLON, KLEV, YDEMIS_DESC%VERTICAL_PROFILE_TYPE, &
             &                           ZPROFILE, &
             &                           PSFCFRAC=YDEMIS_DESC%VERTICAL_SURFACE_FRACTION, &
             &                           PBASEHEIGHT=YDEMIS_DESC%VERTICAL_BASE_HEIGHT, &
             &                           PTOPHEIGHT=YDEMIS_DESC%VERTICAL_TOP_HEIGHT, &
             &                           PTHRESHOLD=YDEMIS_DESC%VERTICAL_THRESHOLD, &
             &                           KBASELEV=YDEMIS_DESC%VERTICAL_BASE_LEVEL, &
             &                           KTOPLEV=YDEMIS_DESC%VERTICAL_TOP_LEVEL, &
             &                           PPARAM=ZPARAM, &
             &                           PBLH=PBLH, PAPHIF=PAPHIF, PDELP=PDELP)

            ! Legacy support: override injection heights for fires, volcanoes and stacks
            IF (IAND(YDEMIS_DESC%LEGACY_CHEM_OVERRIDE,14) /= 0) THEN
              ! LAST override takes precedence: stacks, volcanoes, fires
      
              ! Legacy support: override injection height for stacks
              IF (IAND(YDEMIS_DESC%LEGACY_CHEM_OVERRIDE,8) /= 0) THEN
                IF (.NOT. PRESENT(PLSM)) THEN
                  CALL ABOR1('COMPO_APPLY_EMISSIONS: PLSM REQUIRED FOR LEGACY_CHEM_OVERRIDE=8')
                ENDIF
                SELECT CASE(YDSPECIES(JSPECIES)%CNAME)
                  CASE ('SO2')
                    ZSFCFR_SEA=0.0_JPRB
                    ZSFCFR_LAND=0.1_JPRB
                  CASE ('CO')
                    ZSFCFR_SEA=0.0_JPRB
                    ZSFCFR_LAND=0.8_JPRB
                  CASE ('NO')
                    ZSFCFR_SEA=0.0_JPRB
                    ZSFCFR_LAND=0.95_JPRB
                  CASE DEFAULT
                    ZSFCFR_SEA=1.0_JPRB
                    ZSFCFR_LAND=1.0_JPRB
                END SELECT
      
                ! ships over oceans 30-100 m 
                CALL COMPO_INJECTION_PROFILE( KIDIA, KFDIA, KLON, KLEV, 'HeightRange', ZPROFILE_SEA, &
                                            & PBASEHEIGHT=30._JPRB, PTOPHEIGHT=100._JPRB, &
                                            & PSFCFRAC=ZSFCFR_SEA, PAPHIF=PAPHIF, PDELP=PDELP )
                ! chimneys over land 100-500 m   
                CALL COMPO_INJECTION_PROFILE( KIDIA, KFDIA, KLON, KLEV, 'HeightRange', ZPROFILE_LAND, &
                                            & PBASEHEIGHT=100._JPRB, PTOPHEIGHT=500._JPRB, &
                                            & PSFCFRAC=ZSFCFR_LAND, PAPHIF=PAPHIF, PDELP=PDELP )

                ! "ANT_HIGH" override doesn't touch anywhere that main profile is above the surface.
                DO JK=1,KLEV+1
                  WHERE (ZPROFILE(KIDIA:KFDIA,KLEV+1) == 1._JPRB .AND. PLSM(KIDIA:KFDIA) <= 0.5_JPRB)
                    ZPROFILE(KIDIA:KFDIA,JK) = ZPROFILE_SEA(KIDIA:KFDIA,JK)
                  ELSEWHERE (ZPROFILE(KIDIA:KFDIA,KLEV+1) == 1._JPRB)
                    ZPROFILE(KIDIA:KFDIA,JK) = ZPROFILE_LAND(KIDIA:KFDIA,JK)
                  ENDWHERE
                ENDDO
              ENDIF
      
              ! Legacy support: override injection heights for volcanoes 
              IF (IAND(YDEMIS_DESC%LEGACY_CHEM_OVERRIDE,4) /= 0) THEN
                IF (.NOT. ASSOCIATED(ZVOLCALTI)) THEN
                  CALL ABOR1('COMPO_APPLY_EMISSIONS: VOLCANO ALTITUDE REQUIRED FOR LEGACY_CHEM_OVERRIDE=4')
                ENDIF
                CALL COMPO_INJECTION_PROFILE( KIDIA, KFDIA, KLON, KLEV, 'AltitudeMap', ZPROFILE_LAND, &
                                            & PPARAM=ZVOLCALTI, KBASELEV=-1, KTOPLEV=-4, PAPHIF=PAPHIF, PDELP=PDELP )
                DO JK=1,KLEV+1
                  WHERE(ZVOLCALTI(KIDIA:KFDIA) > 200._JPRB)
                    ZPROFILE(KIDIA:KFDIA,JK) = ZPROFILE_LAND(KIDIA:KFDIA,JK)
                  ENDWHERE
                ENDDO
              ENDIF 
      
              ! Legacy support: override injection height where fire present
              IF (IAND(YDEMIS_DESC%LEGACY_CHEM_OVERRIDE,2) /= 0) THEN
                IF (.NOT. LLHAVEFIREPRESENT) THEN
                  CALL ABOR1('COMPO_APPLY_EMISSIONS: FIRE PRESENCE FLAG REQUIRED FOR LEGACY_CHEM_OVERRIDE=2')
                ENDIF
                IF (.NOT. ASSOCIATED(ZINJF)) THEN
                  CALL ABOR1('COMPO_APPLY_EMISSIONS: FIRE INJECTION HEIGHT FOR LEGACY_CHEM_OVERRIDE=2')
                ENDIF
                CALL COMPO_INJECTION_PROFILE( KIDIA, KFDIA, KLON, KLEV, 'GFAS', ZPROFILE_LAND, &
                                            & PPARAM=ZINJF, PBLH=PBLH, PAPHIF=PAPHIF, PDELP=PDELP )
                DO JK=1,KLEV+1
                  WHERE(LLFIREPRESENT(KIDIA:KFDIA))
                    ZPROFILE(KIDIA:KFDIA,JK) = ZPROFILE_LAND(KIDIA:KFDIA,JK)
                  ENDWHERE
                ENDDO
              ENDIF
            ENDIF
      
            ! Note: ZPROFILE is always positive, for 1:KLEV and for KLEV+1;
            !       ZFLUX is negative for an emission as per usual IFS convention
            DO JK=1,KLEV
              PTENC(KIDIA:KFDIA,JK,ITRAC) = PTENC(KIDIA:KFDIA,JK,ITRAC) - ZPROFILE(KIDIA:KFDIA,JK) * ZFLUX(KIDIA:KFDIA)
            ENDDO
            IF (LDDIAGTEND) THEN
              ! FIXME: this assumes that if we're doing diagnostics here, we're working on CHEM,
              !   which are the leading elements in the EXTRA array (*unlike* the tracer arrays,
              !   where GHG comes first!)
              PEXTRA(KIDIA:KFDIA,IDIAGSLOT,IEXTR_EM) = PEXTRA(KIDIA:KFDIA,IDIAGSLOT,IEXTR_EM) &
                                                     & - (1.0_JPRB - ZPROFILE(KIDIA:KFDIA,KLEV+1)) * ZFLUX(KIDIA:KFDIA) * TSPHY
            ENDIF
      
            ZSFCFLUX(KIDIA:KFDIA) = ZPROFILE(KIDIA:KFDIA,KLEV+1) * ZFLUX(KIDIA:KFDIA)
    
            PCFLX(KIDIA:KFDIA,ITRAC) = PCFLX(KIDIA:KFDIA,ITRAC) + ZSFCFLUX(KIDIA:KFDIA)
    
            IF (LDDIAGFLUX) THEN 
              ! FIXME: this assumes that if we're doing diagnostics here, we're working on CHEM,
              !   which are the leading elements in the EXTRA array (*unlike* the tracer arrays,
              !   where GHG comes first!)
              PEXTRA(KIDIA:KFDIA,IDIAGSLOT,IEXTR_EM) = PEXTRA(KIDIA:KFDIA,IDIAGSLOT,IEXTR_EM) - ZSFCFLUX(KIDIA:KFDIA) * TSPHY
            ENDIF
    
            KFOUND = KFOUND + 1
          ENDIF
        ENDDO
      END ASSOCIATE
    END ASSOCIATE
    IF (LHOOK) CALL DR_HOOK('COMPO_APPLY_EMISSIONS:APPLY_2D_EMISSION',1,ZHOOK_HANDLE)
  END SUBROUTINE APPLY_2D_EMISSION


  SUBROUTINE APPLY_3D_EMISSION(KSPECIES, YDSPECIES, YDEMIS_DESC, LDDIAGTEND, KFOUND, KINDEX)
  
    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN)              :: KSPECIES            ! = NCHEM, NAERO or NGHG
    TYPE(TYPE_GFL_COMP), INTENT(IN)             :: YDSPECIES(KSPECIES) ! = YCHEM, YAERO or YGHG
    TYPE(TCOMPO_EMIS), INTENT(IN)               :: YDEMIS_DESC         ! descriptor for this emission
    LOGICAL, INTENT(IN)                         :: LDDIAGTEND          ! do diags here for tendencies (in 3D, for all)
    INTEGER(KIND=JPIM), INTENT(INOUT)           :: KFOUND              ! incremented when emission applied to a field
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL    :: KINDEX(KSPECIES)    ! = KCHEM, KAERO or nothing for GHG

    INTEGER(KIND=JPIM) :: ITRAC ! index of current tracer in all-tracers arrays (looked up via KINDEX)

    REAL(KIND=JPRB) :: ZEMIS3D(KLON,KLEV) ! emission tendency after diurnal cycle and scaling
    REAL(KIND=JPRB) :: ZTENC0(KLON,KLEV) ! original tendency for diagnostics

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('COMPO_APPLY_EMISSIONS:APPLY_3D_EMISSION',0,ZHOOK_HANDLE)
    ASSOCIATE(YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
            & YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
            & YDCHEM=>YDMODEL%YRML_CHEM%YRCHEM)
      ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
              & IEXTR_FREE=>YDCHEM%IEXTR_FREE)
        DO JSPECIES=1,KSPECIES
          IF (YDEMIS_DESC%LEGACY_CHEM_OVERRIDE /= 0) THEN
            CALL ABOR1('COMPO_APPLY_EMISSIONS: LEGACY_CHEM_OVERRIDE NOT SUPPORTED FOR 3D EMISSIONS')
          ENDIF

          IF (YDSPECIES(JSPECIES)%CNAME == YDEMIS_DESC%SPECIES) THEN
            IF (PRESENT(KINDEX)) THEN
              ITRAC = KINDEX(JSPECIES)
            ELSE
              ITRAC = JSPECIES
            ENDIF

            CALL COMPO_DIURNAL(YDRIP, KIDIA, KFDIA, KLON, YDEMIS_DESC%DIURNAL_CYCLE_TYPE, &
             &                 PGELAM, PGELAT, ZDIURN, &
             &                 PAMPLITUDE=YDEMIS_DESC%DIURNAL_AMPLITUDE, &
             &                 PHOURPEAK=YDEMIS_DESC%DIURNAL_PEAK_HOUR, &
             &                 PBASELINE=YDEMIS_DESC%DIURNAL_BASELINE, &
             &                 PTANDEC=ZTANDEC)

            IF (LDDIAGTEND) ZTENC0(KIDIA:KFDIA,1:KLEV) = PTENC(KIDIA:KFDIA,1:KLEV,ITRAC) 

            ZEMIS3D(KIDIA:KFDIA,1:KLEV) = PEMIS3D(KIDIA:KFDIA,1:KLEV,YDEMIS_DESC%PARAM_INDEX) * YDEMIS_DESC%SCALING
            DO JK=1,KLEV
              ZEMIS3D(KIDIA:KFDIA,JK) = ZEMIS3D(KIDIA:KFDIA,JK) * ZDIURN(KIDIA:KFDIA)
            ENDDO
            PTENC(KIDIA:KFDIA,1:KLEV,ITRAC) = PTENC(KIDIA:KFDIA,1:KLEV,ITRAC) &
                                            & + ZEMIS3D(KIDIA:KFDIA,1:KLEV) * (RG / PDELP(KIDIA:KFDIA,1:KLEV))

            ! ZKFIXME: We can't special-case a fixed IEXTR_FREE field here like was done for hard-coded emissions
            ! We must either output it as a potential category for all tracers or not at all.
            IF (LDDIAGTEND) THEN
              SELECT CASE (YDSPECIES(JSPECIES)%CNAME)
                CASE ('NO2')
                  CALL CHEM_INEXT(KIDIA, KFDIA, KLON, KLEV, 1, 1, PDELP, TSPHY, PTENC(:,:,ITRAC), ZTENC0(:,:), &
                   & PEXTRA(:,IEXTR_FREE(1,1), IEXTR_EM))
                CASE ('NO')
                  CALL CHEM_INEXT(KIDIA, KFDIA, KLON, KLEV, 1, 1, PDELP, TSPHY, PTENC(:,:,ITRAC), ZTENC0(:,:), &
                   & PEXTRA(:,IEXTR_FREE(1,2), IEXTR_EM))
                CASE ('CO2_GHG')
                  CALL CHEM_INEXT(KIDIA, KFDIA, KLON, KLEV, 1, 1, PDELP, TSPHY, PTENC(:,:,ITRAC), ZTENC0(:,:), &
                   & PEXTRA(:,IEXTR_FREE(1,5), IEXTR_EM))
              END SELECT
            ENDIF

            KFOUND = KFOUND + 1
          ENDIF
        ENDDO
      END ASSOCIATE
    END ASSOCIATE
    IF (LHOOK) CALL DR_HOOK('COMPO_APPLY_EMISSIONS:APPLY_3D_EMISSION',1,ZHOOK_HANDLE)
  END SUBROUTINE APPLY_3D_EMISSION

END SUBROUTINE COMPO_APPLY_EMISSIONS

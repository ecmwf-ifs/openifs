! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE COMPO_FLUX_UPDATE(YDCHEM,YDML_GCONF,YDEPHY,KIDIA,KFDIA,KLON,KTRAC, &
 &                            PCVL,PCVH,PCO2FLUX,PCH4FLUX,PAG,PRECO,PCGPP,PCREC,PCFLXIN)


!***

!**   *COMPO_FLUX_UPDATE* - UPDATE TOTAL FLUX OF ATMOSPHERIC COMPOSITION TRACERS WITH MODELLED FLUXES
 
!     PURPOSE.
!     --------


!     INTERFACE.
!     ----------


!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET


!     INPUT PARAMETERS (REAL)

!    *PAG*          GROSS CO2 ASSIMILATION OVER CANOPY            KG_CO2/M2/S
!    *PRECO*        ECOSYSTEM RESPIRATION                         KG_CO2/M2/S
!    *PCGPP*        GPP flux adjustment coefficient               -
!    *PCREC*        REC flux adjustment coefficient               -
!    *PCVL*         LOW VEGETATION COVER                          -  
!    *PCVH*         HIGH VEGETATION COVER                         -  

!     UPDATED PARAMETERS (REAL):

!    *PCFLXIN*     TOTAL TRACER FLUX                             KG/M2/S
!    *PCO2FLUX*     CO2 FLUX                                      KG_CO2/M2/S
!    *PCH4FLUX*     CH4 FLUX					KG/M2/S

!     METHOD.
!     -------


!     REFERENCE.
!     ----------



!     AUTHOR.
!     -------
!     Anna Agusti-Panareda

!     MODIFICATIONS.
!     --------------
!
!     ------------------------------------------------------------------

USE YOMCHEM  , ONLY : TCHEM
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOEPHY   , ONLY : TEPHY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOM_GRIB_CODES, ONLY : NGRBGHG

!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE(TCHEM)       ,INTENT(IN)    :: YDCHEM
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
TYPE(TEPHY)       ,INTENT(IN)    :: YDEPHY

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA


REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCO2FLUX(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCH4FLUX(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAG(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRECO(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGPP(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCREC(KLON)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFLXIN(KLON,KTRAC)


!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL,  JEXT, IFCO2NOBIO, IFCO2NOFLX, IFCO2, ICHEM
REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COMPO_FLUX_UPDATE',0,ZHOOK_HANDLE)

ASSOCIATE(YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(CHEM_SCHEME=>YDCHEM%CHEM_SCHEME, &
 & NGHG=>YGFL%NGHG, YGHG=>YGFL%YGHG, &
 & NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM, NAERO=>YGFL%NAERO , &
 & LBFASCO2=>YDEPHY%LBFASCO2, LNEEONLINE=>YDEPHY%LNEEONLINE, LWETONLINE=>YDEPHY%LWETONLINE)

!     ------------------------------------------------------------------


! Note that PCFLXIN(KLON,KTRAC) array may contain all the tracers (KTRAC=NGHG+NCHEM+NAERO)
! but there is the assumption that the first tracers are GHG, then AERO and then CHEM
! following the aggregation of the prescribed fluxes in gems_init and chem_init


IF (LNEEONLINE) THEN 
  !2.b. UPDATE THE CO2 BIOGENIC FLUXES READ IN gems_init WITH THOSE FROM CTESSEL (PCO2FLUX)
  IF (NGHG>0) THEN  
  !2.2.1. APPLY FLUX ADJUSTMENT TO GPP AND REC FROM CTESSEL
    DO JEXT=1,NGHG 
      IF (YGHG(JEXT)%IGRBCODE == NGRBGHG(1)) THEN !CO2
        DO JL=KIDIA,KFDIA
          IF ((PCO2FLUX(JL) > -9000._JPRB).AND.((PCVL(JL) > 0._JPRB).OR.(PCVH(JL) > 0._JPRB))) THEN
            IF (LBFASCO2) THEN !Apply flux adjustment
                 PCO2FLUX(JL) = PCGPP(JL) * PAG(JL) + PCREC(JL) * PRECO(JL)
            ENDIF
            PCFLXIN(JL,JEXT) = PCFLXIN(JL,JEXT) +  PCO2FLUX(JL)
          ENDIF
        ENDDO
      ENDIF
    ENDDO  
  ENDIF  

  !2.b. UPDATE THE CO2 BIOGENIC FLUXES FOR CO2 BIOGENIC TRACER (AS DONE FOR CO2)
  !     Note that the order of the chemical tracers is harcoded in gems_init (nghg,naero,nchem)
  IF (NCHEM > 0 .AND. TRIM(CHEM_SCHEME) == "carbontracers") THEN
 
   DO JEXT=1,NCHEM
     ICHEM = NGHG + NAERO + JEXT
     IFCO2NOBIO = INDEX(TRIM( YCHEM(JEXT)%CNAME ), 'CO2_NBI')
     IFCO2NOFLX = INDEX(TRIM( YCHEM(JEXT)%CNAME ), 'CO2_NFX')
     IFCO2 = INDEX(TRIM( YCHEM(JEXT)%CNAME ), 'CO2')
     IF ( IFCO2NOBIO == 0 .AND. IFCO2NOFLX == 0 .AND. IFCO2 == 1) THEN
      DO JL=KIDIA,KFDIA
       IF ((PCO2FLUX(JL) > -9000._JPRB).AND.((PCVL(JL) > 0._JPRB).OR.(PCVH(JL) > 0._JPRB))) THEN
             IF (LBFASCO2) THEN !Apply flux adjustment
                PCO2FLUX(JL) = PCGPP(JL) * PAG(JL) + PCREC(JL) * PRECO(JL)
             ENDIF
  !orig           PCFLX(JL,ICHEM) = PCFLX(JL,ICHEM) + PCO2FLUX(JL)
                 PCFLXIN(JL,ICHEM) = PCFLXIN(JL,ICHEM) +  PCO2FLUX(JL)
       ENDIF
      ENDDO
     ENDIF
   ENDDO  
  ENDIF  

ENDIF  

IF (LWETONLINE) THEN
  !2.b. UPDATE THE CH4 WETLAND FLUXES READ IN gems_init WITH THOSE FROM CTESSEL (PCH4FLUX)
  IF (NGHG>0) THEN
  !2.2.1. APPLY FLUX 
    DO JEXT=1,NGHG
      IF (YGHG(JEXT)%IGRBCODE == NGRBGHG(2)) THEN !CH4
        DO JL=KIDIA,KFDIA
          IF (PCH4FLUX(JL) > -9000._JPRB) THEN
            PCFLXIN(JL,JEXT) = PCFLXIN(JL,JEXT) +  PCH4FLUX(JL)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('COMPO_FLUX_UPDATE',1,ZHOOK_HANDLE)
END SUBROUTINE COMPO_FLUX_UPDATE

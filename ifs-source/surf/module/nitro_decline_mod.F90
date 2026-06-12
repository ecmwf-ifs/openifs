MODULE NITRO_DECLINE_MOD
CONTAINS
SUBROUTINE NITRO_DECLINE(KIDIA,KFDIA,KVT,KLON,KSTEP, PTSPHY,&
 & PCV,&
 & PTSKM1M,PTSOIL,PLAT,PLAI,&
 & YDCST,YDAGS,&
 & PANFM,PRESPBSTR,PRESPBSTR2,&
 & PBIOMASS_LAST,PBIOMASSTR_LAST,PBIOMASSTR2_LAST,&
 & PBIOMASS,PBLOSS)

! (C) Copyright 2003- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!**  *NITRO_DECLINE* 

!     PURPOSE
!     -------
!     Calculates the time change in LAI due to senescence.
!     This in turn changes the dry biomass of the canopy.

!     INTERFACE
!     ---------
!     NITRO_DECLINE IS CALLED BY VEGETATION_EVOL

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!     *KVT*          NUMBER OF VEGETATION TYPE
!                    1  DECIDUOUS
!                    2  CONIFEROUS
!                    3  EVERGREEN
!                    4  C3 GRASS
!                    5  C4 GRASS
!                    6  C3 CROPS
!                    7  C4 CROPS 

! MATCHING BATS table with ECOCLIMAP table 

! (1)  ! Crops, Mixed Farming			=>! 7  C4 CROPS
! (2)  ! Short Grass				=>! 4  C3 GRASS
! (3)  ! Evergreen Needleleaf Trees		=>! 2  CONIFEROUS
! (4)  ! Deciduous Needleleaf Trees		=>! 2  CONIFEROUS
! (5)  ! Deciduous Broadleaf Trees		=>! 1  DECIDUOUS
! (6)  ! Evergreen Broadleaf Trees		=>! 3  EVERGREEN
! (7)  ! Tall Grass				=>! 5  C4 GRASS
! (8)  ! Desert					=> ?
! (9)  ! Tundra					=>! 5  C4 GRASS
! (10) ! Irrigated Crops			=>! 6  C3 CROPS
! (11) ! Semidesert				=>! 5  C4 GRASS
! (12) ! Ice Caps and Glaciers			=>?
! (13) ! Bogs and Marshes			=>! 4  C3 GRASS
! (14) ! Inland Water				=>?
! (15) ! Ocean					=>?
! (16) ! Evergreen Shrubs			=>! 5  C4 GRASS
! (17) ! Deciduous Shrubs			=>! 4  C3 GRASS
! (18) ! Mixed Forest/woodland			=>! 1  DECIDUOUS
! (19) ! Interrupted Forest			=>! 1  DECIDUOUS
! (20) ! Water and Land Mixtures		=>! 4  C3 GRASS

!    *KLON*           NUMBER OF GRID POINT
!    *KSTEP*        CURRENT TIME STEP INDEX

!     INPUT PARAMETERS (REAL):

!    *PTSPHY*       TIME STEP FOR THE PHYSICS
!    *PCV*          TILE FRACTION   
!    *PTSKM1M*      SKIN TEMPERATURE                              K
!    *PTSOIL*       SOIL TEMPERATURE OF LAYER X (SEE CALL FROM VEGETATION_EVOL) K
!    *PLAT*         LATITUDE
!    *PLAI*         LEAF AREA INDEX                           M2/M2

!     UPDATED PARAMETERS (REAL)

!    *PANFM*        MAXIMUM LEAF ASSIMILATION              KG_CO2/KG_AIR M/S          
!    *PRESPBSTR*    RESPIRATION OF ABOVE GROUND STRUCTURAL BIOMASS    KG_CO2/M2
!    *PRESPBSTR2*   RESPIRATION OF BELOW GROUND STRUCTURAL BIOMASS    KG_CO2/M2
!    *PBIOMASS_LAST* (ACTIVE) LEAF BIOMASS OF PREVIOUS DAY               KG/M2
!    *PBIOMASSTR_LAST* ABOVE GROUND STRUCTURAL BIOMASS OF PREVIOUS DAY   KG/M2
!    *PBIOMASSTR2_LAST* BELOW GROUND STRUCTURAL BIOMASS OF PREVIOUS DAY  KG/M2
!    *PBIOMASS*     ACTIVE BIOMASS                  KG/M2
!    *PBLOSS* 	    ACTIVE BIOMASS LOSS             KG/M2

!     METHOD
!     ------
!     Calvet and Soussana  (2001)

!     REFERENCE
!     ---------
!     Calvet and Soussana (2001), "Modelling CO2-enrichment effects using an
!     interactive vegetation SVAT scheme", Agricultural and Forest Meteorology, Vol. 108
!     pp. 129-152
      
!     AUTHOR
!     ------
! 	A. Boone           * Meteo-France *
!       V. Rivalland
!       (following Belair)

!     MODIFICATIONS
!     -------------
!     Original    27/01/03 
!     M.H. Voogt (KNMI) "C-Tessel"  09/2005 
!     S. lafont (ECMWF) externalised "C-Tessel" 05/2006

!-------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_AGS  , ONLY : TAGS

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KVT(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCV(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSOIL(:)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAI(:) 
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TAGS)        ,INTENT(IN)    :: YDAGS
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PANFM(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRESPBSTR(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRESPBSTR2(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBIOMASS_LAST(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBIOMASSTR_LAST(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBIOMASSTR2_LAST(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBIOMASS(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBLOSS(:) 

!*         0.     LOCAL VARIABLES.
!                 ----- ----------

REAL(KIND=JPRB) :: ZLAIB_NITRO  ! LAI correction parameter: 
! for high LAI values, the assimilation is underestimated, becasue
! of the shadow effect. The correction increases the e-folding time and as a
! result lowers mortality to avoid an underestimation of the biomass. 
! Based on simulations with ISBA-Ags by Meteo-France.
REAL(KIND=JPRB) :: ZSEFOLD  ! e-folding time for senescence corrected (days)
REAL(KIND=JPRB) :: ZXM  ! Biomass loss (kg m-2)
REAL(KIND=JPRB) :: ZBIOMASSTOTLIM
REAL(KIND=JPRB) :: ZMORTBACT
REAL(KIND=JPRB) :: ZMORTBSTR
REAL(KIND=JPRB) :: ZMORTBSTR2
REAL(KIND=JPRB) :: ZASSIM
REAL(KIND=JPRB) :: ZSTORAGE
REAL(KIND=JPRB) :: ZCC_NITRO
REAL(KIND=JPRB) :: ZCA_NITRO
REAL(KIND=JPRB) :: ZBIOMASSTOT
REAL(KIND=JPRB) :: ZBIOMASSTR
REAL(KIND=JPRB) :: ZBIOMASSTR2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: Z1,Z2,zTENTHLOG2
INTEGER(KIND=JPIM) :: JL,IVT
!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('NITRO_DECLINE_MOD:NITRO_DECLINE',0,ZHOOK_HANDLE)
ASSOCIATE(RCA_1X_CO2_NIT=>YDAGS%RCA_1X_CO2_NIT, &
 & RCA_2X_CO2_NIT=>YDAGS%RCA_2X_CO2_NIT, RCC_NIT=>YDAGS%RCC_NIT, &
 & RCNS_NIT=>YDAGS%RCNS_NIT, RRESPFACTOR_NIT=>YDAGS%RRESPFACTOR_NIT, &
 & RTAULIM=>YDAGS%RTAULIM, RVANMAX=>YDAGS%RVANMAX, RVBSLAI=>YDAGS%RVBSLAI, &
 & RVBSLAI_NITRO=>YDAGS%RVBSLAI_NITRO, RVGMES=>YDAGS%RVGMES, &
 & RVLAIMIN=>YDAGS%RVLAIMIN, RVSEFOLD=>YDAGS%RVSEFOLD, &
 & RDAY=>YDCST%RDAY, RPI=>YDCST%RPI, RTT=>YDCST%RTT)

ZCC_NITRO=0.0_JPRB
ZSEFOLD=0.0_JPRB
ZLAIB_NITRO=0.0_JPRB
ZXM=0.0_JPRB
ZMORTBACT=0.0_JPRB
ZMORTBSTR=0.0_JPRB
ZMORTBSTR2=0.0_JPRB
ZASSIM=0.0_JPRB
ZSTORAGE=0.0_JPRB
ZBIOMASSTOT=0.0_JPRB
ZBIOMASSTR=0.0_JPRB
ZBIOMASSTR2=0.0_JPRB
ZBIOMASSTOTLIM=0.0_JPRB
!--------------------------------------------------------------------


! coef c for biomass in kg/m2 
  ZCC_NITRO=RCC_NIT/(10._JPRB**RCA_1x_CO2_NIT)
  ZCA_NITRO=RCA_1x_CO2_NIT
  
  
! limited total biomass when structure biomass=0
  ZBIOMASSTOTLIM=ZCC_NITRO**(1.0_JPRB/ZCA_NITRO)
  
! Deep structural biomass  
DO JL=KIDIA,KFDIA
   IVT=KVT(JL)
IF (PCV(JL) > 0._JPRB) THEN
  ZBIOMASSTR2=PBIOMASSTR2_LAST(JL)

  ZLAIB_NITRO=MAX(5.76_JPRB-0.64_JPRB*ATAN(ABS(PLAT(JL))*RPI/180._JPRB),3.8_JPRB)

  ZSEFOLD=RVSEFOLD(IVT)*MIN(1.0_JPRB,PANFM(JL)/RVANMAX(IVT))*   &
   & MAX(((RVGMES(IVT)*1000._JPRB)**0.321_JPRB)*PLAI(JL)/ZLAIB_NITRO,1._JPRB)/RDAY


!  ZSEFOLD=SIGN(MAX(1.0E-8_JPRB,ABS(ZSEFOLD)),ZSEFOLD)
  ZSEFOLD=MAX(1.0E-8_JPRB,ZSEFOLD)
! In order to avoid immediate loss of total vegetation when maximum leaf 
! assimilation is small, the minimum e-folding time is set to a percentage
! (the user-defined RTAULIM) of the maximum e-folding time.  
  ZSEFOLD=MAX((RTAULIM/100._JPRB)*RVSEFOLD(IVT)/(RDAY),ZSEFOLD)
  
! senesence of active biomass
  ZXM=MIN(PBIOMASS(JL)-RVLAIMIN(IVT)*RVBSLAI_NITRO(IVT), &
   & PBIOMASS(JL)*(1.0_JPRB-EXP(-1.0_JPRB/ZSEFOLD)))
! avoid negative values due to computation precision
  ZXM=MAX(ZXM,0.0_JPRB)

! daily active biomass assimilation
  ZASSIM=PBIOMASS(JL)-PBIOMASS_LAST(JL)

! new biomass value:
  PBIOMASS(JL)=PBIOMASS(JL)-ZXM
  
! senesence of deep-structural biomass
  ZMORTBSTR2=ZBIOMASSTR2*(1.0_JPRB-EXP(-1.0_JPRB*RDAY/RVSEFOLD(IVT)))
!1
!ENDIF

  IF (PBIOMASS(JL) .GE. PBIOMASS_LAST(JL)) THEN
! growing phase (net assimilation exceeds the B-decline term): the plant N
! decline model can be applied

!   the growth allometric law is applied
!   repartition of total biomass
    ZBIOMASSTOT=MAX(PBIOMASS(JL),(PBIOMASS(JL)/ZCC_NITRO)**(1.0_JPRB/(1.0_JPRB-ZCA_NITRO)))

!   above-ground structure biomass increment and storage
    ZBIOMASSTR=ZBIOMASSTOT-PBIOMASS(JL)
    ZMORTBSTR=ZBIOMASSTR*(1.0_JPRB-EXP(-1.0_JPRB*RDAY/RVSEFOLD(IVT)))
!2
    ZSTORAGE=ZBIOMASSTR-PBIOMASSTR_LAST(JL)+ZMORTBSTR+PRESPBSTR(JL)

  ELSE
! senescent phase (the B-decline term exceeds net assimilation ): the plant N
! decline model cannot be applied  
    
!   the structural biomass dies exponetially at the lowest rate    
    ZSTORAGE=0.0_JPRB
    ZMORTBSTR=PBIOMASSTR_LAST(JL)*(1.0_JPRB-EXP(-1.0_JPRB*RDAY/RVSEFOLD(IVT)))
!3
    ZBIOMASSTR=PBIOMASSTR_LAST(JL)-ZMORTBSTR-PRESPBSTR(JL)

!   Avoid negative values of biomass
!   No test on ZMORTBSTR as it is not used after, or recalculated
!   No test on PRESPBSTR as it should be smaller than PBIOMASSTR_LAST,
!   otherwise there are irrealistic values of temperature 
    ZBIOMASSTR=MAX(ZBIOMASSTR,0.0_JPRB)

    ZBIOMASSTOT=PBIOMASS(JL)+ZBIOMASSTR

  ENDIF

! flow to the deep structural biomass 
  ZBIOMASSTR2=ZBIOMASSTR2-MIN(ZSTORAGE,0.0_JPRB)
  ZSTORAGE=MAX(0.0_JPRB,ZSTORAGE)

  IF (ZSTORAGE .GT. ZXM) THEN
! When storage exceeds the mortality of active biomass, an alternative
! formulation of B-decline is employed  

!   Active biomass mortality = structural storage, so Mb is zero
    ZMORTBSTR=PBIOMASSTR_LAST(JL)*(1.0_JPRB-EXP(-1.0_JPRB*RDAY/RVSEFOLD(IVT)))
!4
    ZBIOMASSTOT=PBIOMASS_LAST(JL)+PBIOMASSTR_LAST(JL)+ZASSIM-ZMORTBSTR-PRESPBSTR(JL)
    PBIOMASS(JL)=ZCC_NITRO*(ZBIOMASSTOT**(1.0_JPRB-ZCA_NITRO))
    ZBIOMASSTR=ZBIOMASSTOT-PBIOMASS(JL)
    ZXM=PBIOMASS_LAST(JL)+ZASSIM-PBIOMASS(JL)
    ZSTORAGE=ZBIOMASSTR-PBIOMASSTR_LAST(JL)+ZMORTBSTR+PRESPBSTR(JL)  

  ENDIF

  ZMORTBACT=MAX(0.0_JPRB,ZXM-ZSTORAGE)
  
  IF (ZBIOMASSTOT .LE. ZBIOMASSTOTLIM .AND. ZSEFOLD .GT. 1.0_JPRB) THEN
! emergency deep structural biomass  
    ZBIOMASSTR2=ZBIOMASSTR2+ZMORTBACT
  ENDIF
  
! Deep structural decline  
  ZBIOMASSTR2=ZBIOMASSTR2-ZMORTBSTR2-PRESPBSTR2(JL)  
  
! Avoid negative values of biomass
! No test on ZMORTBSTR2 as it is not used after
! No test on PRESPBSTR2 as it should be smaller than PBIOMASSTR2_LAST,
! otherwise there are irrealistic values of temperature 
  ZBIOMASSTR2=MAX(ZBIOMASSTR2,0.0_JPRB)  
!  ZBIOMASSTR=MAX(ZBIOMASSTR,0.0_JPRB)  

! re-initialisation of biomass compartments values: X(day) <-- X(day-1)
  PBIOMASS_LAST(JL)=PBIOMASS(JL)
  PBIOMASSTR_LAST(JL)=ZBIOMASSTR
  PBIOMASSTR2_LAST(JL)=ZBIOMASSTR2

! re-initialisation of respiration, mortality and assimilation terms
  PRESPBSTR(JL)=0.0_JPRB
  PRESPBSTR2(JL)=0.0_JPRB
  PANFM(JL)=0.0_JPRB

! loss output
  PBLOSS(JL)=ZXM

  ENDIF !PCV>0

  END DO

!ENDIF ! midnight

! respiration of above-ground and deep structural biomass
! calculated at every time step and accumulated over 1 day

ZTENTHLOG2=LOG(2.0_JPRB)*0.1_JPRB
DO JL=KIDIA,KFDIA
 IF (PCV(JL) > 0._JPRB) THEN
   PRESPBSTR(JL)=PRESPBSTR(JL)+PBIOMASSTR_LAST(JL)*RRESPFACTOR_NIT &
    *EXP(ZTENTHLOG2*(PTSKM1M(JL)-RTT-25._JPRB))*PTSPHY
   PRESPBSTR2(JL)=PRESPBSTR2(JL)+PBIOMASSTR2_LAST(JL)*RRESPFACTOR_NIT &
    *EXP(ZTENTHLOG2*(PTSOIL(JL)-RTT-25._JPRB))*PTSPHY
  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('NITRO_DECLINE_MOD:NITRO_DECLINE',1,ZHOOK_HANDLE)
END SUBROUTINE NITRO_DECLINE
END MODULE NITRO_DECLINE_MOD

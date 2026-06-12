MODULE COTWORESTRESS_MOD
CONTAINS
SUBROUTINE COTWORESTRESS(KIDIA,KFDIA,KLON,KVTYPE,KCO2TYP,PFRTI,&
     & PTM1, PQM1,PCM1, PAPHM1,&
     & PTSKM1M, PTSOIL,&
     & PEVAP, PLAI,  &
     & PSRFD ,PRAQ, PMU0, &
     & PF2, PQS,  &
     & YDCST, YDAGS, YDAGF, YDVEG, YDFLAKE,  &
     & PAN, PAG,PRD, &
     & PWET,PDSP,PDMAXT)

! (C) Copyright 1997- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
  !**   *COTWORESTRESS* - CALCULATES NET ASSIMILATION OF CO2 AND CANOPY CONDUCTANCE  

  !     A. Boone       * Meteo-France *     27/10/97 
  !     (following Belair)
  !     MODIFIED BY
  !     V. Masson and V. Rivailland            12/03 
  !     modification of ISBA routines order
  !     M.H. Voogt (KNMI) "C-Tessel"  09/2005 
  !     S. Lafont (ECMWF) C-TESSEL 18/05
  !     G. Balsamo (ECMWF) 24/3/2014 cleaning and LDLAND protection
  !     A. Agusti-Panareda Nov 2020  couple atm CO2 tracer (new input) with photosynthesis 
  !     V.Bastrikov,F.Maignan,P.Peylin,A.Agusti-Panareda/S.Boussetta Feb 2021 Add Farquhar photosynthesis model
   
  !     PURPOSE
  !     -------
  !     Calculates net assimilation of CO2 and leaf/canopy conductance.

  !     INTERFACE
  !     ---------
  !     COTWORESTRESS IS CALLED BY *VSURF_MOD* 

  !     PARAMETER     DESCRIPTION                                   UNITS
  !     ---------     -----------                                   -----
  !     INPUT PARAMETERS (INTEGER):

  !     *KVTYPE*       VEGETATION TYPE CORRESPONDING TO TILE 

  !     INPUT PARAMETERS (REAL)
  !     *PFRTI*      TILE FRACTIONS                                   (0-1)
  !            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
  !            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
  !            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
  !            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
  !     *PTM1*         TEMPERATURE AT T-1                            K
  !     *PQM1*         SPECIFIC HUMIDITY AT T-1                      KG/KG 
  !     *PCM1*         ATMOSPHERIC CO2 AT T-1                      KG/KG 
  !     *PAPHM1*       PRESSURE AT T-1				   PA
  !     *PTSKM1M*      SURFACE TEMPERATURE                           K
  !     *PTSOIL*       SOIL TEMPERATURE LEVEL 3 (28 - 100cm)         K
  !     *PEVAP*        PRELIMINARY MOISTURE FLUX                     KG/M2/S
  !     *PLAI*         LEAF AREA INDEX                               M2/M2
  !     *PSRFD*        DOWNWARD SHORT WAVE RADIATION FLUX AT SURFACE W/M2
  !     *PRAQ*         PRELIMINARY AERODYNAMIC RESISTANCE            S/M
  !     *PMU0*        LOCAL COSINE OF INSTANTANEOUS MEAN SOLAR ZENITH ANGLE
  !     *PF2*	     SOIL MOISTURE STRESS FUNCTION 	           -
  !     *PQS*          SATURATION Q AT SURFACE			   KG/KG

  !     OUTPUT PARAMETERS (REAL):

  !     *PAN*          NET CO2 ASSIMILATION OVER CANOPY          KG_CO2/M2/S
  !                    positive downwards, to be changed for diagnostic output
  !     *PAG*          GROSS CO2 ASSIMILATION OVER CANOPY        KG_CO2/M2/S
  !                    positive downwards, to be changed for diagnostic output
  !     *PRD*          DARK RESPIRATION                          KG_CO2/M2/S
  !                    positive upwards
  !     *PWET*         "BULK" STOMATAL RESISTANCE = canopy resistance  S/M  
  !     *PDSP*         specific humidity deficit for PDHVEGS
  !     *PDMAXT*       maximum specific humidity deficit for PDHVEGS

  !     METHOD
  !     ------
  !     Calvet et al. 1998 Forr. Agri. Met. [from model of Jacobs(1994)]

  !     REFERENCE
  !     ---------
  !     Calvet et al. 1998 Forr. Agri. Met. 

  !     ------------------------------------------------------------------------

  USE PARKIND1   ,ONLY : JPIM, JPRB
  USE YOMHOOK    ,ONLY : LHOOK, DR_HOOK, JPHOOK
  USE YOS_CST	 ,ONLY : TCST
  USE YOS_AGS    ,ONLY : TAGS
  USE YOS_AGF    ,ONLY : TAGF
  USE YOS_VEG    ,ONLY : TVEG
  USE YOS_FLAKE  ,ONLY : TFLAKE
  USE CCETR_MOD
  USE COTWO_MOD
  USE FARQUHAR_MOD
  USE ABORT_SURF_MOD

  IMPLICIT NONE
  INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
  INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
  INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
  INTEGER(KIND=JPIM),INTENT(IN)    :: KVTYPE(:)
  INTEGER(KIND=JPIM),INTENT(IN)    :: KCO2TYP(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:) 
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM1(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSOIL(:) 
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAP(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAI(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PSRFD(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAQ(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PF2(:)
  REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(:)
  TYPE(TCST)        ,INTENT(IN)    :: YDCST
  TYPE(TAGS)        ,INTENT(IN)    :: YDAGS
  TYPE(TAGF)        ,INTENT(IN)    :: YDAGF
  TYPE(TVEG)        ,INTENT(IN)    :: YDVEG
  TYPE(TFLAKE)      ,INTENT(IN)    :: YDFLAKE
  REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAN(:)
  REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAG(:)
  REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRD(:)
  REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWET(:)
  REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDSP(:)
  REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMAXT(:)

  !*         0.     LOCAL VARIABLES.
  !                 ----- ----------

  INTEGER(KIND=JPIM) :: JL,JINT   ! index for loops

  REAL(KIND=JPRB)    :: ZTSK(KLON), ZIA(KLON), ZGAMMT(KLON), ZRHO(KLON), ZQSURF(KLON), &
       & ZDSP(KLON), ZCO2(KLON), ZMU0(KLON), ZAMMAX(KLON), ZGMEST(KLON), &
       & ZTAN(KLON), ZTAG(KLON), ZTRD(KLON), ZTGS(KLON), ZXIA(KLON), &
       & ZAN0(KLON), ZAG0(KLON), ZRD0(KLON), ZGS0(KLON), ZXTGS(KLON), &
       & ZDMAXSTAR(KLON), ZFZEROT(KLON), ZDMAXT(KLON), ZFZEROSTAR(KLON), &
       & ZFZERON(KLON), ZGMESTN(KLON),ZRVGC(KLON),ZRVEPSO(KLON),Z1(KLON)  
  LOGICAL            :: LDLAND(KLON)

! VEGETATION TYPES ARE:
! 1  DECIDUOUS
! 2  CONIFEROUS
! 3  EVERGREEN
! 4  C3 GRASS
! 5  C4 GRASS
! 6  C3 CROPS
! 7  C4 CROPS 

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

  !  ZTSK    = surface/skin temperature (degrees C) 
  !  ZIA     = absorbed PAR at the top of the canopy (W m-2) 
  !  ZGAMMT  = CO2 compensation point at T=Tskin (kgCO2 kgAir-1)
  !  ZRHO    = air density (kg m-3)
  !  ZQSURF  = specific humidity near the surface
  !  ZDSP    = corrected specific humidity deficit (kg kg-1)
  !  ZCO2    = corrected atmospheric concentration of CO2 at canopy level (kgCO2 kgAir-1)
  !  ZMU0   = cosine of solar zenith angle
  !  ZAMMAX  = leaf photosynthetic capacity at T=Tskin (kgCO2 kgAir-1 m s-1)
  !  ZGMEST  = mesophyll conductance at T=Tskin corrected for moisture stress (m s-1)
  !  ZTAN    = sum for integrated net assimilation (kgCO2 kgAir-1 m s-1)
  !  ZTAG    = sum for integrated gross assimilation (kgCO2 kgAir-1 m s-1)
  !  ZTRD    = sum for integrated dark respiration (kgCO2 kgAir-1 m s-1) 
  !  ZTGS    = sum for integrated leaf conductance (m s-1)
  !  ZXIA    = absorbed PAR at different Gaussian levels in the canopy (W m-2) 
  !  ZAN0    = net assimilation at each interval in the canopy (kgCO2 kgAir-1 m s-1)
  !  ZAG0    = gross assimilation at each interval in the canopy (kgCO2 kgAir-1 m s-1)
  !  ZRD0    = dark respiration at each interval in the canopy (kgCO2 kgAir-1 m s-1) 
  !  ZGS0    = leaf conductance at each interval in the canopy (m s-1)
  !  ZXTGS   = total canopy conductance (m s-1)    

  !  ZFZEROT = fzero: ideal value of f, no photorespiration or specific humidity deficit
  !            (-)
  !  ZDMAXT  = maximum specific humidity deficit of atmosphere tolerated by 
  !            vegetation (kg kg-1)
  !  ZDMAXSTAR = maximum specific humidity deficit of atmosphere tolerated
  !              by vegetation without soil water stress (kg kg-1)
  !  ZFZEROSTAR= fzero without soil water stress for woody vegetation
  !              (-)
  !  ZFZERON = minimum value for fzero in defensive woody strategy
  !  ZGMESTN = gmest value at pf2=rvf2i in offensive woody strategy

  !     -------------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('COTWORESTRESS_MOD:COTWORESTRESS',0,ZHOOK_HANDLE)


ASSOCIATE(LVSTRESS=>YDAGS%LVSTRESS, RABC=>YDAGS%RABC, RASW=>YDAGS%RASW, &
 & RAW=>YDAGS%RAW, RBW=>YDAGS%RBW, RCO2=>YDAGS%RCO2, &
 & RCONDCTMIN=>YDAGS%RCONDCTMIN, RDMAXN=>YDAGS%RDMAXN, RDMAXX=>YDAGS%RDMAXX, &
 & RPARCF=>YDAGS%RPARCF, RPOI=>YDAGS%RPOI, RRACCF=>YDAGS%RRACCF, &
 & RVAH=>YDAGS%RVAH, RVAMMAX=>YDAGS%RVAMMAX, RVBH=>YDAGS%RVBH, &
 & RVDMAX=>YDAGS%RVDMAX, RVEPSO=>YDAGS%RVEPSO, RVF2I=>YDAGS%RVF2I, &
 & RVFZEROST=>YDAGS%RVFZEROST, RVGAMM=>YDAGS%RVGAMM, RVGC=>YDAGS%RVGC, &
 & RVGMES=>YDAGS%RVGMES, RVQDAMMAX=>YDAGS%RVQDAMMAX, RVQDGAMM=>YDAGS%RVQDGAMM, &
 & RVQDGMES=>YDAGS%RVQDGMES, RVT1AMMAX=>YDAGS%RVT1AMMAX, &
 & RVT1GMES=>YDAGS%RVT1GMES, RVT2AMMAX=>YDAGS%RVT2AMMAX, &
 & RVT2GMES=>YDAGS%RVT2GMES, &
 & RD=>YDCST%RD, RETV=>YDCST%RETV, RTT=>YDCST%RTT, &
 & LEFLAKE=>YDFLAKE%LEFLAKE,LEAIRCO2COUP=>YDVEG%LEAIRCO2COUP,&
 & LEFARQUHAR=>YDVEG%LEFARQUHAR)

IF (LEFLAKE) THEN
  DO JL=KIDIA,KFDIA
    LDLAND(JL)=((PFRTI(JL,1)+PFRTI(JL,2)+PFRTI(JL,9))<=0.5_JPRB)
  ENDDO
ELSE
  DO JL=KIDIA,KFDIA
    LDLAND(JL)=((PFRTI(JL,1)+PFRTI(JL,2))<=0.5_JPRB)
  ENDDO
ENDIF
ZIA(:)=0.0_JPRB  ! To secure unconditioned use in CCETR

! initialization of local variables
  DO JL=KIDIA,KFDIA
    IF (LDLAND(JL)) THEN

     ZGAMMT(JL)=0.0_JPRB
     ZRHO(JL)=0.0_JPRB
     ZQSURF(JL)=0.0_JPRB
     ZDSP(JL)=0.0_JPRB
     ZMU0(JL)=0.0_JPRB
     ZAMMAX(JL)=0.0_JPRB
     ZGMEST(JL)=0.0_JPRB
     ZTAN(JL)=0.0_JPRB
     ZTAG(JL)=0.0_JPRB
     ZTRD(JL)=0.0_JPRB
     ZTGS(JL)=0.0_JPRB
     ZXIA(JL)=0.0_JPRB
     ZAN0(JL)=0.0_JPRB
     ZAG0(JL)=0.0_JPRB
     ZRD0(JL)=0.0_JPRB
     ZGS0(JL)=0.0_JPRB
     ZXTGS(JL)=0.0_JPRB
     ZDMAXSTAR(JL)=0.0_JPRB
     ZFZEROT(JL)=0.0_JPRB
     ZDMAXT(JL)=0.0_JPRB
     ZFZEROSTAR(JL)=0.0_JPRB
     ZFZERON(JL)=0.0_JPRB
     ZGMESTN(JL)=0.0_JPRB
     ZRVGC(JL)=0.0_JPRB
     ZRVEPSO(JL)=0.0_JPRB

 ! Skin temperature : convert from K to C
     ! Skin temperature : convert from K to C
     ZTSK(JL)=PTSKM1M(JL)-RTT  
     Z1(JL)=0.1_JPRB*(ZTSK(JL)-25._JPRB)
     ! Absorbed PAR
     ZIA(JL)=PSRFD(JL)*RPARCF               
  !-------------------------------------
  ! moisture stress response (offensive/defensive strategy) for woody (trees) and 
  ! herbaceaus (grass,crops) vegetation: 
  ! Add soil moisture stress effect to leaf conductance:
     IF (KVTYPE(JL)>= 0._JPRB) THEN
       ZGMEST(JL)=RVGMES(KVTYPE(JL))
       ZFZEROT(JL)=RVFZEROST(KVTYPE(JL))
       ZDMAXT(JL)=RVDMAX(KVTYPE(JL))
       ZRVGC(JL)=RVGC(KVTYPE(JL))   
       ZRVEPSO(JL)=RVEPSO(KVTYPE(JL))
     ELSE
       ZGMEST(JL)=0._JPRB
       ZFZEROT(JL)=0._JPRB
       ZDMAXT(JL)=0._JPRB
       ZRVGC(JL)=0._JPRB
       ZRVEPSO(JL)=0._JPRB
     ENDIF
     IF (KVTYPE(JL) == 1 .OR. KVTYPE(JL) == 2 .OR. KVTYPE(JL) == 7 .OR. KVTYPE(JL) == 9 .OR. &
       & KVTYPE(JL) == 10 .OR. KVTYPE(JL) == 11 .OR. KVTYPE(JL) == 13 .OR. KVTYPE(JL) == 16 .OR. &
       & KVTYPE(JL) == 17 .OR. KVTYPE(JL) == 20 .OR. KVTYPE(JL) == 0) THEN  !BATS classification matching
 
       ZDMAXSTAR(JL)=EXP((LOG(ZGMEST(JL)*1000._JPRB)-RVAH(KVTYPE(JL)))/RVBH(KVTYPE(JL)))/1000._JPRB
    ! compensation point: temperature response
       ZGAMMT(JL)=RVGAMM(KVTYPE(JL))*EXP(Z1(JL)*LOG(RVQDGAMM(KVTYPE(JL))))

        ! defensive soil water stress response AND soil moist above critical value 
       IF (LVSTRESS(KVTYPE(JL)) .AND. (PF2(JL) >= RVF2I(KVTYPE(JL)))) THEN
           ZDMAXT(JL)=RDMAXN+(ZDMAXSTAR(JL)-RDMAXN)*(PF2(JL)-RVF2I(KVTYPE(JL)))/(1._JPRB-RVF2I(KVTYPE(JL)))
           ZGMEST(JL)=EXP(RVAH(KVTYPE(JL))+RVBH(KVTYPE(JL))*LOG(ZDMAXT(JL)*1000._JPRB))/1000._JPRB

           ! defensive soil water stress response AND soil moist below critical value 
       ELSEIF (LVSTRESS(KVTYPE(JL)) .AND. (PF2(JL) < RVF2I(KVTYPE(JL)))) THEN
           ZDMAXT(JL)=RDMAXN
           ZGMEST(JL)=(EXP(RVAH(KVTYPE(JL))+RVBH(KVTYPE(JL))*LOG(RDMAXN*1000._JPRB))/1000._JPRB)*&
                & PF2(JL)/RVF2I(KVTYPE(JL))

           ! offensive soil water stress response AND soil moist above critical value 
       ELSEIF (.NOT. LVSTRESS(KVTYPE(JL)) .AND. (PF2(JL) >= RVF2I(KVTYPE(JL)))) THEN
           ZDMAXT(JL)=RDMAXX+(ZDMAXSTAR(JL)-RDMAXX)*(PF2(JL)-RVF2I(KVTYPE(JL)))/(1._JPRB-RVF2I(KVTYPE(JL)))
           ZGMEST(JL)=EXP(RVAH(KVTYPE(JL))+RVBH(KVTYPE(JL))*LOG(ZDMAXT(JL)*1000._JPRB))/1000._JPRB

           ! offensive soil water stress response AND soil moist below critical value
       ELSEIF (.NOT. LVSTRESS(KVTYPE(JL)) .AND. (PF2(JL) < RVF2I(KVTYPE(JL)))) THEN
           ZDMAXT(JL)=RDMAXX*PF2(JL)/RVF2I(KVTYPE(JL))
           ZGMEST(JL)=EXP(RVAH(KVTYPE(JL))+RVBH(KVTYPE(JL))*LOG(RDMAXX*1000._JPRB))/1000._JPRB
       ENDIF

       ! Modification Seb
       !to  limit photosynthesis below wilting point 
       IF (.NOT. LVSTRESS(KVTYPE(JL)) .AND. ((ZDMAXT(JL) <= RDMAXN)) ) THEN

            ZDMAXT(JL) = RDMAXN
            ZGMEST(JL) = (EXP(RVAH(KVTYPE(JL))+RVBH(KVTYPE(JL))*LOG(RDMAXN*1000._JPRB))/1000._JPRB)*PF2(JL)/RVF2I(KVTYPE(JL))
       ENDIF
     ELSEIF (KVTYPE(JL) == 3 .OR. KVTYPE(JL) == 4 .OR. KVTYPE(JL) == 5 .OR. KVTYPE(JL) == 6 .OR. &
           & KVTYPE(JL) == 18 .OR. KVTYPE(JL) == 19) THEN  !BATS classification matching
     ! compensation point: temperature response
       ZGAMMT(JL)=RVGAMM(KVTYPE(JL))*EXP(Z1(JL)*LOG(RVQDGAMM(KVTYPE(JL))))
        ZFZEROSTAR(JL)=(RAW-LOG(ZGMEST(JL)*1000._JPRB))/RBW
        ZFZERON(JL)=(RASW-LOG(ZGMEST(JL)*1000._JPRB))/RBW
        ZGMESTN(JL)=EXP(RASW-RBW*ZFZEROSTAR(JL))/1000._JPRB

        ! defensive soil water stress response AND soil moist above critical value 
       IF (LVSTRESS(KVTYPE(JL)) .AND. (PF2(JL) >= RVF2I(KVTYPE(JL)))) THEN
           ZFZEROT(JL)=ZFZERON(JL)+(ZFZEROSTAR(JL)-ZFZERON(JL))*(PF2(JL)-RVF2I(KVTYPE(JL)))/&
                &(1._JPRB-RVF2I(KVTYPE(JL)))
!           ZGMEST(JL)=EXP(RASW-RBW*ZFZEROT(JL))

           ! defensive soil water stress response AND soil moist below critical value 
       ELSEIF (LVSTRESS(KVTYPE(JL)) .AND. (PF2(JL) < RVF2I(KVTYPE(JL)))) THEN
           ZGMEST(JL)=MAX(1.0E-6_JPRB,ZGMEST(JL)*PF2(JL)/RVF2I(KVTYPE(JL)))
           ZFZEROT(JL)=MIN(.95_JPRB,(RASW-LOG(ZGMEST(JL)*1000._JPRB))/RBW)

           ! offensive soil water stress response AND soil moist above critical value 
       ELSEIF (.NOT. LVSTRESS(KVTYPE(JL)) .AND. (PF2(JL) >= RVF2I(KVTYPE(JL)))) THEN    
           ZGMEST(JL)=ZGMESTN(JL)+(ZGMEST(JL)-ZGMESTN(JL))*(PF2(JL)-RVF2I(KVTYPE(JL)))/(1._JPRB-RVF2I(KVTYPE(JL)))
           ZFZEROT(JL)=ZFZEROSTAR(JL)

           ! offensive soil water stress response AND soil moist below critical value
       ELSEIF (.NOT. LVSTRESS(KVTYPE(JL)) .AND. (PF2(JL) < RVF2I(KVTYPE(JL)))) THEN    
           ZGMEST(JL)=ZGMESTN(JL)*PF2(JL)/RVF2I(KVTYPE(JL))
           ZFZEROT(JL)=MIN(.95_JPRB,(RASW-LOG(ZGMEST(JL)*1000._JPRB))/RBW)
       ENDIF
     ENDIF

     PDMAXT(JL)=ZDMAXT(JL)
     !-------------------------

     ! specific humidity deficit (kg kg-1) corrected to canopy level
     ZRHO(JL)=PAPHM1(JL)/( RD*PTM1(JL)*(1._JPRB+RETV*PQM1(JL)))
     ! specific humidity at the surface (corrected with evaporation flux at the surface)
     ZQSURF(JL)=PQM1(JL)-MIN(0._JPRB,PEVAP(JL))*PRAQ(JL)/ZRHO(JL)
     !ZDSP=MAX(0._JPRB,PQS-PQM1-PEVAP*PRAQ/ZRHO)
     ! bug correction S.lafont  PEVAP is <0 so it is +PEVAP 
     !PDSP(JL)=MAX(0._JPRB,PQS(JL)-PQM1(JL)+PEVAP(JL)*PRAQ(JL)/ZRHO(JL))
     ! Limit PEVAP to evaporation only - avoid deposition, very small impact but with failures in single precision 
     PDSP(JL)=MAX(0._JPRB,PQS(JL)-PQM1(JL)+MIN(0._JPRB,PEVAP(JL))*PRAQ(JL)/ZRHO(JL))     


     ! Cosine of solar zenith angle 
     ! the minimum value of 0.1 is chosen, in order to prevent a division by zero
     ! in CCETR. (SWdown is interpolated, so it can have values higher than zero, 
     ! while the sun has not risen yet). 
     ZMU0(JL)=MAX(PMU0(JL),0.1_JPRB)

     ! Compute temperature response functions for Am,max and gmes:
     IF (KVTYPE(JL)>=0._JPRB) THEN
       ZAMMAX(JL)=(RVAMMAX(KVTYPE(JL))*EXP(Z1(JL)*LOG(RVQDAMMAX(KVTYPE(JL))))) &
          & /((1._JPRB+EXP(0.3_JPRB*(RVT1AMMAX(KVTYPE(JL))-ZTSK(JL))))* &
          & (1._JPRB+EXP(0.3_JPRB*(ZTSK(JL)-RVT2AMMAX(KVTYPE(JL))))))
       ZGMEST(JL)=(ZGMEST(JL)*EXP(Z1(JL)*LOG(RVQDGMES(KVTYPE(JL))))) &
          & /((1._JPRB+EXP(0.3_JPRB*(RVT1GMES(KVTYPE(JL))-ZTSK(JL))))* &
          & (1._JPRB+EXP(0.3_JPRB*(ZTSK(JL)-RVT2GMES(KVTYPE(JL))))))
     ENDIF
    ENDIF !LDLAND

    ! Atmospheric concentration for CO2, not taking into account the aerodynamic
    ! resistance for CO2, since we work with a respiration calibration  
!old    IF (LEAIRCO2COUP .AND. PCM1(JL) > 0.0_JPRB) THEN
    IF (PCM1(JL) > 0.0_JPRB) THEN
    ! Variable atmospheric CO2
       ZCO2(JL)=PCM1(JL)
    ELSE 
       CALL ABORT_SURF('COTWORESTRESS_MOD: ATM CO2 SHOULD BE POSITIVE ')
    ENDIF


!old    ! Constanst CO2
!old       ZCO2(JL)=RCO2 
!old !        write (*,*) 'COTWORESTRESS. Farquhar_mode:',FARQUHAR_MODE,'Fixed CO2' 
!old    ! originally it was:
!old    !     ZCO2=RCO2+PRAQ*RRACCF*PCO2FLUX/ZRHO
!old    ENDIF

  ! Integration over the canopy: 
  ! SIZE(RABC) increments are used to approximate the integral.
     ZTAN(JL)=0._JPRB
     ZTAG(JL)=0._JPRB
     ZTRD(JL)=0._JPRB
     ZTGS(JL)=0._JPRB
  ENDDO

IF (LEFARQUHAR) THEN 

  ! Farquhar photosynthesis model
  !
  CALL FARQUHAR(KIDIA, KFDIA, KLON, LDLAND, KVTYPE, KCO2TYP, YDAGS, YDAGF, &
 &               YDCST, PSRFD, &
 &               PAPHM1, PQS, ZQSURF, PQM1, PTM1, PTSKM1M, PTSOIL, PF2, &
 &               ZCO2, PLAI, ZXTGS, PAG, PRD, PAN)

  ! Net assimilation, gross assimilation and dark respiration over canopy
  ! (kgCO2 m-2 s-1) An,Ag positive downwards (to be changed in co2.F90 for
  ! diagnostic output), Rd positive upwards
  DO JL=KIDIA,KFDIA
     IF (LDLAND(JL) .AND. KVTYPE(JL)>0._JPRB) THEN
        ! Canopy resistance from Ags (s m-1)
        PWET(JL)=1._JPRB/MAX(RCONDCTMIN,ZXTGS(JL))
     ELSE
        PWET(JL) =0._JPRB
     ENDIF
  ENDDO

ELSE

  ! A-gs photosynthesis model

  DO JINT=1,SIZE(RABC)
     ! jint=1 is the lowest part of the canopy, jint=size(rabc) is the highest part. 
     !  Diffusion of incident radiation:

     CALL CCETR(KIDIA,KFDIA,KLON,KVTYPE,LDLAND,ZIA,ZMU0,RABC(JINT),PLAI,YDAGS,ZXIA)

     !  Compute conductance and assimilation of CO2: 
     CALL COTWO(KIDIA,KFDIA,KLON,LDLAND,ZAN0,ZAG0,ZRD0,ZGS0,ZRVGC,ZCO2, &
          & PDSP,ZDMAXT,ZXIA,ZGAMMT,YDAGS, &
          & ZFZEROT,ZGMEST,ZRVEPSO,ZAMMAX)

     DO JL=KIDIA,KFDIA
      IF (LDLAND(JL)) THEN
        ZTAN(JL)=ZTAN(JL)+ZAN0(JL)*RPOI(JINT) 
        ZTAG(JL)=ZTAG(JL)+ZAG0(JL)*RPOI(JINT) 
        ZTRD(JL)=ZTRD(JL)+ZRD0(JL)*RPOI(JINT) 
        ZTGS(JL)=ZTGS(JL)+ZGS0(JL)*RPOI(JINT) 
      ENDIF
     ENDDO
  ENDDO

  ! Net assimilation, gross assimilation and dark respiration over canopy
  ! (kgCO2 m-2 s-1) An,Ag positive downwards (to be changed in co2.F90 for
  ! diagnostic output), Rd positive upwards
  DO JL=KIDIA,KFDIA
     IF (LDLAND(JL) .AND. KVTYPE(JL)>0._JPRB) THEN
        PAN(JL)=ZTAN(JL)*PLAI(JL)*ZRHO(JL)
        PAG(JL)=ZTAG(JL)*PLAI(JL)*ZRHO(JL) 
        PRD(JL)=ZTRD(JL)*PLAI(JL)*ZRHO(JL)
        
        ! Total conductance over canopy (m s-1)
        ZXTGS(JL)=ZTGS(JL)*PLAI(JL)

        ! Canopy resistance from Ags (s m-1)
        PWET(JL)=1._JPRB/MAX(RCONDCTMIN,ZXTGS(JL))
     ELSE
        PAN(JL)  =0._JPRB
        PAG(JL)  =0._JPRB
        PRD(JL)  =0._JPRB
        PWET(JL) =0._JPRB
     ENDIF
  ENDDO


   
ENDIF !photosynthesis model

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('COTWORESTRESS_MOD:COTWORESTRESS',1,ZHOOK_HANDLE)
END SUBROUTINE COTWORESTRESS
END MODULE COTWORESTRESS_MOD

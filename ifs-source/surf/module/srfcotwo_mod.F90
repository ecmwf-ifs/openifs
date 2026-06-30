MODULE SRFCOTWO_MOD
CONTAINS
SUBROUTINE SRFCOTWO(KIDIA,KFDIA,KLON,KLEVS,KTILES,&
 & KVEG,KVTTL,KCO2TYP,KSOTY,&
 & PTSTEP, &
 & PTMLEV,PQMLEV,PAPHMS,&
 & PCVT,PFRTI,PLAIVT,PWLIQ ,&
 & PWSOIL,PTSOIL,&
 & PDSN,PFWET,PLAT, &
 & PANTI,PAGTI,PRDTI,&
 & YDCST,YDVEG,YDSOIL,YDAGS,&
 & PANDAYVT,PANFMVT,&
 & PAG,PRD,PAN,PRSOIL_STR,PRECO,PCO2FLUX,PCH4FLUX,&
 & PDHCO2S)

 ! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

 !     ------------------------------------------------------------------

!**   *SRFCOTWO* - DOES THE POST-PROCESSING OF CO2 VALUES

!     Marita Voogt (KNMI)         16/09/2005
!     Sebastien Lafont (ECMWF)    18/05/2006
!     S. Boussetta/G.Balsamo June 2010 Add soil moisture scaling factor for Reco
!     S. Boussetta/G.Balsamo October 2010 Add snowpack and cold season temeprature effects on respiration

!     PURPOSE
!     -------

!     POST-PROCESSING OF THE CO2 VALUES FROM VDFSURF

!     INTERFACE
!     ---------

!     *SRFCOTWO* IS CALLED BY *SURFEXCDRIVER*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KTILES*       NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                    OF SURFACE BOUNDARY CONDITION)
!     *KVTTL*        CROSS-REFERENCE BETWEEN TILES AND VEGETATION TYPES 
!     *KCO2TYP*      TYPE OF PHOTOSYNTHETIC PATHWAY FOR LOW VEGETATION(C3/C4) (INDEX 1/2)

!     INPUT PARAMETERS (REAL):

!     *PTSTEP*       TIMESTEP
!     *PTLEV*         TEMPERATURE AT T-1                            K
!     *PQLEV*         SPECIFIC HUMIDITY AT T-1                      KG/KG 
!     *PAPHMS*       PRESSURE AT T-1				   PA
!     *PCVT*         VEGETATION TYPE FRACTION                    (0-1)
!            1: LOW  : 0-1
!            2: HIGH : 0-1
!            3: needed for cross-reference with non-vegetation tiles: 0

!     *PFRTI*        TILE FRACTION                              (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!     *PLAIVT*       LEAF AREA INDEX                                (-)
!     *PWSOIL*       SOIL MOISTURE OF LAYER X (SEE CALL FROM VDFMAIN)
!     *PTSOIL*       SOIL TEMPERATURE OF LAYER X (SEE CALL FROM VDFMAIN) 

!     *PDSN*         Total Snow depth (m) 

!     *PANTI*        NET CO2 ASSIMILATION OVER CANOPY          KG_CO2/M2/S
!                    positive downwards, to be changed for diagnostic output
!     *PAGTI*        GROSS CO2 ASSIMILATION OVER CANOPY        KG_CO2/M2/S
!                    positive downwards, to be changed for diagnostic output
!     *PRDTI*        DARK RESPIRATION                          KG_CO2/M2/S
!                    positive upwards


!     UPDATED PARAMETERS (REAL):

!     *PANDAYVT*     DAILY NET CO2 ASSIM.OVER CANOPY PER VEGTYPE   KG_CO2/M2
!     *PANFMVT*      MAXIMUM LEAF ASSIMILATION PER VEGTYPE    KG_CO2/KG_AIR M/S
!     OUTPUT PARAMETERS (REAL):

!     all positive upwards:
!     *PAN*          NET CO2 ASSIMILATION OVER CANOPY          KG_CO2/M2/S
!     *PAG*          GROSS CO2 ASSIMILATION OVER CANOPY        KG_CO2/M2/S
!     *PRD*          DARK RESPIRATION                          KG_CO2/M2/S
!     *PRSOIL_STR*   RESPIRATION FROM SOIL AND STRUCTURAL BIOMASS KG_CO2/M2/S
!     *PRECO*        ECOSYSTEM RESPIRATION                     KG_CO2/M2/S

!     *PCO2FLUX*     CO2 FLUX                                  KG_CO2/M2/S
!     *PCH4FLUX*     CH4 FLUX                                  KG_CO2/M2/S

!     *PDHCO2S*      Diagnostic array for CO2 (see module yomcdh) 
!                      (kgCO2 m-2 s-1 for fluxes)


!     METHOD
!     ------

!     TO BE DONE

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_CST  , ONLY : TCST
USE YOS_VEG  , ONLY : TVEG
USE YOS_SOIL , ONLY : TSOIL
USE YOS_AGS  , ONLY : TAGS


IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KVTTL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVT(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIVT(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSOIL(:,:)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSN(:)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSOIL(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFWET(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PANTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAGTI(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDTI(:,:)

TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TVEG)        ,INTENT(IN)    :: YDVEG
TYPE(TSOIL)       ,INTENT(IN)    :: YDSOIL
TYPE(TAGS)        ,INTENT(IN)    :: YDAGS
INTEGER(KIND=JPIM),INTENT(IN)    :: KCO2TYP(:)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PANDAYVT(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PANFMVT(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAN(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAG(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRD(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSOIL_STR(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRECO(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCO2FLUX(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCH4FLUX(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDHCO2S(:,:,:)

REAL(KIND=JPRB)   ,INTENT(IN) :: PWLIQ(KLON,KLEVS,KTILES)

!*         0.     LOCAL VARIABLES.
!                 ----- ----------
!ZLIQ is passed to compute soil moisture scaling factor in Reco (CO2 routine) CTESSEL

REAL(KIND=JPRB) ::     ZFWCO2(KLON),ZWROOT(KLON)
INTEGER(KIND=JPIM) :: KVEG(KLON,KTILES)
REAL(KIND=JPRB) ::ZRWCAP

INTEGER(KIND=JPIM) :: JL, JTILE, JVT, JS, JK ,IVT, ICTYPE
REAL(KIND=JPRB)    :: ZRHO,ZCVS
REAL(KIND=JPRB)    :: ZAGVT(KLON,YDVEG%NVTILES),ZRDVT(KLON,YDVEG%NVTILES),&
 &                    ZANVT(KLON,YDVEG%NVTILES),ZRSOIL_STRVT(KLON,YDVEG%NVTILES),ZRNOQ10VT(KLON,YDVEG%NVTILES),&
 &                    ZRECOVT(KLON,YDVEG%NVTILES),ZCO2FLUXVT(KLON,YDVEG%NVTILES),ZCH4FLUXVT(KLON,YDVEG%NVTILES)
REAL(KIND=JPRB)    :: ZTSOIL,ZWSOIL,ZRSWC,ZS
REAL(KIND=JPRB)    :: ZFSN(KLON)
REAL(KIND=JPRB)    :: ZOMEGASN(KLON)
REAL(KIND=JPRB)    :: ZFITAUSN(KLON)

REAL(KIND=JPRB)    :: ZBARE(KLON)
REAL(KIND=JPRB)    :: ZRQ10,ZCH4S
REAL(KIND=JPRB)    :: ZEPSR  !limit value for ZRQ10 when temperature is very high (>56deg)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
 

!  ZTSOIL  = 2nd soil layer temperature (degrees C) 
!  ZWSOIL  = 2nd soil layer mooisture (m3 m-3)
!  ZRSWC   = relative soil water content used in respiration parameterization
!  ZS      = soil moist stress function used in respiration parameterization


! It is assumed that wet leaves assimilate CO2 in the same way as dry leaves.
! Stomata are mostly located at the bottom side of leaves, so they are assumed 
! not to be influenced by the interception water. 
! Therefore, to get the grid box averaged Ag, Rd and An, the tiled values must
! not be multiplied by the tile fraction (PFRTI = (1-Ci)(1-Csn)*PCVT) but by 
! (1-Csn)*PCVT. Since there is only one shaded snow tile, this needs special
! treatment too. 
! For convenience, first Ag, Rd and An are derived per vegetation type. From
! those values the grid box averages are calculated. 


!*         1.1     Ag, Rd AND An PER VEGETATION TYPE 
!                  ---------------------------------

! For interactive vegetation, assimilation from the shaded snow tile needs to be
! addresses to the corresponding vegetation type. 
! For "tiled" output validation, it is also more convenient to work with values
! per vegetation type than values per tile.

! The shaded snow tile only contains the dominant high vegetation type. The 
! assimilation of the snow-covered part of the non-dominant high vegetation 
! types is assumed to be zero.

! ASSIMILATION DIAGRAM, valid for Agross, Rdark and Anet NOT normalised to the
! grid square
! sn = snow
! un = uncovered part
! _____________________________________________________________________________
!|               |           |       |             |                           |
!|  Vegetation   | Intercep. | Snow  | Uncov. part | Corrected assimilation    |
!|_______________|___________|_______|_____________|___________________________|
!|               |           |       |             |  	                       |
!| low           |   Aun     | Asn=0 |     Aun     |    (1-Csn)*Aun            |
!|_______________|___________|_______|_____________|___________________________|
!|               |           |       |             |  	                       |
!| dominant high |   Aun     |  Asn  |     Aun     | Csn*Asn+(1-Csn)*Aun       |
!|               |           |       |             |                           |
!|               |           |       |             |                           |
!|_______________|___________|_______|_____________|___________________________| 
!|               |           |       |             |  	                       |
!| not dom. high |   Aun     | Asn=0 |     Aun     |    (1-Csn)*Aun            |
!|_______________|___________|_______|_____________|___________________________|

! Asn is the assimilation from the shaded snow tile. Aun is the assimilation 
! from the vegetation tiles.
! Csn can either be derived from the exposed or shaded snow tile. If no
! vegetation is present at all, there is no CO2 assimilation. In that case, Csn 
! doesn't play a role and is set to 0.Here, the
! shaded snow tile is used.
IF (LHOOK) CALL DR_HOOK('SRFCOTWO_MOD:SRFCOTWO',0,ZHOOK_HANDLE)
ASSOCIATE(RQ10=>YDAGS%RQ10, RVR0VT=>YDAGS%RVR0VT, RVCH4QVT=>YDAGS%RVCH4QVT, RVCH4S=>YDAGS%RVCH4S, &
 & RD=>YDCST%RD, RETV=>YDCST%RETV, RTT=>YDCST%RTT, RPI=>YDCST%RPI, &
 & LEVGEN=>YDSOIL%LEVGEN, RWCAP=>YDSOIL%RWCAP, RWCAPM=>YDSOIL%RWCAPM, &
 & RWPWP=>YDSOIL%RWPWP, &
 & NVTILES=>YDVEG%NVTILES, RVROOTSA=>YDVEG%RVROOTSA)

! initialization of local variables

    ZFSN(KIDIA:KFDIA)=1._JPRB
    ZFWCO2(KIDIA:KFDIA)=0._JPRB
    ZWROOT(KIDIA:KFDIA)=0._JPRB
    ZAGVT(KIDIA:KFDIA,:)=0._JPRB
    ZRDVT(KIDIA:KFDIA,:)=0._JPRB
    ZANVT(KIDIA:KFDIA,:)=0._JPRB
    ZRSOIL_STRVT(KIDIA:KFDIA,:)=0._JPRB
    ZRNOQ10VT(KIDIA:KFDIA,:)=0._JPRB
    ZRECOVT(KIDIA:KFDIA,:)=0._JPRB
    ZCO2FLUXVT(KIDIA:KFDIA,:)=0._JPRB
    ZCH4FLUXVT(KIDIA:KFDIA,:)=0._JPRB
    ZBARE(KIDIA:KFDIA)=0._JPRB

    ZRWCAP=0._JPRB
    ZRHO=0._JPRB
    ZCVS=0._JPRB
    ZTSOIL=0._JPRB
    ZWSOIL=0._JPRB
    ZRSWC=0._JPRB
    ZS=0._JPRB
    ZRQ10=2.0_JPRB
    ZCH4S=1.0_JPRB

ZEPSR=0.0001_JPRB !limit value for ZRQ10 when soil temperature is above 56deg 

! End initialization

DO JL=KIDIA,KFDIA
   ZBARE(JL)=1._JPRB-PCVT(JL,KVTTL(4))-PCVT(JL,KVTTL(6))
   ! Csn if there is some vegetation  
   ! high veg
   IF (PCVT(JL,KVTTL(6)) .GT. 0._JPRB) THEN 
      ZCVS=PFRTI(JL,7)/PCVT(JL,KVTTL(6))
   ! low veg
   ELSEIF  (PCVT(JL,KVTTL(4)) .GT. 0._JPRB) THEN
      ZCVS=PFRTI(JL,5)/(PCVT(JL,KVTTL(4))+ZBARE(JL))
   ELSE
      ZCVS=0._JPRB
   ENDIF
     ZCVS=MAX(ZCVS,0._JPRB)
     ZCVS=MIN(ZCVS,1._JPRB)

   ! MIND: values are positive downwards

   ! low vegetation
   IF (PCVT(JL,KVTTL(4)) == 0._JPRB) THEN
      ZAGVT(JL,KVTTL(4))=0._JPRB
      ZRDVT(JL,KVTTL(4))=0._JPRB
      ZANVT(JL,KVTTL(4))=0._JPRB      
   ELSE
      ZAGVT(JL,KVTTL(4))=(1._JPRB-ZCVS)*PAGTI(JL,4)
      ZAGVT(JL,KVTTL(4))=MAX(ZAGVT(JL,KVTTL(4)),0.0_JPRB)
      ZRDVT(JL,KVTTL(4))=(1._JPRB-ZCVS)*PRDTI(JL,4)
      ZANVT(JL,KVTTL(4))=(1._JPRB-ZCVS)*PANTI(JL,4)
   ENDIF

   ! high vegetation
   ! dominant
   ! fraction shaded snow cover  added to fraction high veg
   IF (PCVT(JL,KVTTL(6)) > 0._JPRB) THEN
      ZAGVT(JL,KVTTL(6))=ZCVS*PAGTI(JL,7)+(1._JPRB-ZCVS)*PAGTI(JL,6)
      ZAGVT(JL,KVTTL(6))=MAX(ZAGVT(JL,KVTTL(6)),0.0_JPRB)
      ZRDVT(JL,KVTTL(6))=ZCVS*PRDTI(JL,7)+(1._JPRB-ZCVS)*PRDTI(JL,6)
      ZANVT(JL,KVTTL(6))=ZCVS*PANTI(JL,7)+(1._JPRB-ZCVS)*PANTI(JL,6)
   ELSE 
      ZAGVT(JL,KVTTL(6))=0._JPRB
      ZRDVT(JL,KVTTL(6))=0._JPRB
      ZANVT(JL,KVTTL(6))=0._JPRB  
   ENDIF

ENDDO

!*         1.2     RESPIRATION PER VEGETATION TYPE 
!                  -----------------------------------------

!Introducing the Soil moisture effect on Reco through a scaling factor
!after Alberget et al. 2010 Biogeosciences
! ZFWCO2=wg/wfc
DO JL=KIDIA,KFDIA  
    ZWROOT(JL)=0._JPRB
  DO JTILE=1,KTILES
    DO JK=1,KLEVS
      ZWROOT(JL)=ZWROOT(JL)+PWLIQ(JL,JK,JTILE)*RVROOTSA(JK,KVEG(JL,JTILE))
    ENDDO
  ENDDO
  IF (LEVGEN) THEN
    JS=KSOTY(JL)
    IF (JS.GE.1_JPIM) ZRWCAP=1.0_JPRB/RWCAPM(JS)
  ELSE
      ZRWCAP=1.0_JPRB/RWCAP
  ENDIF
  ZFWCO2(JL)=MIN(1.0_JPRB,(ZWROOT(JL)*ZRWCAP))
ENDDO

! Reco is the sum of Rd and Rsoil_str. Rsoil_str is calculated via a Q10 
! function for each vegetation type in each grid box.
! The reference respiration R0 is calculated via the respiration
! calibration (outside the model code: multi-annual accumulated net assimilation 
! equals Rsoil_str + "harvest").
! The temperature is from the second soil layer.

DO JL=KIDIA,KFDIA

!averaging temperature over the 2 first layers
   ZTSOIL=0._JPRB
    DO JK=1,KLEVS-2_JPIM
       ZTSOIL=ZTSOIL+PTSOIL(JL,JK)
    ENDDO
  ZTSOIL=ZTSOIL/(KLEVS-2_JPIM)
  ZTSOIL=ZTSOIL-RTT
!  ZTSOIL=PTSOIL(JL,2)-RTT

   ZBARE(JL)=1._JPRB-PCVT(JL,KVTTL(4))-PCVT(JL,KVTTL(6))

   IF (PCVT(JL,KVTTL(6)) .GT. 0._JPRB) THEN 
      ZCVS=PFRTI(JL,7)/PCVT(JL,KVTTL(6))
   ! low veg
   ELSEIF  (PCVT(JL,KVTTL(4)) .GT. 0._JPRB) THEN
      ZCVS=PFRTI(JL,5)/(PCVT(JL,KVTTL(4))+ZBARE(JL))
   ELSE
      ZCVS=0._JPRB
   ENDIF
     ZCVS=MAX(ZCVS,0._JPRB)
     ZCVS=MIN(ZCVS,1._JPRB)

!Cold season respiration

!* Respiration attenuation function for Snow pack
!FSN=1-fsn*(1-exp(ws*z))
!z=snow depth, ws=snow attenuation factor=2 10^4*(1-Fi*tau)*Dco2*rho(Mco2,STP)
!Fi=snow porosity = 1-snd/973 , tau=tortuosity=Fi^(1/3) , Dco2= 0.1381 10^-4 m2/s Co2 diffusivity in the air
!rho(Mco2,STP)= 44.613 molecular density of CO2 at STP, T0 standard temperature = 273.16K, 
! ==> ws= 12.322*(1-Fi*tau)

!ZFITAUSN(JL)=(1._JPRB-PDSN(JL)/973._JPRB)**1.3333_JPRB
!ZOMEGASN(JL)=-12.322_JPRB*(1._JPRB-ZFITAUSN(JL))

!ZFSN(JL)=1._JPRB-ZCVS*(1._JPRB-exp(ZOMEGASN(JL)*ZSNDEP(JL)))
!for simplicity ZOMEGASN(JL)=-2
ZFSN(JL)=1._JPRB-ZCVS*(1._JPRB-exp(-2._JPRB*MAX(PDSN(JL),0._JPRB)))


!For respiration Q10 is computed after (McGuire et al., Global Biochemical cycles vol6 n 2, 101-124, June 1992)
! to take into account Q10 temperature dependency esspecially for low values  
ZRQ10=2.5665_JPRB-0.05308_JPRB*ZTSOIL+0.00238_JPRB*ZTSOIL*ZTSOIL-0.00004_JPRB*ZTSOIL*ZTSOIL*ZTSOIL
!ZRQ10=RQ10 previously fixed to 2.
ZRQ10=MAX(ZRQ10,ZEPSR)
  IF ((PLAT(JL)*180._JPRB/RPI) .GT. 25._JPRB) THEN
     ZCH4S=RVCH4S(1)
  ELSEIF ((PLAT(JL)*180._JPRB/RPI) .LT. -25._JPRB) THEN
     ZCH4S=RVCH4S(3)
  ELSE
     ZCH4S=RVCH4S(2)
  ENDIF

  ICTYPE=1
  IF (KCO2TYP(JL).EQ.4) THEN
      ICTYPE = 2
  ENDIF 

! Correct scaling factor to fit climatology
 ZCH4S=ZCH4S*2.0_JPRB

  DO JVT=1,NVTILES
! ! For low temperature season respiration should be lower (Q10 higher) (McDowell et al.Tree Physiology 20, 2000) 

!  correspondancy between tiles and high/low types
     IF  (JVT  ==  1 ) THEN
        IVT=4 !type low veg
     ELSEIF ((JVT .EQ.  2) .AND. (PFRTI(JL,7) .EQ. 0._JPRB)) THEN  ! (KVEG(JL,6)=KVEG(JL,7)=KTVH )
        IVT=6 !type high veg KVEG
     ELSEIF ((JVT .EQ.  2) .AND. (PFRTI(JL,7) .GT. 0._JPRB)) THEN  !(KVEG(JL,6)=KVEG(JL,7)=KTVH )
        IVT=7 !shaded snow same type high veg
     ENDIF

    ZRSOIL_STRVT(JL,JVT)=ZFWCO2(JL)*RVR0VT(KVEG(JL,IVT),ICTYPE)*ZRQ10**(0.1_JPRB*(ZTSOIL-25._JPRB))
    ZRSOIL_STRVT(JL,JVT)=ZFSN(JL)*ZRSOIL_STRVT(JL,JVT)
 
    ZRNOQ10VT(JL,JVT)=ZFWCO2(JL)*RVR0VT(KVEG(JL,IVT),ICTYPE)
    ZRNOQ10VT(JL,JVT)=ZFSN(JL)*ZRNOQ10VT(JL,JVT)
    IF (RVCH4QVT(KVEG(JL,IVT)) .EQ. 0._JPRB .OR. PFWET(JL) .EQ. 0._JPRB .OR. ZRNOQ10VT(JL,JVT) .EQ. 0._JPRB .OR. PTSOIL(JL,1) .LT. 273.15) THEN
        ZCH4FLUXVT(JL,JVT)=0._JPRB
    ELSE
        ZCH4FLUXVT(JL,JVT)=ZCH4S*PFWET(JL)*ZRNOQ10VT(JL,JVT)*RVCH4QVT(KVEG(JL,IVT))**(0.1_JPRB*(PTSOIL(JL,1)-RTT-25._JPRB))
    ENDIF

    ZRECOVT(JL,JVT)=ZRDVT(JL,JVT)+ZRSOIL_STRVT(JL,JVT)
  ENDDO
ENDDO

!*         1.3     NET CO2-FLUX PER VEGETATION TYPE 
!                  --------------------------------

! The net CO2 flux is the sum of Ag and Reco
DO JL=KIDIA,KFDIA
  DO JVT=1,NVTILES
    ZCO2FLUXVT(JL,JVT)=ZAGVT(JL,JVT)*(-1._JPRB)+ZRECOVT(JL,JVT)
  ENDDO
ENDDO


!*         2.1     GRID AVERAGED VALUES 
!                  --------------------

PAG(KIDIA:KFDIA)=0._JPRB
PRD(KIDIA:KFDIA)=0._JPRB
PAN(KIDIA:KFDIA)=0._JPRB
PRSOIL_STR(KIDIA:KFDIA)=0._JPRB
PRECO(KIDIA:KFDIA)=0._JPRB
PCO2FLUX(KIDIA:KFDIA)=0._JPRB
PCH4FLUX(KIDIA:KFDIA)=0._JPRB
DO JL=KIDIA,KFDIA
  DO JVT=1,NVTILES
!sign changed to ECMWF convention (CO2 flux positif downward)
    PAG(JL)=PAG(JL)+PCVT(JL,JVT)*ZAGVT(JL,JVT)
    PRD(JL)=PRD(JL)+PCVT(JL,JVT)*ZRDVT(JL,JVT)*(-1._JPRB)     
    PAN(JL)=PAN(JL)+PCVT(JL,JVT)*ZANVT(JL,JVT)
    PRSOIL_STR(JL)=PRSOIL_STR(JL)+PCVT(JL,JVT)*ZRSOIL_STRVT(JL,JVT)*(-1._JPRB)
    PRECO(JL)=PRECO(JL)+PCVT(JL,JVT)*ZRECOVT(JL,JVT)*(-1._JPRB)
    PCO2FLUX(JL)=PCO2FLUX(JL)+PCVT(JL,JVT)*ZCO2FLUXVT(JL,JVT)*(-1._JPRB)
    PCH4FLUX(JL)=PCH4FLUX(JL)+PCVT(JL,JVT)*ZCH4FLUXVT(JL,JVT)*(-1._JPRB)
  ENDDO

ENDDO



!*         4.1     DIAGNOSTICS
!                  -----------

! Diagnostis in PDH-arrays are NOT normalised to grid box quantities, so values
! are NOT multiplied by the vegetation type fraction (PCVT).
DO JL=KIDIA,KFDIA
  DO JVT=1,NVTILES
    ! not normalised to grid square
!sign changed to ECMWF convention (CO2 flux positif downward)

    PDHCO2S(JL,JVT,1)=ZAGVT(JL,JVT) 
    PDHCO2S(JL,JVT,2)=ZRDVT(JL,JVT)*(-1._JPRB)
    PDHCO2S(JL,JVT,3)=ZANVT(JL,JVT)
    PDHCO2S(JL,JVT,4)=ZRSOIL_STRVT(JL,JVT)*(-1._JPRB)
    PDHCO2S(JL,JVT,5)=ZRECOVT(JL,JVT)*(-1._JPRB)
    PDHCO2S(JL,JVT,6)=ZCO2FLUXVT(JL,JVT)*(-1._JPRB)
!    PDHCO2S(JL,JVT,7)=ZRNOQ10VT(JL,JVT)*(-1._JPRB)
  ENDDO
ENDDO

!*         4.2     INTERACTIVE VEGETATION
!                  ----------------------
DO JL=KIDIA,KFDIA
  ZRHO=PAPHMS(JL)/(RD*PTMLEV(JL)*(1._JPRB+RETV*PQMLEV(JL)))
  DO JVT=1,NVTILES
    PANDAYVT(JL,JVT)=PANDAYVT(JL,JVT)+ZANVT(JL,JVT)*PTSTEP
    IF ((PCVT(JL,JVT).GT. 0._JPRB).AND.(PLAIVT(JL,JVT).GT. 0._JPRB)) THEN
      PANFMVT(JL,JVT)=MAX(PANFMVT(JL,JVT),ZANVT(JL,JVT)/(ZRHO*PLAIVT(JL,JVT)))
    ELSE
      PANFMVT(JL,JVT)=0._JPRB
    ENDIF
  ENDDO
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFCOTWO_MOD:SRFCOTWO',1,ZHOOK_HANDLE)
END SUBROUTINE SRFCOTWO
END MODULE SRFCOTWO_MOD

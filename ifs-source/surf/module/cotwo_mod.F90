MODULE COTWO_MOD

CONTAINS
SUBROUTINE COTWO(KIDIA,KFDIA,KLON,LDLAND, PAN, PAG, PRD, PGS, PGC, PCSP, &
 & PDS, PDMAX, PIA, PGAMMT, YDAGS, &
 & PFZERO, PGMEST, PEPSO, PAMMAX  )

!***
! (C) Copyright 1997- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**   *COTWO* - CALCULATES NET ASSIMILATION OF CO2 AND LEAF CONDUCTANCE

!     A. Boone       * Meteo-France *     27/10/97 
!     (following Belair)
!     MODIFIED BY
!     M.H. Voogt (KNMI) "C-Tessel"  09/2005 
!     S. Lafont (ECMWF) vectorisation 04/2006 
!     S. Lafont (ECMWF) optimization  10/2006
!     G. Balsamo (ECMWF) 24/3/2014 cleaning and LDLAND protection
!      F. Vana  05-Mar-2015  Support for single precision
!     PURPOSE
!     -------
!     Calculates net assimilation of CO2 and leaf conductance.
              
!     INTERFACE
!     ---------
!     COTWO IS CALLED BY SUCOTWO AND COTWORESTRESS 

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (REAL):

!     *PGC*          CUTICULAR CONDUCTANCE                        M/S
!     *PCSP*         ATMOSPHERIC CONCENTRATION OF CO2      KG_CO2 KG_AIR-1
!     *PDS*          SPECIFIC HUMIDITY DEFICIT                    KG/KG
!     *PDMAX*	     MAXIMUM SPECIFIC HUMIDITY DEFICIT 
!                      TOLERATED BY VEGETATION                    KG/KG
!     *PIA*          ABSORBED PAR                                 W/M
!     *PGAMMT*       CO2 COMPENSATION POINT AT T=TSKIN     KG_CO2 KG_AIR-1
!     *PFZERO*       IDEAL VALUE OF F, NO PHOTORESPIRATION OR 
!                      SATURATION DEFICIT                         -
!     *PGMEST*       MESOPHYLL CONDUCTANCE AT T=TSKIN             M/S
!     *PEPSO*        MAXIMUM INITIAL QUANTUM USE EFFICIENCY KG_CO2 J-1 PAR M3 KG_AIR-1
!     *PAMMAX*       LEAF PHOTOSYNTHETIC CAPACITY           KG_CO2 KG_AIR-1 M S-1

!     OUTPUT PARAMETERS (REAL):      

!     *PAN*          NET ASSIMILATION OF CO2                KG_CO2 KG_AIR-1 M S-1
!     *PAG*          GROSS ASSIMILATION OF CO2              KG_CO2 KG_AIR-1 M S-1
!     *PRD*          DARK RESPIRATION OF CO2                KG_CO2 KG_AIR-1 M S-1

!     *PGS*          LEAF CONDUCTANCE TO H20 (CUTICULAR AND STOMATAL) M/S

!     METHOD
!     ------
!     Calvet et al. 1998 Forr. Agri. Met. [from model of Jacobs(1994)]

!     REFERENCE
!     ---------
!     Calvet et al. 1998 Forr. Agri. Met. 
      
!     --------------------------------------------------------------------------
USE PARKIND1, ONLY : JPIM, JPRB,JPRD
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_AGS,  ONLY : TAGS

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
LOGICAL           ,INTENT(IN)    :: LDLAND(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGC(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSP(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDS(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDMAX(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIA(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAMMT(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFZERO(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMEST(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEPSO(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAMMAX(:)
TYPE(TAGS)        ,INTENT(IN)    :: YDAGS
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAN(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAG(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRD(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGS(:)

!*         0.     LOCAL VARIABLES.
!                 ----- ----------

REAL(KIND=JPRD) :: ZEVAP(KLON), ZAMIN(KLON), ZGSC(KLON), ZGS(KLON), ZAM(KLON), ZCI(KLON), ZCSP(KLON) 
INTEGER(KIND=JPIM) JL, ITERCSTOM, INTER
!  ZEVAP   = leaf transpiration (kgH20 kgAir-1 m s-1)
!  ZAMIN   = Amin; minimum net assimilation (kgCO2 kgAir-1 m s-1)
!  ZGSC    = stomatal conductance to CO2 (m s-1)
!  ZGS     = stomatal conductance to H2O (m s-1)
!  ZAM     = Am; net assimilation as a function of CO2 deficit (kgCO2 kgAir-1 m s-1)
!  ZCI     = Leaf internal concentration of CO2 (kgCO2 kgAir-1)
!  ZCSP    = variable for max(PCSP,PGAMMT)

!      -------------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRD) :: ZF(KLON), ZFMIN(KLON), ZDRAP(KLON)
   !  ZF    = factor related to diffusion (-)
   !  ZFMIN = minimum f factor (-)
   !  ZDRAP = ratio Ds/Dmax (-)
REAL(KIND=JPRD) :: ZEPS(KLON)
REAL(KIND=JPRD) :: ZEXP1
   !  ZEPS     = initial quantum use efficiency (kgCO2 J-1 PAR m3 kgAir-1)
REAL(KIND=JPRD) :: ZCMIN(KLON)
   !  ZCMIN    = minimim internal leaf CO2 concentration (kgCO2 kgAir-1)

REAL(KIND=JPRD) :: ZAGR(KLON), ZAG(KLON), Z1
   !  ZAGR    = assimilation rate ratio (-)
   !  ZAG     = modified gross assimilation rate (kgCO2 kgAir-1 m s-1)

REAL(KIND=JPRD) :: ZEPSILON

IF (LHOOK) CALL DR_HOOK('COTWO_MOD:COTWO',0,ZHOOK_HANDLE)

ASSOCIATE(RAIRTOH2O=>YDAGS%RAIRTOH2O, RCO2TOH2O=>YDAGS%RCO2TOH2O, &
 & RCONDSTMIN=>YDAGS%RCONDSTMIN, RRDCF=>YDAGS%RRDCF)

ZEPSILON=100._JPRD*EPSILON(ZEPSILON)

!*       1.     COMPUTE PRELIMINARY QUANITIES NEEDED FOR CO2 MODEL
!               --------------------------------------------------
! initialization of local variable
DO JL=KIDIA,KFDIA
  ZEVAP(JL)=0._JPRD    ! initialize leaf transpiration
  ZAMIN(JL)=0._JPRD
  ZGSC(JL)=0._JPRD
  ZGS(JL)=0._JPRD
  ZAM(JL)=0._JPRD
  ZCI(JL)=0._JPRD
  ZF(JL)=0._JPRD
  ZFMIN(JL)=0._JPRD
  ZDRAP(JL)=0._JPRD
  ZEPS(JL)=0._JPRD
  ZCMIN(JL)=0._JPRD
  ZAGR(JL)=0._JPRD
  ZAG(JL)=0._JPRD
END DO
ZEXP1=0._JPRD
Z1=0._JPRD
!*       2.     BEGIN CO2 MODEL:
!               ----------------
DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN
    ZCSP(JL)=MAX(PCSP(JL),PGAMMT(JL)+1.E-6_JPRD)

! Equations from Jacobs (1994) Ph.D. Thesis:

    ZFMIN(JL)=PGC(JL)/(PGC(JL)+PGMEST(JL))                                 
! fmin <= f0, and so f <= f0
    ZFMIN(JL)=MIN(ZFMIN(JL),PFZERO(JL))
! fmin > 0 so ZCI > PGAMMT   
    ZFMIN(JL)=MAX(ZFMIN(JL),1.E-7_JPRD)

! f from specific humidity deficit ds
   
    ZDRAP(JL)=MIN(1._JPRD,PDS(JL)/PDMAX(JL))
    ZF(JL)=PFZERO(JL)*(1._JPRD-ZDRAP(JL))+ZFMIN(JL)*ZDRAP(JL)
   
! ci/cs ratio = f+(1.-f)*gammt/cs ; internal leaf CO2 concentration
    ZCI(JL)=ZCSP(JL)*(ZF(JL)+(1._JPRD-ZF(JL))*PGAMMT(JL)/ZCSP(JL))
! fmin > 0 so ZCI > PGAMMT   
    ZCI(JL)=MAX(ZCI(JL),PGAMMT(JL))
! f0 <= 1 so f <= 1
    ZCI(JL)=MIN(ZCI(JL),ZCSP(JL))

    ZCMIN(JL)=(PGC(JL)*ZCSP(JL)+PGMEST(JL)*PGAMMT(JL))/(PGMEST(JL)+PGC(JL))
    ZAMIN(JL) = PGMEST(JL)*(ZCMIN(JL)-PGAMMT(JL))
! Initial quantum use efficiency 
    ZEPS(JL)=PEPSO(JL)*(ZCI(JL)-PGAMMT(JL))/(ZCI(JL)+2._JPRD*PGAMMT(JL))

! Light response curve 
    ZAM(JL)=PGMEST(JL)*(ZCI(JL)-PGAMMT(JL))
    ZAM(JL)=PAMMAX(JL)*(1._JPRD-EXP(-ZAM(JL)/PAMMAX(JL)))
    ZAM(JL)=MAX(ZAM(JL),ZAMIN(JL))

    PRD(JL)=ZAM(JL)*RRDCF
    ZEXP1=(1._JPRD-EXP(-ZEPS(JL)*PIA(JL)/ &
    & SIGN(MAX(ABS(ZAM(JL)+PRD(JL)),ZEPSILON),ZAM(JL)+PRD(JL))))

    PAN(JL)=(ZAM(JL)+PRD(JL))*ZEXP1-PRD(JL)
    PAN(JL)=MAX(-PRD(JL),PAN(JL))

    PAG(JL)=(ZAM(JL)+PRD(JL))*ZEXP1
  ELSE
    PAN(JL)=0.0_JPRD
    PAG(JL)=0.0_JPRD
    PRD(JL)=0.0_JPRD
    PGS(JL)=0.0_JPRD
  ENDIF
END DO

! Iterations are for stomatal conductance and stomatal evaporation only
INTER=3
DO ITERCSTOM=1,INTER

  DO JL=KIDIA,KFDIA
    IF (LDLAND(JL)) THEN
      ZAGR(JL)=(PAN(JL)+PRD(JL))/SIGN(MAX(ABS(ZAM(JL)+PRD(JL)),ZEPSILON),ZAM(JL)+PRD(JL))
      ZDRAP(JL)=MIN(1._JPRD,PDS(JL)/PDMAX(JL))
      ZAG(JL)=PAN(JL)-ZAMIN(JL)*ZDRAP(JL)*ZAGR(JL)+PRD(JL)*(1._JPRD-ZAGR(JL))
      Z1=1._JPRD/(PCSP(JL)-ZCI(JL))
      ZGSC(JL)=ZAG(JL)*Z1
      ZGSC(JL)=MAX(RCONDSTMIN,ZGSC(JL))
      ZGSC(JL)=ZGSC(JL)+RAIRTOH2O*ZEVAP(JL)*((PCSP(JL)+ZCI(JL))/2._JPRD*Z1)
      ZGS(JL)=RCO2TOH2O*ZGSC(JL)
      ZEVAP(JL)=ZGS(JL)*PDS(JL)
    ENDIF
  END DO

ENDDO
! End of iterations

! Final calculation of leaf conductance (stomatal AND cuticular)
DO JL=KIDIA,KFDIA  
  IF (LDLAND(JL)) THEN
    PGS(JL)=RCO2TOH2O*ZGSC(JL)+PGC(JL)
  ENDIF
END DO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('COTWO_MOD:COTWO',1,ZHOOK_HANDLE)
END SUBROUTINE COTWO
END MODULE COTWO_MOD

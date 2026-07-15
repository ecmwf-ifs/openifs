
! (C) Copyright 1990- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

!**   *VUPDZ0* - COMPUTES Z0M,Z0H,Z0Q OVER SEA; SETS Z0H,Z0Q OVER LAND

!     Original   A.C.M. BELJAARS       E.C.M.W.F.    26/03/90.
!     Modified   A.C.M. BELJAARS  26/03/99   Surface tiling
!     Modified   P. Viterbo  ECMWF 12/05/2005 Externalize SURF
!     Modified   A. Beljaars ECMWF 03/12/2005 Roughness tables + TOFD
!     Modified   A. Beljaars ECMWF 17/05/2007 Clean-up of z0 initialization
!     Modified   E. Dutra/G. Balsamo 01/05/2008 Lake tile
!     Modified   J. Bidlot 21/10/2013 make sea ice roughness a function of sea ice cover.
!     Modified   I. Sandu and G. Balsamo 21/01/2015 gradual change of roughness for exposed snow
!     Modified   J. Bidlot 15/12/2018 to use PZ0WN to initialise Z0 over the oceans for step 0
!     Modified   M. Kelbling and S. Thober (UFZ) 11/6/2020 implemented spatially distributed parameters and
!                                               use of parameter values defined in namelist
!     Modified   J. Bidlot 15/02/2021 Sea state effect in Z0H and Z0Q over the oceans (under LWCOU2W and LWCOUHMF switches).
!     Modified   J. McNorton 24/08/2022 urban tile
!     Modified   I. Ayan-Miguez (BSC) Sep 2023 Added PSSDP2 object for surface spatially distributed parameters
!     Modified   J. Bidlot 19/12/2023 introduce the option of sea ice roughness for heat and momentum variable (not yet active!).

!     PURPOSE
!     -------

!     DERIVE Z0M,Z0H AND Z0Q FROM SURFACE FLUXES OVER SEA, SET Z0H AND
!     Z0Q OVER LAND AND DERIVE THE BUOYANCY FLUX.
!     (THE T-1 VALUES ARE UPDATED WITH FLUXES FROM THE PREVIOUS TIME
!      STEP)

!     INTERFACE
!     ---------

!     *VUPDZ0* IS CALLED BY *SURFEXCDRIVER_CTL*

!     Integer (In):
!     *KIDIA*        START OF LOOPS
!     *KFDIA*        END OF LOOPS
!     *KLON*         NUMBER OF POINTS IN PACKET
!     *KTILES*       NUMBER OF TILES
!     *KSTEP*        Time step index

!    Characters (In):
!     *CDCONF*       IFS Configuration


!    Integer (in):
!     *KTVL*         LOW VEGETATION TYPE 
!     *KTVH*         HIGH VEGETATION TYPE 

!    Reals (In):
!     *PCVL*         LOW VEGETATION COVER (CLIMATOLOGICAL)
!     *PCVH*         HIGH VEGETATION COVER (CLIMATOLOGICAL)
!     *PCUR*         URBAN COVER
!     *PUMLEV*       WIND X-COMPONENT AT T-1, lowest model level
!     *PVMLEV*       WIND Y-COMPONENT AT T-1, lowest model level
!     *PTMLEV*       TEMPERATURE AT T-1, lowest model level
!     *PQMLEV*       SPECIFIC HUMUDITY AT T-1, lowest model level
!     *PAPHMS*       PRESSURE AT T-1, surface
!     *PGEOMLEV*     GEOPOTENTIAL T-1, lowest model level
!     *PDSN*         Total snow depth (m) 
!     *PUSTRTI*      X-STRESS
!     *PVSTRTI*      Y-STRESS
!     *PAHFSTI*      SENSIBLE HEAT FLUX
!     *PEVAPTI*      MOISTURE FLUX
!     *PHLICE*       LAKE ICE THICKNESS
!     *PTSKTI*       SURFACE TEMPERATURE
!     *PCHAR*        CHARNOCK PARAMETER
!     *PCHARHQ*      EQUIVALENT CHARNOCK PARAMETER FOR HEAT AND MOISTURE
!     *PUCURR*       OCEAN CURRENT U-COMPONENT
!     *PVCURR*       OCEAN CURRENT V-COMPONENT
!     *PFRTI*        TILE FRACTION

!    Logicals (In):
!    *LDSICE*     SEA ICE MASK (.T. OVER SEA ICE)

!    Reals (Out):
!     *PZ0MTI*       NEW AERODYNAMIC ROUGHNESS LENGTH
!     *PZ0HTI*       NEW ROUGHNESS LENGTH FOR HEAT
!     *PZ0QTI*       NEW ROUGHNESS LENGTH FOR MOISTURE
!     *PBUOMTI*      BUOYANCY FLUX
!     *PZDLTI*       Z/L AT LOWEST MODEL LEVEL
!     *PRAQTI*       PRELIMINARY AERODYNAMIC RESISTANCE FOR MOISTURE 

!    Additional parameters for boundary condition (in SCM model):

!    *LROUGH*       If .TRUE. surface roughness length is externally specified 
!    *REXTZ0M*      Roughness length for momentum [m]
!    *REXTZ0H*      Roughness length for heat [m]

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

MODULE VUPDZ0_MOD

USE PARKIND1  , ONLY : JPIM, JPRB, JPRD

IMPLICIT NONE

PRIVATE PZ0WN
PUBLIC VUPDZ0

CONTAINS
SUBROUTINE VUPDZ0(KIDIA,KFDIA,KLON,KTILES,KSTEP,CDCONF,&
 & KTVL,KTVH,PCVL,PCVH,PCUR,PUMLEV,PVMLEV,&
 & PTMLEV,PQMLEV,PAPHMS,PGEOMLEV,PDSN,&
 & PUSTRTI,PVSTRTI,PAHFSTI,PEVAPTI,&
 & PHLICE, &
 & PTSKTI,PCHAR,PCHARHQ,PFRTI,PUCURR,PVCURR,&
 & YDCST,YDEXC,YDVEG,YDFLAKE,YDURB,&
 & PZ0MTI,PZ0HTI,PZ0QTI,PBUOMTI,PZDLTI,PRAQTI)

USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_EXCS  , ONLY : RCHBCD, DRITBL, RCHBBCD, RCHBB, RCHB23A, RITBL, &
 & RCHBA, RCHBHDL, RCDHALF, RCHETB, RCHBD, RCHETA, RCDHPI2, JPRITBL
USE YOS_EXC   , ONLY : TEXC
USE YOS_CST   , ONLY : TCST
USE YOS_VEG   , ONLY : TVEG
USE YOS_FLAKE , ONLY : TFLAKE
USE YOS_URB   , ONLY : TURB

!     ------------------------------------------------------------------

!**   *VUPDZ0* - COMPUTES Z0M,Z0H,Z0Q OVER SEA; SETS Z0H,Z0Q OVER LAND

!     Original   A.C.M. BELJAARS       E.C.M.W.F.    26/03/90.
!     Modified   A.C.M. BELJAARS  26/03/99   Surface tiling
!     Modified   P. Viterbo  ECMWF 12/05/2005 Externalize SURF
!     Modified   A. Beljaars ECMWF 03/12/2005 Roughness tables + TOFD
!     Modified   A. Beljaars ECMWF 17/05/2007 Clean-up of z0 initialization
!     Modified   E. Dutra/G. Balsamo 01/05/2008 Lake tile
!     Modified   J. Bidlot 21/10/2013 make sea ice roughness a function of sea ice cover.
!     Modified   I. Sandu and G. Balsamo 21/01/2015 gradual change of roughness for exposed snow
!     Modified   J. Bidlot 15/12/2018 to use PZ0WN to initialise Z0 over the oceans for step 0
!     Modified   J. Bidlot 15/02/2021 Sea state effect in Z0H and Z0Q over the oceans (under LWCOU2W and LWCOUHMF switches).
!     Modified   J. McNorton 24/08/2022 urban tile

!     PURPOSE
!     -------

!     DERIVE Z0M,Z0H AND Z0Q FROM SURFACE FLUXES OVER SEA, SET Z0H AND
!     Z0Q OVER LAND AND DERIVE THE BUOYANCY FLUX.
!     (THE T-1 VALUES ARE UPDATED WITH FLUXES FROM THE PREVIOUS TIME
!      STEP)

!     INTERFACE
!     ---------

!     *VUPDZ0* IS CALLED BY *SURFEXCDRIVER_CTL*

!     Integer (In):
!     *KIDIA*        START OF LOOPS
!     *KFDIA*        END OF LOOPS
!     *KLON*         NUMBER OF POINTS IN PACKET
!     *KTILES*       NUMBER OF TILES
!     *KSTEP*        Time step index

!    Characters (In):
!     *CDCONF*       IFS Configuration


!    Integer (in):
!     *KTVL*         LOW VEGETATION TYPE 
!     *KTVH*         HIGH VEGETATION TYPE 

!    Reals (In):
!     *PCVL*         LOW VEGETATION COVER (CLIMATOLOGICAL)
!     *PCVH*         HIGH VEGETATION COVER (CLIMATOLOGICAL)
!     *PCUR*         URBAN COVER
!     *PUMLEV*       WIND X-COMPONENT AT T-1, lowest model level
!     *PVMLEV*       WIND Y-COMPONENT AT T-1, lowest model level
!     *PTMLEV*       TEMPERATURE AT T-1, lowest model level
!     *PQMLEV*       SPECIFIC HUMUDITY AT T-1, lowest model level
!     *PAPHMS*       PRESSURE AT T-1, surface
!     *PGEOMLEV*     GEOPOTENTIAL T-1, lowest model level
!     *PDSN*         Total snow depth (m) 
!     *PUSTRTI*      X-STRESS
!     *PVSTRTI*      Y-STRESS
!     *PAHFSTI*      SENSIBLE HEAT FLUX
!     *PEVAPTI*      MOISTURE FLUX
!     *PHLICE*       LAKE ICE THICKNESS
!     *PTSKTI*       SURFACE TEMPERATURE
!     *PCHAR*        CHARNOCK PARAMETER
!     *PCHARHQ*      EQUIVALENT CHARNOCK PARAMETER FOR HEAT AND MOISTURE
!     *PUCURR*       OCEAN CURRENT U-COMPONENT
!     *PVCURR*       OCEAN CURRENT V-COMPONENT
!     *PFRTI*        TILE FRACTION

!    Reals (Out):
!     *PZ0MTI*       NEW AERODYNAMIC ROUGHNESS LENGTH
!     *PZ0HTI*       NEW ROUGHNESS LENGTH FOR HEAT
!     *PZ0QTI*       NEW ROUGHNESS LENGTH FOR MOISTURE
!     *PBUOMTI*      BUOYANCY FLUX
!     *PZDLTI*       Z/L AT LOWEST MODEL LEVEL
!     *PRAQTI*       PRELIMINARY AERODYNAMIC RESISTANCE FOR MOISTURE 

!    Additional parameters for boundary condition (in SCM model):

!    *LROUGH*       If .TRUE. surface roughness length is externally specified 
!    *REXTZ0M*      Roughness length for momentum [m]
!    *REXTZ0H*      Roughness length for heat [m]

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP

CHARACTER(LEN=1)  ,INTENT(IN)   ,OPTIONAL :: CDCONF 
 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCUR(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSN(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHARHQ(:) 
REAL(KIND=JPRB)   ,INTENT(IN)   ,OPTIONAL :: PUCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)   ,OPTIONAL :: PVCURR(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0MTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0HTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0QTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBUOMTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZDLTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAQTI(:,:) 
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TEXC)        ,INTENT(IN)    :: YDEXC
TYPE(TVEG)        ,INTENT(IN)    :: YDVEG
TYPE(TFLAKE)      ,INTENT(IN)    :: YDFLAKE
TYPE(TURB)        ,INTENT(IN)    :: YDURB

!*    LOCAL STORAGE
!     ----- -------

INTEGER(KIND=JPIM) :: JL,JTILE

REAL(KIND=JPRB) :: Z1DZ0Q, ZCON2, ZIPBL, ZNLEV, ZPRH1,&
 & ZPRQ, ZPRQ0, ZROWQ, ZROWT, &
 & ZWST2, &
 & ZXLNQ,ZZCDN,ZCDFC  
REAL(KIND=JPRB) :: Z0M, Z0H, Z0Q, Z0WHQ
REAL(KIND=JPRB) :: Z0W(KLON), ZUSTM1(KLON)
REAL(KIND=JPRB) :: ZUST(KLON,KTILES),ZUST2(KLON,KTILES)
REAL(KIND=JPRB) :: ZDUA(KLON), ZDU2(KLON),ZRHO(KLON), ZSNWGHT(KLON), ZMLOW, ZHLOW
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZLICE(KLON),ZLWAT(KLON)
REAL(KIND=JPRD) :: ZDUMMY
REAL(KIND=JPRD) :: zat_sn, zbt_sn,  zct_sn , zdt_sn, zaq_sn, zbq_sn, zcq_sn, zdq_sn, Rstar
REAL(KIND=JPRB) :: ZURBF,ZPCVL, ZPCVH, ZPCVB
REAL(KIND=JPRB) :: ZZ0HSNOW, ZZ0QSNOW


LOGICAL :: LLCURR,LLINIT

!             INCLUDE STABILITY FUNCTIONS
!             ------- --------- ---------

#include "fcsvdfs.h"

!             INCLUDE ROUGNESS LENGTH FUNCTIONS
!             ------- -------- ------ ---------
#include "fcz0.h"

!     ------------------------------------------------------------------

!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------
IF (LHOOK) CALL DR_HOOK('VUPDZ0_MOD:VUPDZ0',0,ZHOOK_HANDLE)
ASSOCIATE(RCPD=>YDCST%RCPD, RD=>YDCST%RD, RETV=>YDCST%RETV, RG=>YDCST%RG, &
 & LROUGH=>YDEXC%LROUGH, LSCMEC=>YDEXC%LSCMEC, REPDU2=>YDEXC%REPDU2, &
 & REPUST=>YDEXC%REPUST, REXTZ0H=>YDEXC%REXTZ0H, REXTZ0M=>YDEXC%REXTZ0M, &
 & RKAP=>YDEXC%RKAP, RNUH=>YDEXC%RNUH, RNUM=>YDEXC%RNUM, RNUQ=>YDEXC%RNUQ, &
 & RPARZI=>YDEXC%RPARZI, RZ0ICE=>YDEXC%RZ0ICE, &
 & LEFLAKE=>YDFLAKE%LEFLAKE, LEURBAN=>YDURB%LEURBAN, RH_ICE_MIN_FLK=>YDFLAKE%RH_ICE_MIN_FLK, &
 & LWCOU2W=>YDEXC%LWCOU2W, LWCOUHMF=>YDEXC%LWCOUHMF, &
 & RVZ0H=>YDVEG%RVZ0H, RVZ0M=>YDVEG%RVZ0M, &
 & RURBZTM=>YDURB%RURBZTM,RURBZTH=>YDURB%RURBZTH)

IF (KTILES.LT.8) THEN
  STOP "Wrong number of tiles in VDFUPDZ0"
ENDIF
IF (LEFLAKE .OR. KTILES .GT. 8) THEN
!* FIND LAKE POINTS WITH ICE COVER 
  DO JL=KIDIA,KFDIA
    IF(PHLICE(JL) > RH_ICE_MIN_FLK ) THEN
      ZLICE(JL)=1._JPRB
      ZLWAT(JL)=0._JPRB
    ELSE
      ZLICE(JL)=0._JPRB
      ZLWAT(JL)=1._JPRB
    ENDIF
  ENDDO
ENDIF 

ZURBF = 0._JPRB

ZCON2  =2.0_JPRB/3._JPRB

! Constants for Andreas z0h computation for snow
zat_sn=-0.0061586_JPRD 
zbt_sn=-0.12756_JPRD
zct_sn=-0.66267_JPRD
zdt_sn=0.25344_JPRD
zaq_sn=-0.0054869_JPRB
zbq_sn=-0.12027_JPRB
zcq_sn=-0.68407_JPRB
zdq_sn=0.48260_JPRB


!     PBL HEIGHT FOR W* - EFFECT

ZIPBL=RPARZI

IF(PRESENT(PUCURR).AND.PRESENT(PVCURR)) THEN
  LLCURR=.TRUE.
ELSE
  LLCURR=.FALSE.
ENDIF

LLINIT= ( KSTEP == 0)


!     ------------------------------------------------------------------

!*         2.      PRE-COMPUTATION OF TILE INDEPENDENT ARRAYS
!                  
DO JL=KIDIA,KFDIA
  ZRHO(JL)=PAPHMS(JL)/( RD*PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)) )
  IF(LLCURR) THEN
    ZDU2(JL)=MAX(REPDU2,(PUMLEV(JL)-PUCURR(JL))**2+&
     & (PVMLEV(JL)-PVCURR(JL))**2)
  ELSE
    ZDU2(JL)=MAX(REPDU2,PUMLEV(JL)**2+PVMLEV(JL)**2)
  ENDIF
  ZSNWGHT(JL)=MIN(MAX(PDSN(JL),0.0_JPRB)/0.25_JPRB,1.0_JPRB)
  ZDUA(JL)=SQRT(ZDU2(JL))
ENDDO
!*         3.   ESTIMATE SURF.FL. FOR STEP 0
!*              (ASSUME NEUTRAL STRATIFICATION)
IF (LLINIT) THEN
  DO JL=KIDIA,KFDIA
    ZDUA(JL)=SQRT(ZDU2(JL))
!       - Preliminary value for water
    PZ0MTI(JL,1)=PZ0WN(ZDUA(JL),PGEOMLEV(JL), PCHAR(JL), RG, RNUM, RKAP)
!       - Sea ice
    PZ0MTI(JL,2)=PZ0ICE(RZ0ICE,PFRTI(JL,2))
!       - Wet skin
    PZ0MTI(JL,3)=PCVL(JL)*RVZ0M(KTVL(JL))+PCVH(JL)*RVZ0M(KTVH(JL))&
      & +(1.-PCVL(JL)-PCVH(JL))*RVZ0M(0)
!       - Low Vegetation
    PZ0MTI(JL,4)=RVZ0M(KTVL(JL))
!       - Exposed snow
! blend between roughness of low vegetation and soil and roughness of ice caps, when snow depth below 25cm 
    ZMLOW=PCVL(JL)/(1.0_JPRB-PCVH(JL))*RVZ0M(KTVL(JL))+(1.0_JPRB-PCVL(JL)-PCVH(JL))/(1.0_JPRB-PCVH(JL))*RVZ0M(0)
    PZ0MTI(JL,5)=ZSNWGHT(JL)*RVZ0M(12)+(1.0_JPRB-ZSNWGHT(JL))*ZMLOW
!       - High vegetation
    PZ0MTI(JL,6)=RVZ0M(KTVH(JL))
!       - Sheltered snow
    PZ0MTI(JL,7)=RVZ0M(KTVH(JL))
!       - Bare soil
    PZ0MTI(JL,8)=RVZ0M(0)
    IF (LEFLAKE .OR. KTILES .GT. 8 ) THEN
!       - LAKES   
      PZ0MTI(JL,9)=PZ0MTI(JL,1)*ZLWAT(JL)+PZ0MTI(JL,2)*ZLICE(JL)   
    ENDIF
!       - URBAN
    IF (LEURBAN) THEN
      ZURBF=PCUR(JL)
      ZPCVL=MAX(PCVL(JL)*(1.0_JPRB-ZURBF),0.0_JPRB)
      ZPCVH=MAX(PCVH(JL)*(1.0_JPRB-ZURBF),0.0_JPRB)
      ZPCVB=MAX(1.0_JPRB-ZPCVL-ZPCVH-PCUR(JL),0.0_JPRB)
      IF (ZURBF.GT.0.0_JPRB) THEN
      PZ0MTI(JL,3)=(ZPCVL*RVZ0M(KTVL(JL))+ZPCVH*RVZ0M(KTVH(JL))&
      & +ZPCVB*RVZ0M(0)+RURBZTM*ZURBF)
      ZMLOW=ZPCVL/(1.0_JPRB-ZPCVH)*RVZ0M(KTVL(JL))+ZPCVB/(1.0_JPRB-ZPCVH)*RVZ0M(0)
      PZ0MTI(JL,5)=ZSNWGHT(JL)*RVZ0M(12)+(1.0_JPRB-ZSNWGHT(JL))*ZMLOW+ZURBF/(1.0_JPRB-ZPCVH)*RURBZTM
      ENDIF
      PZ0MTI(JL,10)=RURBZTM
    ENDIF
  ENDDO

  DO JTILE=1,KTILES
    DO JL=KIDIA,KFDIA
      ZZCDN=(RKAP/LOG(1.0_JPRB+PGEOMLEV(JL)/(RG*PZ0MTI(JL,JTILE))))**2
      PUSTRTI(JL,JTILE)=ZRHO(JL)*PUMLEV(JL)*ZDUA(JL)*ZZCDN
      PVSTRTI(JL,JTILE)=ZRHO(JL)*PVMLEV(JL)*ZDUA(JL)*ZZCDN
      PAHFSTI(JL,JTILE)=0.0_JPRB
      PEVAPTI(JL,JTILE)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF

!*         4.      STABILITY PARAMETERS AND FREE CONVECTION
!                  VELOCITY SCALE
!
!                  ZCDFC is a rough estimate of the drag coefficient used  
!                  to convert w*-gustiness at the lowest model level into u*. 
!                  The value is choosen from Fig. 1 on page 37 of the ECMWF 
!                  seminar proceedings on "Atmopshere-surface interaction", 
!                  ie charqacteristic for 1 m/s in unstable situations. 

!* DETERMINE the friction velocity u*

ZCDFC=2.E-3_JPRB
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZROWQ=PEVAPTI(JL,JTILE)
    ZROWT=PAHFSTI(JL,JTILE)/RCPD
    PBUOMTI(JL,JTILE)=RG*(-RETV*ZROWQ-ZROWT/PTSKTI(JL,JTILE))/ZRHO(JL)
    ZUST2(JL,JTILE)=SQRT(PUSTRTI(JL,JTILE)**2+PVSTRTI(JL,JTILE)**2)/ZRHO(JL)

!     APPLY W* CORRECTION

    IF (PBUOMTI(JL,JTILE)  >  0.0_JPRB) THEN
      ZWST2=(PBUOMTI(JL,JTILE)*ZIPBL)**ZCON2
      ZUST2(JL,JTILE)=ZUST2(JL,JTILE)+ZCDFC*ZWST2
    ENDIF

    ZUST(JL,JTILE)=MAX(SQRT(ZUST2(JL,JTILE)),REPUST)
    PZDLTI(JL,JTILE)=-PGEOMLEV(JL)*RKAP*PBUOMTI(JL,JTILE)/(RG*ZUST(JL,JTILE)**3)
  ENDDO
ENDDO

!*         5.    SETTING OF ROUGHNESS LENGTHS 
!                ----------------------------

JTILE=1
!   - Ocean open water
! Momentum:
DO JL=KIDIA,KFDIA
  ZUSTM1(JL)=1.0_JPRB/ZUST(JL,JTILE)
  Z0M=RNUM*ZUSTM1(JL)
  Z0W(JL)=PZ0SEA(RG,PCHAR(JL),ZUST2(JL,JTILE))
  PZ0MTI(JL,JTILE)=Z0M+Z0W(JL)
ENDDO

! Heat and moisture:
IF( LWCOU2W .AND. LWCOUHMF) THEN
  ! With sea state effect:
  DO JL=KIDIA,KFDIA
    Z0WHQ=PZ0SEA(RG,PCHARHQ(JL),ZUST2(JL,JTILE))
    Z0H=RNUH*ZUSTM1(JL)
    PZ0HTI(JL,JTILE)=PZNSEA(Z0WHQ,Z0H,Z0W(JL),ZUST2(JL,JTILE))
    Z0Q=RNUQ*ZUSTM1(JL)
    PZ0QTI(JL,JTILE)=PZNSEA(Z0WHQ,Z0Q,Z0W(JL),ZUST2(JL,JTILE))
  ENDDO
ELSE
  ! Purely diffusive
  DO JL=KIDIA,KFDIA
    PZ0HTI(JL,JTILE)=RNUH*ZUSTM1(JL)
    PZ0QTI(JL,JTILE)=RNUQ*ZUSTM1(JL)
  ENDDO
ENDIF

JTILE = 2
!   - Sea ice
DO JL=KIDIA,KFDIA
!!!! I use the tile fraction because it should be equivalent to sea ice cover but it should actually be programmed with the sea ice cover !!!
  PZ0MTI(JL,JTILE)=PZ0ICE(RZ0ICE,PFRTI(JL,JTILE))
  PZ0HTI(JL,JTILE)=RZ0ICE
  PZ0QTI(JL,JTILE)=RZ0ICE
ENDDO

JTILE = 3
!   - Wet skin
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=PCVL(JL)*RVZ0M(KTVL(JL))+PCVH(JL)*RVZ0M(KTVH(JL))+(1.-PCVL(JL)-PCVH(JL))*RVZ0M(0)
  PZ0HTI(JL,JTILE)=PCVL(JL)*RVZ0H(KTVL(JL))+PCVH(JL)*RVZ0H(KTVH(JL))+(1.-PCVL(JL)-PCVH(JL))*RVZ0H(0)
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
  IF (LEURBAN) THEN
    ZURBF=PCUR(JL)
    PZ0MTI(JL,JTILE)=PZ0MTI(JL,JTILE)*(1.0_JPRB-ZURBF)+RURBZTM*ZURBF
    PZ0HTI(JL,JTILE)=PZ0HTI(JL,JTILE)*(1.0_JPRB-ZURBF)+RURBZTH*ZURBF
    PZ0QTI(JL,JTILE)=PZ0QTI(JL,JTILE)*(1.0_JPRB-ZURBF)+RURBZTH*ZURBF
  ENDIF
ENDDO

JTILE = 4
!   - Low Vegetation
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RVZ0M(KTVL(JL))
  PZ0HTI(JL,JTILE)=RVZ0H(KTVL(JL))
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

JTILE = 5
!   - Exposed snow
! blend between roughness of low vegetation and soil and roughness of ice caps, when snow depth below 25cm 
DO JL=KIDIA,KFDIA
  ZMLOW=PCVL(JL)/(1.0_JPRB-PCVH(JL))*RVZ0M(KTVL(JL))+(1.0_JPRB-PCVL(JL)-PCVH(JL))/(1.0_JPRB-PCVH(JL))*RVZ0M(0)
  ZHLOW=PCVL(JL)/(1.0_JPRB-PCVH(JL))*RVZ0H(KTVL(JL))+(1.0_JPRB-PCVL(JL)-PCVH(JL))/(1.0_JPRB-PCVH(JL))*RVZ0H(0)
  PZ0MTI(JL,JTILE)=ZSNWGHT(JL)*RVZ0M(12)+(1.0_JPRB-ZSNWGHT(JL))*ZMLOW
!*  PZ0HTI(JL,JTILE)=ZSNWGHT(JL)*RVZ0H(12)+(1.0_JPRB-ZSNWGHT(JL))*ZHLOW
!*  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)

  ! Use Andreas (2002) formula for z0h for exposed snow, tile 5
  ! doi: 10.1175/1525-7541(2002)003<0417:PSTOSA>2.0.CO;2
  ! The actual z0h is a log-average of vegetation and snow, dependending on ZSNWGHT
  ! Kinematic viscosity of air at 0C: 1.35E-5;
  Rstar=(REAL(RVZ0M(12)*ZUST(JL,JTILE),KIND=JPRD))/(1.35E-5_JPRD)
  ZZ0HSNOW=REAL(RVZ0M(12),KIND=JPRD)*exp( zat_sn*log(Rstar)**3._JPRB + zbt_sn*log(Rstar)**2._JPRB + zct_sn*log(Rstar) + zdt_sn )

  ZZ0QSNOW=REAL(RVZ0M(12),KIND=JPRD)*exp( zaq_sn*log(Rstar)**3._JPRB + zbq_sn*log(Rstar)**2._JPRB + zcq_sn*log(Rstar) + zdq_sn )

  PZ0HTI(JL,JTILE)=ZSNWGHT(JL)*ZZ0HSNOW+(1.0_JPRB-ZSNWGHT(JL))*ZHLOW
  PZ0QTI(JL,JTILE)=ZSNWGHT(JL)*ZZ0QSNOW+(1.0_JPRB-ZSNWGHT(JL))*ZHLOW
  IF (LEURBAN) THEN
    ZURBF=PCUR(JL)
    ZPCVL=MAX(PCVL(JL)*(1.0_JPRB-ZURBF),0.0_JPRB)
    ZPCVH=MAX(PCVH(JL)*(1.0_JPRB-ZURBF),0.0_JPRB)
    ZPCVB=MAX(1.0_JPRB-ZPCVL-ZPCVH-ZURBF,0.0_JPRB)
    IF (ZURBF.GT.0.0_JPRB) THEN
      ZMLOW=ZPCVL/(1.0_JPRB-ZPCVH)*RVZ0M(KTVL(JL))+ZPCVB/(1.0_JPRB-ZPCVH)*RVZ0M(0)
      ZHLOW=ZPCVL/(1.0_JPRB-ZPCVH)*RVZ0H(KTVL(JL))+ZPCVB/(1.0_JPRB-ZPCVH)*RVZ0H(0)
      PZ0MTI(JL,JTILE)=ZSNWGHT(JL)*RVZ0M(12)+(1.0_JPRB-ZSNWGHT(JL))*ZMLOW+ZURBF/(1.0_JPRB-ZPCVH)*RURBZTM
      PZ0HTI(JL,JTILE)=ZSNWGHT(JL)*RVZ0H(12)+(1.0_JPRB-ZSNWGHT(JL))*ZHLOW+ZURBF/(1.0_JPRB-ZPCVH)*RURBZTH
      PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
    ENDIF
  ENDIF

  PZ0HTI(JL,JTILE)=ZSNWGHT(JL)*RVZ0H(12)+(1.0_JPRB-ZSNWGHT(JL))*ZHLOW
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)


ENDDO

JTILE = 6
!   - High vegetation
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RVZ0M(KTVH(JL))
  PZ0HTI(JL,JTILE)=RVZ0H(KTVH(JL))
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

JTILE = 7
!   - Sheltered snow
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RVZ0M(KTVH(JL))
  PZ0HTI(JL,JTILE)=RVZ0H(KTVH(JL))
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

JTILE = 8
!   - Bare soil
DO JL=KIDIA,KFDIA
  PZ0MTI(JL,JTILE)=RVZ0M(0)
  PZ0HTI(JL,JTILE)=RVZ0H(0)
  PZ0QTI(JL,JTILE)=PZ0HTI(JL,JTILE)
ENDDO

IF (LEFLAKE .OR. KTILES .GT. 8) THEN
  JTILE = 9 
  DO JL=KIDIA,KFDIA
    PZ0MTI(JL,JTILE)=PZ0MTI(JL,1)*ZLWAT(JL)+PZ0MTI(JL,2)*ZLICE(JL)
    PZ0HTI(JL,JTILE)=PZ0HTI(JL,1)*ZLWAT(JL)+PZ0HTI(JL,2)*ZLICE(JL)
    PZ0QTI(JL,JTILE)=PZ0QTI(JL,1)*ZLWAT(JL)+PZ0QTI(JL,2)*ZLICE(JL)
  ENDDO 
ENDIF
IF (LEURBAN) THEN
 JTILE = 10
 !   - Urban
 DO JL=KIDIA,KFDIA
   PZ0MTI(JL,JTILE)=RURBZTM
   PZ0HTI(JL,JTILE)=RURBZTH
   PZ0QTI(JL,JTILE)=RURBZTH
 ENDDO
ENDIF

!*        6.   SCM: Fixed roughness lengths

IF (LSCMEC .AND. LROUGH) THEN
  PZ0MTI(:,:) = REXTZ0M   ! scm namelist parameters
  PZ0HTI(:,:) = REXTZ0H
ENDIF


!*        7.   COMPUTE PRELIMINARY AERODYNAMIC RESISTANCE FOR COMPUTATION 
!*             OF APPARENT SURFACE HUMIDITY IN VDFSURF
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZNLEV=PGEOMLEV(JL)/RG+PZ0MTI(JL,JTILE)
    ZXLNQ=LOG(ZNLEV/PZ0QTI(JL,JTILE))
    Z1DZ0Q=ZNLEV/PZ0QTI(JL,JTILE)
    IF (PZDLTI(JL,JTILE)  >  0.0_JPRB) THEN
       ZDUMMY = PZDLTI(JL,JTILE)
       ZPRH1=PSIHS(ZDUMMY)
       ZDUMMY = PZDLTI(JL,JTILE)/Z1DZ0Q
       ZPRQ0=PSIHS(ZDUMMY)
    ELSE
       ZDUMMY = PZDLTI(JL,JTILE)
       ZPRH1=PSIHU(ZDUMMY)
       ZDUMMY = PZDLTI(JL,JTILE)/Z1DZ0Q
       ZPRQ0=PSIHU(ZDUMMY)
    ENDIF
    ZPRQ   =ZXLNQ-ZPRH1+ZPRQ0
    PRAQTI(JL,JTILE)=(ZXLNQ-ZPRH1+ZPRQ0)/(ZUST(JL,JTILE)*RKAP)
  ENDDO
ENDDO
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VUPDZ0_MOD:VUPDZ0',1,ZHOOK_HANDLE)

END SUBROUTINE VUPDZ0

! Function to compute Z0M for neutral wind conditions
#include "fcz0wn.h"

END MODULE VUPDZ0_MOD

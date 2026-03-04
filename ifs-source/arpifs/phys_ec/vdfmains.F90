! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFMAINS (YDMODEL,KIDIA,KFDIA,KLON,KLEV,KLEVS,KLEVSN,KTILES,&
 & KTRAC,PTSPHY,KTVL,KTVH,PCVL,PCVH,PCUR, &
 & PLAIL, PLAIH, PFWET, PSNS  ,PRSN  ,PSIGFLT, &
 & PUM1  ,PVM1  ,PTM1  ,PQM1  ,PCM1  ,PAPHM1  ,PAPM1, PGEOM1, PGEOH , &
 & PIM1  ,PLM1, &
 & PTSKM1M,PTSAM1M,PWSAM1M, &
 & PSSRFL,PSLRFL ,PEMIS ,&
 & PTHKICE, PSNTICE, &
 & PTSNOW, PTICE,&
 ! & PHLICE, PTLICE, PTLWML, &  !lake passive
 & PSST  ,KSOTY  ,PFRTI ,PALBTI ,PWLMX ,&
 & PCHAR ,PTSKRAD, PCFLX,&
 ! OUTPUT (TRAJECTORY)           
 & PZ0M   ,PZ0H , &
 & PU10M  , PV10M  , PT2M   , PD2M   , PQ2M ,&
 & P10NU  , P10NV  ,&
 & PZINV,&
 & PSSRFLTI, PEVAPSNW,&
 ! INPUT TENDENCIES (TRAJECTORY)           
 & PTE0  ,PQE0  ,PVOM0 ,PVOL0,&
 ! OUTPUT TENDENCIES (TRAJECTORY)           
 & PTE   ,PQE   ,PVOM  ,PVOL ,PTENC  ,PTSKE1,&
 !-UPDATED FIELDS FOR TILES (TRAJECTORY)
 & PUSTRTI,PVSTRTI,PAHFSTI,PEVAPTI,PTSKTI,&
 ! OUTPUT FLUXES (TRAJECTORY)           
 & PDIFTS,PDIFTQ,PSTRTU,PSTRTV,&
 ! CO2 VEGETATION FLUX ADJUSTMENT COEFFICIENTS FROM BFAS
 & PCGPP  , PCREC , PAG, PRECO)
!-----------------------------------------------------------------------     

!***

!**   *VDFMAINS* - DOES THE VERTICAL EXCHANGE OF U,V,T,Q BY TURBULENCE.
!                   (Nonlinear version for trajectory in adjoint)

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE FOUR
!     PROGNOSTIC VARIABLES U,V,T AND Q DUE TO THE VERTICAL EXCHANGE BY
!     TURBULENT (= NON-MOIST CONVECTIVE) PROCESSES. THESE TENDENCIES ARE
!     OBTAINED AS THE DIFFERENCE BETWEEN THE RESULT OF AN IMPLICIT
!     TIME-STEP STARTING FROM VALUES AT T-1 AND THESE T-1 VALUES. ALL
!     THE DIAGNOSTIC COMPUTATIONS (EXCHANGE COEFFICIENTS, ...) ARE DONE
!      FROM THE T-1 VALUES. AS A BY-PRODUCT THE ROUGHNESS LENGTH OVER SEA
!     IS UPDATED ACCORDINGLY TO THE *CHARNOCK FORMULA. HEAT AND MOISTURE
!     SURFACE FLUXES AND THEIR DERIVATIVES AGAINST TS, WS AND WL
!     (THE LATTER WILL BE LATER WEIGHTED WITH THE SNOW FACTOR IN
!     *VDIFF*), LATER TO BE USED FOR SOIL PROCESSES TREATMENT, ARE ALSO
!     COMPUTED AS WELL AS A STABILITY VALUE TO BE USED AS A DIAGNOSTIC
!     OF THE DEPTH OF THE WELL MIXED LAYER IN CONVECTIVE COMPUTATIONS.

!     INTERFACE.
!     ----------
!          *VDIFF* TAKES THE MODEL VARIABLES AT T-1 AND RETURNS THE VALUES
!     FOR THE PROGNOSTIC TIME T+1 DUE TO VERTICAL DIFFUSION.
!     THE MODEL VARIABLES, THE MODEL DIMENSIONS AND THE DIAGNOSTICS DATA
!     ARE PASSED AS SUBROUTINE ARGUMENTS. CONSTANTS THAT DO NOT CHANGE
!     DURING A MODEL RUN (E.G. PHYSICAL CONSTANTS, SWITCHES ETC.) ARE
!     STORED IN A SINGLE COMMON BLOCK *YOMVDF*, WHICH IS INITIALIZED
!     BY SET-UP ROUTINE *SUVDF*.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLEV*         NUMBER OF LEVELS
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*         NUMBER OF SOIL LAYERS
!    *KTILES*       NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT
!                   OF SURFACE BOUNDARY CONDITION)
!    *KTRAC*        NUMBER OF TRACERS

!    *KTVL*         VEGETATION TYPE FOR LOW VEGETATION FRACTION
!    *KTVH*         VEGETATION TYPE FOR HIGH VEGETATION FRACTION
!    *KSOTY*        SOIL TYPE                                   (1-7)

!     INPUT PARAMETERS (LOGICAL)

!     INPUT PARAMETERS (REAL)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):

!    *PCVL*        LOW VEGETATION COVER                          -  
!    *PCVH*        HIGH VEGETATION COVER                         -  
!    *PCUR*        URBAN (PASSIVE)                               (0-1)
!    *PLAIL*        LOW VEGETATION LAI                           m2/m2
!    *PLAIH*        HIGH VEGETATION LAI                          m2/m2
!    *PFWET*       WETLAND FRACTION                              - 
!    *PSNS*         SNOW MASS                                    kg/m2
!    *PRSN*         SNOW DENSITY                                 kg/m3
!    *PSIGFLT*     STANDARD DEVIATION OF FILTERED OROGRAPHY      M
!    *PUM1*        X-VELOCITY COMPONENT                          M/S
!    *PVM1*        Y-VELOCITY COMPONENT                          M/S
!    *PTM1*        TEMPERATURE                                   K
!    *PQM1*        SPECIFIC HUMIDITY                             KG/KG
!    *PCM1*        TRACER CONCENTRATION                          KG/KG
!    *PAPHM1*      PRESSURE ON HALF-LEVELS                       PA
!    *PAPM1*       PRESSURE ON FULL-LEVELS                       PA
!    *PGEOM1*      GEOPOTENTIAL                                  M2/S2
!    *PGEOH*       GEOPOTENTIAL AT HALF LEVELS                   M2/S2
!    *PQI1*        CLOUD ICE WATER CONTENT                       KG/KG
!    *PQL1*        CLOUD LIQUID WATER CONTENT                    KG/KG
!    *PTSKM1M*     SKIN TEMPERATURE                              K
!    *PTSKRAD*     SKIN TEMPERATURE at last full radiation
!                    time step                                   K
!    *PTSAM1M*     SURFACE TEMPERATURE                           K
!    *PWSAM1M*     SOIL MOISTURE ALL LAYERS                      M**3/M**3
!    *PSSRFL*      NET SHORTWAVE RADIATION FLUX AT SURFACE       W/M2
!    *PSLRFL*      NET LONGWAVE RADIATION FLUX AT SURFACE        W/M2
!    *PEMIS*       MODEL SURFACE LONGWAVE EMISSIVITY            
!    *PTSNOW*      SNOW TEMPERATURE                              K
!    *PTICE*       ICE TEMPERATURE (TOP SLAB)                    K
!    *PHLICE*      LAKE ICE THICKNESS                            M 
!    *PTLICE*      LAKE ICE TEMPERATURE                          K
!    *PTLWML*      LAKE MIXED LAYER TEMPERATURE                  K
!    *PSST*         (OPEN) SEA SURFACE TEMPERATURE               K 

!    *PFRTI*        TILE FRACTIONS                              (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!    *PALBTI*      BROADBAND ALBEDO FOR TILE FRACTIONS        
!    *PWLMX*       MAXIMUM SKIN RESERVOIR CAPACITY                KG/M**2
!    *PCHAR*        "EQUIVALENT" CHARNOCK PARAMETER            
!    *PCFLX*       TRACER SURFACE FLUX                          KG/(M2 S)

!    UPDATED PARAMETERS (REAL):

!    *PTE*         TEMPERATURE TENDENCY                           K/S
!    *PQE*         MOISTURE TENDENCY                              KG/(KG S)
!    *PVOM*        MERIODINAL VELOCITY TENDENCY (DU/DT)           M/S2
!    *PVOL*        LATITUDE TENDENCY            (DV/DT)           M/S2
!    *PTENC*       TRACER TENDENCY                                KG/(KG S)
!    *PTSKE1*      SKIN TEMPERATURE TENDENCY                      K/S
!    *PZ0M*        AERODYNAMIC ROUGHNESS LENGTH                   M
!    *PZ0H*        ROUGHNESS LENGTH FOR HEAT                      M


!    UDATED PARAMETERS FOR TILES (REAL):

!    *PU10M*        U-COMPONENT WIND AT 10 M                      M/S
!    *PV10M*        V-COMPONENT WIND AT 10 M                      M/S
!    *P10NV*        U-COMPONENT NEUTRAL WIND AT 10 M              M/S
!    *P10NV*        V-COMPONENT NEUTRAL WIND AT 10 M              M/S
!    *PT2M*         TEMPERATURE AT 2M                             K
!    *PD2M*         DEW POINT TEMPERATURE AT 2M                   K
!    *PQ2M*         SPECIFIC HUMIDITY AT 2M                       KG/KG
!    *PUSTRTI*      SURFACE U-STRESS                              N/M2    
!    *PVSTRTI*      SURFACE V-STRESS                              N/M2   
!    *PAHFSTI*      SURFACE SENSIBLE HEAT FLUX                    W/M2 
!    *PEVAPTI*      SURFACE MOISTURE FLUX                         KG/M2/S
!    *PTSKTI*       SKIN TEMPERATURE                              K      

!    OUTPUT PARAMETERS (REAL):

!    *PSSRFLTI*     NET SHORTWAVE RADIATION FLUX AT SURFACE, FOR
!                      EACH TILE                                  W/M2
!    *PEVAPSNW*     EVAPORATION FROM SNOW UNDER FOREST            KG/(M2*S)
!    *PSTRTU*       TURBULENT FLUX OF U-MOMEMTUM                 KG*(M/S)/(M2*S)
!    *PSTRTV*       TURBULENT FLUX OF V-MOMEMTUM                 KG*(M/S)/(M2*S)
!    *PDIFTS*       TURBULENT FLUX OF HEAT                        J/(M2*S)
!    *PDIFTQ*       TRUBULENT FLUX OF SPECIFIC HUMIDITY           KG/(M2*S)

!     OPTIONAL PARAMETERS:

!    *PCGPP*        GPP flux adjustment coefficient               -
!    *PCREC*        REC flux adjustment coefficient               -
!    *PAG*          GPP flux                -
!    *PRECO*        REC flux                -

!     METHOD.
!     -------

!     SEE VDFMAIN

!     --------------------------------------------------------
!     THIS TANGENT LINEAR VERSION IS ONLY COMPATIBLE WITH THE
!     OPTION LPHYLIN=.TRUE. IN VDFMAIN
!     --------------------------------------------------------

!     EXTERNALS.
!     ----------

!     *VDFMAINS* CALLS SUCESSIVELY:
!         *SURFEXCDRIVERS*
!         *VDFEXCUS*
!         *VDFTOFDCS*
!         *VDFDIFMS*
!         *VDFDIFHS*
!         *VDFINCRS*
!         *VDFDIFCS*

!     REFERENCE.
!     ----------

!     SEE VERTICAL DIFFUSION'S PART OF THE MODEL'S DOCUMENTATION
!     FOR DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     AUTHOR.
!     -------
!      J.F. MAHFOUF     E.C.M.W.F.         05/10/95

!     MODIFICATIONS.
!     --------------
!      M. Janiskova              11/06/N9 (TL for surface tile)
!      M. Janiskova              03/09/03 (TL of new tile coupling)
!      P. Viterbo                24/05/2004 (Change surface units)
!      P. Lopez                  24/02/2005 (nonlinear version for trajectory in adjoint)
!      M. Janiskova              27/06/2005 (removed option for MASS
!                                           vector functions noy used
!                                           in corresponding TL/AD)
!      A. Beljaars               10/12/2005 TOFD
!      P. Lopez                  13/02/2006 (Added tracer diffusion)
!      G. Balsamo                15/01/2007 (Added soil type)
!      G. Balsamo                21/10/2008 (Added lake tile)
!      H. Hersbach               04/12/2009 (10-m neutral wind)
!      S. Boussetta/G.Balsamo      May 2009 (Add variable LAI)
!      L. Magnusson              28-09-2010 For NEMO-LIM
!      M. Janiskova               Sept 2011 (PSSRFLTI,PEVAPSNW as outputs)
!      P. Lopez                  18/01/2021 (New inversion height computation)
!      S. Massart                 March 2021 Adding BFAS parameters 
!     ------------------------------------------------------------------

USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RG, RCPD, RETV, RLVTT, RLSTT, RTT  
USE YOETHF    , ONLY : RVTMP2, R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                     R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
 &                     RTWAT_RTICE_R, RTWAT_RTICECU_R  
USE YOMCT3    , ONLY : NSTEP

IMPLICIT NONE

TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCUR(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIL(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAIH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFWET(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNS(KLON,KLEVSN)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSN(KLON,KLEVSN)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIGFLT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM1(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSAM1M(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSAM1M(KLON,KLEVS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTHKICE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNTICE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNOW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(KLON) 
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(KLON) !lake passive
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICE(KLON) !lake passive
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLWML(KLON) !lake passive
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWLMX(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHAR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKRAD(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFLX(KLON,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0M(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0H(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P10NU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P10NV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2M(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZINV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSRFLTI(KLON,KTILES)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEVAPSNW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTE0(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQE0(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOM0(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVOL0(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVOM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVOL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKE1(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTSKTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQ(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGPP(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCREC(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAG(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRECO(KLON)


REAL(KIND=JPRB) ::    ZVDIS(KLON)
REAL(KIND=JPRB) ::    ZCPTGZ(KLON,KLEV), ZCFM(KLON,KLEV) , ZCFH(KLON,KLEV),&
                    & ZUDIF(KLON,KLEV) , ZVDIF(KLON,KLEV) ,&
                    & ZTDIF(KLON,KLEV) , ZQDIF(KLON,KLEV) ,ZTOFDC(KLON,KLEV)  
REAL(KIND=JPRB) ::    ZKHFL(KLON)       , ZKQFL(KLON)       , ZKMFL(KLON)  
REAL(KIND=JPRB) ::    ZQEA(KLON,KLEV), ZTEA(KLON,KLEV), ZTEWODIS(KLON,KLEV),&
                    & ZUEA(KLON,KLEV), ZVEA(KLON,KLEV)
REAL(KIND=JPRB) ::    ZZ0MW(KLON),  ZZ0HW(KLON), ZZ0QW(KLON)
REAL(KIND=JPRB) ::    ZZCPTS(KLON), ZZQSA(KLON), ZZBUOM(KLON)
          
REAL(KIND=JPRB) ::    ZCPTSTI(KLON,KTILES),&
                    & ZQSTI(KLON,KTILES)  , ZDQSTI(KLON,KTILES),&
                    & ZCSATTI(KLON,KTILES), ZCAIRTI(KLON,KTILES), &
                    & ZCFHTI(KLON,KTILES) , ZCFQTI(KLON,KTILES),&
                    & ZAHFLTI(KLON,KTILES),&
                    & ZTSKTIP1(KLON,KTILES),ZQSTIP1(KLON,KTILES),&
                    & ZCPTSTIP1(KLON,KTILES)  
REAL(KIND=JPRB) ::    ZSTR(KLON,KTILES),ZG0(KLON,KTILES)
REAL(KIND=JPRB) ::    ZZINV(KLON), ZKHVFL(KLON)
!LOP REAL(KIND=JPRB) ::    ZPBLH(KLON)
REAL(KIND=JPRB) ::    ZCFLXIN(KLON,KTRAC), ZCO2FLUX(KLON), ZCH4FLUX(KLON)

REAL(KIND=JPRB) :: ZQTM1, ZSLGM1
REAL(KIND=JPRB) :: ZCLDBASE(KLON), ZEIS(KLON)
REAL(KIND=JPRB) :: ZSLGSM1(KLON,KLEV), ZRI(KLON,KLEV)
INTEGER(KIND=JPIM) :: IINV(KLON), ICBASE(KLON), ICTOP(KLON), IPBLTYPE(KLON)

INTEGER(KIND=JPIM) :: ITOP, JK, JL, JEXT

REAL(KIND=JPRB) :: ZCP, ZGDPH, ZTMST
REAL(KIND=JPRB) :: ZRG

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "surfexcdrivers.h"
#include "surfpps.h"

#include "vdfdifhs.intfb.h"
#include "vdfdifms.intfb.h"
#include "vdfexcus.intfb.h"
#include "vdfincrs.intfb.h"
#include "vdftofdcs.intfb.h"
#include "vdfdifcs.intfb.h"
!LOP #include "vdfdpbls.intfb.h"
#include "compo_flux_update.intfb.h"
#include "vdfhghtns.intfb.h"

!     ------------------------------------------------------------------

#include "fcttre.func.h"

!     ------------------------------------------------------------------

!*         1.     INITIALIZE CONSTANTS
!                 ---------- ---------

IF (LHOOK) CALL DR_HOOK('VDFMAINS',0,ZHOOK_HANDLE)
ASSOCIATE(YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YDEPHLI=>YDMODEL%YRML_PHY_SLIN%YREPHLI, &
 & YDPHNC=>YDMODEL%YRML_PHY_SLIN%YRPHNC,YDVDF=>YDMODEL%YRML_PHY_G%YRVDF, &
 & YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP, YDECUMF2=>YDMODEL%YRML_PHY_SLIN%YRECUMF2 &
 & )
ASSOCIATE(LEPPCFLS=>YDEPHLI%LEPPCFLS, &
 & LVDFTRAC=>YDEPHY%LVDFTRAC, YSURF=>YDEPHY%YSURF, &
 & LESURF2=>YDPHNC%LESURF2,PVDIFTS=>YDVDF%RVDIFTS )
 
ZTMST=PTSPHY
ZRG=-1.0_JPRB/RG

! Setup tendencies
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PQE(JL,JK)=PQE0(JL,JK)
    PTE(JL,JK)=PTE0(JL,JK)
    PVOM(JL,JK)=PVOM0(JL,JK)
    PVOL(JL,JK)=PVOL0(JL,JK)
  ENDDO
ENDDO

!          1.1  GET PROVISONAL T+DT WINDS

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQEA(JL,JK)=PQE(JL,JK)
    ZTEA(JL,JK)=PTE(JL,JK)
    ZUEA(JL,JK)=PVOM(JL,JK)
    ZVEA(JL,JK)=PVOL(JL,JK)
  ENDDO
ENDDO

!*         1.2  dry static energy cp(q)*T + gz

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZCPTGZ(JL,JK) = PGEOM1(JL,JK)+PTM1(JL,JK)*RCPD &
     &            * (1.0_JPRB+RVTMP2*PQM1(JL,JK))
  ENDDO
ENDDO

!  Update total flux of compo tracers with fluxes computed in land surface model

IF (LVDFTRAC .AND. KTRAC > 0) THEN 
  ! make local copy 
  DO JEXT=1,KTRAC
    DO JL=KIDIA,KFDIA
      ZCFLXIN(JL,JEXT) = PCFLX(JL,JEXT)
      ZCO2FLUX(JL) = 0.0_JPRB
      ZCH4FLUX(JL) = 0.0_JPRB
    ENDDO
  ENDDO 
  CALL COMPO_FLUX_UPDATE(YDMODEL%YRML_CHEM%YRCHEM,YDMODEL%YRML_GCONF,YDEPHY,&
         &               KIDIA,KFDIA,KLON,KTRAC,&
         &               PCVL,PCVH,ZCO2FLUX,ZCH4FLUX,PAG,PRECO,PCGPP,PCREC,ZCFLXIN)
ENDIF


!*         2.  Compute all surface related quantities
!          ------------------------------------------

CALL SURFEXCDRIVERS( YDSURF=YSURF,&
 & KIDIA=KIDIA, KFDIA=KFDIA, KLON=KLON, KLEVS=KLEVS, KTILES=KTILES,&
 & KSTEP=NSTEP,&
 & PTSTEP=PTSPHY, PRVDIFTS=PVDIFTS,&
 & LDSURF2=LESURF2,&
! input data, non-tiled
 & KTVL=KTVL, KTVH=KTVH, PCVL=PCVL, PCVH=PCVH,&
 & PLAIL=PLAIL, PLAIH=PLAIH, &
 & PUMLEV=PUM1(:,KLEV), PVMLEV=PVM1(:,KLEV), PTMLEV=PTM1(:,KLEV), &
 & PQMLEV=PQM1(:,KLEV), PAPHMS=PAPHM1(:,KLEV), PGEOMLEV=PGEOM1(:,KLEV), &
 & PCPTGZLEV=ZCPTGZ(:,KLEV), PSST=PSST, PTSKM1M=PTSKM1M, PCHAR=PCHAR, &
 & PSSRFL=PSSRFL, PTICE=PTICE, PTSNOW=PTSNOW, PWLMX=PWLMX, &
! input data, soil
 & PTSAM1M=PTSAM1M, PWSAM1M=PWSAM1M, KSOTY=KSOTY,&
! input data, tiled
 & PFRTI=PFRTI, PALBTI=PALBTI, &
! updated data, tiled
 & PUSTRTI=PUSTRTI, PVSTRTI=PVSTRTI, PAHFSTI=PAHFSTI, PEVAPTI=PEVAPTI, &
 & PTSKTI=PTSKTI, &
! updated data, non-tiled
 & PZ0M=PZ0M, PZ0H=PZ0H, &
! output data, tiled
 & PSSRFLTI=PSSRFLTI, PQSTI=ZQSTI, PDQSTI=ZDQSTI, PCPTSTI=ZCPTSTI, &
 & PCFHTI=ZCFHTI, PCFQTI=ZCFQTI, PCSATTI=ZCSATTI, PCAIRTI=ZCAIRTI, &
! output data, non-tiled
 & PCFMLEV=ZCFM(:,KLEV), PKMFL=ZKMFL, PKHFL=ZKHFL, PKQFL=ZKQFL, &
 & PEVAPSNW=PEVAPSNW, &
 & PZ0MW=ZZ0MW, PZ0HW=ZZ0HW, PZ0QW=ZZ0QW, PCPTSPP=ZZCPTS, PQSAPP=ZZQSA, &
 & PBUOMPP=ZZBUOM &
 & )

!     ------------------------------------------------------------------

!*         3.   NEW VARIABLES S, SLG, QT
!*              (at initial time level)
!               ------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    !* total water and generalized liquid water static energy:
    !*       slg = cp*T + gz - Lcond*ql
    ZSLGM1 = ZCPTGZ(JL,JK) - RLVTT * PLM1(JL,JK) - RLSTT * PIM1(JL,JK) 
    ZQTM1  = PQM1(JL,JK) + PLM1(JL,JK) + PIM1(JL,JK)
    ZSLGSM1(JL,JK) = ZSLGM1 + RCPD * PTM1(JL,JK) * 5.87_JPRB * ZQTM1
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         5.     EXCHANGE COEFFICIENTS
!                 ---------------------

!LOP !*         5.5  BOUNDARY LAYER HEIGHT FOR DIANOSTICS ONLY

!LOP CALL VDFDPBL(KIDIA,KFDIA,KLON,KLEV,&
!LOP  & PUM1,PVM1,PTM1,PQM1,PGEOM1,&
!LOP  & ZKMFL,ZKHFL,ZKQFL,ZPBLH)

!LOP CALL VDFDPBLS(YDEPHY,KIDIA,KFDIA,KLON,KLEV,&
!LOP  & PUM1,PVM1,PTM1,PQM1,PGEOM1,ZKHFL,ZKQFL,ZPBLH)
!LOP 
!LOP DO JL=KIDIA,KFDIA
!LOP   ZKHVFL(JL)  = ZKHFL(JL) + RETV * PTM1(JL,KLEV) * ZKQFL(JL) !w'theta,v'
!LOP   IF ( ZKHVFL(JL) >= 0.0_JPRB ) THEN
!LOP     ZZINV(JL) = 0.0_JPRB   ! stable, therefore 0 PBL depth
!LOP   ELSE
!LOP     ZZINV(JL) = ZPBLH(JL)
!LOP   ENDIF
!LOP ENDDO

!*         5.6  PARCEL UPDRAFT

CALL VDFHGHTNS (YDEPHLI, YDECLDP, YDECUMF2, YDVDF,&
             &  KIDIA    , KFDIA   , KLON    , KLEV    , ZTMST ,&
             &  PTM1     , PQM1    , ZSLGSM1 , PUM1    , PVM1  , ZCPTGZ ,&
             &  PAPHM1   , PAPM1   , PGEOM1  , PGEOH   ,&
             &  ZKMFL    , ZKHFL   , ZKQFL   , ZKHVFL  ,&
             &  PZINV    , IINV    , ICBASE  , ICTOP   ,&
             &  ZCLDBASE , ZRI     , ZEIS    , IPBLTYPE)

!*         5.7  EXCHANGE COEFFICIENTS ABOVE THE SURFACE LAYER

CALL VDFEXCUS(YDEPHLI, YDEPHY, &
 &            KIDIA ,  KFDIA ,  KLON , KLEV,&
 &            IINV  ,  ICBASE,  ICTOP, IPBLTYPE, &
 &            ZTMST ,  PVDIFTS, PZ0M, &
 &            PUM1  ,  PVM1  ,  PTM1 ,  PQM1  , &
 &            PAPHM1,  PGEOM1,  PGEOH,  ZCPTGZ, &
 &            ZKMFL ,  ZKHFL ,  ZKQFL,  ZZINV ,  ZCFM,  ZCFH)
 

!*         5.8     TURBULENT OROGRAPHIC DRAG COEFFICIENTS 

CALL VDFTOFDCS(YDVDF,KIDIA,KFDIA,KLON,KLEV,ZTMST,PVDIFTS,&
 & PUM1,PVM1,PGEOM1,PSIGFLT,&
 & ZTOFDC)  

!     ------------------------------------------------------------------
!*         6.     SOLVE DIFFUSION EQUATION.
!                 ----- --------- ---------

!          6.1    MOMENTUM

ITOP=1

CALL VDFDIFMS (KIDIA, KFDIA, KLON, KLEV, ITOP,&
 & ZTMST ,PVDIFTS, PUM1 ,PVM1  ,PAPHM1 , ZCFM, ZTOFDC,&
 & PVOM ,PVOL ,ZUDIF ,ZVDIF) 

!          6.2    DRY STATIC ENERGY AND MOISTURE

CALL VDFDIFHS (YDMODEL%YRML_AOC%YRMCC,YDEPHY,YDMODEL%YRML_GCONF%YRRIP, &
 & KIDIA ,KFDIA   ,KLON     ,KLEV  , KLEVSN,ITOP  ,KTILES,&
 & KTVL,KTVH,ZTMST  , PVDIFTS, PFRTI   ,PSSRFLTI,PSLRFL ,PEMIS  ,PEVAPSNW,&
! & PHLICE , PTLICE , PTLWML , & !Lake passive
 & ZCPTGZ  ,PTM1     ,PQM1     ,PAPHM1 ,&
 & ZCFH    ,ZCFHTI   ,ZCFQTI   ,&
 & PTHKICE, PSNTICE, &
 & ZTDIF   ,ZQDIF    ,ZCPTSTI  ,ZQSTI  ,ZCAIRTI  ,ZCSATTI,&
 & ZDQSTI  ,PTSKTI   ,PTSKRAD  ,&
 & PTSAM1M(1,1),PTSNOW,PSNS    ,PRSN   ,PTICE    ,PSST,&
 & ZTSKTIP1,ZQSTIP1  ,ZCPTSTIP1,PTE ,PQE  ,&
 & PEVAPTI ,PAHFSTI  ,ZAHFLTI  ,ZSTR ,ZG0)

!          6.3     INCREMENTATION OF U AND V TENDENCIES, STORAGE OF
!                  THE DISSIPATION, COMPUTATION OF MULTILEVEL FLUXES.

CALL VDFINCRS (KIDIA ,KFDIA   ,KLON   ,KLEV   ,ITOP   ,ZTMST ,PVDIFTS,&
 & PUM1   ,PVM1  , PTM1   ,PQM1  ,PAPHM1, &
 & PGEOM1 ,ZCFM  , ZTOFDC ,ZCPTGZ ,&
 & ZUDIF  ,ZVDIF , ZTDIF  ,ZQDIF , &
 & PVOM   ,PVOL  , PTE    ,PQE   ,&
 & ZTEWODIS, ZVDIS ,PSTRTU, PSTRTV)

!          6.4  Solve for tracers

IF (LVDFTRAC .AND. KTRAC > 0) THEN 
  CALL VDFDIFCS(KIDIA,KFDIA,KLON,KLEV,ITOP,KTRAC,&
              & ZTMST,PCM1,PTENC,PAPHM1,ZCFH,ZCFLXIN)
ENDIF

!     -----------------------------------------------------------------

!*         7.     SURFACE FLUXES - TILES
!                 ----------------------

! Wrap-up computations for the surface,
! compute weather parameters, and skin T tendency.

CALL SURFPPS( YSURF,KIDIA=KIDIA,KFDIA=KFDIA,KLON=KLON,KTILES=KTILES, &
 & LDPPCFLS=LEPPCFLS,PTSTEP=PTSPHY, &
! input
 & PFRTI=PFRTI, PTSKTIP1=ZTSKTIP1, PAHFLTI=ZAHFLTI, PG0TI=ZG0, &
 & PSTRTULEV=PSTRTU(:,KLEV), PSTRTVLEV=PSTRTV(:,KLEV), PTSKM1M=PTSKM1M, &
 & PUMLEV=PUM1(:,KLEV), PVMLEV=PVM1(:,KLEV), PQMLEV=PQM1(:,KLEV), &
 & PGEOMLEV=PGEOM1(:,KLEV), PCPTSPP=ZZCPTS, PCPTGZLEV=ZCPTGZ(:,KLEV), &
 & PAPHMS=PAPHM1(:,KLEV), PZ0MW=ZZ0MW, PZ0HW=ZZ0HW, PZ0QW=ZZ0QW, &
 & PQSAPP=ZZQSA, PBUOM=ZZBUOM, &
! updated
 & PAHFSTI=PAHFSTI, PEVAPTI=PEVAPTI, PTSKE1=PTSKE1, &
! output
 & PDIFTSLEV=PDIFTS(:,KLEV), PDIFTQLEV=PDIFTQ(:,KLEV), PUSTRTI=PUSTRTI, &
 & PVSTRTI=PVSTRTI,  PTSKTI=PTSKTI, &
 & PU10M=PU10M, PV10M=PV10M, PT2M=PT2M, PD2M=PD2M, PQ2M=PQ2M, &
 & P10NU=P10NU, P10NV=P10NV &
 & )

!     ------------------------------------------------------------------

!*         9.     FLUX COMPUTATIONS
!                 -----------------

DO JL = KIDIA, KFDIA
  PDIFTQ(JL,0) = 0.0_JPRB
  PDIFTS(JL,0) = 0.0_JPRB
  PSTRTU(JL,0) = 0.0_JPRB
  PSTRTV(JL,0) = 0.0_JPRB
ENDDO

DO JK = KLEV, 2, -1
  DO JL = KIDIA, KFDIA
    ZGDPH=(PAPHM1(JL,JK-1)-PAPHM1(JL,JK))*ZRG
!...change of dry static energy is converted to flux of d.s.e.
    ZCP=RCPD*(1.0_JPRB+RVTMP2*PQM1(JL,JK))
    PDIFTS(JL,JK-1)=ZCP*(ZTEWODIS(JL,JK)-ZTEA(JL,JK))*&
     & ZGDPH+PDIFTS(JL,JK)  
!...changes in Q,U,V tendencies are converted to fluxes
    PDIFTQ(JL,JK-1)=(PQE (JL,JK)-ZQEA(JL,JK))*ZGDPH+PDIFTQ(JL,JK)
    PSTRTU(JL,JK-1)=(PVOM(JL,JK)-ZUEA(JL,JK))*ZGDPH+PSTRTU(JL,JK)
    PSTRTV(JL,JK-1)=(PVOL(JL,JK)-ZVEA(JL,JK))*ZGDPH+PSTRTV(JL,JK)
  ENDDO
ENDDO

! Extract tendencies
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PQE(JL,JK)=PQE(JL,JK)-PQE0(JL,JK)
    PTE(JL,JK)=PTE(JL,JK)-PTE0(JL,JK)
    PVOM(JL,JK)=PVOM(JL,JK)-PVOM0(JL,JK)
    PVOL(JL,JK)=PVOL(JL,JK)-PVOL0(JL,JK)
  ENDDO
ENDDO

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VDFMAINS',1,ZHOOK_HANDLE)
END SUBROUTINE VDFMAINS

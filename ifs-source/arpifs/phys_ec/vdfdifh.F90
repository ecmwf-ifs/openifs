! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE VDFDIFH(&
 & YDMCC  , YDEPHY , YDRIP,&
 & KIDIA  , KFDIA  , KLON   , KLEV   , KTOP   , KTILES, &
 & KTVL   , KTVH   , PTMST  , PFSH1D , PFLH1D , LDFLUX1D, &
 & PFRTI  , PSSRFLTI,PSLRFL , PEMIS  , PEVAPSNW, &
 & PHLICE , PTLICE , PTLWML , &
 & PSLM1  , PQTM1  , PAPHM1 , &
 & PCFH   , PCFHTI , PCFQTI , PMFLX  , PSLUH  , PQTUH  , &
 & PTDIF  , PQDIF  , PCPTSTI, PQSTI  , PCAIRTI, PCSATTI, &
 & PCPTSTIU,PCAIRTIU,PCSATTIU, PTSRF ,PLAMSK , &
 & PDQSTI , PTSKTI , PTSKRAD, PTSM1M , &
 & PTSNOW , PSNM, PRSN, PTICE  , PSST, &
 & PTHKICE, PSNTICE, &
 & PTSKTIH ,PSLGE  , PQTE   , PSLGEVDF, PQTEVDF, PTSKTITE,&
 & PJQ    , PSSH   , PSLH   , PSTR   , PG0, PJQU)  
!     ------------------------------------------------------------------

!**   *VDFDIFH* - DOES THE IMPLICIT CALCULATION FOR DIFFUSION OF S. L.

!     DERIVED FROM VDIFF (CY34) BY
!     A.C.M. BELJAARS       E.C.M.W.F.    10-11-89

!     OBUKHOV-L UPDATE      ACMB          26/03/90.
!     SKIN T CLEANING       P. VITERBO    15-11-96.
!     TILE BOUNDARY COND.   ACMB          20-11-98.
!     Surface DDH for TILES P. Viterbo    17-05-2000.
!     New tile coupling, 
!     DDH moved to VDFMAIN  A. Beljaars   2-05-2003.
!     Mass flux terms,
!     Flux b.c. for SCM,
!     Moist generalization  A. Beljaars/M. Ko"hler 3-12-2004.    
!     Multiple mass fluxes  M. Ko"hler               05-2005.
!     Removed option for linearized    P. Lopez   02/06/2005
!     physics (now called separately)   
!     Add lake tile                    G. Balsamo 18/04/2008
!     Offline Jacobians EKF P.de Rosnay/G.Balsamo 07/03/2009

!     PURPOSE
!     -------

!     SOLVE TRIDIAGONAL MATRICES FOR DIFFUSION OF DRY STATIC ENERGY
!     AND MOISTURE; IN SO DOING, IT ALSO SOLVES THE SKIN TEMPERATURE
!     EQUATION AND THE SURFACE ENERGY BALANCE. 
!     The time step can contain the over implicit factor for 
!     extrapolation in time, or the factor for the predictor corrector 
!     method. The output is the transformed variable (e.g. the extrapolated 
!     in time or the hat parameters in the documentation). Output also includes 
!     tendencies and fluxes.    

!     INTERFACE
!     ---------

!     *VDFDIFH* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KLEV*         NUMBER OF LEVELS
!     *KTOP*         INDEX FOR BOUNDARY LAYER TOP

!     INPUT PARAMETERS (REAL):

!     *PTMST*        TIME STEP including possible extrapolation in time
!     *PFRTI*        FRACTION OF SURFACE AREA COVERED BY TILES
!     *PSLM1*        GENERALIZED LIQUID WATER STATIC ENERGY    AT T-1
!                    (NOTE: In lin/adj physics = Dry static energy)
!     *PQTM1*        SPECIFIC TOTAL WATER   AT T-1
!     *PAPHM1*       PRESSURE AT T-1
!     *PCFH*         Rho*K/dz (C,K-STAR IN DOC.)
!     *PCFHTI*       IDEM FOR HEAT (SURFACE LAYER ONLY)
!     *PCFQTI*       IDEM FOR MOISTURE (SURFACE LAYER ONLY)
!     *PMFLX*        MASSFLUX (kgm-2s-1)
!     *PSLUH*        UPDRAFT GENERALIZED LIQUID WATER STATIC ENERGY AT HALF LEVEL
!     *PQTUH*        UPDRAFT SPECIFIC TOTAL WATER AT HALF LEVEL
!     *PCPTSTI*      DRY STATIC ENRGY AT SURFACE
!     *PQSTI*        SATURATION Q AT SURFACE
!     *PCAIRTI*      MULTIPLICATION FACTOR FOR Q AT LOWEST MODEL LEVEL
!                    FOR SURFACE FLUX COMPUTATION
!     *PCSATTI*      MULTIPLICATION FACTOR FOR QS AT SURFACE
!                    FOR SURFACE FLUX COMPUTATION
!     *PCPTSTIU*     AS PCPTSTI FOR UNSTRESSED EVAPORARTION FROM LOW VEGET.
!     *PCAIRTIU*     AS PCAIRTI FOR UNSTRESSED EVAPORARTION FROM LOW VEGET.
!     *PCSATTIU*     AS PCSATTI FOR UNSTRESSED EVAPORARTION FROM LOW VEGET.
!     *PDQSTI*       D/DT (PQS)
!     *PSSRFLTI*     NET SOLAR RADIATION AT THE SURFACE, FOR EACH TILE
!     *PSLRFL*       NET THERMAL RADIATION AT THE SURFACE
!     *PEMIS*        MODEL SURFACE LONGWAVE EMISSIVITY
!     *PEVAPSNW*     EVAPORATION FROM SNOW UNDER FOREST
!     *PTLICE*       LAKE ICE TEMPERATURE [K]
!     *PHLICE*       LAKE ICE THICKNESS   [m]
!     *PTLWML*       LAKE MEAN WATER TEMPERATURE [K]
!     *PTSKTI*       SKIN TEMPERATURE AT T-1
!     *PTSKRAD*      SKIN TEMPERATURE OF LATEST FULL RADIATION TIMESTEP
!     *PTSM1M*       TOP SOIL LAYER TEMPERATURE
!     *PTSNOW*       SNOW TEMPERATURE 
!     *PSNS*         SNOW MASS                                     
!     *PRSN*         SNOW DENSITY
!     *PTICE*        ICE TEMPERATURE (TOP SLAB)
!     *PSST*         (OPEN) SEA SURFACE TEMPERATURE
!     *PTHKICE*      SEA-ICE THICKNESS
!     *PSNTICE*      SNOW THICKNESS ON SEA-ICE
!     *PSLGE*        GENERALIZED DRY STATIC ENERGY TENDENCY
!                    (NOTE: In lin/adj physics = Temperature tendency)
!     *PQTE*         TOTAL WATER TENDENCY
!                    (NOTE: In lin/adj physics = Humidity tendency)

!     OUTPUT PARAMETERS (REAL):

!     *PTDIF*        SLG-hat (SLG at new time level)
!     *PQDIF*        QT-hat  (SLG at new time level)
!     *PTSKTIH*      Tsk-hat (Tsk at new time level)
!     *PSLGEVDF*     Generalized dry static energy tendency from VDF
!     *PQTEVDF*      Total water tendency from VDF
!      PTSKTITE      Skin temperature tendency
!     *PJQ*          Surface moisture flux                      (kg/m2s)
!     *PSSH*         Surface sensible heat flux                 (W/m2)
!     *PSLH*         Surface latent heat flux                   (W/m2)
!     *PSTR*         Surface net thermal radiation              (W/m2)
!     *PG0*          Surface ground heat flux (solar radiation  (W/m2)
!                    leakage is not included in this term)
!     *PJQU*         Surface moisture flux unstressed low veg   (kg/m2s)

!     Additional parameters for flux boundary condtion (in 1D model):

!     *LDFLUX1D*     If .TRUE. flux boundary condtion is used 
!     *PFSH1D*       Specified sensible heat flux (W/m2)
!     *PFLH1D*       Specified latent heat flux (W/m2)

!     METHOD
!     ------

!     *LU*-DECOMPOSITION (DOWNWARD SCAN), FOLLOWED BY SKIN-TEMPERATURE
!     SOLVER, AND BACK SUBSTITUTION (UPWARD SCAN).

!     EXTERNALS.
!     ----------

!     *VDFDIFH* CALLS:
!         *SURFSEB*

!     AUTHOR.
!     -------
!      A. BELJAARS           10-11-1989  DERIVED FROM VDIFF (CY34)

!     MODIFICATIONS.
!     --------------
!      A. Beljaars            2-05-2003  DDH moved to VDFMAIN
!      M. Ko"hler                        Mass flux terms,
!      M. Ko"hler                        Flux b.c. for SCM,
!      A. Beljaars/M. Ko"hler 3-12-2004  Moist generalization    
!      M. Ko"hler              May-2005  Multiple mass fluxes  M. Ko"hler
!      P. Lopez              02/06/2005  Removed option for linearized physics
!                                        (now called separately)
!      G. Balsamo            18/04/2008  Add lake tile
!      P.de Rosnay/G.Balsamo 07/03/2009  Offline Jacobians EKF
!      P.de Rosnay/G.Balsamo   Oct-2009  Offline Jacobians EKF commented
!      N. Semane/P.Bechtold  04-10-2012  Add RPLDARE factor for small planet
!      Snow scheme stability fix  08-08-2013 G. Balsamo/A. Beljaars 
!      I. Sandu    24-02-2014  Lambda skin values by vegetation type instead of tile
!      A. Beljaars             Jan 2014  Clean-up for flexible numerics
!      Compute unstressed evaporation  A. Beljaars 26/02/2010
!      E. Dutra    October 2016  Lambda skin / PTSRF as input array 
!      E. Dutra    July 2017: 2nd call to surface energy balance to limit skinT of snow for multi-layer
!      P. Bechtold           09-01-2019  Cleaning and simplifying of dry mass flux
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RCPD, RG, RLSTT, RLVTT, RTT
USE YOETHF    , ONLY : RVTMP2
USE YOMDYNCORE, ONLY : RPLDARE
USE YOEPHY    , ONLY : TEPHY
USE YOMMCC    , ONLY : TMCC
USE YOMRIP    , ONLY : TRIP

IMPLICIT NONE

#include "surfseb.h"

TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
TYPE(TMCC)        ,INTENT(INOUT) :: YDMCC
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFSH1D(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFLH1D(KLON) 
LOGICAL           ,INTENT(IN)    :: LDFLUX1D 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFLTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVAPSNW(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLWML(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFHTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFQTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFLX(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLUH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTUH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQDIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCAIRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSATTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTSTIU(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCAIRTIU(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSATTIU(KLON,KTILES)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSRF(KLON,KTILES)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAMSK(KLON,KTILES)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDQSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKRAD(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSM1M(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSNOW(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSN(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTICE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSST(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTHKICE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNTICE(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSKTIH(KLON,KTILES)  
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLGE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQTE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLGEVDF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQTEVDF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSKTITE(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PJQ(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSH(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLH(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTR(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PG0(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PJQU(KLON,KTILES)


REAL(KIND=JPRB) ::    ZMSL(KLON,0:KLEV), ZMQT(KLON,0:KLEV)
REAL(KIND=JPRB) ::    ZAA(KLON,KLEV) ,ZBB(KLON,KLEV) ,ZCC(KLON,KLEV) ,&
                    & ZTYY(KLON,KLEV),ZQYY(KLON,KLEV),ZGAM(KLON,KLEV),&
                    & Z1DP(KLON,KLEV)
REAL(KIND=JPRB) ::    Z1BET(KLON)    ,&
                    & ZAQL(KLON)     ,ZBQL(KLON)     ,ZASL(KLON)     ,&
                    & ZBSL(KLON)     ,ZSL(KLON)      ,ZQL(KLON)
REAL(KIND=JPRB) ::    ZTSRF(KLON,KTILES)  ,& 
                    & ZJS(KLON,KTILES)      ,&
                    & ZSSK(KLON,KTILES) , ZLAMSK(KLON,KTILES) 

INTEGER(KIND=JPIM) :: JK, JL, JT
LOGICAL         :: LREPEAT

REAL(KIND=JPRB) ::    ZQDP, ZCSNQ, ZCSNS, ZCONS1
REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('VDFDIFH',0,ZHOOK_HANDLE)
ASSOCIATE(LEFLAKE=>YDEPHY%LEFLAKE, LEOCCO=>YDEPHY%LEOCCO, LEOCWA=>YDEPHY%LEOCWA, &
 & YSURF=>YDEPHY%YSURF, &
 & LNEMOLIMTHK=>YDMCC%LNEMOLIMTHK, &
 & TSTEP=>YDRIP%TSTEP,&
 & LESNML =>YDEPHY%LESNML)

ZCONS1=RG*PTMST

!     ------------------------------------------------------------------

!*         1.     FULL MODEL PHYSICS WITH MOIST MASS FLUX PBL
!                 -------------------------------------------

!*         1.0    PRECALCULATION OF MULTIPLE MASS-FLUX TERMS

DO JK=0,KLEV
  DO JL=KIDIA,KFDIA
    ZMSL (JL,JK) = PMFLX(JL,JK)*PSLUH(JL,JK)
    ZMQT (JL,JK) = PMFLX(JL,JK)*PQTUH(JL,JK)
  ENDDO
ENDDO

!*         1.1    SETTING OF THE MATRIX A, B AND C.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    Z1DP(JL,JK)=ZCONS1/(PAPHM1(JL,JK)-PAPHM1(JL,JK-1))
    ZAA(JL,JK) =(-PCFH(JL,JK-1)-PMFLX(JL,JK-1))*Z1DP(JL,JK)
    ZCC(JL,JK) =(-PCFH(JL,JK)                 )*Z1DP(JL,JK)
    ZBB(JL,JK) =1.0_JPRB+(PCFH(JL,JK-1)+PCFH(JL,JK)&
     & +PMFLX(JL,JK))*Z1DP(JL,JK)  
  ENDDO
ENDDO

!          1.1a   THE SURFACE BOUNDARY CONDITION

DO JL=KIDIA,KFDIA
  Z1DP(JL,KLEV)=ZCONS1/(PAPHM1(JL,KLEV)-PAPHM1(JL,KLEV-1))
  ZCC(JL,KLEV) =0.0_JPRB
  ZAA(JL,KLEV) =        (-PCFH(JL,KLEV-1)-PMFLX(JL,KLEV-1))*Z1DP(JL,KLEV)
  ZBB(JL,KLEV) =1.0_JPRB+(PCFH(JL,KLEV-1)                 )*Z1DP(JL,KLEV)  
ENDDO

!          1.1b   THE TOP BOUNDARY CONDITION    

DO JL=KIDIA,KFDIA
  Z1DP(JL,KTOP)=ZCONS1/(PAPHM1(JL,KTOP)-PAPHM1(JL,KTOP-1))
  ZAA(JL,KTOP) =0.0_JPRB
  ZCC(JL,KTOP) =         (-PCFH(JL,KTOP))*Z1DP(JL,KTOP)
  ZBB(JL,KTOP) =1.0_JPRB+( PCFH(JL,KTOP)+PMFLX(JL,KTOP))*Z1DP(JL,KTOP)
ENDDO

!*         1.2    SETTING OF RIGHT HAND SIDES.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZTYY(JL,JK) = PSLM1(JL,JK) &
     & + PTMST * PSLGE(JL,JK) &
     & + (ZMSL(JL,JK)-ZMSL(JL,JK-1)) * Z1DP(JL,JK)
    ZQYY(JL,JK) = PQTM1(JL,JK) &
     & + PTMST * PQTE(JL,JK) &
     & + (ZMQT(JL,JK)-ZMQT(JL,JK-1)) * Z1DP(JL,JK)
  ENDDO
ENDDO

!          1.2a   SURFACE

JK=KLEV
DO JL=KIDIA,KFDIA
  ZTYY(JL,JK) = PSLM1(JL,JK) &
   & + PTMST * PSLGE(JL,JK) &
   & - ZMSL(JL,JK-1) * Z1DP(JL,JK)
  ZQYY(JL,JK) = PQTM1(JL,JK) &
   & + PTMST * PQTE(JL,JK) &
   & - ZMQT(JL,JK-1) * Z1DP(JL,JK)
ENDDO

!          1.2b   TOP

JK=KTOP
DO JL=KIDIA,KFDIA
  ZTYY(JL,JK) = PSLM1(JL,JK) &
   & + PTMST * PSLGE(JL,JK) &
   & + ZMSL(JL,JK) * Z1DP(JL,JK)
  ZQYY(JL,JK) = PQTM1(JL,JK) &
   & + PTMST * PQTE(JL,JK) &
   & + ZMQT(JL,JK) * Z1DP(JL,JK)
ENDDO

!*         1.3    ADD MOISTURE FLUX FROM SNOW FROM TILE 7 AS EXPLICIT TERM

JK=KLEV
DO JL=KIDIA,KFDIA
  ZCSNQ=PFRTI(JL,7)*PEVAPSNW(JL)*Z1DP(JL,JK)
  ZCSNS=RCPD*RVTMP2*PTSKTI(JL,7)*ZCSNQ
  ZTYY(JL,JK)=ZTYY(JL,JK)-ZCSNS
  ZQYY(JL,JK)=ZQYY(JL,JK)-ZCSNQ
ENDDO

!*         1.4    TOP LAYER ELIMINATION.

DO JL=KIDIA,KFDIA
  Z1BET(JL)=1.0_JPRB/ZBB(JL,KTOP)
  PTDIF(JL,KTOP)=ZTYY(JL,KTOP)*Z1BET(JL)
  PQDIF(JL,KTOP)=ZQYY(JL,KTOP)*Z1BET(JL)
ENDDO

!*         1.5    ELIMINATION FOR MIDDLE LAYERS.

DO JK=KTOP+1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZGAM(JL,JK)=ZCC(JL,JK-1)*Z1BET(JL)
    Z1BET(JL)=1.0_JPRB/(ZBB(JL,JK)-ZAA(JL,JK)*ZGAM(JL,JK))
    PTDIF(JL,JK)=(ZTYY(JL,JK)-ZAA(JL,JK)*PTDIF(JL,JK-1))*Z1BET(JL)
    PQDIF(JL,JK)=(ZQYY(JL,JK)-ZAA(JL,JK)*PQDIF(JL,JK-1))*Z1BET(JL)
  ENDDO
ENDDO

!*         1.6    BOTTOM LAYER, LINEAR RELATION BETWEEN LOWEST
!                 MODEL LEVEL S AND Q AND FLUXES.

DO JL=KIDIA,KFDIA
  ZGAM(JL,KLEV)=ZCC(JL,KLEV-1)*Z1BET(JL)
  Z1BET(JL)=1.0_JPRB/(ZBB(JL,KLEV)-ZAA(JL,KLEV)*ZGAM(JL,KLEV))
  ZBSL(JL)=(ZTYY(JL,KLEV)-ZAA(JL,KLEV)*PTDIF(JL,KLEV-1))*Z1BET(JL)
  ZBQL(JL)=(ZQYY(JL,KLEV)-ZAA(JL,KLEV)*PQDIF(JL,KLEV-1))*Z1BET(JL)
  ZQDP=1.0_JPRB/(PAPHM1(JL,KLEV)-PAPHM1(JL,KLEV-1))
  ZASL(JL)=-RG*RPLDARE*PTMST*ZQDP*Z1BET(JL)
  ZAQL(JL)=ZASL(JL)
ENDDO

!*         1.7    PREPARE ARRAY'S FOR CALL TO SURFACE ENERGY
!                 BALANCE ROUTINE


! E. Dutra:Oct 2016
! We keep ZTSRF association here but it's no longer used
! in surfseb, instead we use PTSRF that comes from surfexcdriver. 
! These lines can be removed latter 
IF (LEOCWA .OR. LEOCCO) THEN
  ZTSRF(KIDIA:KFDIA,1)=PTSKTI(KIDIA:KFDIA,1)
ELSE
  ZTSRF(KIDIA:KFDIA,1)=PSST(KIDIA:KFDIA)
ENDIF
ZTSRF(KIDIA:KFDIA,2)=PTICE(KIDIA:KFDIA)
ZTSRF(KIDIA:KFDIA,3)=PTSM1M(KIDIA:KFDIA)
ZTSRF(KIDIA:KFDIA,4)=PTSM1M(KIDIA:KFDIA)
ZTSRF(KIDIA:KFDIA,5)=PTSNOW(KIDIA:KFDIA)
ZTSRF(KIDIA:KFDIA,6)=PTSM1M(KIDIA:KFDIA)
ZTSRF(KIDIA:KFDIA,7)=PTSNOW(KIDIA:KFDIA)
ZTSRF(KIDIA:KFDIA,8)=PTSM1M(KIDIA:KFDIA)
IF (LEFLAKE) THEN
  DO JL=KIDIA,KFDIA
    IF(PHLICE(JL) > 1.E-9_JPRB) THEN ! 1.E-9 or H_ICE_MIN_FLK present
      ZTSRF(JL,9)=PTLICE(JL)
    ELSE
      ZTSRF(JL,9)=PTLWML(JL)
    ENDIF
  ENDDO
ENDIF


!*         1.8    CALL TO SURFACE ENERGY BALANCE ROUTINE
!                 REMEMBER: OUTPUT IS EXTRAPOLATED IN TIME

CALL SURFSEB(YDSURF=YSURF,KIDIA=KIDIA,KFDIA=KFDIA,KLON=KLON,&
 & KTILES=KTILES,KTVL=KTVL,KTVH=KTVH,&
 & PTMST=TSTEP,PSSKM1M=PCPTSTIU,PTSKM1M=PTSKTI,PQSKM1M=PQSTI,&
 & PDQSDT=PDQSTI,PRHOCHU=PCFHTI,PRHOCQU=PCFQTI,&
 & PALPHAL=PCAIRTIU,PALPHAS=PCSATTIU,&
 & PSSRFL=PSSRFLTI,PFRTI=PFRTI,PTSRF=PTSRF,PLAMSK=PLAMSK,&
 & PSNM=PSNM,PRSN=PRSN,PHLICE=PHLICE,&
 & PSLRFL=PSLRFL,PTSKRAD=PTSKRAD,PEMIS=PEMIS,&
 & PASL=ZASL,PBSL=ZBSL,PAQL=ZAQL,PBQL=ZBQL,&
 & PTHKICE=PTHKICE,PSNTICE=PSNTICE,&
!out
 & PJS=ZJS,PJQ=PJQU,PSSK=ZSSK,PTSK=PTSKTIH,&
 & PSSH=PSSH,PSLH=PSLH,PSTR=PSTR,PG0=PG0,&
 & PSL=ZSL,PQL=ZQL,LNEMOLIMTHK=LNEMOLIMTHK) 

CALL SURFSEB(YDSURF=YSURF,KIDIA=KIDIA,KFDIA=KFDIA,KLON=KLON,&
 & KTILES=KTILES,PTMST=TSTEP,&
 & KTVL=KTVL,KTVH=KTVH,PSSKM1M=PCPTSTI,PTSKM1M=PTSKTI,PQSKM1M=PQSTI,&
 & PDQSDT=PDQSTI,PRHOCHU=PCFHTI,PRHOCQU=PCFQTI,&
 & PALPHAL=PCAIRTI,PALPHAS=PCSATTI,&
 & PSSRFL=PSSRFLTI,PFRTI=PFRTI,PTSRF=PTSRF,PLAMSK=PLAMSK,&
 & PSNM=PSNM,PRSN=PRSN,PHLICE=PHLICE,&
 & PSLRFL=PSLRFL,PTSKRAD=PTSKRAD,PEMIS=PEMIS,&
 & PASL=ZASL,PBSL=ZBSL,PAQL=ZAQL,PBQL=ZBQL,&
 & PTHKICE=PTHKICE,PSNTICE=PSNTICE,&
 !out
 & PJS=ZJS,PJQ=PJQ,PSSK=ZSSK,PTSK=PTSKTIH,&
 & PSSH=PSSH,PSLH=PSLH,PSTR=PSTR,PG0=PG0,&
 & PSL=ZSL,PQL=ZQL,LNEMOLIMTHK=LNEMOLIMTHK)  

! 2nd call to surface energy balance 
! in case of skkin temperature above freezing point 
! only active with the new multi-layer snow 
IF (LESNML) THEN
  !! Check stuff 
  LREPEAT=.FALSE.
  ZTSRF(KIDIA:KFDIA,:) = PTSRF(KIDIA:KFDIA,:)
  ZLAMSK(KIDIA:KFDIA,:) = PLAMSK(KIDIA:KFDIA,:)
  DO JL=KIDIA,KFDIA
    IF (PFRTI(JL,5) >= 1.E-6_JPRB ) THEN
      IF (PTSKTIH(JL,5) >= RTT-1.E-6_JPRB ) THEN
        ZTSRF(JL,5)=RTT
        ZLAMSK(JL,5)= 100._JPRB 
        LREPEAT=.TRUE.
      ENDIF
    ENDIF
  ENDDO
  IF (LREPEAT) THEN
    ! call again surfseb replacing ZTSRF and ZLAMSK to avoid overshoot of tile 5
    ! skin temperatures
    CALL SURFSEB(YDSURF=YSURF,KIDIA=KIDIA,KFDIA=KFDIA,KLON=KLON,&
     & KTILES=KTILES,PTMST=TSTEP,&
     & KTVL=KTVL,KTVH=KTVH,PSSKM1M=PCPTSTI,PTSKM1M=PTSKTI,PQSKM1M=PQSTI,&
     & PDQSDT=PDQSTI,PRHOCHU=PCFHTI,PRHOCQU=PCFQTI,&
     & PALPHAL=PCAIRTI,PALPHAS=PCSATTI,&
     & PSSRFL=PSSRFLTI,PFRTI=PFRTI,PTSRF=ZTSRF,PLAMSK=ZLAMSK,&
     & PSNM=PSNM,PRSN=PRSN,PHLICE=PHLICE,&
     & PSLRFL=PSLRFL,PTSKRAD=PTSKRAD,PEMIS=PEMIS,&
     & PASL=ZASL,PBSL=ZBSL,PAQL=ZAQL,PBQL=ZBQL,&
     & PTHKICE=PTHKICE,PSNTICE=PSNTICE,&
     !out
     & PJS=ZJS,PJQ=PJQ,PSSK=ZSSK,PTSK=PTSKTIH,&
     & PSSH=PSSH,PSLH=PSLH,PSTR=PSTR,PG0=PG0,&
     & PSL=ZSL,PQL=ZQL,LNEMOLIMTHK=LNEMOLIMTHK)  
  ENDIF
ENDIF 

!*         1.9    ADD SNOW EVAPORATION TO FLUXES

PJQ (KIDIA:KFDIA,7)=PJQ (KIDIA:KFDIA,7)+PEVAPSNW(KIDIA:KFDIA)
PSLH(KIDIA:KFDIA,7)=PSLH(KIDIA:KFDIA,7)+PEVAPSNW(KIDIA:KFDIA)*RLSTT

!*         1.10   Flux boundary condition for 1D model (fluxes in W/m2)
!                 (Over-write output of SURFSEB)

IF (LDFLUX1D) THEN
  DO JT=1,KTILES
    DO JL=KIDIA,KFDIA
      ZJS(JL,JT)=PFSH1D(JL)+RCPD*PTSKTI(JL,JT)*RVTMP2*PFLH1D(JL)/RLVTT
      PJQ(JL,JT)=PFLH1D(JL)/RLVTT

      ZSSK(JL,JT)=ZBSL(JL)-ZJS(JL,JT)*(ZASL(JL)-1.0_JPRB/PCFHTI(JL,JT)) 
      PTSKTIH(JL,JT)=ZSSK(JL,JT)/(RCPD*(1.+RVTMP2*PQSTI(JL,JT)))
      PSSH(JL,JT)=PFSH1D(JL)
      PSLH(JL,JT)=PFLH1D(JL)
      PSTR(JL,JT)=PSLRFL(JL)
      PG0 (JL,JT)=PFSH1D(JL)+PFLH1D(JL)+PSLRFL(JL)+PSSRFLTI(JL,JT)
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    ZSL(JL)=ZJS(JL,1)*ZASL(JL)+ZBSL(JL)
    ZQL(JL)=PJQ(JL,1)*ZAQL(JL)+ZBQL(JL)
  ENDDO
ENDIF

!*         1.11   COMPUTE PARAMETERS AT NEW TIME LEVEL 

!*         1.12   COPY LOWEST MODEL SOLUTION FROM SURFSEB

DO JL=KIDIA,KFDIA
  PTDIF(JL,KLEV)=ZSL(JL)
  PQDIF(JL,KLEV)=ZQL(JL)
ENDDO

!*         1.13   BACK-SUBSTITUTION.

DO JK=KLEV-1,KTOP,-1
  DO JL=KIDIA,KFDIA
    PTDIF(JL,JK)=PTDIF(JL,JK)-ZGAM(JL,JK+1)*PTDIF(JL,JK+1)
    PQDIF(JL,JK)=PQDIF(JL,JK)-ZGAM(JL,JK+1)*PQDIF(JL,JK+1)
  ENDDO
ENDDO

!*         1.6     COMPUTATION OF TENDENCIES.

!-----------------------------------------------------------------------------
! Input: PSLM1    liq. static energy at the start of the time step
!        PSLGE    liq. static energy tendency dynamics
!        PQTM1    total water at the start of the time step
!        PQTE     total water tendency dynamics
!
!Output: PTDIF    liq. static energy at the end the time step including dynamics
!        PQDIF    total water at the end the time step including dynamics
!        PSLGEVDF liq. static energy tendency from diffusion only
!        PQTEVDF  total water tendency from diffusion only
!-----------------------------------------------------------------------------

PSLGEVDF(KIDIA:KFDIA,1:KTOP-1)=0.0_JPRB
PQTEVDF(KIDIA:KFDIA,1:KTOP-1)=0.0_JPRB

DO JK=KTOP,KLEV
  DO JL=KIDIA,KFDIA
    PSLGEVDF(JL,JK)=(PTDIF(JL,JK)-PSLM1(JL,JK))/PTMST - PSLGE(JL,JK)
    PQTEVDF(JL,JK) =(PQDIF(JL,JK)-PQTM1(JL,JK))/PTMST - PQTE(JL,JK)
  ENDDO
ENDDO

DO JT=1,KTILES
  DO JL=KIDIA,KFDIA
    PTSKTITE(JL,JT)=(PTSKTIH(JL,JT)-PTSKTI(JL,JT))/PTMST 
  ENDDO
ENDDO
  
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VDFDIFH',1,ZHOOK_HANDLE)
END SUBROUTINE VDFDIFH


! (C) Copyright 1993- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!     ------------------------------------------------------------------

!**   *SPPCFL* - COMPUTES THE SURFACE (2 M) TEMPERATURE AND HUMIDITY
!                WITH STABILITY AS FUNCTION OF OBUKHOV-L

!     P. VITERBO         E.C.M.W.F.    08/10/93. (MODIFIED FROM VDFT2M)
!     Modified   P. Viterbo ECMWF 12/05/2005 Externalize SURF
!                                   (based on vdfppcfl)
!                H. Hersbach ECMWF 04/12/2009 Add 10m neutral wind and friction velocity
!     N.Semane+P.Bechtold 04-10-2012 Add PRPLRG factor for small planet
!     F. Vana  05-Mar-2015  Support for single precision
!     A. van Niekerk 05/2023 Remove limiter for 10m wind

!     PURPOSE
!     -------

!     COMPUTE WIND COMPONENTS, TEMPERATURE AND DEWPOINT TEMPERATURE
!     AT SCREEN LEVEL HEIGHT

!     INTERFACE
!     ---------

!     *SPPCFL* IS CALLED BY *SURFPP*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PUMLEV*       U-COMPONENT WIND AT T-1, lowest atmospheric level
!     *PVMLEV*       V-COMPONENT WIND AT T-1, lowest atmospheric level
!     *PQMLEV*       SPECIFIC HUMIDITY AT T-1, lowest atmospheric level
!     *PAPHMS*       surface pressure
!     *PGEOMLEV*     GEOPOTENTIAL AT T-1, lowest atmospheric level
!     *PCPTS*        DRY STATIC ENRGY AT SURFACE
!     *PCPTGZLEV*    DRY STATIC ENERGY    AT T-1, lowest atmospheric level
!     *PZ0MM*        AERODYNAMIC ROUGHNESS LENGTH
!     *PZ0MSM*       AERODYNAMIC ROUGHNESS LENGTH FOR SCALARS
!     *PZ0HM*        ROUGHNESS LENGTH FOR TEMPERATURE
!     *PZ0QM*        ROUGHNESS LENGTH FOR MOISTURE
!     *PZDL*         ZNLEV+Z0M DEVIDED BY OBUKHOV LENGTH
!     *PZDL_WIND*    ZNLEV+Z0M DEVIDED BY OBUKHOV LENGTH (for wind calculation)
!     *PQSA*         APPARENT SURFACE HUMIDITY
!     *PBLEND*       HEIGHT FROM WHICH WIND SPEED IS INTERPOLATED TO 10 M
!     *PFBLEND*      WIND SPEED AT HEIGHT PBLEND
!     *PUCURR*       U COMPONENT OF OCEAN SURFACE CURRENT
!     *PVCURR*       V COMPONENT OF OCEAN SURFACE CURRENT

!     OUTPUT PARAMETERS (REAL):

!     *PU10*         U-COMPONENT WIND AT 10 M
!     *PV10*         V-COMPONENT WIND AT 10 M
!     *P10NU*        U-COMPONENT NEUTRAL WIND AT 10 M
!     *P10NV*        V-COMPONENT NEUTRAL WIND AT 10 M
!     *PUST*         FRICTION VELOCITY
!     *PT2*          TEMPERATURE AT 2 M
!     *PD2*          DEW POINT TEMPERATURE AT 2 M
!     *PQ2*          SPECIFIC HUMIDITY AT 2 M

!     METHOD
!     ------

!     UNIVERSAL FUNCTIONS ARE USED TO INTERPOLATE BETWEEN SURF. AND NLEV

!     ------------------------------------------------------------------

MODULE SPPCFL_MOD
CONTAINS
SUBROUTINE SPPCFL(KIDIA, KFDIA, KLON,KTILE &
 & , PUMLEV, PVMLEV, PQMLEV, PGEOMLEV, PCPTS, PCPTGZLEV &
 & , PAPHMS, PZ0MM, PZDL_WIND &
 & , PZ0MSM, PZ0HM, PZ0QM, PZDL, PQSA &
 & , PBLEND, PFBLEND, PUCURR, PVCURR &
 & , YDCST, YDEXC &
 ! OUTPUTS
 & , PU10, PV10, P10NU, P10NV, PUST, PT2, PD2, PQ2, PRPLRG &
 & , LWIND )  

USE PARKIND1  , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF   , ONLY : R4LES, R2ES, R3LES, RVTMP2
USE YOS_EXCS  , ONLY : RCHBCD, RCHBBCD, RCHBB, RCHBD, RCHBA, RCHBHDL, &
 & RCDHALF, RCHETB, RCHB23A, RCHETA, RCDHPI2
USE YOS_CST   , ONLY : TCST
USE YOS_EXC    ,ONLY : TEXC

!     ------------------------------------------------------------------

!**   *SPPCFL* - COMPUTES THE SURFACE (2 M) TEMPERATURE AND HUMIDITY
!                WITH STABILITY AS FUNCTION OF OBUKHOV-L

!     P. VITERBO         E.C.M.W.F.    08/10/93. (MODIFIED FROM VDFT2M)
!     Modified   P. Viterbo ECMWF 12/05/2005 Externalize SURF
!                                   (based on vdfppcfl)
!                H. Hersbach ECMWF 04/12/2009 Add 10m neutral wind and friction velocity
!     N.Semane+P.Bechtold 04-10-2012 Add PRPLRG factor for small planet
!     F. Vana  05-Mar-2015  Support for single precision
!     A. van Niekerk 05/2023 Remove limiter for 10m wind

!     PURPOSE
!     -------

!     COMPUTE WIND COMPONENTS, TEMPERATURE AND DEWPOINT TEMPERATURE
!     AT SCREEN LEVEL HEIGHT

!     INTERFACE
!     ---------

!     *SPPCFL* IS CALLED BY *SURFPP*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PUMLEV*       U-COMPONENT WIND AT T-1, lowest atmospheric level
!     *PVMLEV*       V-COMPONENT WIND AT T-1, lowest atmospheric level
!     *PQMLEV*       SPECIFIC HUMIDITY AT T-1, lowest atmospheric level
!     *PAPHMS*       surface pressure
!     *PGEOMLEV*     GEOPOTENTIAL AT T-1, lowest atmospheric level
!     *PCPTS*        DRY STATIC ENRGY AT SURFACE
!     *PCPTGZLEV*    DRY STATIC ENERGY    AT T-1, lowest atmospheric level
!     *PZ0MM*        AERODYNAMIC ROUGHNESS LENGTH
!     *PZ0MSM*       AERODYNAMIC ROUGHNESS LENGTH FOR SCALARS
!     *PZ0HM*        ROUGHNESS LENGTH FOR TEMPERATURE
!     *PZ0QM*        ROUGHNESS LENGTH FOR MOISTURE
!     *PZDL*         ZNLEV+Z0M DEVIDED BY OBUKHOV LENGTH
!     *PZDL_WIND*    ZNLEV+Z0M DEVIDED BY OBUKHOV LENGTH (for wind calculation)
!     *PQSA*         APPARENT SURFACE HUMIDITY
!     *PBLEND*       HEIGHT FROM WHICH WIND SPEED IS INTERPOLATED TO 10 M
!     *PFBLEND*      WIND SPEED AT HEIGHT PBLEND
!     *PUCURR*       U COMPONENT OF OCEAN SURFACE CURRENT
!     *PVCURR*       V COMPONENT OF OCEAN SURFACE CURRENT

!     OUTPUT PARAMETERS (REAL):

!     *PU10*         U-COMPONENT WIND AT 10 M
!     *PV10*         V-COMPONENT WIND AT 10 M
!     *P10NU*        U-COMPONENT NEUTRAL WIND AT 10 M
!     *P10NV*        V-COMPONENT NEUTRAL WIND AT 10 M
!     *PUST*         FRICTION VELOCITY
!     *PT2*          TEMPERATURE AT 2 M
!     *PD2*          DEW POINT TEMPERATURE AT 2 M
!     *PQ2*          SPECIFIC HUMIDITY AT 2 M

!     METHOD
!     ------

!     UNIVERSAL FUNCTIONS ARE USED TO INTERPOLATE BETWEEN SURF. AND NLEV

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILE
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPTGZLEV(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZDL_WIND(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0MSM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0HM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0QM(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZDL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSA(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLEND(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFBLEND(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUCURR(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCURR(:) 
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TEXC)        ,INTENT(IN)    :: YDEXC
LOGICAL           ,INTENT(IN)    :: LWIND
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU10(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV10(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P10NU(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P10NV(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUST(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT2(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD2(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ2(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPLRG

!*    LOCAL STORAGE
!     ----- -------

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRD) :: Z10DNL,Z10DNN, Z10M, Z10MP, Z2M, Z2MP,&
 & ZAPH2M, ZCPT2M, ZCVM3, ZCVM4, ZDL, &
 & ZF1, ZFRAC, ZHTQ, ZL, ZNLEV, ZPRH0, ZPRH1, &
 & ZPRH2, ZPRM0, ZPRM1, ZPRM10, ZPRM2, ZPRQ0, &
 & ZWIND, ZZQM1, ZDUMMY,&
 & ZNLEV_S, ZDL_S, ZL_S
REAL(KIND=JPRB) :: ZEPSILON
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! New t2m pp in stable conditions
REAL(KIND=JPRD)   :: ZFRAC_MOISTGRAD, ZFRAC_TEMPGRAD, ZLIMFUNC
REAL(KIND=JPRB)   :: ZSCALE = 0.3_JPRB, ZLIMVAL = 0.75_JPRB
REAL(KIND=JPRD)   :: ZPI
!
!CGL       ZLIMVAL determines limit behavior for very stable
!          -> stable ==> (t2m-ts)/(tnlev-ts) -> ZLIMFUNC
!CGL       ZSCALE  determines how fast this limit behavior is approached

#include "fcsvdfs.h"

!     ------------------------------------------------------------------

!*       1.   INITIALIZE CONSTANTS
!             ---------- ----------

ZEPSILON=100._JPRB*EPSILON(ZEPSILON)

!     THIS SPECIFIES THE HEIGHT FOR U,V (10 M) AND T,Q (2 M)

IF (LHOOK) CALL DR_HOOK('SPPCFL_MOD:SPPCFL',0,ZHOOK_HANDLE)
ASSOCIATE(RCPD=>YDCST%RCPD, RD=>YDCST%RD, RETV=>YDCST%RETV, RG=>YDCST%RG, &
 & RLSTT=>YDCST%RLSTT, RLVTT=>YDCST%RLVTT, RTT=>YDCST%RTT, &
 & RKAP=>YDEXC%RKAP)

Z10M=10._JPRD/PRPLRG
Z2M=2.0_JPRD/PRPLRG
ZHTQ=Z2M*RG

!     ------------------------------------------------------------------

!        2.   COMPUTE  WIND COMPONENTS AND TEMPERATURE
!             ----------------------------------------

DO JL=KIDIA,KFDIA
! z/L limiters are removed for T2m and Q2m, and substitute with a new
! formula. We keep it for 10-metre wind and ustar.
  ZNLEV=PGEOMLEV(JL)/RG+PZ0MM(JL)
  ZDL=SIGN(MAX(ABS(PZDL_WIND(JL)),ZEPSILON),PZDL_WIND(JL))
  ZL=ZNLEV/ZDL

  IF (ZDL  >  RCHBHDL) THEN
    !*ZNLEV=RCHBHDL*ZL+MAX(PZ0MM(JL),PZ0HM(JL))
    ZNLEV=RCHBHDL*ZL+PZ0MM(JL) !remove max
    ZDL=ZNLEV*ZDL/(PGEOMLEV(JL)/RG+PZ0MM(JL))
  ENDIF

! We use different z0m for scalars and wind,
! so different definition of heights etc.
  ZNLEV_S=PGEOMLEV(JL)/RG+PZ0MSM(JL)
  ZDL_S=SIGN(MAX(ABS(PZDL(JL)),ZEPSILON),PZDL(JL))
  ZL_S=ZNLEV_S/ZDL_S

  Z2MP =Z2M+PZ0MSM(JL)

  Z10MP =Z10M+PZ0MM(JL)

  ZWIND =PBLEND(JL)+PZ0MM(JL)

! We use different z0m and z/L for scalars and wind,
! Scalars
  IF (PZDL(JL)  >  0.0_JPRB) THEN
     ZDUMMY = ZDL_S*Z2MP/ZNLEV_S
     ZPRH2=PSIHS(ZDUMMY)
     ZDUMMY = ZDL_S
     ZPRH1=PSIHS(ZDUMMY)
     ZDUMMY = ZDL_S*PZ0HM(JL)/ZNLEV_S
     ZPRH0=PSIHS(ZDUMMY)
     
     ZDUMMY = ZDL_S*PZ0QM(JL)/ZNLEV_S
     ZPRQ0=PSIHS(ZDUMMY)
  ELSE
     ZDUMMY = ZDL_S*Z2MP/ZNLEV_S
     ZPRH2=PSIHU(ZDUMMY)
     ZDUMMY = ZDL_S
     ZPRH1=PSIHU(ZDUMMY)
     ZDUMMY = ZDL_S*PZ0HM(JL)/ZNLEV_S
     ZPRH0=PSIHU(ZDUMMY)
     
     ZDUMMY = ZDL_S*PZ0QM(JL)/ZNLEV_S
     ZPRQ0=PSIHU(ZDUMMY)
  ENDIF
! Wind 
  IF (LWIND) THEN    
    IF (PZDL_WIND(JL)  >  0.0_JPRB) THEN
       ZDUMMY = ZDL*ZWIND/ZNLEV
       ZPRM2 =PSIMS(ZDUMMY)
       ZDUMMY = ZDL*Z10MP/ZNLEV
       ZPRM10=PSIMS(ZDUMMY)
       ZDUMMY = ZDL
       ZPRM1 =PSIMS(ZDUMMY)
       ZDUMMY = ZDL*PZ0MM(JL)/ZNLEV
       ZPRM0 =PSIMS(ZDUMMY)
     ELSE
       ZDUMMY = ZDL*ZWIND/ZNLEV
       ZPRM2 =PSIMU(ZDUMMY)
       ZDUMMY = ZDL*Z10MP/ZNLEV
       ZPRM10=PSIMU(ZDUMMY)
       ZDUMMY = ZDL
       ZPRM1 =PSIMU(ZDUMMY)
       ZDUMMY = ZDL*PZ0MM(JL)/ZNLEV
       ZPRM0 =PSIMU(ZDUMMY)
     ENDIF
  ENDIF

  IF (LWIND) THEN    
    !     10 M WIND COMPONENT
    !      REAL
    Z10DNL= (LOG(Z10MP/PZ0MM(JL))-ZPRM10+ZPRM0)&
     & /(LOG(ZWIND/PZ0MM(JL))-ZPRM2+ZPRM0)  
    ZF1=MAX(0.1_JPRB,SQRT((PUMLEV(JL)-PUCURR(JL))**2&
     & +(PVMLEV(JL)-PVCURR(JL))**2))
    PU10(JL)=(PUMLEV(JL)-PUCURR(JL))*Z10DNL*PFBLEND(JL)/ZF1+PUCURR(JL)
    PV10(JL)=(PVMLEV(JL)-PVCURR(JL))*Z10DNL*PFBLEND(JL)/ZF1+PVCURR(JL)

    !     10M NEUTRAL WIND COMPONENT
    !     is related to stress. Therefore, no ocean current is to be added
    Z10DNN= (LOG(Z10MP/PZ0MM(JL)))&
     & /(LOG(ZWIND/PZ0MM(JL))-ZPRM2+ZPRM0)  
    P10NU(JL)=(PUMLEV(JL)-PUCURR(JL))*Z10DNN*PFBLEND(JL)/ZF1
    P10NV(JL)=(PVMLEV(JL)-PVCURR(JL))*Z10DNN*PFBLEND(JL)/ZF1

    !     FRICTION VELOCITY
      PUST(JL)= PFBLEND(JL)*RKAP&
       & /(LOG(ZWIND/PZ0MM(JL))-ZPRM2+ZPRM0)  
  ELSE
     ! If skipping wind calculation, set it to 0.001
     PU10(JL) = 0.001_JPRB
     PV10(JL) = 0.001_JPRB
     P10NU(JL)= 0.001_JPRB
     P10NV(JL)= 0.001_JPRB
     PUST(JL) = 0.001_JPRB
  ENDIF

!     2M TEMPERATURE AND 2M SPECIFIC HUMIDITY
! New formula based on knmi, Meijgaard et al. Technical report TR - 302
  ZZQM1=MAX(ZEPSILON,PQMLEV(JL))
  ZPI=4.0_JPRD*ATAN(1.0_JPRD)
  ZFRAC_TEMPGRAD = (LOG(Z2MP /PZ0HM(JL))-ZPRH2+ZPRH0)&
                 &/(LOG(ZNLEV_S/PZ0HM(JL))-ZPRH1+ZPRH0)
  ZLIMFUNC = 2.0_JPRD*ATAN(ZSCALE*MAX(PZDL(JL),0.0_JPRB))/ZPI
  ZFRAC_TEMPGRAD = MIN(1._JPRD,ZFRAC_TEMPGRAD*(1.0_JPRD - ZLIMFUNC )  +  ZLIMVAL * ZLIMFUNC)
  ZCPT2M=PCPTS(JL)+ (PCPTGZLEV(JL)-PCPTS(JL))*ZFRAC_TEMPGRAD
  
  ZFRAC_MOISTGRAD = (LOG(Z2MP /PZ0QM(JL))-ZPRH2+ZPRQ0)&
                  &/(LOG(ZNLEV_S/PZ0QM(JL))-ZPRH1+ZPRQ0)
  ZFRAC_MOISTGRAD = MIN(1._JPRD, ZFRAC_MOISTGRAD*(1.0_JPRD - ZLIMFUNC)  + ZLIMVAL * ZLIMFUNC)
  PQ2(JL)=PQSA(JL) + (ZZQM1-PQSA(JL))*ZFRAC_MOISTGRAD

! Standard formula pre-49r1:
!*  ZZQM1=MAX(ZEPSILON,PQMLEV(JL))
!*  ZCPT2M=PCPTS(JL)+ (PCPTGZLEV(JL)-PCPTS(JL))&
!*   & *(LOG(Z2MP /PZ0HM(JL))-ZPRH2+ZPRH0)&
!*   & /(LOG(ZNLEV_S/PZ0HM(JL))-ZPRH1+ZPRH0)  
!*  PQ2(JL)=PQSA(JL) + (ZZQM1-PQSA(JL))&
!*   & *(LOG(Z2MP /PZ0QM(JL))-ZPRH2+ZPRQ0)&
!*   & /(LOG(ZNLEV_S/PZ0QM(JL))-ZPRH1+ZPRQ0)  

!     APPROXIMATE MOISTURE CORRECTION IN CP WITH QNLEV

  PT2(JL)=(ZCPT2M-RG*Z2M)/( RCPD*(1.0_JPRB+RVTMP2*ZZQM1) )
ENDDO

!     ------------------------------------------------------------------

!        3.   COMPUTE  DEW POINT TEMPERATURE
!             ------------------------------

DO JL=KIDIA,KFDIA
  ZZQM1=MAX(ZEPSILON,PQMLEV(JL))
  ZAPH2M=PAPHMS(JL)*(1.0_JPRB-ZHTQ/(RD*PT2(JL)*(1.0_JPRB+RETV*ZZQM1)))
! Note the use ofs saturation with respect to WATER only, in conformity to
!   WMO observation reporting standards
  ZCVM3=R3LES
  ZCVM4=R4LES
  ZFRAC=LOG(ZAPH2M*PQ2(JL)/(R2ES*(1.0_JPRB+RETV*PQ2(JL))))/ZCVM3
  PD2(JL)=(RTT-ZFRAC*ZCVM4)/(1.0_JPRB-ZFRAC)

!     LIMIT DEW POINT TEMPERATURE < TEMPERATURE

  PD2(JL)=MIN(PT2(JL),PD2(JL))
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPPCFL_MOD:SPPCFL',1,ZHOOK_HANDLE)
END SUBROUTINE SPPCFL
END MODULE SPPCFL_MOD

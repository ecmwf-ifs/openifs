! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUPDRA &
 & ( YDECUMF2,YDEPHLI,YDECLDP,&
 & KIDIA,    KFDIA,    KLON,     KLEV,    KINDEX,&
 & K_JKSTART,  LD_LLGO_ON, LDRAIN1D, LDMIXS,&
 & PQENH,    PGEO,     PGEOH,    PAPH,    PQSEN,&
 & PTU,      PQU,      PLU,&
 & PSENH,    PWU2H,    PSUH,&
 & PBUOH,    PLGLAC,   PQPRCV,   PCAPE,&
 & KLAB,     LDCUM,    LDSC,     KCBOT,&
 & KBOTSC,   KCTOP,    KSTUP  )  

!          THIS ROUTINE CALCULATES CHARACTERISTICS OF 
!          THE CUMULUS UPDRAFTS, CLOUD BASE HEIGHT AND 
!          CLOUD TOP HEIGHT

!          P. LOPEZ     ECMWF   (01/2002)
!          inspired from A. Pier Siebesma (KNMI)
!          and C. Jakob (ECMWF) (01/2001) 

!          PURPOSE.
!          --------
!          TO PRODUCE UPDRAFT CALCULATION, CLOUD BASE AND CLOUD
!          TOP VALUES FOR SIMPLIFIED CONVECTIVE PARAMETRIZATION

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUBASEN2*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=0 FOR STABLE LAYERS
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CLOUD LEVELS LEVEL

!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP.
!          ENTRAINING PLUME WITH ENTRAINMENT PROPORTIONAL TO 1/Z
!          FOR SHALLOW CONVECTION (STARTING FROM LOWEST MODEL LEVEL)
!          AND FIXED OTHERWISE.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KINDEX*       TOP INDEX FOR OUTER VERTICAL LOOP
!    *K_JKSTART*    STARTING LEVEL OF CURRENT ASCENDING PLUME

!    INPUT PARAMETERS (LOGICAL):

!    *LD_LLGO_ON*    WHETHER TO CONTINUE ASCENT COMPUTATION (FOR EACH GRID POINT).
!    *LDRAIN1D*      ROUTINE IS CALLED FROM 1D-VAR.
!    *LDMIXS*        WEAK (FALSE) OR STRONG (TRUE) CLOUD MIXING FOR SURFACE PARCEL ONLY

!    INPUT PARAMETERS (REAL):

!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PGEO*         GEOPOTENTIAL ON FULL LEVELS                   M2/S2 
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PQSEN*        PROVISIONAL ENVIRONMENT SATU. HUMIDITY (T+1)  KG/KG

!    UPDATED PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PBUOH*        BUOYANCY ON HALF LEVELS                       M/S2

!    UPDATED PARAMETERS (INTEGER):

!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST

!    OUTPUT PARAMETERS (INTEGER):

!    *KSTUP*       LEVEL AT WHICH THE UPDRAFT STARTS 
!    *KCBOT*       CLOUD BASE LEVEL     
!    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL 
!                  WITH A NON-ZERO CLOUD UPDRAFT.
!    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS

!    OUTPUT PARAMETERS (REAL):

!    *PLGLAC*      FROZEN CLOUD WATER CONTENT                   KG/KG 
!    *PQPRCV*      CONVECTIVE PRECIPITATION CONTENT             KG/KG 
!    *PCAPE*       CONVECTIVE AVAILABLE POTENTIAL ENERGY        J/KG 

!    EXTERNALS
!    ---------
!    *CUADJTQ*  FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
!    *CUADJTQS* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT

!    MODIFICATIONS
!    -------------
!      P.Lopez    01-Oct-2005  New version of linearized convection
!      P.Lopez    11-Jan-2007  Added option LDRAIN (1D-Var rain)
!      P.Lopez    04-Sep-2007  Revised version to improve match to Tiedtke scheme
!      A.Geer     01-Oct-2008  LDRAIN1D name change to reflect usage
!      P.Lopez    09-Nov-2012  Revision to improve match to non-linear version
!      P.Lopez    15-Oct-2015  Added CAPE computation (excl. liquid water loading)
!      P.Lopez    30-Jan-2018  Changed in buoyancy computation to reflect full NL change
!      P.Lopez    28-Feb-2019  Changed entrainment and buoyancy computation to match ref NL
!      P.Lopez    01-Dec-2020  Added entrainment scaling (ENTSTPC32) for PBL height computation

!----------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM,  JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RETV, RD, RG, RCPD, RLVTT, RLSTT, RTT  
USE YOECUMF2  , ONLY : TECUMF2
USE YOECLDP   , ONLY : TECLDP
USE YOETHF    , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                     R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 &                     RTWAT_RTICE_R, RTWAT_RTICECU_R
USE YOEPHLI   , ONLY : TEPHLI
USE YOMDYNCORE, ONLY : RPLRG

IMPLICIT NONE

TYPE(TECUMF2)     ,INTENT(INOUT) :: YDECUMF2
TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
TYPE(TECLDP)      ,INTENT(INOUT) :: YDECLDP
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDEX
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JKSTART 
LOGICAL           ,INTENT(INOUT) :: LD_LLGO_ON(KLON) 
LOGICAL           ,INTENT(IN)    :: LDRAIN1D 
LOGICAL           ,INTENT(IN)    :: LDMIXS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQENH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSENH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWU2H(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSUH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBUOH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLGLAC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQPRCV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCAPE(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLAB(KLON,KLEV) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON) 
LOGICAL           ,INTENT(INOUT) :: LDSC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KBOTSC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KSTUP(KLON) 

!             LOCAL STORAGE
!             ----- -------

INTEGER(KIND=JPIM) :: ICALL, IK, JOK, JK, JL, JKT2

REAL(KIND=JPRB) :: ZQOLD(KLON),ZPH(KLON)
REAL(KIND=JPRB) :: ZMIX(KLON),ZDMIX(KLON),ZDMIX2
REAL(KIND=JPRB) :: ZDZ(KLON),ZCBASE(KLON)

REAL(KIND=JPRB) :: ZBUOF, ZLUOLD, ZCOND, ZLNEW, ZLMAX
REAL(KIND=JPRB) :: ZEPS      ! FRACTIONAL ENTRAINMENT RATE   [M^-1]
REAL(KIND=JPRB) :: ZTVENH    ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)  
REAL(KIND=JPRB) :: ZTVUH     ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
REAL(KIND=JPRB) :: ZQSU, ZCOR, ZDQ, ZALFAW, ZFACW, ZFACI, ZFAC,&
                 & ZESDP,ZDQSDT,ZDTDP,ZDP,ZPDIFFTOP,ZPDIFFBOT,&
                 & ZSF,ZQF,ZAW,ZBW,ZFACT3 
REAL(KIND=JPRB) :: Z1DPH
REAL(KIND=JPRB) :: ZRG, ZRCPD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cuadjtq.intfb.h"
#include "cuadjtqs.intfb.h"

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"

!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUPDRA',0,ZHOOK_HANDLE)
ASSOCIATE(RLMIN=>YDECLDP%RLMIN, &
 & ENTRORG2=>YDECUMF2%ENTRORG2, ENTSTPC12=>YDECUMF2%ENTSTPC12, &
 & ENTSTPC22=>YDECUMF2%ENTSTPC22, NJKT22=>YDECUMF2%NJKT22, &
 & ENTSTPC32=>YDECUMF2%ENTSTPC32,&
 & LPHYLIN=>YDEPHLI%LPHYLIN)
ZAW=1.0_JPRB
ZBW=1.0_JPRB
ZRG=1.0_JPRB/RG
ZRCPD=1.0_JPRB/RCPD

JKT2=NJKT22

DO JL=KIDIA,KFDIA
  PCAPE(JL) = 0.0_JPRB
ENDDO

!----------------------------------------------------------------------

!     1.0          DO ASCENT IN SUBCLOUD AND LAYER,
!                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  CHECK FOR BUOYANCY AND SET FLAGS
!                  -------------------------------------
!       ------------------------------------------------------------
!        1.1  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
!       ------------------------------------------------------------

DO JK=K_JKSTART,JKT2,-1
  JOK=0
  DO JL=KIDIA,KFDIA
    IF (LD_LLGO_ON(JL)) THEN
      JOK=JOK+1
      ZDZ(JL) = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
      IF (K_JKSTART == KLEV-1) THEN
        ZEPS = ENTSTPC12/((PGEO(JL,JK)-PGEOH(JL,KLEV+1))*ZRG*RPLRG) + ENTSTPC22
        ! The following line is only active in the computation of PBL height (no TL/AD).
        IF (LDMIXS .AND. KINDEX == KLEV-1 .AND. PLU(JL,JK+1) > 0.0_JPRB) ZEPS = ZEPS * ENTSTPC32
        ZMIX(JL)= 0.5_JPRB*RPLRG*ZDZ(JL)*ZEPS
        ZMIX(JL) = MIN(ZMIX(JL),1.0_JPRB)
        ZDMIX(JL)= 1.0_JPRB/(1.0_JPRB+ZMIX(JL))
        ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5_JPRB
        ZSF = (PSENH(JL,JK+1) + PSENH(JL,JK))*0.5_JPRB
        PQU(JL,JK)= (PQU(JL,JK+1)*(1.0_JPRB-ZMIX(JL))&
         & + 2.0_JPRB*ZQF*ZMIX(JL))*ZDMIX(JL)  
        PSUH (JL,JK)= (PSUH(JL,JK+1)*(1.0_JPRB-ZMIX(JL))&
         & + 2.0_JPRB*ZSF*ZMIX(JL))*ZDMIX(JL)  
      ELSE
        ZEPS = 0.4_JPRB*ENTRORG2*MIN(1.0_JPRB,(PQSEN(JL,JK)/PQSEN(JL,KLEV))**3)
        ZMIX(JL)= ZDZ(JL)*ZEPS
        ZMIX(JL) = MIN(ZMIX(JL),1.0_JPRB)
        ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5_JPRB
        ZSF = (PSENH(JL,JK+1) + PSENH(JL,JK))*0.5_JPRB
        PQU(JL,JK) = PQU(JL,JK+1) *(1.0_JPRB-ZMIX(JL)) + ZQF*ZMIX(JL)
        PSUH(JL,JK)= PSUH(JL,JK+1)*(1.0_JPRB-ZMIX(JL)) + ZSF*ZMIX(JL) 
      ENDIF
      ZQOLD(JL)   = PQU(JL,JK)
      PTU (JL,JK) = (PSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
      ZPH  (JL)   = PAPH(JL,JK)
    ENDIF
  ENDDO
  
  IK=JK
  IF (LPHYLIN .OR. LDRAIN1D) THEN
    ICALL=1
    CALL CUADJTQS &
     & ( KIDIA,    KFDIA,    KLON,    KLEV,        IK,&
     &   ZPH,      PTU,      PQU,     LD_LLGO_ON,  ICALL)  
  ELSE
    ICALL=1
    CALL CUADJTQ &
     & ( YDEPHLI,  KIDIA,    KFDIA,    KLON,    KLEV,        IK,&
     &   ZPH,      PTU,      PQU,     LD_LLGO_ON,  ICALL)  
  ENDIF

!DIR$ IVDEP
!OCL NOVREC

  DO JL=KIDIA,KFDIA
    IF(LD_LLGO_ON(JL)) THEN

! condensation 

      ZCOND=MAX(ZQOLD(JL)-PQU(JL,JK),0.0_JPRB)

! microphysics

      ZLUOLD=PLU(JL,JK+1)
      ZLMAX=MAX(ZLUOLD+ZCOND,0.0_JPRB)

! Simple microphysics at this stage (to mimic Tiedtke scheme)

      IF (K_JKSTART == KLEV-1) THEN  ! No precip for shallow convection
        ZLNEW=MIN(ZLMAX,5.E-3_JPRB)
      ELSE 
        ZLNEW=0.5_JPRB*ZLMAX
      ENDIF
      PQPRCV(JL,JK)=MAX(0.0_JPRB,ZLMAX-ZLNEW)
      PLU(JL,JK)=ZLNEW

! frozen water

      PLGLAC(JL,JK)=ZCOND*((1.0_JPRB-FOEALFCU(PTU(JL,JK)))-&
       & (1.0_JPRB-FOEALFCU(PTU(JL,JK+1))))   
                
! update dry static energy after condensation + freezing

      PSUH(JL,JK)    = RCPD*(PTU(JL,JK)+RALFDCP*PLGLAC(JL,JK))+PGEOH(JL,JK)
      
! Buoyancy on half and full levels
         
      ZTVUH           = (1.0_JPRB+RETV*PQU(JL,JK)-PLU(JL,JK))*PTU(JL,JK)&
       & +RALFDCP*PLGLAC(JL,JK)  
      ZTVENH          = (1.0_JPRB+RETV*PQENH(JL,JK)) &
       & *(PSENH(JL,JK)-PGEOH(JL,JK))*ZRCPD  
      PBUOH(JL,JK)   = (ZTVUH-ZTVENH)*RG/ZTVENH
      ZBUOF = (PBUOH(JL,JK) + PBUOH(JL,JK+1))*0.5_JPRB

      IF (ZBUOF > 0.0_JPRB) THEN
        PCAPE(JL) = PCAPE(JL) + ZBUOF*ZDZ(JL)
      ENDIF

! solve kinetic energy equation

      ZDMIX2=1.0_JPRB/(1.0_JPRB+2.0_JPRB*ZBW*ZMIX(JL))
      PWU2H(JL,JK) = (PWU2H(JL,JK+1)*(1.0_JPRB-2.0_JPRB*ZBW*ZMIX(JL))&
       & + 2.0_JPRB*ZAW*ZBUOF*ZDZ(JL))*ZDMIX2  

! first layer with liquid water - find exact cloud base

      IF(PLU(JL,JK) >0.0_JPRB.AND.KLAB(JL,JK+1)==1) THEN
        
        IK=JK+1

        ZALFAW=FOEALFA(PTU(JL,IK))
        Z1DPH=1.0_JPRB/PAPH(JL,IK)
        ZQSU=FOEEWM(PTU(JL,IK))*Z1DPH
        ZQSU=MIN(0.5_JPRB,ZQSU)
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZQSU)
        ZQSU=ZQSU*ZCOR
        ZDQ=PQU(JL,IK)-ZQSU
        ZFACW=R5LES/((PTU(JL,IK)-R4LES)**2)
        ZFACI=R5IES/((PTU(JL,IK)-R4IES)**2)
        ZFAC=ZALFAW*ZFACW+(1.-ZALFAW)*ZFACI
        ZESDP=FOEEWM(PTU(JL,IK))*Z1DPH
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZESDP)
        ZDQSDT=ZFAC*ZCOR*ZQSU
        ZDTDP=RD*PTU(JL,IK)*ZRCPD*Z1DPH
        ZFACT3=1.0_JPRB/(ZDQSDT*ZDTDP)
        ZDP=ZDQ*ZFACT3
        ZCBASE(JL)=PAPH(JL,IK)+ZDP
        
! chose nearest half level as cloud base

        ZPDIFFTOP=ZCBASE(JL)-PAPH(JL,JK)
        ZPDIFFBOT=PAPH(JL,JK+1)-ZCBASE(JL)
        
        IF(ZPDIFFTOP > ZPDIFFBOT.AND.PWU2H(JL,JK+1)>0.0_JPRB.AND.&
           & JK < KLEV-1) THEN  
          KLAB(JL,JK)=2
          KLAB(JL,JK+1)=2
          LDSC(JL)   =.TRUE.
          KBOTSC(JL) =JK+1
          KCBOT(JL)  =JK+1
          PLU(JL,JK+1) = RLMIN
        ELSEIF((ZPDIFFTOP <= ZPDIFFBOT .OR. &
           & (ZPDIFFTOP > ZPDIFFBOT .AND. JK == KLEV-1)) .AND. &
           & PWU2H(JL,JK)>0.0_JPRB) THEN  
          KLAB(JL,JK)=2
          LDSC(JL)   =.TRUE.
          KBOTSC(JL) =JK
          KCBOT(JL)  =JK
        ENDIF
      ENDIF

! decide on presence of convection, cloud base and cloud top based on
! kinetic energy

      IF (PWU2H(JL,JK) < 0.0_JPRB) THEN
        LD_LLGO_ON(JL) = .FALSE.
        IF (PLU(JL,JK+1) > 0.0_JPRB .AND. JK < K_JKSTART) THEN
          KCTOP(JL)   = JK
          KSTUP(JL)   = K_JKSTART+1
          LDCUM(JL)   = .TRUE.
        ELSE
          LDCUM(JL)   = .FALSE.
          KCBOT(JL)   = -1
          KCTOP(JL)   = -1
        ENDIF
      ELSE
        IF (PLU(JL,JK)>0.0_JPRB) THEN
          KLAB(JL,JK) = 2
        ELSE
          KLAB(JL,JK) = 1
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  IF (JOK == 0) EXIT

ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUPDRA',1,ZHOOK_HANDLE)
END SUBROUTINE CUPDRA

! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUASCN2 &
 & (YDECUMF2, YDEPHLI, KIDIA,   KFDIA,    KLON,    KLEV,&
 & LDRAIN1D, LDLAND,   PTEN,    PTENH,    PQENH,&
 & PGEOH,    PQEN,     PQSEN,&
 & PAPH,&
 & PLGLAC,   PQPRCV,&
 & PBUOH,    PWUBASE,  LDCUM,    KTYPE,&
 & PTU,      PQU,      PLU,&
 & PMFU,     PMFUB,&
 & PMFUS,    PMFUQ,    PMFUL,    PLUDE,    PDMFUP,&
 & KCBOT,    KCTOP,    KSTUP, &
 & PMFUDE_RATE, PWMEAN, PERATEU )  

!          THIS ROUTINE DOES THE CALCULATIONS OF THE MASS-FLUX PROFILE
!          FOR THE SIMPLIFIED CUMULUS PARAMETERIZATION

!          P.LOPEZ           E.C.M.W.F.     Jan 2001

!          PURPOSE.
!          --------
!          TO PRODUCE MASS-FLUX PROFILE FOR SIMPLIFIED 
!          CONVECTIVE PARAMETRIZATION

!          INTERFACE
!          ---------

!          THIS ROUTINE IS CALLED FROM *CUMASTRN2*.

!          METHOD.
!          --------

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (LOGICAL)

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)

!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KSTUP*        LEVEL AT WHICH THE UPDRAFT STARTS 

!    INPUT PARAMETERS (REAL):

!    *PTEN*         TEMPERATURE ON FULL LEVELS                   K 
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
!    *PQEN*         SPECIFIC HUMIDITY ON FULL LEVELS             KG/KG 
!    *PQSEN*        SATURATION SPECIFIC HUMIDITY ON FULL LEVELS  KG/KG 
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS          PA

!    UPDATED PARAMETERS (REAL):

!    *PTENH*        TEMPERATURE ON HALF LEVELS                   K 
!    *PQENH*        SPECIFIC HUMIDITY ON HALF LEVELS             KG/KG 
!    *PQPRCV*       CONVECTIVE PRECIPITATION CONTENT             KG/KG
!    *PBUOH*        BUOYANCY ON HALF-LEVELS                      M/S2
!    *PTU*          TEMPERATURE IN UPDRAFTS                      K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                   KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG

!    UPDATED PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 

!    UPDATED PARAMETERS (INTEGER):

!    *KCTOP*        CLOUD TOP LEVEL

!    OUTPUT PARAMETERS (REAL):

!    *PMFU*         MASSFLUX IN UPDRAFTS                         KG/(M2*S)
!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE           KG/(M2*S)
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS         J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS           KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS             KG/(M2*S)
!    *PLUDE*        DETRAINED LIQUID WATER                       KG/(M2*S)
!    *PLGLAC*       FROZEN CLOUD WATER/RAIN CONTENT              KG/KG
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS       KG/(M2*S)
!    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                     KG/(M2*S)
!    *PWMEAN*       MEAN UPDRAFT VELOCITY                         M/S
!    *PERATEU*      UPDRAFT ENTRAINMENT RATE                      

!          EXTERNALS
!          ---------

!          REFERENCE
!          ---------

!          MODIFICATIONS
!          -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.Lopez       01-Oct-2005 Changed entrainment/detrainment 
!        P.Lopez       19-Oct-2007 Revised version to improve match to Tiedtke scheme
!        A.Geer        01-Oct-2008 LDRAIN1D name change to reflect usage
!        P.Lopez       10-May-2016 Revised precipitation glaciation

!----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RETV,   RG,     RCPD,&
 &                    RLVTT,  RLSTT, RTT  
USE YOECUMF2 , ONLY : TECUMF2
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                    R4IES    ,R5LES    ,R5IES   ,R5ALVCP  ,R5ALSCP  ,&
 &                    RALVDCP  ,RALSDCP  ,RALFDCP ,RTWAT    ,RTICE    ,RTICECU  ,&
 &                    RTWAT_RTICE_R      ,RTWAT_RTICECU_R    ,RTBER    ,RTBERCU  
USE YOEPHLI  , ONLY : TEPHLI

IMPLICIT NONE

TYPE(TECUMF2)     ,INTENT(INOUT) :: YDECUMF2
TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
LOGICAL           ,INTENT(IN)    :: LDRAIN1D 
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQENH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQEN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLGLAC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQPRCV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBUOH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWUBASE(KLON) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KTYPE(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFUB(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFUS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFUQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFUL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLUDE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDMFUP(KLON,KLEV) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KSTUP(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFUDE_RATE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWMEAN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PERATEU(KLON,KLEV) 

!     LOCAL SCALAR ARRAYS
REAL(KIND=JPRB) ::  ZWU2H(KLON,KLEV), ZSUH(KLON,KLEV) 
REAL(KIND=JPRB) ::  ZSENH(KLON,KLEV)
REAL(KIND=JPRB) ::  ZERATE(KLON,KLEV), ZDERATE(KLON,KLEV)
REAL(KIND=JPRB) ::  ZPMID(KLON)
REAL(KIND=JPRB) ::  ZDPMEAN(KLON)
REAL(KIND=JPRB) ::  ZQOLD(KLON), ZPH(KLON)
REAL(KIND=JPRB) ::  ZEPSORG(KLON)
REAL(KIND=JPRB) ::  ZENT(KLON), ZDET(KLON), ZDZ(KLON)

LOGICAL :: LLTEST(KLON)

INTEGER(KIND=JPIM) :: JK, JL, IKB, JKT2, IK, ICALL, IKS

REAL(KIND=JPRB) :: ZRH(KLON), ZFACTQS(KLON)
REAL(KIND=JPRB) :: ZRATQS, ZDNOPRC
REAL(KIND=JPRB) :: ZBUOF, ZAW, ZBW
REAL(KIND=JPRB) :: ZRG, ZCONS1
REAL(KIND=JPRB) :: ZRCPD, ZLUOLD
REAL(KIND=JPRB) :: ZPRCDGW, Z_CLCRIT, Z_CPRC2, ZPRCON, ZDT, ZZCO, ZC
REAL(KIND=JPRB) :: ZLCRIT, ZDFI, ZDT1, ZLCRITI, ZCBFI, ZEXP, ZSQRT
REAL(KIND=JPRB) :: ZFACT2, ZFACT22, ZFACT3, ZFACT4, ZFACT7
REAL(KIND=JPRB) :: ZD, ZCOND, ZINT, ZLNEW, ZCBF, ZWU, ZWU0, Z1ZD
REAL(KIND=JPRB) :: ZQCMAX, ZLMAX, ZTVENH, ZTVUH, ZFAC, ZRHFACT
REAL(KIND=JPRB) :: ZOCUDET, ZRATW, ZSQRATW, ZDETNEW, ZGLAC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cuadjtq.intfb.h"
#include "cuadjtqs.intfb.h"

#include "fcttre.func.h"
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUASCN2',0,ZHOOK_HANDLE)
ASSOCIATE(DETRPEN2=>YDECUMF2%DETRPEN2, ENTRORG2=>YDECUMF2%ENTRORG2, &
 & ENTSHALP2=>YDECUMF2%ENTSHALP2, NJKT22=>YDECUMF2%NJKT22, &
 & RPRCON2=>YDECUMF2%RPRCON2, LMFGLAC2=>YDECUMF2%LMFGLAC2, &
 & LPHYLIN=>YDEPHLI%LPHYLIN)
!*    1.           SPECIFY PARAMETERS
!                  ------------------

ZAW=0.5_JPRB/(1.0_JPRB+0.5_JPRB)
ZBW=1.0_JPRB+(3._JPRB/8._JPRB)*0.506_JPRB/0.2_JPRB
ZQCMAX=5.E-03_JPRB
ZPRCDGW=RPRCON2/RG
Z_CLCRIT=0.5E-3_JPRB
Z_CPRC2=0.5_JPRB
ZRCPD=1.0_JPRB/RCPD
ZRG=1.0_JPRB/RG
ZCONS1=3.0_JPRB/1.5E+04_JPRB

JKT2=NJKT22

IF(LMFGLAC2) THEN
  ZGLAC=1.0_JPRB
ELSE
  ZGLAC=0.0_JPRB
ENDIF

!----------------------------------------------------------------------

!     2.           SET DEFAULT VALUES
!                  ------------------

DO JL=KIDIA,KFDIA
  IF(.NOT.LDCUM(JL)) THEN
    KCBOT(JL)=-1
    KCTOP(JL)=-1
    KTYPE(JL)=0
    KSTUP(JL)=-1
  ELSE
    ZPMID(JL)=0.5_JPRB*(PAPH(JL,KCTOP(JL))+PAPH(JL,KCBOT(JL)))
  ENDIF
  ZEPSORG(JL) = 0.0_JPRB
ENDDO

! initialize various quantities

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PMFU(JL,JK)=0.0_JPRB
    PMFUS(JL,JK)=0.0_JPRB
    PMFUQ(JL,JK)=0.0_JPRB
    PMFUL(JL,JK)=0.0_JPRB
    PLUDE(JL,JK)=0.0_JPRB
    PDMFUP(JL,JK)=0.0_JPRB
    ZWU2H(JL,JK)=0.0_JPRB
    PMFUDE_RATE(JL,JK)=0.0_JPRB
    ZSENH(JL,JK)=RCPD*PTENH(JL,JK)+PGEOH(JL,JK)
    ZSUH(JL,JK)=0.0_JPRB
    PERATEU(JL,JK)=0.0_JPRB
    ZERATE(JL,JK)=0.0_JPRB
    ZDERATE(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!----------------------------------------------------------------------

!     3.0          INITIALIZE VALUES AT cloud base LEVEL
!                  -------------------------------------

DO JL=KIDIA,KFDIA
  KCTOP(JL)=KCBOT(JL)
  PWMEAN(JL)=0.0_JPRB
  ZDPMEAN(JL)=0.0_JPRB
  IF (LDCUM(JL)) THEN
    IKB=KCBOT(JL)
    IKS=KSTUP(JL)-1
    ZWU2H(JL,IKB)=PWUBASE(JL)**2
    PMFU(JL,IKB)=PMFUB(JL)
    PMFUS(JL,IKB)=PTU(JL,IKB)
    PMFUQ(JL,IKB)=PQU(JL,IKB)
    PMFUL(JL,IKB)=PLU(JL,IKB)
    PDMFUP(JL,IKB)=PQPRCV(JL,IKB)*PMFU(JL,IKB)
    ZSUH(JL,IKB)=RCPD*PTU(JL,IKB)+PGEOH(JL,IKB)
    DO JK=1,IKB-1
      PLU(JL,JK)=0.0_JPRB
      PLGLAC(JL,JK)=0.0_JPRB 
      PQPRCV(JL,JK)=0.0_JPRB 
      PBUOH(JL,JK)=0.0_JPRB 
    ENDDO 
!LOP Initialise buoyancy to zero at cloud base to mimic Tietdke scheme
    PBUOH(JL,IKB)=0.0_JPRB 
  ENDIF
ENDDO


!----------------------------------------------------------------------

!     4.0          COMPUTE MASS-FLUX PROFILE
!                  -------------------------

DO JK=KLEV-1,JKT2,-1
  DO JL=KIDIA,KFDIA   
    IF (LDCUM(JL) .AND. JK <= KCBOT(JL)) THEN
      IKB=KCBOT(JL)
      ZDZ(JL)=(PGEOH(JL,JK)-PGEOH(JL,JK+1))*ZRG
      ZFACT7=1.0_JPRB/PQSEN(JL,IKB)
      ZRATQS=PQSEN(JL,JK)*ZFACT7
      ZFACTQS(JL)=MIN(1.0_JPRB,ZRATQS)**3
      ZFACT4=1.0_JPRB/PQSEN(JL,JK)
      ZRH(JL)=PQEN(JL,JK)*ZFACT4
    ENDIF
    LLTEST(JL)=(LDCUM(JL) .AND. JK < KCBOT(JL) .AND.&
             & PMFU(JL,JK+1) > 0.0_JPRB .AND. ZWU2H(JL,JK+1) > 0.0_JPRB)
    IF (LLTEST(JL)) THEN  

      ZENT(JL) = ZEPSORG(JL)
      IF (KTYPE(JL) >= 2) THEN
        ZENT(JL) = ENTSHALP2 * ZENT(JL)
        ZDET(JL) = ZENT(JL)
      ELSE
        ZDET(JL) = DETRPEN2 * ZDZ(JL)
      ENDIF 
      ZDET(JL) = ZDET(JL) * (1.6_JPRB - MIN(1.0_JPRB,ZRH(JL)))

      PQU(JL,JK) = PQU(JL,JK+1) *(1.0_JPRB-ZENT(JL)) + PQENH(JL,JK+1)*ZENT(JL)
      ZSUH(JL,JK)= ZSUH(JL,JK+1)*(1.0_JPRB-ZENT(JL)) + ZSENH(JL,JK+1)*ZENT(JL) 
      ZQOLD(JL)  = PQU(JL,JK)
      PTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
      ZPH (JL) = PAPH(JL,JK)

    ENDIF
  ENDDO

  IK=JK
  IF (LPHYLIN .OR. LDRAIN1D) THEN
    ICALL=1
    CALL CUADJTQS &
     & ( KIDIA,    KFDIA,    KLON,    KLEV,     IK,&
     &   ZPH,      PTU,      PQU,     LLTEST,   ICALL)  
  ELSE
    ICALL=1
    CALL CUADJTQ &
     & ( YDEPHLI,  KIDIA,  KFDIA,  KLON,   KLEV,  IK,&
     &   ZPH,      PTU,    PQU,    LLTEST, ICALL)
  ENDIF

  DO JL=KIDIA,KFDIA
    IF (LLTEST(JL)) THEN  

      IF(LDLAND(JL)) THEN
        ZDNOPRC=5.E-4_JPRB
      ELSE
        ZDNOPRC=3.E-4_JPRB
      ENDIF

! condensation 

      ZCOND=MAX(ZQOLD(JL)-PQU(JL,JK),0.0_JPRB)

! microphysics

      ZLUOLD=PLU(JL,JK+1)
      ZLMAX=MAX(ZLUOLD+ZCOND,0.0_JPRB)

      IF(ZLMAX > ZDNOPRC) THEN
        ZWU0=SQRT(ZWU2H(JL,JK+1))
! Note: 0.447 corresponds to threshold of 0.1 for kinetic energy (m2/s2)
        ZWU=MIN(15._JPRB,MAX(0.447_JPRB,ZWU0))
        ZPRCON=ZPRCDGW/(0.75_JPRB*ZWU)
        ZDT1=RTBER-PTU(JL,JK)
        ZDT=MIN(RTBERCU-RTICECU,MAX(ZDT1,0.0_JPRB))
        ZSQRT=SQRT(ZDT)
        ZCBF=1.0_JPRB+Z_CPRC2*ZSQRT
        ZZCO=ZPRCON*ZCBF
        ZCBFI=1.0_JPRB/ZCBF
        ZLCRIT=ZDNOPRC*ZCBFI
        ZDFI=PGEOH(JL,JK)-PGEOH(JL,JK+1)
        ZC=ZLMAX-ZLUOLD
        ZLCRITI=1.0_JPRB/ZLCRIT
        ZEXP=EXP(-(ZLMAX*ZLCRITI)**2)
        ZD=ZZCO*ZDFI*(1.0_JPRB-ZEXP)
        ZINT=EXP(-ZD)
        Z1ZD=1.0_JPRB/ZD
        ZLNEW=ZLUOLD*ZINT+ZC*Z1ZD*(1.0_JPRB-ZINT)
        ZLNEW=MAX(0.0_JPRB,MIN(ZLMAX,ZLNEW))
      ELSE
        ZLNEW=ZLMAX
      ENDIF

      ZLNEW=MIN(ZQCMAX,ZLNEW)
      PQPRCV(JL,JK)=MAX(0.0_JPRB,ZLMAX-ZLNEW)
      PLU(JL,JK)=ZLNEW

! frozen water

      PLGLAC(JL,JK)=ZCOND*((1.0_JPRB-FOEALFCU(PTU(JL,JK)))&
                 & -(1.0_JPRB-FOEALFCU(PTU(JL,JK+1))))   

! add glaciation of rain

      IF (RTT > PTEN(JL,JK)) THEN
        PLGLAC(JL,JK) = PLGLAC(JL,JK) + ZGLAC*FOEALFCU(PTEN(JL,JK))*PQPRCV(JL,JK)
      ENDIF
                
! update dry static energy after condensation + freezing

      ZSUH(JL,JK) = RCPD*(PTU(JL,JK)+RALFDCP*PLGLAC(JL,JK))+PGEOH(JL,JK)
      
! Buoyancy on half and full levels
         
      ZTVUH        = PTU(JL,JK)*(1.0_JPRB+RETV*PQU(JL,JK)-PLU(JL,JK+1) &
                 & - PQPRCV(JL,JK+1)) + RALFDCP*PLGLAC(JL,JK)
      ZTVENH       = PTENH(JL,JK)*(1.0_JPRB+RETV*PQENH(JL,JK)) 
      ZFACT2 = 1.0_JPRB/ZTVENH
      ZFACT22 = ZFACT2
      PBUOH(JL,JK) = RG*(ZTVUH-ZTVENH)*ZFACT2

      IF (PBUOH(JL,JK) < 0.0_JPRB .AND. PLU(JL,JK) > 0.0_JPRB) THEN
        PTENH(JL,JK) = 0.5_JPRB*(PTEN(JL,JK)+PTEN(JL,JK-1))
        PQENH(JL,JK) = 0.5_JPRB*(PQEN(JL,JK)+PQEN(JL,JK-1))
        ZSENH(JL,JK) = RCPD*PTENH(JL,JK)+PGEOH(JL,JK)
        ZTVENH = PTENH(JL,JK)*(1.0_JPRB+RETV*PQENH(JL,JK))
        ZFACT22 = 1.0_JPRB/ZTVENH
        PBUOH(JL,JK) = RG*(ZTVUH-ZTVENH)*ZFACT22
      ENDIF 

      ZBUOF = (PBUOH(JL,JK) + PBUOH(JL,JK+1))*0.5_JPRB

! Solve kinetic energy equation

      ZWU2H(JL,JK) = (ZWU2H(JL,JK+1)*(1.0_JPRB-ZBW*ZENT(JL))&
         & + 2.0_JPRB*ZAW*ZBUOF*ZDZ(JL))/(1.0_JPRB+ZBW*ZENT(JL))  

! Accumulate updraught kinetic energy in the vertical

      PWMEAN(JL)  = PWMEAN(JL)  + ZWU2H(JL,JK)*(PAPH(JL,JK+1) - PAPH(JL,JK))
      ZDPMEAN(JL) = ZDPMEAN(JL) + PAPH(JL,JK+1) - PAPH(JL,JK)

! Entrainment and detrainment according to variation of kinetic energy
! on the vertical.

      IF (ZWU2H(JL,JK) >= ZWU2H(JL,JK+1)) THEN
        PMFU(JL,JK) = PMFU(JL,JK+1)*(1.0_JPRB+ZENT(JL)-ZDET(JL))
        PLUDE(JL,JK) = ZDET(JL)*PLU(JL,JK+1)
        PMFUDE_RATE(JL,JK) = ZDET(JL)*PMFU(JL,JK+1)
      ELSE
! Organised detrainment
        IF (PBUOH(JL,JK) < 0.0_JPRB .AND. PLU(JL,JK) > 0.0_JPRB) THEN
          ZFACT3=1.0_JPRB/ZWU2H(JL,JK+1)
          ZRATW=ZWU2H(JL,JK)*ZFACT3
          IF (ZRATW > 0.0_JPRB) THEN
            ZSQRATW=SQRT(ZRATW)
            IF (ZRH(JL) < 1.0_JPRB) THEN
              ZOCUDET=1.6_JPRB-ZRH(JL)
            ELSE
              ZOCUDET=0.6_JPRB
            ENDIF
            ZDETNEW=1.0_JPRB-ZOCUDET*ZSQRATW
          ELSE
            ZDETNEW=1.0_JPRB
          ENDIF
          ZDET(JL)=MAX(ZDET(JL),ZDETNEW)
        ENDIF
        IF (ZWU2H(JL,JK) > 0.0_JPRB) THEN
          PMFU(JL,JK) = PMFU(JL,JK+1)*(1.0_JPRB+ZENT(JL)-ZDET(JL))
          PLUDE(JL,JK) = ZDET(JL)*PLU(JL,JK+1)
          PMFUDE_RATE(JL,JK) = ZDET(JL)*PMFU(JL,JK+1)
        ENDIF
      ENDIF
 
      IF (PMFU(JL,JK) <= 0.0_JPRB .OR. ZWU2H(JL,JK) <= 0.0_JPRB &
        & .OR. ZBUOF == 0.0_JPRB) THEN
        KCTOP(JL)=JK+1
        PMFU(JL,JK)=0.0_JPRB
        PLUDE(JL,JK)=PLU(JL,JK+1)
        ZDET(JL)=1.0_JPRB
        PMFUDE_RATE(JL,JK)=PMFU(JL,JK+1)
      ENDIF

! Store entrainment and detrainment rates for momentum computations
      ZERATE (JL,JK)=ZENT(JL)
      ZDERATE(JL,JK)=ZDET(JL)

! Organised entrainment for next level
      IF (ZBUOF > -0.2_JPRB*RG*ZFACT22 .AND. &
        & PLU(JL,JK+1) > 0.0_JPRB) THEN
        ZRHFACT=0.30_JPRB-(MIN(1.0_JPRB,ZRH(JL))-1.0_JPRB)
        ZEPSORG(JL)=ENTRORG2*ZRHFACT*ZDZ(JL)*ZFACTQS(JL)
        ZEPSORG(JL)=MIN(0.4_JPRB,ZEPSORG(JL))
      ELSE
        ZEPSORG(JL)=0.0_JPRB
      ENDIF

      PMFUS (JL,JK)=PTU(JL,JK)
      PMFUQ (JL,JK)=PQU(JL,JK)
      PMFUL (JL,JK)=PLU(JL,JK)
      PDMFUP(JL,JK)=PQPRCV(JL,JK)*PMFU(JL,JK)
    ENDIF

  ENDDO
ENDDO

! Average updraught velocity
DO JL=KIDIA,KFDIA
  IF(LDCUM(JL)) THEN
    IF (ZDPMEAN(JL) > 1.0_JPRB) THEN
      ZFACT2=1.0_JPRB/ZDPMEAN(JL)
      PWMEAN(JL)=PWMEAN(JL)*ZFACT2
    ENDIF
    PWMEAN(JL)=MAX(1.E-2_JPRB,PWMEAN(JL))
    PWMEAN(JL)=SQRT(PWMEAN(JL))
  ENDIF
ENDDO

! Disregard single layer convective clouds
! to avoid problem in CAPE closure.
DO JL=KIDIA,KFDIA
  IF (LDCUM(JL) .AND. KCTOP(JL) == KCBOT(JL)) THEN
    LDCUM(JL)=.FALSE.
    KCTOP(JL)=-2
    KCBOT(JL)=-2
    KTYPE(JL)=0
    PWMEAN(JL)=0.0_JPRB
    DO JK=1,KLEV
      ZERATE(JL,JK)=0.0_JPRB
      ZDERATE(JL,JK)=0.0_JPRB
    ENDDO
  ENDIF
ENDDO

! Updraught modified entrainment rates for momentum computations
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF (LDCUM(JL)) THEN
      ZFAC=0.0_JPRB
      IF (KTYPE(JL) == 1 .OR. KTYPE(JL) == 3) ZFAC=2.0_JPRB
      IF (KTYPE(JL) == 1 .AND. JK <= KCTOP(JL)+2) ZFAC=3.0_JPRB
      PERATEU(JL,JK)=ZERATE(JL,JK)+ZFAC*ZDERATE(JL,JK)
    ENDIF
  ENDDO
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUASCN2',1,ZHOOK_HANDLE)
END SUBROUTINE CUASCN2

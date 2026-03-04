! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SWR &
 & ( YDEOVLP,YDERAD,YDPHNC,YDECLD,KIDIA , KFDIA , KLON , KLEV  , KNU,&
 & PALBD , PCG   , PCLD , POMEGA, PSEC , PTAU,&
 & PCGAZ , PPIZAZ, PRAY1, PRAY2 , PREFZ, PRJ  , PRK , PRMUE,&
 & PTAUAZ, PTRA1 , PTRA2, PTRCLD, &
! for adjoint computation
 & PRE1  , PRE2  , PTR1 , PTR2  &
 & )  

!**** *SWR* - CONTINUUM SCATTERING COMPUTATIONS

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
!     CONTINUUM SCATTERING

!**   INTERFACE.
!     ----------

!          *SWR* IS CALLED EITHER FROM *SW1S*
!                              OR FROM *SWNI*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
!     OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)

!     EXTERNALS.
!     ----------

!          *SWDE*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!      JEAN-JACQUES MORCRETTE  *ECMWF*
!      ORIGINAL : 89-07-14

!     MODIFICATIONS.
!     --------------
!      03-10-10 Deborah Salmond and Marta Janiskova Optimisation
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Y.Seity 06-09-09 : add modset from O.Thouron (MesoNH) under NOVLP tests
!      M.Janiskova 08-03-27 : modifications for optimized TL/AD code and cleaning
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERAD   , ONLY : TERAD
USE YOECLD   , ONLY : TECLD
USE YOEOVLP  , ONLY : TEOVLP
USE YOPHNC   , ONLY : TPHNC

IMPLICIT NONE

TYPE(TECLD)       ,INTENT(IN)    :: YDECLD
TYPE(TEOVLP)      ,INTENT(IN)    :: YDEOVLP
TYPE(TERAD)       ,INTENT(IN)    :: YDERAD
TYPE(TPHNC)       ,INTENT(IN)    :: YDPHNC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,YDERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCG(KLON,YDERAD%NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMEGA(KLON,YDERAD%NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEC(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU(KLON,YDERAD%NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAY1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAY2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREFZ(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRJ(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRK(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRMUE(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRA1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRA2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRCLD(KLON) 
! for adjoint computation from swde
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE2(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTR1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTR2(KLON,KLEV)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZC1I(KLON,KLEV+1) &
 & ,  ZCLEAR(KLON)         , ZCLOUD(KLON) &
 & ,  ZGG(KLON)            , ZRE1(KLON)   , ZRE2(KLON)&
 & ,  ZR21(KLON)           , ZR22(KLON)&
 & ,  ZSS1(KLON)           , ZTO1(KLON)   , ZTR(KLON,2,KLEV+1)&
 & ,  ZTR1(KLON)           , ZTR2(KLON)&
 & ,  ZW(KLON)  

INTEGER(KIND=JPIM) :: IKL, IKLP1, JA, JAJ, JK, JKM1, JL, INU1

REAL(KIND=JPRB) :: ZBMU0, ZBMU1, ZCORAE, ZCORCD, ZDEN, ZDEN1,&
 & ZFACOA, ZFACOC, ZMUE, ZRE11, &
 & ZALPHA1, ZCHKAE, ZCHKCD  
REAL(KIND=JPRB) :: ZICLEAR, ZDIV1, ZDIV2, ZDIV3, ZDIV4
REAL(KIND=JPRB) :: ZRR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "swde.intfb.h"

!     ------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

IF (LHOOK) CALL DR_HOOK('SWR',0,ZHOOK_HANDLE)
ASSOCIATE(REPSEC=>YDECLD%REPSEC, &
 & RA1OVLP=>YDEOVLP%RA1OVLP, &
 & NOVLP=>YDERAD%NOVLP, NSW=>YDERAD%NSW, &
 & LWSOPT=>YDPHNC%LWSOPT)
DO JK = 1 , KLEV+1
  DO JA = 1 , 6
    DO JL = KIDIA,KFDIA
      PRJ(JL,JA,JK) = 0.0_JPRB
      PRK(JL,JA,JK) = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
!                ----------------------------------------------

DO JL = KIDIA,KFDIA
  ZC1I(JL,KLEV+1) = 0.0_JPRB
  ZCLEAR(JL) = 1.0_JPRB
  ZCLOUD(JL) = 0.0_JPRB
ENDDO

JK = 1
IKL = KLEV+1 - JK
IKLP1 = IKL + 1
ZALPHA1=RA1OVLP( IKL )
DO JL = KIDIA,KFDIA
!++MODIFCODE
  IF (NOVLP >= 5) THEN !MESONH VERSION
    ZFACOA =PTAUAZ(JL,IKL) 
    ZFACOC = 1.0_JPRB - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
    ZCORAE = ZFACOA * PSEC(JL)
    ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
  ELSE !ECMWF VERSION
    ZFACOA = 1.0_JPRB - PPIZAZ(JL,IKL)*PCGAZ(JL,IKL)*PCGAZ(JL,IKL)
    ZFACOC = 1.0_JPRB - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
    ZCORAE = ZFACOA * PTAUAZ(JL,IKL) * PSEC(JL)
    ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
  ENDIF
!--MODIFCODE
  ZCHKAE = MIN( 200._JPRB, ZCORAE )
  ZCHKCD = MIN( 200._JPRB, ZCORCD )
  ZR21(JL) = EXP( - ZCHKAE )
  ZR22(JL) = EXP( - ZCHKCD )
  
  ZSS1(JL) = PCLD(JL,IKL)*(1.0_JPRB-ZR21(JL)*ZR22(JL))&
   & + (1.0_JPRB-PCLD(JL,IKL))*(1.0_JPRB-ZR21(JL))  

!++MODIFCODE
  IF ((NOVLP == 1).OR.(NOVLP == 8)) THEN
!--MODIFCODE
!* maximum-random  
    ZICLEAR = 1.0_JPRB/(1.0_JPRB-MIN(ZCLOUD(JL),1.0_JPRB-REPSEC))
    ZCLEAR(JL) = ZCLEAR(JL)&
     & *(1.0_JPRB-MAX(ZSS1(JL),ZCLOUD(JL))) * ZICLEAR
    ZC1I(JL,IKL) = 1.0_JPRB - ZCLEAR(JL)
    ZCLOUD(JL) = ZSS1(JL)
  ELSEIF (NOVLP == 2) THEN
!* maximum
    ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
    ZC1I(JL,IKL) = ZCLOUD(JL)
!++MODIFCODE
  ELSEIF ((NOVLP == 3).OR.((NOVLP  >=  5).AND.(NOVLP /= 8))) THEN
!--MODIFCODE
!* random
    ZCLEAR(JL) = ZCLEAR(JL)*(1.0_JPRB - ZSS1(JL))
    ZCLOUD(JL) = 1.0_JPRB - ZCLEAR(JL)
    ZC1I(JL,IKL) = ZCLOUD(JL)
  ELSEIF (NOVLP == 4) THEN
!* Hogan & Illingworth, 2001
    ZICLEAR = 1.0_JPRB/(1.0_JPRB-MIN(ZCLOUD(JL),1.0_JPRB-REPSEC))
    ZCLEAR(JL)=ZCLEAR(JL)*(&
     & ZALPHA1*(1.0_JPRB-MAX(ZSS1(JL),ZCLOUD(JL))) * ZICLEAR&
     & +(1.0_JPRB-ZALPHA1)*(1.0_JPRB-ZSS1(JL)) )  
    ZC1I(JL,IKL) = 1.0_JPRB - ZCLEAR(JL) 
    ZCLOUD(JL) = ZSS1(JL)
  ENDIF
ENDDO

DO JK = 2 , KLEV
  IKL = KLEV+1 - JK
  IKLP1 = IKL + 1
  ZALPHA1=RA1OVLP( IKL )
  DO JL = KIDIA,KFDIA
!++MODIFCODE
    IF (NOVLP >= 5) THEN !MESONH VERSION
      ZFACOA =PTAUAZ(JL,IKL) 
      ZFACOC = 1.0_JPRB - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
      ZCORAE = ZFACOA * PSEC(JL)
      ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
    ELSE !ECMWF VERSION
      ZFACOA = 1.0_JPRB - PPIZAZ(JL,IKL)*PCGAZ(JL,IKL)*PCGAZ(JL,IKL)
      ZFACOC = 1.0_JPRB - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,IKL) * PSEC(JL)
      ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
    ENDIF
!--MODIFCODE
!    ZR21(JL) = EXP(-ZCORAE   )
!    ZR22(JL) = EXP(-ZCORCD   )

    ZCHKAE = MIN( 200._JPRB, ZCORAE )
    ZCHKCD = MIN( 200._JPRB, ZCORCD )
    ZR21(JL) = EXP( - ZCHKAE )
    ZR22(JL) = EXP( - ZCHKCD )

    ZSS1(JL) = PCLD(JL,IKL)*(1.0_JPRB-ZR21(JL)*ZR22(JL))&
     & + (1.0_JPRB-PCLD(JL,IKL))*(1.0_JPRB-ZR21(JL))  

!++MODIFCODE
    IF ((NOVLP == 1).OR.(NOVLP == 8)) THEN
!--MODIFCODE
!* maximum-random 
      ZICLEAR = 1.0_JPRB/(1.0_JPRB-MIN(ZCLOUD(JL),1.0_JPRB-REPSEC)) 
      ZCLEAR(JL) = ZCLEAR(JL)&
       & *(1.0_JPRB-MAX(ZSS1(JL),ZCLOUD(JL))) * ZICLEAR
      ZC1I(JL,IKL) = 1.0_JPRB - ZCLEAR(JL)
      ZCLOUD(JL) = ZSS1(JL)
    ELSEIF (NOVLP == 2) THEN
!* maximum
      ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
      ZC1I(JL,IKL) = ZCLOUD(JL)
!++MODIFCODE
  ELSEIF ((NOVLP == 3).OR.((NOVLP  >=  5).AND.(NOVLP /= 8))) THEN
!--MODIFCODE
!* random
      ZCLEAR(JL) = ZCLEAR(JL)*(1.0_JPRB - ZSS1(JL))
      ZCLOUD(JL) = 1.0_JPRB - ZCLEAR(JL)
      ZC1I(JL,IKL) = ZCLOUD(JL)
    ELSEIF (NOVLP == 4) THEN
!* Hogan & Illingworth, 2001
      ZICLEAR = 1.0_JPRB/(1.0_JPRB-MIN(ZCLOUD(JL),1.0_JPRB-REPSEC))
      ZCLEAR(JL)=ZCLEAR(JL)*(&
       & ZALPHA1*(1.0_JPRB-MAX(ZSS1(JL),ZCLOUD(JL))) * ZICLEAR&
       & +(1.0_JPRB-ZALPHA1)*(1.0_JPRB-ZSS1(JL)) )  
      ZC1I(JL,IKL) = 1.0_JPRB - ZCLEAR(JL) 
      ZCLOUD(JL) = ZSS1(JL)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
!                -----------------------------------------------

DO JL = KIDIA,KFDIA
  PRAY1(JL,KLEV+1) = 0.0_JPRB
  PRAY2(JL,KLEV+1) = 0.0_JPRB
  PREFZ(JL,2,1) = PALBD(JL,KNU)
  PREFZ(JL,1,1) = PALBD(JL,KNU)
  PTRA1(JL,KLEV+1) = 1.0_JPRB
  PTRA2(JL,KLEV+1) = 1.0_JPRB
ENDDO

DO JK = 2 , KLEV+1
  JKM1 = JK-1
  DO JL = KIDIA,KFDIA

!     ------------------------------------------------------------------

!*         3.1  EQUIVALENT ZENITH ANGLE
!               -----------------------

    ZMUE = (1.0_JPRB-ZC1I(JL,JK)) * PSEC(JL)+ ZC1I(JL,JK) * 1.66_JPRB
    PRMUE(JL,JK) = 1.0_JPRB/ZMUE

!     ------------------------------------------------------------------

!*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
!               ----------------------------------------------------

    ZBMU0 = 0.5_JPRB - 0.75_JPRB * PCGAZ(JL,JKM1) * PRMUE(JL,JK)
    ZDEN = 1.0_JPRB + (1.0_JPRB - PPIZAZ(JL,JKM1)&
     & + ZBMU0 * PPIZAZ(JL,JKM1)) * PTAUAZ(JL,JKM1) * ZMUE&
     & + (1-PPIZAZ(JL,JKM1)) * (1.0_JPRB - PPIZAZ(JL,JKM1)&
     & + 2.0_JPRB*ZBMU0*PPIZAZ(JL,JKM1))*PTAUAZ(JL,JKM1)*PTAUAZ(JL,JKM1)&
     & * ZMUE*ZMUE  
    PTRA1(JL,JKM1) = 1.0_JPRB/ZDEN
    PRAY1(JL,JKM1) = ZBMU0 * PPIZAZ(JL,JKM1) * PTAUAZ(JL,JKM1)&
     & * ZMUE * PTRA1(JL,JKM1)

    ZBMU1 = 0.5_JPRB - 0.75_JPRB * PCGAZ(JL,JKM1) * 0.5_JPRB
    ZDEN1= 1.0_JPRB + (1.0_JPRB - PPIZAZ(JL,JKM1)&
     & + ZBMU1 * PPIZAZ(JL,JKM1)) * PTAUAZ(JL,JKM1) * 2.0_JPRB&
     & + (1-PPIZAZ(JL,JKM1)) * (1.0_JPRB - PPIZAZ(JL,JKM1)&
     & + 2.0_JPRB*ZBMU1*PPIZAZ(JL,JKM1))*PTAUAZ(JL,JKM1)&
     & * PTAUAZ(JL,JKM1)*4.0_JPRB 
    PTRA2(JL,JKM1) = 1.0_JPRB/ZDEN1
    PRAY2(JL,JKM1) = ZBMU1 * PPIZAZ(JL,JKM1) * PTAUAZ(JL,JKM1)&
     & * 2.0_JPRB * PTRA2(JL,JKM1)

!     ------------------------------------------------------------------

!*         3.3  EFFECT OF CLOUD LAYER
!               ---------------------


!++MODIFCODE
    IF (NOVLP >= 5)THEN !MESONH VERSION
      ZW(JL) =PCG(JL,KNU,JKM1)*PCG(JL,KNU,JKM1)
      ZTO1(JL) = PTAU(JL,KNU,JKM1)*(1-(POMEGA(JL,KNU,JKM1)*ZW(JL)))
      ZW(JL) =POMEGA(JL,KNU,JKM1)*(1-ZW(JL))/(1-(POMEGA(JL,KNU,JKM1)*ZW(JL)))
      ZGG(JL) = PCG(JL,KNU,JKM1)/(1+PCG(JL,KNU,JKM1))
      ZGG(JL)=ZTO1(JL)*ZW(JL)*ZGG(JL)+PTAUAZ(JL,JKM1)*PPIZAZ(JL,JKM1)&
       & *PCGAZ(JL,JKM1)
      ZW(JL) =ZTO1(JL)*ZW(JL)+PTAUAZ(JL,JKM1)*PPIZAZ(JL,JKM1)
      ZTO1(JL) = ZTO1(JL) +  PTAUAZ(JL,JKM1)
      ZGG(JL)=ZGG(JL)/ZW(JL)
      ZW(JL) =ZW(JL)/ZTO1(JL)
    ELSE !ECMWF VERSION
      ZDIV1 = 1.0_JPRB/POMEGA(JL,KNU,JKM1)
      ZDIV2 = 1.0_JPRB/PPIZAZ(JL,JKM1)
      ZTO1(JL) = PTAU(JL,KNU,JKM1)*ZDIV1 + PTAUAZ(JL,JKM1)*ZDIV2
      ZR21(JL) = PTAU(JL,KNU,JKM1) + PTAUAZ(JL,JKM1)
      ZDIV3 = 1.0_JPRB/ZR21(JL) 
      ZR22(JL) = PTAU(JL,KNU,JKM1)*ZDIV3
      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)&
      & + (1.0_JPRB - ZR22(JL)) * PCGAZ(JL,JKM1)  
      IF (POMEGA(JL,KNU,JKM1)==1.0_JPRB .AND. PPIZAZ(JL,JKM1)==1.0_JPRB) THEN
        ZW(JL)=1.0_JPRB
      ELSE
        ZDIV4 = 1.0_JPRB/ZTO1(JL)
        ZW(JL) = ZR21(JL) * ZDIV4
      ENDIF
    ENDIF
!--MODIFCODE

  ENDDO

  CALL SWDE(YDERAD,KIDIA,KFDIA,KLON,ZGG,PREFZ(1,1,JKM1),PRMUE(1,JK),ZTO1,ZW,ZRE1,ZRE2,ZTR1,ZTR2)

  DO JL = KIDIA,KFDIA

    IF (LWSOPT) THEN
! saved for adjoint computation
      PRE1(JL,JKM1) = ZRE1(JL)
      PRE2(JL,JKM1) = ZRE2(JL)
      PTR1(JL,JKM1) = ZTR1(JL)
      PTR2(JL,JKM1) = ZTR2(JL)
    ENDIF

    ZRR=1.0_JPRB/(1.0_JPRB-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1))

    PREFZ(JL,1,JK) = (1.0_JPRB-PCLD(JL,JKM1)) * (PRAY1(JL,JKM1)&
     & + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1) * PTRA2(JL,JKM1) * ZRR )&
     & + PCLD(JL,JKM1) * ZRE2(JL)  

    ZTR(JL,1,JKM1) = PCLD(JL,JKM1) * ZTR2(JL) + (PTRA1(JL,JKM1)&
     & * ZRR )&
     & * (1.0_JPRB-PCLD(JL,JKM1))  

    PREFZ(JL,2,JK) = (1.0_JPRB-PCLD(JL,JKM1)) * (PRAY1(JL,JKM1)&
     & + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
     & * PTRA2(JL,JKM1) )&
     & + PCLD(JL,JKM1) * ZRE1(JL)  

    ZTR(JL,2,JKM1) = PCLD(JL,JKM1)*ZTR1(JL)&
     & + PTRA1(JL,JKM1)*(1.0_JPRB-PCLD(JL,JKM1))

  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZMUE = (1.0_JPRB-ZC1I(JL,1))*PSEC(JL)+ZC1I(JL,1)*1.66_JPRB
  PRMUE(JL,1)=1.0_JPRB/ZMUE
  PTRCLD(JL)=1.0_JPRB-ZC1I(JL,1)
ENDDO

!     ------------------------------------------------------------------

!*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!                 -------------------------------------------------

IF (NSW <= 4) THEN
  INU1=1
ELSEIF (NSW == 6) THEN
  INU1=3
ENDIF    

IF (KNU <= INU1) THEN
  JAJ = 2
  DO JL = KIDIA,KFDIA
    PRJ(JL,JAJ,KLEV+1) = 1.0_JPRB
    PRK(JL,JAJ,KLEV+1) = PREFZ(JL, 1,KLEV+1)
  ENDDO

  DO JK = 1 , KLEV
    IKL = KLEV+1 - JK
    IKLP1 = IKL + 1
    DO JL = KIDIA,KFDIA
      ZRE11= PRJ(JL,JAJ,IKLP1) * ZTR(JL,  1,IKL)
      PRJ(JL,JAJ,IKL) = ZRE11
      PRK(JL,JAJ,IKL) = ZRE11 * PREFZ(JL,  1,IKL)
    ENDDO
  ENDDO

ELSE

  DO JAJ = 1 , 2
    DO JL = KIDIA,KFDIA
      PRJ(JL,JAJ,KLEV+1) = 1.0_JPRB
      PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      IKL = KLEV+1 - JK
      IKLP1 = IKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11= PRJ(JL,JAJ,IKLP1) * ZTR(JL,JAJ,IKL)
        PRJ(JL,JAJ,IKL) = ZRE11
        PRK(JL,JAJ,IKL) = ZRE11 * PREFZ(JL,JAJ,IKL)
      ENDDO
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SWR',1,ZHOOK_HANDLE)
END SUBROUTINE SWR

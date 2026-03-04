! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE UVR &
 &( YDEOVLP,YDERAD,YDECLD,KIDIA , KFDIA , KLON , KLEV , &
 &  PALBD , PCG   , PCLD , POMEGA, PSEC , PTAU,&
 &  PCGAZ , PPIZAZ, PRAY1, PRAY2 , PREFZ, PRJ  , PRK , PRMUE,&
 &  PTAUAZ, PTRA1 , PTRA2, PTRCLD &
 &)  

!**** *UVR* - CONTINUUM SCATTERING COMPUTATIONS

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
!     CONTINUUM SCATTERING

!**   INTERFACE.
!     ----------

!          *UVR* IS CALLED FROM *UVFLX*

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

!          *UVDE*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2005-10-04 from SWR (but top to bottom)
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERAD   , ONLY : TERAD
USE YOECLD   , ONLY : TECLD
USE YOEOVLP  , ONLY : TEOVLP

IMPLICIT NONE

TYPE(TECLD)       ,INTENT(INOUT):: YDECLD
TYPE(TEOVLP)      ,INTENT(INOUT):: YDEOVLP
TYPE(TERAD)       ,INTENT(INOUT):: YDERAD
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PALBD(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCLD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: POMEGA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PSEC(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PTAU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCGAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PPIZAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PTAUAZ(KLON,KLEV) 

REAL(KIND=JPRB)   ,INTENT(OUT) :: PRAY1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PRAY2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PREFZ(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PRJ(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PRK(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PRMUE(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PTRA1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PTRA2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PTRCLD(KLON) 
!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!     ------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZC1I(KLON,KLEV+1)    , ZCLEQ(KLON,KLEV)&
 & ,  ZCLEAR(KLON)         , ZCLOUD(KLON) &
 & ,  ZGG(KLON)            , ZREF(KLON)&
 & ,  ZRE1(KLON)           , ZRE2(KLON)&
 & ,  ZRMUZ(KLON)          , ZRNEB(KLON)&
 & ,  ZR21(KLON)           , ZR22(KLON)&
 & ,  ZR23(KLON)           , ZSS1(KLON)&
 & ,  ZTO1(KLON)           , ZTR(KLON,2,KLEV+1)&
 & ,  ZTR1(KLON)           , ZTR2(KLON)&
 & ,  ZW(KLON)  

INTEGER(KIND=JPIM) :: JA, JAJ, JK, JKP1, JL

REAL(KIND=JPRB) :: ZBMU0, ZBMU1, ZCORAE, ZCORCD, ZDEN, ZDEN1,&
 & ZFACOA, ZFACOC, ZGAP, ZMU1, ZMUE, ZRE11, &
 & ZTO, ZWW, ZALPHA1, ZCHKAE, ZCHKCD  
REAL(KIND=JPRB) :: ZRR,ZIMU1,ZI2MU1,ZIDEN,ZIDEN1

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "uvde.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UVR',0,ZHOOK_HANDLE)
ASSOCIATE(REPSEC=>YDECLD%REPSEC, &
 & RA1OVLP=>YDEOVLP%RA1OVLP, &
 & NOVLP=>YDERAD%NOVLP)
!*         1.    INITIALIZATION
!                --------------

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
  ZR23(JL) = 0.0_JPRB
  ZC1I(JL,1) = 0.0_JPRB
  ZCLEAR(JL) = 1.0_JPRB
  ZCLOUD(JL) = 0.0_JPRB
ENDDO

DO JK = 1 , KLEV
  ZALPHA1=RA1OVLP( JK )
  DO JL = KIDIA,KFDIA
    ZFACOA = 1.0_JPRB - PPIZAZ(JL,JK)*PCGAZ(JL,JK)*PCGAZ(JL,JK)
    ZFACOC = 1.0_JPRB - POMEGA(JL,JK) * PCG(JL,JK)* PCG(JL,JK)
    ZCORAE = ZFACOA * PTAUAZ(JL,JK) * PSEC(JL)
    ZCORCD = ZFACOC * PTAU(JL,JK) * PSEC(JL)

    ZCHKAE = MIN( 200._JPRB, ZCORAE )
    ZCHKCD = MIN( 200._JPRB, ZCORCD )
    ZR21(JL) = EXP( - ZCHKAE )
    ZR22(JL) = EXP( - ZCHKCD )
  
    ZSS1(JL) = PCLD(JL,JK)*(1.0_JPRB-ZR21(JL)*ZR22(JL))&
     & + (1.0_JPRB-PCLD(JL,JK))*(1.0_JPRB-ZR21(JL))  
    ZCLEQ(JL,JK) = ZSS1(JL)

    IF (NOVLP == 1) THEN
!* maximum-random      
      ZCLEAR(JL) = ZCLEAR(JL)&
       & *(1.0_JPRB-MAX(ZSS1(JL),ZCLOUD(JL)))&
       & /(1.0_JPRB-MIN(ZCLOUD(JL),1.0_JPRB-REPSEC))  
      ZC1I(JL,JK+1) = 1.0_JPRB - ZCLEAR(JL)
      ZCLOUD(JL) = ZSS1(JL)
    ELSEIF (NOVLP == 2) THEN
!* maximum
      ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
      ZC1I(JL,JK+1) = ZCLOUD(JL)
    ELSEIF (NOVLP == 3) THEN
!* random
      ZCLEAR(JL) = ZCLEAR(JL)*(1.0_JPRB - ZSS1(JL))
      ZCLOUD(JL) = 1.0_JPRB - ZCLEAR(JL)
      ZC1I(JL,JK+1) = ZCLOUD(JL)
    ELSEIF (NOVLP == 4) THEN
!* Hogan & Illingworth, 2001  
      ZCLEAR(JL)=ZCLEAR(JL)*(&
       & ZALPHA1*(1.0_JPRB-MAX(ZSS1(JL),ZCLOUD(JL)))&
       & /(1.0_JPRB-MIN(ZCLOUD(JL),1.0_JPRB-REPSEC))&
       & +(1.0_JPRB-ZALPHA1)*(1.0_JPRB-ZSS1(JL)) )  
      ZC1I(JL,JK+1) = 1.0_JPRB - ZCLEAR(JL) 
      ZCLOUD(JL) = ZSS1(JL)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
!                -----------------------------------------------

DO JL = KIDIA,KFDIA
  PRAY1(JL,1) = 0.0_JPRB
  PRAY2(JL,1) = 0.0_JPRB
  PREFZ(JL,2,KLEV+1) = PALBD(JL)
  PREFZ(JL,1,KLEV+1) = PALBD(JL)
  PTRA1(JL,1) = 1.0_JPRB
  PTRA2(JL,1) = 1.0_JPRB
ENDDO

DO JK = KLEV, 1 , -1
  JKP1 = JK+1
  DO JL = KIDIA,KFDIA
    ZRNEB(JL)= PCLD(JL,JK)
    ZRE1(JL)=0.0_JPRB
    ZTR1(JL)=0.0_JPRB
    ZRE2(JL)=0.0_JPRB
    ZTR2(JL)=0.0_JPRB

!     ------------------------------------------------------------------

!*         3.1  EQUIVALENT ZENITH ANGLE
!               -----------------------

    ZMUE = (1.0_JPRB-ZC1I(JL,JK)) * PSEC(JL)+ ZC1I(JL,JK) * 1.66_JPRB
!-- no fiddling
!    ZMUE = PSEC(JL)
    PRMUE(JL,JK) = 1.0_JPRB/ZMUE

!     ------------------------------------------------------------------

!*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
!               ----------------------------------------------------

    ZGAP = PCGAZ(JL,JK)
    ZBMU0 = 0.5_JPRB - 0.75_JPRB * ZGAP / ZMUE
    ZWW = PPIZAZ(JL,JK)
    ZTO = PTAUAZ(JL,JK)
    ZDEN = 1.0_JPRB + (1.0_JPRB - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE&
     & + (1-ZWW) * (1.0_JPRB - ZWW +2.0_JPRB*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE  
    ZIDEN=1.0_JPRB/ZDEN
    PRAY1(JL,JKP1) = ZBMU0 * ZWW * ZTO * ZMUE * ZIDEN
    PTRA1(JL,JKP1) = ZIDEN

    ZMU1  = 0.5_JPRB
    ZIMU1 = 2.0_JPRB
    ZI2MU1= 4.0_JPRB
    ZBMU1 = 0.5_JPRB - 0.75_JPRB * ZGAP * ZMU1
    ZDEN1= 1.0_JPRB + (1.0_JPRB - ZWW + ZBMU1 * ZWW) * ZTO * ZIMU1&
     & + (1-ZWW) * (1.0_JPRB - ZWW +2.0_JPRB*ZBMU1*ZWW)*ZTO*ZTO*ZI2MU1  
    ZIDEN1=1.0_JPRB/ZDEN1
    PRAY2(JL,JKP1) = ZBMU1 * ZWW * ZTO * ZIMU1 * ZIDEN1
    PTRA2(JL,JKP1) = ZIDEN1

!     ------------------------------------------------------------------

!*         3.3  EFFECT OF CLOUD LAYER
!               ---------------------

    ZW(JL) = POMEGA(JL,JK)
    ZTO1(JL) = PTAU(JL,JK)/ZW(JL)+ PTAUAZ(JL,JK)/PPIZAZ(JL,JK)
    ZR21(JL) = PTAU(JL,JK) + PTAUAZ(JL,JK)
    ZR22(JL) = PTAU(JL,JK) / ZR21(JL)
    ZGG(JL) = ZR22(JL) * PCG(JL,JK)&
     & + (1.0_JPRB - ZR22(JL)) * PCGAZ(JL,JK)  
    IF (ZW(JL) == 1.0_JPRB .AND. PPIZAZ(JL,JK) == 1.0_JPRB) THEN
      ZW(JL)=1.0_JPRB
    ELSE
      ZW(JL) = ZR21(JL) / ZTO1(JL)
    ENDIF
    ZREF(JL) = PREFZ(JL,1,JKP1)
    ZRMUZ(JL) = PRMUE(JL,JK)
  ENDDO
!  JL=KIDIA
!  print 9001,ZGG(JL),ZREF(JL),ZRMUZ(JL),ZTO1(JL),ZW(JL)
9001 FORMAT(1X,'UVDE G Re Mu To W ',5E12.5)

  CALL UVDE ( KIDIA, KFDIA , KLON,&
   & ZGG  , ZREF  , ZRMUZ , ZTO1 , ZW,&
   & ZRE1 , ZRE2  , ZTR1  , ZTR2      )
     
!  JL=KIDIA
!  print 9002,ZRE1(JL),ZRE2(JL),ZTR1(JL),ZTR2(JL)
9002 FORMAT(1X,'UVDE Re1 Re2 Tr1 Tr2 ',5E12.5)
  DO JL = KIDIA,KFDIA

    ZRR=1.0_JPRB/(1.0_JPRB-PRAY2(JL,JKP1)*PREFZ(JL,1,JKP1))

    PREFZ(JL,1,JK) = (1.0_JPRB-ZRNEB(JL)) * (PRAY1(JL,JKP1)&
     & + PREFZ(JL,1,JKP1) * PTRA1(JL,JKP1)&
     & * PTRA2(JL,JKP1)&
     & * ZRR )&
     & + ZRNEB(JL) * ZRE2(JL)  

    ZTR(JL,1,JKP1) = ZRNEB(JL) * ZTR2(JL) + (PTRA1(JL,JKP1) * ZRR )&
     & * (1.0_JPRB-ZRNEB(JL))  

    PREFZ(JL,2,JK) = (1.0_JPRB-ZRNEB(JL)) * (PRAY1(JL,JKP1)&
     & + PREFZ(JL,2,JKP1) * PTRA1(JL,JKP1)&
     & * PTRA2(JL,JKP1) )&
     & + ZRNEB(JL) * ZRE1(JL)  

    ZTR(JL,2,JKP1) = ZRNEB(JL) * ZTR1(JL)+ PTRA1(JL,JKP1) * (1.0_JPRB-ZRNEB(JL))

!    print 9003,JK,PREFZ(JL,1,JK),ZTR(JL,1,JKP1),PREFZ(JL,2,JK),ZTR(JL,2,JKP1),&
!      &PRAY1(JL,JKP1),PREFZ(JL,1,JKP1),PTRA1(JL,JKP1),PTRA2(JL,JKP1),ZRR
9003 FORMAT(1X,'UVR:',10F10.6)
  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZMUE = (1.0_JPRB-ZC1I(JL,KLEV+1))*PSEC(JL)+ZC1I(JL,KLEV+1)*1.66_JPRB
!-- no fiddling
!  ZMUE = PSEC(JL)
  PRMUE(JL,KLEV+1)=1.0_JPRB/ZMUE
  PTRCLD(JL)=1.0_JPRB-ZC1I(JL,KLEV+1)
ENDDO

!     ------------------------------------------------------------------

!*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!                 -------------------------------------------------

JAJ = 2
DO JL = KIDIA,KFDIA
  PRJ(JL,JAJ,1) = 1.0_JPRB
  PRK(JL,JAJ,1) = PREFZ(JL,1,1)
ENDDO

DO JK = 1 , KLEV
  DO JL = KIDIA,KFDIA
    ZRE11= PRJ(JL,JAJ,JK) * ZTR(JL,1,JK+1)
    PRJ(JL,JAJ,JK+1) = ZRE11
    PRK(JL,JAJ,JK+1) = ZRE11 * PREFZ(JL,1,JK+1)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UVR',1,ZHOOK_HANDLE)
END SUBROUTINE UVR


! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE UVFLXA &
 &( YDML_PHY_RAD,YDECLD,KIDIA, KFDIA , KLON  , KLEV  , KUV , &
 &  PAKI , PALBD , PALBP , PCG   , PCLD, PCLEAR, &
 &  PDSIG, POMEGA, PRAYL , PRMU  , PSEC, &
 &  PTAU , PCGAZ , PPIZAZ, PTAUAZ, PUD , &
 &  PFDOWN,PCDOWN &
 & )  

!**** *UVFLXA* - SHORTWAVE RADIATION, UV SPECTRAL INTERVAL (v.2)

!     PURPOSE.
!     --------

!          COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE NEAR-INFRARED 
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *UVFLXA* IS CALLED FROM *UVRADI*.

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
!     CONTINUUM SCATTERING
!          2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
!     A GREY MOLECULAR ABSORPTION
!          3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
!     OF ABSORBERS
!          4. APPLY OZONE TRANSMISSION FUNCTIONS

!     EXTERNALS.
!     ----------

!          *SWCLR*, *SWR*, *SWDE*, *SWTT*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2005-11-14
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_RADIATION_MOD , ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE YOECLD   , ONLY : TECLD
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERDU   , ONLY : REPLOG 

IMPLICIT NONE

TYPE(TECLD)       ,INTENT(INOUT):: YDECLD
TYPE(MODEL_PHYSICS_RADIATION_TYPE),INTENT(INOUT):: YDML_PHY_RAD
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KUV 

REAL(KIND=JPRB)   ,INTENT(IN)  :: PAKI(KLON,2) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PALBD(KLON), PALBP(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCLD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PCLEAR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDSIG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: POMEGA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PRAYL(KLON), PRMU(KLON), PSEC(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PTAU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PTAUAZ(KLON,KLEV),PCGAZ(KLON,KLEV),PPIZAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PUD(KLON,2,KLEV+1) 

REAL(KIND=JPRB)   ,INTENT(OUT) :: PFDOWN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PCDOWN(KLON,KLEV+1) 

!#include "yoeaer.h"
!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!     ------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: &
! &    ZFD(KLON,KLEV+1)  , ZFU(KLON,KLEV+1) &
 &    ZFD(KLON,KLEV+1)   &
! & ,  ZFDOWN(KLON,KLEV+1),ZFUP(KLON,KLEV+1)&
 & ,  ZFDOWN(KLON,KLEV+1)&
 & ,  ZG(KLON)          , ZGG(KLON)  
REAL(KIND=JPRB) :: &
 &    ZRAY1(KLON,KLEV+1), ZRAY2(KLON,KLEV+1)&
 & ,  ZREF(KLON)        , ZREFZ(KLON,2,KLEV+1)&
 & ,  ZRE1(KLON)        , ZRE2(KLON)&
 & ,  ZRJ(KLON,6,KLEV+1), ZRJ0(KLON,6,KLEV+1)&
 & ,  ZRK(KLON,6,KLEV+1), ZRK0(KLON,6,KLEV+1)&
 & ,  ZRL(KLON,8)&
 & ,  ZRMUE(KLON,KLEV+1), ZRMU0(KLON,KLEV+1)  , ZRMUZ(KLON)&
 & ,  ZRNEB(KLON)       , ZRUEF(KLON,8)       &
 & ,  ZR2(KLON,2)       , ZR3(KLON,6)         &
 & ,  ZR21(KLON)        , ZR22(KLON)  
REAL(KIND=JPRB) :: ZS(KLON)&
 & ,  ZTO1(KLON)        , ZTR(KLON,2,KLEV+1)&
 & ,  ZTRA1(KLON,KLEV+1), ZTRA2(KLON,KLEV+1)&
 & ,  ZTRCLD(KLON)      , ZTRCLR(KLON)&
 & ,  ZTR1(KLON)        , ZTR2(KLON)&
 & ,  ZW(KLON)          , ZW2(KLON,2)&
 & ,  ZW3(KLON,6)  

INTEGER(KIND=JPIM) :: JABS, JAJ, JAJP, JK, JKKI,&
 & JKKP4, JKM1, JKP1, JL, JN, JN2J, JREF  

REAL(KIND=JPRB) :: ZAA, ZRE11, ZRKI, ZCHKG, ZCHKS
REAL(KIND=JPRB) :: ZRR, ZRRJ 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "uvclr.intfb.h"
#include "uvde.intfb.h"
#include "uvr.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UVFLXA',0,ZHOOK_HANDLE)
ASSOCIATE(RK250=>YDML_PHY_RAD%YREUVRAD%RK250)
!*         1.     VALID FOR SPECTRAL INTERVAL WITHIN 280-400 NM
!                 ---------------------------------------------

!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING  (in *UVRADI*)
!                 -------------------------------------------------------

!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------

!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------

CALL UVCLR &
 &( YDML_PHY_RAD%YRERAD,YDML_PHY_RAD%YRERDI, &
 &  KIDIA , KFDIA , KLON  , KLEV , &
 &  PALBP , PDSIG , PRAYL , PSEC , &
 &  PCGAZ , PPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 , &
 &  ZRK0  , ZRMU0 , PTAUAZ, ZTRA1, ZTRA2, ZTRCLR &
 &)  

!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------

CALL UVR &
 &( YDML_PHY_RAD%YREOVLP,YDML_PHY_RAD%YRERAD,YDECLD, &
 &  KIDIA , KFDIA , KLON , KLEV  , &
 &  PALBD , PCG   , PCLD , POMEGA, PSEC , PTAU , &
 &  PCGAZ , PPIZAZ, ZRAY1, ZRAY2 , ZREFZ, ZRJ  , ZRK, ZRMUE , &
 &  PTAUAZ, ZTRA1 , ZTRA2, ZTRCLD &
 & )  

!     ------------------------------------------------------------------

!*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
!                ------------------------------------------------------
!                 O3 : JABS=1 ; NO2 : JABS=2

ZRL(:,:)=0._JPRB

JN = 2

DO JABS=1,1

!*         3.1  SURFACE CONDITIONS
!               ------------------

  DO JL = KIDIA,KFDIA
    ZREFZ(JL,2,KLEV+1) = PALBD(JL)
    ZREFZ(JL,1,KLEV+1) = PALBD(JL)
  ENDDO

!*         3.2  INTRODUCING CLOUD EFFECTS
!               -------------------------

  DO JK = KLEV,1,-1
    JKP1 = JK+1
    DO JL = KIDIA,KFDIA
      ZRNEB(JL) = PCLD(JL,JK)
      ZAA=PUD(JL,JABS,JK)
      ZRKI = PAKI(JL,JABS)
      
      ZCHKS = MIN( 200._JPRB, ZRKI * ZAA * 1.66_JPRB )
      ZCHKG = MIN( 200._JPRB, ZRKI * ZAA / ZRMUE(JL,JK) )
      ZS(JL) = EXP( - ZCHKS )
      ZG(JL) = EXP( - ZCHKG )
      
      ZTR1(JL) = 0.0_JPRB
      ZRE1(JL) = 0.0_JPRB
      ZTR2(JL) = 0.0_JPRB
      ZRE2(JL) = 0.0_JPRB

      ZW(JL)= POMEGA(JL,JK)
      ZTO1(JL) = PTAU(JL,JK) / ZW(JL)&
       & + PTAUAZ(JL,JK) / PPIZAZ(JL,JK)&
       & + ZAA * ZRKI  

      ZR21(JL) = PTAU(JL,JK) + PTAUAZ(JL,JK)
      ZR22(JL) = PTAU(JL,JK) / ZR21(JL)
      ZGG(JL) = ZR22(JL) * PCG(JL,JK)&
       & + (1.0_JPRB - ZR22(JL)) * PCGAZ(JL,JK)  
      ZW(JL) = ZR21(JL) / ZTO1(JL)
      ZREF(JL) = ZREFZ(JL,1,JK+1)
      ZRMUZ(JL) = ZRMUE(JL,JK)
    ENDDO

!    JL=KIDIA
!    print 9001,ZGG(JL),ZREF(JL),ZRMUZ(JL),ZTO1(JL),ZW(JL),ZRNEB(JL), &
!& ZRKI,PTAU(JL,JK),POMEGA(JL,JK),PCG(JL,JK),PTAUAZ(JL,JK),PPIZAZ(JL,JK), &
!& PCGAZ(JL,JK)
9001 FORMAT(1X,'FLXA G Re Mu To W ',13F10.5)

    CALL UVDE ( KIDIA, KFDIA, KLON,&
     & ZGG  , ZREF , ZRMUZ, ZTO1, ZW,&
     & ZRE1 , ZRE2 , ZTR1 , ZTR2     )  

!    JL=KIDIA
!    print 9002,ZRE1(JL),ZRE2(JL),ZTR1(JL),ZTR2(JL)
9002 FORMAT(1X,'FLXA Re1 Re2 Tr1 Tr2 ',5E12.5)

    DO JL = KIDIA,KFDIA

      ZRR=1.0_JPRB/(1.0_JPRB-ZRAY2(JL,JKP1)*ZREFZ(JL,1,JKP1))

      ZREFZ(JL,2,JK) = (1.0_JPRB-ZRNEB(JL)) * (ZRAY1(JL,JKP1)&
       & + ZREFZ(JL,2,JKP1) * ZTRA1(JL,JKP1)&
       & * ZTRA2(JL,JKP1) ) * ZG(JL) * ZS(JL)&
       & + ZRNEB(JL) * ZRE1(JL)  

      ZTR(JL,2,JKP1)=ZRNEB(JL) * ZTR1(JL)&
       & + (ZTRA1(JL,JKP1)) * ZG(JL) * (1.0_JPRB-ZRNEB(JL))  

      ZREFZ(JL,1,JK)=(1.0_JPRB-ZRNEB(JL))*(ZRAY1(JL,JKP1)&
       & +ZREFZ(JL,1,JKP1)*ZTRA1(JL,JKP1)*ZTRA2(JL,JKP1)&
       & *ZRR ) &
       & *ZG(JL)*ZS(JL)&
       & + ZRNEB(JL) * ZRE2(JL)  

      ZTR(JL,1,JKP1)= ZRNEB(JL) * ZTR2(JL)&
       & + (ZTRA1(JL,JKP1) *ZRR ) * ZG(JL) * (1.0_JPRB -ZRNEB(JL))  

    ENDDO
!    JL=KIDIA
!    print 9003,JK,ZREFZ(JL,1,JK),ZTR(JL,1,JKP1),ZREFZ(JL,2,JK),ZTR(JL,2,JKP1),ZRNEB(JL),ZRAY1(JL,JKP1),&
!      & ZRAY2(JL,JKP1),ZTRA1(JL,JKP1),ZTRA2(JL,JKP1),ZS(JL),ZG(JL)
9003 FORMAT(1X,'FLXA',I3,12F10.7)

  ENDDO

!*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!               -------------------------------------------------

  DO JREF=1,2

    JN = JN + 1

    DO JL = KIDIA,KFDIA
      ZRJ(JL,JN,1) = 1.0_JPRB
      ZRK(JL,JN,1) = ZREFZ(JL,JREF,1)
    ENDDO

    DO JK = 2 , KLEV+1
      JKM1 = JK-1
      DO JL = KIDIA,KFDIA
        ZRE11 = ZRJ(JL,JN,JKM1) * ZTR(JL,JREF,JK)
        ZRJ(JL,JN,JK) = ZRE11
        ZRK(JL,JN,JK) = ZRE11 * ZREFZ(JL,JREF,JK)
      ENDDO
!      JL=KIDIA
!      print 9004,JABS,JREF,JK,ZTR(JL,JREF,JK),ZRJ(JL,JN,JKM1),ZRJ(JL,JN,JK),ZRK(JL,JN,JK),ZREFZ(JL,JREF,JK)
9004  FORMAT(1X,'FLXR',3I3,12F10.7)
    ENDDO
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         4.    INVERT GREY AND CONTINUUM FLUXES
!                --------------------------------

!*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
!                ---------------------------------------------

DO JK = 1 , KLEV+1
  DO JAJ = 1 , 5 , 2
    JAJP = JAJ + 1
    DO JL = KIDIA,KFDIA
      ZRJ(JL,JAJ,JK)=        ZRJ(JL,JAJ,JK) - ZRJ(JL,JAJP,JK)
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , REPLOG )
    ENDDO
!    JL=KIDIA
!    print 9005,JAJ,JK,ZRJ(JL,JAJ,JKM1)
9005 FORMAT(1X,'FLXS',2I3,12F10.7)
  ENDDO
ENDDO

DO JK = 1 , KLEV+1
  DO JAJ = 2 , 6 , 2
    DO JL = KIDIA,KFDIA
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , REPLOG )
    ENDDO
!    JL=KIDIA
!    print 9006,JAJ,JK,ZRJ(JL,JAJ,JKM1)
9006 FORMAT(1X,'FLXT',2I3,12F10.7)
  ENDDO
ENDDO

!*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
!                 ---------------------------------------------

DO JK = 1 , KLEV+1
  JKKI = 1
  DO JAJ = 1 , 1
    DO JN = 1 , 2
      JN2J = JN + 2 * JAJ
      JKKP4 = JKKI + 4

!*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
!                 --------------------------

      DO JL = KIDIA,KFDIA
        ZRR=1.0_JPRB/PAKI(JL,JAJ)
        ZRRJ=ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK)
        ZW2(JL,1) = LOG( ZRRJ ) * ZRR

!*         4.2.2  TRANSMISSION FUNCTION
!                 ---------------------

        ZR2(JL,1)=EXP(-RK250(KUV)*ZW2(JL,1))

        ZRL(JL,JKKI) = ZR2(JL,1)
        ZRUEF(JL,JKKI) = ZW2(JL,1)
      ENDDO
!      JL=KIDIA
!      print 9007,JKKI,JAJ,JN,JK,ZW2(JL,1),ZR2(JL,1),ZRL(JL,JKKI),ZRUEF(JL,JKKI)
9007  FORMAT(1X,'FLXU',4I3,12E12.5)

      JKKI=JKKI+1
    ENDDO
  ENDDO

!*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
!                 ------------------------------------------------------

  DO JL = KIDIA,KFDIA
!!!    ZFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)
    ZFDOWN(JL,JK) = ZRJ(JL,2,JK) * ZRL(JL,2)   
  ENDDO
!  JL=KIDIA
!  print 9008,JK,ZRJ(JL,1,JK),ZRL(JL,1),ZRL(JL,3),ZRJ(JL,2,JK),ZRL(JL,2),ZRL(JL,4),ZFDOWN(JL,JK)
9008 FORMAT(1X,'FLXV',I3,6E12.5)
ENDDO

!     ------------------------------------------------------------------

!*         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
!                ----------------------------------------

!*         5.1   DOWNWARD FLUXES
!                ---------------

JAJ = 2

DO JL = KIDIA,KFDIA
  ZW3(JL,1)=0.0_JPRB
  ZFD(JL,1)= ZRJ0(JL,JAJ,1)
ENDDO
DO JK = 1 , KLEV
  JKP1=JK+1
  DO JL = KIDIA,KFDIA
    ZRR=1.0_JPRB/ZRMU0(JL,JK)
    ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,JK)*ZRR

    ZR3(JL,1)=EXP(-RK250(KUV)*ZW3(JL,1))
 
    ZFD(JL,JKP1) = ZR3(JL,1) * ZRJ0(JL,JAJ,JKP1)
  ENDDO
ENDDO


!     ------------------------------------------------------------------

!*         6.     FINAL FLUXES
!                 ------------

DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    PFDOWN(JL,JK)=(1.0_JPRB-PCLEAR(JL))*ZFDOWN(JL,JK)+PCLEAR(JL)*ZFD(JL,JK)
    PCDOWN(JL,JK)=ZFD(JL,JK)
  ENDDO
ENDDO
!JL=KIDIA
!print 9011,KUV,PCLEAR(JL),PCDOWN(JL,KLEV+1),PFDOWN(JL,KLEV+1),ZFDOWN(JL,KLEV+1)
9011 FORMAT(1X,'UVFLXA',I5,F10.7,3E12.5)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UVFLXA',1,ZHOOK_HANDLE)
END SUBROUTINE UVFLXA


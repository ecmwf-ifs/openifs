! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SWNI &
 & ( YDML_PHY_RAD,YDPHNC,YDECLD,KIDIA , KFDIA , KLON  , KLEV , KAER  , KNU,&
 & PAER  , PAKI  , PALBD , PALBP, PCG   , PCLD, PCLEAR,&
 & PDSIG , POMEGA, POZ   , PRMU , PSEC  , PTAU,&
 & PUD   , PWV   , PQS,&
 & PFDOWN, PFUP  , PCDOWN, PCUP , PSUDU2, PDIFFS, PDIRFS, &
& LDDUST,PPIZA_DST,PCGA_DST,PTAUREL_DST )

!**** *SWNI* - SHORTWAVE RADIATION, NEAR-INFRARED SPECTRAL INTERVALS

!     PURPOSE.
!     --------

!          COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE NEAR-INFRARED 
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SWNI* IS CALLED FROM *SW*.

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
!          4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
!          5. MULTIPLY BY OZONE TRANSMISSION FUNCTION

!     EXTERNALS.
!     ----------

!          *SWCLR*, *SWR*, *SWDE*, *SWTT*

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
!      03-03-17   JJMorcrette        Sunshine duration (correction)
!      03-10-10 Deborah Salmond and Marta Janiskova Optimisation
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      04-11-18   Y.Seity : add 2 arguments for AROME extern. surface
!      Y.Seity  05-10-10 : add add 3 optional arg. for dust SW properties
!      Y.Seity 06-09-09 : add modset from O.Thouron (MesoNH) under NOVLP tests
!      M.Janiskova 08-03-26 : code optimization and cleaning
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_RADIATION_MOD , ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE YOECLD   , ONLY : TECLD
USE YOPHNC   , ONLY : TPHNC
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOERDU   , ONLY : REPLOG, REPSCQ, REPSC
USE YOMCT0   , ONLY : LECMWF

IMPLICIT NONE

TYPE(MODEL_PHYSICS_RADIATION_TYPE),INTENT(IN):: YDML_PHY_RAD
TYPE(TPHNC)       ,INTENT(IN)    :: YDPHNC
TYPE(TECLD)       ,INTENT(IN)    :: YDECLD
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAER 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAKI(KLON,2,YDML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBD(KLON,YDML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,YDML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCG(KLON,YDML_PHY_RAD%YRERAD%NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLEAR(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSIG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMEGA(KLON,YDML_PHY_RAD%YRERAD%NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEC(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAU(KLON,YDML_PHY_RAD%YRERAD%NSW,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUD(KLON,5,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(KLON,KLEV) 
LOGICAL           ,INTENT(IN)    :: LDDUST
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUREL_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFDOWN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFUP(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCDOWN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCUP(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSUDU2(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFFS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIRFS(KLON) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IIND2(2), IIND3(6)
REAL(KIND=JPRB) :: ZCGAZ(KLON,KLEV) &
 & ,  ZRTMP(KLON)       , ZCLRTMP(KLON) &
 & ,  ZFD(KLON,KLEV+1)  , ZFU(KLON,KLEV+1) &
 & ,  ZG(KLON)          , ZGG(KLON)  
REAL(KIND=JPRB) :: ZPIZAZ(KLON,KLEV)&
 & ,  ZRAYL(KLON)       , ZRAY1(KLON,KLEV+1)  , ZRAY2(KLON,KLEV+1)&
 & ,  ZREFZ(KLON,2,KLEV+1)&
 & ,  ZRE1(KLON)        , ZRE2(KLON)&
 & ,  ZRJ(KLON,6,KLEV+1), ZRJ0(KLON,6,KLEV+1)&
 & ,  ZRK(KLON,6,KLEV+1), ZRK0(KLON,6,KLEV+1)&
 & ,  ZRL(KLON,8)&
 & ,  ZRMUE(KLON,KLEV+1), ZRMU0(KLON,KLEV+1) &
 & ,  ZRUEF(KLON,8)     , ZR1(KLON) &
 & ,  ZR2(KLON,2)       , ZR3(KLON,6)         , ZR4(KLON,2)&
 & ,  ZR21(KLON)        , ZR22(KLON)  
REAL(KIND=JPRB) :: ZS(KLON)&
 & ,  ZTAUAZ(KLON,KLEV) , ZTA1(KLON)          , ZTO1(KLON)&
 & ,  ZTR(KLON,2,KLEV+1), ZTRA1(KLON,KLEV+1)  , ZTRA2(KLON,KLEV+1)&
 & ,  ZTRCLD(KLON)      , ZTRCLR(KLON)&
 & ,  ZTR1(KLON)        , ZTR2(KLON)&
 & ,  ZW(KLON)          , ZW1(KLON)           , ZW2(KLON,2)&
 & ,  ZW3(KLON,6)       , ZW4(KLON,2)         , ZW5(KLON,2)  
REAL(KIND=JPRB) :: ZRR2(KLON,2)

INTEGER(KIND=JPIM) :: IABSL, IKL, IKM1, JABS, JAJ, JAJP, JK, &
 & JKKI, JKKP4, JKL, JKLP1, JKM1, JL, JN, JN2J, JREF  

REAL(KIND=JPRB) :: ZAA, ZBB, ZCNEB, ZRE11, ZRMUM1, ZWH2O, ZCHKG, ZCHKS
REAL(KIND=JPRB) :: ZRR,ZRRJ,ZRRK
REAL(KIND=JPRB) :: ZCLDIR

! fields for adjoint computation
REAL(KIND=JPRB) :: ZRESWR1(KLON,KLEV), ZRESWR2(KLON,KLEV)
REAL(KIND=JPRB) :: ZTRSWR1(KLON,KLEV), ZTRSWR2(KLON,KLEV)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZB_ODI(KLON)

#include "swclr.intfb.h"
#include "swde.intfb.h"
#include "swr.intfb.h"
#include "swtt.intfb.h"
#include "swtt1.intfb.h"

!     ------------------------------------------------------------------

!*         1.     NEAR-INFRARED SPECTRAL INTERVAL (0.68-4.00 MICRON)
!                 --------------------------------------------------

!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------

IF (LHOOK) CALL DR_HOOK('SWNI',0,ZHOOK_HANDLE)
ASSOCIATE(NOVLP=>YDML_PHY_RAD%YRERAD%NOVLP, NSW=>YDML_PHY_RAD%YRERAD%NSW, YDESWRT=>YDML_PHY_RAD%YRESWRT)
DO JL = KIDIA,KFDIA
  ZRMUM1 = 1.0_JPRB - PRMU(JL)
  ZRAYL(JL) =  YDESWRT%RRAY(KNU,1) + ZRMUM1   * (YDESWRT%RRAY(KNU,2) + ZRMUM1 &
   & * (YDESWRT%RRAY(KNU,3) + ZRMUM1   * (YDESWRT%RRAY(KNU,4) + ZRMUM1 &
   & * (YDESWRT%RRAY(KNU,5) + ZRMUM1   *  YDESWRT%RRAY(KNU,6)     ))))  
  ZRAYL(JL) = MAX (ZRAYL(JL), 0.0_JPRB) 
ENDDO

!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------

!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------

   CALL SWCLR &
        &( YDML_PHY_RAD, &
        &  KIDIA , KFDIA , KLON ,  KLEV , KAER , KNU , &
        &  PAER  , PALBP , PDSIG , ZRAYL, PSEC , &
        &  ZCGAZ , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 , &
        &  ZRK0  , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2, ZTRCLR , &
        &  LDDUST,PPIZA_DST,PCGA_DST,PTAUREL_DST &
        &)

!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------

CALL SWR &
 & ( YDML_PHY_RAD%YREOVLP,YDML_PHY_RAD%YRERAD,YDPHNC,YDECLD, &
 & KIDIA , KFDIA , KLON , KLEV  , KNU,&
 & PALBD , PCG   , PCLD , POMEGA, PSEC , PTAU,&
 & ZCGAZ , ZPIZAZ, ZRAY1, ZRAY2 , ZREFZ, ZRJ  , ZRK, ZRMUE,&
 & ZTAUAZ, ZTRA1 , ZTRA2, ZTRCLD,&
! for adjoint computation
 & ZRESWR1,ZRESWR2,ZTRSWR1,ZTRSWR2 &
 & )  

!     ------------------------------------------------------------------

!*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
!                ------------------------------------------------------

JN = 2

DO JABS=1,2

!*         3.1  SURFACE CONDITIONS
!               ------------------

  DO JL = KIDIA,KFDIA
    ZREFZ(JL,2,1) = PALBD(JL,KNU)
    ZREFZ(JL,1,1) = PALBD(JL,KNU)
  ENDDO

!*         3.2  INTRODUCING CLOUD EFFECTS
!               -------------------------

  DO JK = 2 , KLEV+1
    JKM1 = JK - 1
    IKL=KLEV+1-JKM1
    DO JL = KIDIA,KFDIA
      IF (JABS == 1.AND. PCLD(JL,JKM1) > REPSC ) THEN
        ZWH2O=MAX(PWV(JL,IKL),REPSCQ)
        ZCNEB=MAX(REPSC ,MIN(PCLD(JL,JKM1),1.0_JPRB-REPSC ))
        ZBB=PUD(JL,JABS,JKM1)*PQS(JL,IKL)/ZWH2O
        ZAA=MAX((PUD(JL,JABS,JKM1)-ZCNEB*ZBB)/(1.0_JPRB-ZCNEB),REPSCQ)
      ELSE
        ZAA=PUD(JL,JABS,JKM1)
        ZBB=ZAA
      ENDIF
            
!      ZS(JL) = EXP(-PAKI(JL,JABS,KNU) * ZAA * 1.66_JPRB)
!      ZG(JL) = EXP(-PAKI(JL,JABS,KNU) * ZAA / ZRMUE(JL,JK) )
      
      ZCHKS = MIN( 200._JPRB, PAKI(JL,JABS,KNU) * ZAA * 1.66_JPRB )
      ZCHKG = MIN( 200._JPRB, PAKI(JL,JABS,KNU) * ZAA / ZRMUE(JL,JK))
      ZS(JL) = EXP( - ZCHKS )
      ZG(JL) = EXP( - ZCHKG )
      
      IF (NOVLP >= 5)THEN !MESONH VERSION
        ZW(JL) =PCG(JL,KNU,JKM1)*PCG(JL,KNU,JKM1)
        ZTO1(JL) = PTAU(JL,KNU,JKM1)*(1-(POMEGA(JL,KNU,JKM1)*ZW(JL)))
        ZW(JL) =POMEGA(JL,KNU,JKM1)*(1-ZW(JL))/(1-(POMEGA(JL,KNU,JKM1)*ZW(JL)))
        ZGG(JL) =PCG(JL,KNU,JKM1)/(1+PCG(JL,KNU,JKM1))
        ZGG(JL)=ZTO1(JL)*ZW(JL)*ZGG(JL) &
         & +ZTAUAZ(JL,JKM1)*ZPIZAZ(JL,JKM1)*ZCGAZ(JL,JKM1)
        ZGG(JL)=ZGG(JL)/(ZTO1(JL)*ZW(JL)+ZTAUAZ(JL,JKM1)*ZPIZAZ(JL,JKM1))
        ZB_ODI(JL)=ZTO1(JL) / ZW(JL)&
         &+ ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)&
     !if g=0 tau/w=tau'/w'
         &+ ZBB * PAKI(JL,JABS,KNU)
        ZB_ODI(JL)=(1/( (ZTO1(JL) / ZW(JL))&
         &+ (ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)) ))-(1/ZB_ODI(JL))
        ZB_ODI(JL)=((ZTO1(JL) +  ZTAUAZ(JL,JKM1))**2)*ZB_ODI(JL)
        ZW(JL)=ZTO1(JL)*ZW(JL)+ZTAUAZ(JL,JKM1)*ZPIZAZ(JL,JKM1)-ZB_ODI(JL)
        ZTO1(JL) = ZTO1(JL) +  ZTAUAZ(JL,JKM1)
        ZW(JL)=ZW(JL)/ZTO1(JL)
      ELSE !ECMWF VERSION
        ZW(JL)= POMEGA(JL,KNU,JKM1)
        ZTO1(JL) = PTAU(JL,KNU,JKM1) / ZW(JL)&
         & + ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)&
         & + ZBB * PAKI(JL,JABS,KNU)  
        ZR21(JL) = PTAU(JL,KNU,JKM1) + ZTAUAZ(JL,JKM1)
        ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
        ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)&
         & + (1.0_JPRB - ZR22(JL)) * ZCGAZ(JL,JKM1)  
        ZW(JL) = ZR21(JL) / ZTO1(JL)
      ENDIF
    ENDDO

    CALL SWDE ( YDML_PHY_RAD%YRERAD, KIDIA, KFDIA, KLON,&
     & ZGG  , ZREFZ(1,1,JKM1), ZRMUE(1,JK), ZTO1, ZW,&
     & ZRE1 , ZRE2 , ZTR1 , ZTR2   )  

     DO JL = KIDIA,KFDIA

      ZRR=1.0_JPRB/(1.0_JPRB-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1))
      ZREFZ(JL,2,JK) = (1.0_JPRB-PCLD(JL,JKM1)) * (ZRAY1(JL,JKM1)&
       & + ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)&
       & * ZTRA2(JL,JKM1) ) * ZG(JL) * ZS(JL)&
       & + PCLD(JL,JKM1) * ZRE1(JL)  

      ZTR(JL,2,JKM1)=PCLD(JL,JKM1)*ZTR1(JL)&
       & + (ZTRA1(JL,JKM1)) * ZG(JL) * (1.0_JPRB-PCLD(JL,JKM1))  

      ZREFZ(JL,1,JK)=(1.0_JPRB-PCLD(JL,JKM1))*(ZRAY1(JL,JKM1)&
       & +ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)&
       & *ZRR ) &
       & *ZG(JL)*ZS(JL)&
       & + PCLD(JL,JKM1) * ZRE2(JL)  

      ZTR(JL,1,JKM1)= PCLD(JL,JKM1) * ZTR2(JL)&
       & + (ZTRA1(JL,JKM1) &
       & *ZRR ) &
       & * ZG(JL) * (1.0_JPRB -PCLD(JL,JKM1))  

    ENDDO
  ENDDO

!*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!               -------------------------------------------------

  DO JREF=1,2

    JN = JN + 1

    DO JL = KIDIA,KFDIA
      ZRJ(JL,JN,KLEV+1) = 1.0_JPRB
      ZRK(JL,JN,KLEV+1) = ZREFZ(JL,JREF,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      JKL = KLEV+1 - JK
      JKLP1 = JKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11 = ZRJ(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
        ZRJ(JL,JN,JKL) = ZRE11
        ZRK(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
      ENDDO
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
      ZRK(JL,JAJ,JK)=        ZRK(JL,JAJ,JK) - ZRK(JL,JAJP,JK)
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , REPLOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , REPLOG )
    ENDDO
  ENDDO
ENDDO

DO JK = 1 , KLEV+1
  DO JAJ = 2 , 6 , 2
    DO JL = KIDIA,KFDIA
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , REPLOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , REPLOG )
    ENDDO
  ENDDO
ENDDO

!*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
!                 ---------------------------------------------

DO JAJ = 1 , 2
  DO JL = KIDIA,KFDIA
    ZRR2(JL,JAJ)=1.0_JPRB/PAKI(JL,JAJ,KNU)
  ENDDO
ENDDO

DO JK = 1 , KLEV+1
  JKKI = 1
  DO JAJ = 1 , 2
    IIND2(1)=JAJ
    IIND2(2)=JAJ
    DO JN = 1 , 2
      JN2J = JN + 2 * JAJ
      JKKP4 = JKKI + 4

!*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
!                 --------------------------

      DO JL = KIDIA,KFDIA
!        ZRR=1.0_JPRB/PAKI(JL,JAJ,KNU)
        ZRRJ=ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK)
        ZRRK=ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK)
        IF (LECMWF) THEN
          ZW2(JL,1) = LOG( ZRRJ ) * ZRR2(JL,JAJ)
          ZW2(JL,2) = LOG( ZRRK ) * ZRR2(JL,JAJ)
        ELSE
          ZW2(JL,1) = LOG( MAX(1.0_JPRB,ZRRJ) ) * ZRR2(JL,JAJ)
          ZW2(JL,2) = LOG( MAX(1.0_JPRB,ZRRK) ) * ZRR2(JL,JAJ)
        ENDIF
      ENDDO

!*         4.2.2  TRANSMISSION FUNCTION
!                 ---------------------

      CALL SWTT1 ( YDESWRT, KIDIA,KFDIA,KLON, KNU, 2, IIND2,&
       & ZW2,&
       & ZR2                              )  

      DO JL = KIDIA,KFDIA
        ZRL(JL,JKKI) = ZR2(JL,1)
        ZRUEF(JL,JKKI) = ZW2(JL,1)
        ZRL(JL,JKKP4) = ZR2(JL,2)
        ZRUEF(JL,JKKP4) = ZW2(JL,2)
      ENDDO

      JKKI=JKKI+1
    ENDDO
  ENDDO

!*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
!                 ------------------------------------------------------

  DO JL = KIDIA,KFDIA
    PFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)&
     & + ZRJ(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)  
    PFUP(JL,JK)   = ZRK(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)&
     & + ZRK(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)  
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
!                ----------------------------------------

!*         5.1   DOWNWARD FLUXES
!                ---------------

JAJ = 2
IIND3(1)=1
IIND3(2)=2
IIND3(3)=3
IIND3(4)=1
IIND3(5)=2
IIND3(6)=3

DO JL = KIDIA,KFDIA
  ZW3(JL,1)=0.0_JPRB
  ZW3(JL,2)=0.0_JPRB
  ZW3(JL,3)=0.0_JPRB
  ZW3(JL,4)=0.0_JPRB
  ZW3(JL,5)=0.0_JPRB
  ZW3(JL,6)=0.0_JPRB

  ZW4(JL,1)=0.0_JPRB
  ZW5(JL,1)=0.0_JPRB
  ZR4(JL,1)=1.0_JPRB
  ZW4(JL,2)=0.0_JPRB
  ZW5(JL,2)=0.0_JPRB
  ZR4(JL,2)=1.0_JPRB
  ZFD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1)
ENDDO
DO JK = 1 , KLEV
  IKL = KLEV+1-JK
  DO JL = KIDIA,KFDIA
    ZRR=1.0_JPRB/ZRMU0(JL,IKL)
    ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKL)*ZRR
    ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKL)*ZRR
    ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKL)*ZRR
    ZW4(JL,1)=ZW4(JL,1)+PUD(JL,4,IKL)*ZRR
    ZW5(JL,1)=ZW5(JL,1)+PUD(JL,5,IKL)*ZRR

    ZRR=1.0_JPRB/ZRMUE(JL,IKL)
    ZW3(JL,4)=ZW3(JL,4)+PUD(JL,1,IKL)*ZRR
    ZW3(JL,5)=ZW3(JL,5)+PUD(JL,2,IKL)*ZRR
    ZW3(JL,6)=ZW3(JL,6)+POZ(JL,  IKL)*ZRR
    ZW4(JL,2)=ZW4(JL,2)+PUD(JL,4,IKL)*ZRR
    ZW5(JL,2)=ZW5(JL,2)+PUD(JL,5,IKL)*ZRR
  ENDDO

  CALL SWTT1 ( YDESWRT, KIDIA,KFDIA,KLON, KNU, 6, IIND3,&
   & ZW3,&
   & ZR3                              )  

  DO JL = KIDIA,KFDIA
    ZR4(JL,1) = EXP(-YDESWRT%RSWCE(KNU)*ZW4(JL,1)-YDESWRT%RSWCP(KNU)*ZW5(JL,1))
    ZR4(JL,2) = EXP(-YDESWRT%RSWCE(KNU)*ZW4(JL,2)-YDESWRT%RSWCP(KNU)*ZW5(JL,2))
    ZFD(JL,IKL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL,1)* ZRJ0(JL,JAJ,IKL)
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  ZRTMP(JL) = ZR3(JL,4)*ZR3(JL,5)*ZR3(JL,6)*ZR4(JL,2)*ZTRCLD(JL)
  ZCLRTMP(JL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL,1)*ZTRCLR(JL)
  PSUDU2(JL) = ((1.0_JPRB-PCLEAR(JL)) * ZRTMP(JL) &
   & +PCLEAR(JL) * ZCLRTMP(JL)) * YDESWRT%RSUN(KNU)  
ENDDO

!*         5.2   UPWARD FLUXES
!                -------------

DO JL = KIDIA,KFDIA
  ZFU(JL,1) = ZFD(JL,1)*PALBP(JL,KNU)
ENDDO

DO JK = 2 , KLEV+1
  IKM1=JK-1
  DO JL = KIDIA,KFDIA
    ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
    ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
    ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKM1)*1.66_JPRB
    ZW4(JL,1)=ZW4(JL,1)+PUD(JL,4,IKM1)*1.66_JPRB
    ZW5(JL,1)=ZW5(JL,1)+PUD(JL,5,IKM1)*1.66_JPRB
  ENDDO

  CALL SWTT1 ( YDESWRT, KIDIA,KFDIA,KLON, KNU, 3, IIND3,&
   & ZW3,&
   & ZR3                              )  

  DO JL = KIDIA,KFDIA
    ZR4(JL,1) = EXP(-YDESWRT%RSWCE(KNU)*ZW4(JL,1)-YDESWRT%RSWCP(KNU)*ZW5(JL,1))
    ZFU(JL,JK) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL,1)* ZRK0(JL,JAJ,JK)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
!                 --------------------------------------------------

IABSL=3

!*         6.1    DOWNWARD FLUXES
!                 ---------------

DO JL = KIDIA,KFDIA
  ZW1(JL)=0.0_JPRB
  ZW4(JL,1)=0.0_JPRB
  ZW5(JL,1)=0.0_JPRB
  ZR1(JL)=0.0_JPRB
  PFDOWN(JL,KLEV+1) = ((1.0_JPRB-PCLEAR(JL))*PFDOWN(JL,KLEV+1)&
   & + PCLEAR(JL) * ZFD(JL,KLEV+1)) * YDESWRT%RSUN(KNU)  
  PCDOWN(JL,KLEV+1) = ZFD(JL,KLEV+1) * YDESWRT%RSUN(KNU)
ENDDO

DO JL = KIDIA,KFDIA
  ZTA1(JL)=0.0_JPRB
  ZTO1(JL)=0.0_JPRB
ENDDO
DO JK = 1 , KLEV
  IKL  = KLEV+1-JK
  DO JL = KIDIA,KFDIA
    ZRR=1.0_JPRB/ZRMUE(JL,IKL)
    ZW1(JL) = ZW1(JL)+POZ(JL,  IKL) * ZRR
    ZW4(JL,1) = ZW4(JL,1)+PUD(JL,4,IKL) * ZRR
    ZW5(JL,1) = ZW5(JL,1)+PUD(JL,5,IKL) * ZRR
    ZR4(JL,1) = EXP(-YDESWRT%RSWCE(KNU)*ZW4(JL,1)-YDESWRT%RSWCP(KNU)*ZW5(JL,1))
  ENDDO

  CALL SWTT ( YDESWRT, KIDIA, KFDIA, KLON, KNU, IABSL, ZW1, ZR1 )

  DO JL = KIDIA,KFDIA
    PFDOWN(JL,IKL) = ((1.0_JPRB-PCLEAR(JL))*ZR1(JL)*ZR4(JL,1)*PFDOWN(JL,IKL)&
     & +PCLEAR(JL)*ZFD(JL,IKL)) * YDESWRT%RSUN(KNU)  
    PCDOWN(JL,IKL) = ZFD(JL,IKL) * YDESWRT%RSUN(KNU)
    ZTA1(JL)=ZTA1(JL)+ZTAUAZ(JL,IKL)
    ZTO1(JL) = PTAU(JL,KNU,IKL)*(1.-(POMEGA(JL,KNU,IKL)* &
   &           PCG(JL,KNU,IKL)*PCG(JL,KNU,IKL))) + ZTO1(JL)
  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZCLDIR     = ZFD(JL,IKL)/ZRJ0(JL,JAJ,IKL)*EXP(-ZTA1(JL)/PRMU(JL))
  PDIRFS(JL) = ((1.0_JPRB-PCLEAR(JL))*ZCLDIR*EXP(-ZTO1(JL)/PRMU(JL)) + &
 &              PCLEAR(JL)*ZCLDIR) * YDESWRT%RSUN(KNU)
  PDIRFS(JL) = MIN(PFDOWN(JL,1),PDIRFS(JL))
  PDIFFS(JL) = PFDOWN(JL,1) - PDIRFS(JL)
ENDDO

!*         6.2    UPWARD FLUXES
!                 -------------

DO JL = KIDIA,KFDIA
  PFUP(JL,1) = ((1.0_JPRB-PCLEAR(JL))*ZR1(JL)*ZR4(JL,1) * PFUP(JL,1)&
   & +PCLEAR(JL)*ZFU(JL,1)) * YDESWRT%RSUN(KNU)  
  PCUP(JL,1) = ZFU(JL,1) * YDESWRT%RSUN(KNU)
ENDDO

DO JK = 2 , KLEV+1
  IKM1=JK-1
  DO JL = KIDIA,KFDIA
    ZW1(JL) = ZW1(JL)+POZ(JL  ,IKM1)*1.66_JPRB
    ZW4(JL,1) = ZW4(JL,1)+PUD(JL,4,IKM1)*1.66_JPRB
    ZW5(JL,1) = ZW5(JL,1)+PUD(JL,5,IKM1)*1.66_JPRB
    ZR4(JL,1) = EXP(-YDESWRT%RSWCE(KNU)*ZW4(JL,1)-YDESWRT%RSWCP(KNU)*ZW5(JL,1))
  ENDDO

  CALL SWTT ( YDESWRT, KIDIA, KFDIA, KLON, KNU, IABSL, ZW1, ZR1 )

  DO JL = KIDIA,KFDIA
    PFUP(JL,JK) = ((1.0_JPRB-PCLEAR(JL))*ZR1(JL)*ZR4(JL,1) * PFUP(JL,JK)&
     & +PCLEAR(JL)*ZFU(JL,JK)) * YDESWRT%RSUN(KNU)  
    PCUP(JL,JK) = ZFU(JL,JK) * YDESWRT%RSUN(KNU)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SWNI',1,ZHOOK_HANDLE)
END SUBROUTINE SWNI

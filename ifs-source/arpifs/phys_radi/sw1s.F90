! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SW1S &
 & ( YDML_PHY_RAD,YDPHNC,YDECLD,KIDIA , KFDIA , KLON , KLEV , KAER , KNU,&
 & PAER  , PALBD , PALBP, PCG  , PCLD , PCLEAR,&
 & PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD,&
 & PFD   , PFU   , PCD  , PCU  , PSUDU1,PDIFFS, PDIRFS, &
!++MODIFCODE
 & LDDUST,PPIZA_DST,PCGA_DST,PTAUREL_DST  &
!--MODIFCODE
 &)

!**** *SW1S* - SHORTWAVE RADIATION, FIRST SPECTRAL INTERVAL

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SW1S* IS CALLED FROM *SW*.

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES QUANTITIES FOR THE CLEAR-SKY FRACTION OF THE
!     COLUMN
!          2. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
!     CONTINUUM SCATTERING
!          3. MULTIPLY BY OZONE TRANSMISSION FUNCTION

!     EXTERNALS.
!     ----------

!          *SWCLR*, *SWR*, *SWTT*, *SWUVO3*

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
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Y.Seity  04-11-19 : add two arguments for AROME externalized surface
!      Y.Seity  05-10-10 : add 3 optional arg. for dust SW properties
!      Y.Seity 06-09-09 : add modset from O.Thouron (MesoNH) under NOVLP tests
!      M.Janiskova 08-03-18 modifications for optimized TL/AD code
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_RADIATION_MOD , ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE YOECLD   , ONLY : TECLD
USE YOPHNC   , ONLY : TPHNC
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(TECLD)       ,INTENT(IN)    :: YDECLD
TYPE(MODEL_PHYSICS_RADIATION_TYPE),INTENT(IN):: YDML_PHY_RAD
TYPE(TPHNC)       ,INTENT(IN)    :: YDPHNC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAER 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) 
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
!++MODIFCODE
LOGICAL           ,INTENT(IN)    :: LDDUST          ! flag for DUST
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUREL_DST(KLON,KLEV)
!--MODIFCODE
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFD(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCD(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSUDU1(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFFS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIRFS(KLON) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IIND(6)

REAL(KIND=JPRB) :: ZCGAZ(KLON,KLEV)&
 & ,  ZRTMP(KLON)        , ZCLRTMP(KLON)        &
 & ,  ZDIFT(KLON)        , ZDIRT(KLON)        &
 & ,  ZPIZAZ(KLON,KLEV)&
 & ,  ZRAYL(KLON), ZRAY1(KLON,KLEV+1), ZRAY2(KLON,KLEV+1)&
 & ,  ZREFZ(KLON,2,KLEV+1)&
 & ,  ZRJ(KLON,6,KLEV+1), ZRJ0(KLON,6,KLEV+1)&
 & ,  ZRK(KLON,6,KLEV+1), ZRK0(KLON,6,KLEV+1)&
 & ,  ZRMUE(KLON,KLEV+1), ZRMU0(KLON,KLEV+1)&
 & ,  ZR(KLON,6)&
 & ,  ZTAUAZ(KLON,KLEV)&
 & ,  ZTRA1(KLON,KLEV+1), ZTRA2(KLON,KLEV+1)&
 & ,  ZTRCLD(KLON)      , ZTRCLR(KLON)&
 & ,  ZW(KLON,6)        , ZO(KLON,2) ,ZT(KLON,2)   
REAL(KIND=JPRB) :: ZRE, ZR0
REAL(KIND=JPRB) :: ZTA1(KLON), ZTO1(KLON)
REAL(KIND=JPRB) :: ZCLDIR

! fields for adjoint computation
REAL(KIND=JPRB) :: ZRESWR1(KLON,KLEV), ZRESWR2(KLON,KLEV)
REAL(KIND=JPRB) :: ZTRSWR1(KLON,KLEV), ZTRSWR2(KLON,KLEV)

INTEGER(KIND=JPIM) :: IKL, IKM1, JAJ, JK, JL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "swclr.intfb.h"
#include "swr.intfb.h"
#include "swtt1.intfb.h"
#include "swuvo3.intfb.h"

!     ------------------------------------------------------------------

!*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
!                 ----------------------- ------------------

!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------

IF (LHOOK) CALL DR_HOOK('SW1S',0,ZHOOK_HANDLE)
ASSOCIATE(NSW=>YDML_PHY_RAD%YRERAD%NSW, &
 & RRAY=>YDML_PHY_RAD%YRESWRT%RRAY, RSUN=>YDML_PHY_RAD%YRESWRT%RSUN)
DO JL = KIDIA,KFDIA
  ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)&
   & * (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)&
   & * (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))  
ENDDO
!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------

!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------

!++MODIFCODE
CALL SWCLR(YDML_PHY_RAD,KIDIA,KFDIA,KLON,KLEV,KAER,KNU,PAER,PALBP,PDSIG, &
 &         ZRAYL,PSEC,ZCGAZ,ZPIZAZ,ZRAY1,ZRAY2,ZREFZ,ZRJ0,ZRK0,ZRMU0,ZTAUAZ,ZTRA1,ZTRA2,ZTRCLR,LDDUST,&
 &         PPIZA_DST,PCGA_DST,PTAUREL_DST)
!--MODIFCODE

!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------

CALL SWR(YDML_PHY_RAD%YREOVLP,YDML_PHY_RAD%YRERAD,YDPHNC,YDECLD,KIDIA,KFDIA,KLON,KLEV,KNU,PALBD,PCG,PCLD, &
 &       POMEGA,PSEC,PTAU,ZCGAZ,ZPIZAZ,ZRAY1,ZRAY2,ZREFZ,ZRJ,ZRK,ZRMUE,ZTAUAZ,ZTRA1,ZTRA2,ZTRCLD,ZRESWR1, &
 &       ZRESWR2,ZTRSWR1,ZTRSWR2)

!     ------------------------------------------------------------------

!*         3.    OZONE ABSORPTION
!                ----------------

IF (NSW <= 4) THEN

!*         3.1   TWO OR FOUR SPECTRAL INTERVALS
!                ------------------------------

  IIND(1)=1
  IIND(2)=2
  IIND(3)=3
  IIND(4)=1
  IIND(5)=2
  IIND(6)=3

!*         3.1.1  DOWNWARD FLUXES
!                 ---------------

  JAJ = 2

  DO JL = KIDIA,KFDIA
    ZW(JL,1)=0.0_JPRB
    ZW(JL,2)=0.0_JPRB
    ZW(JL,3)=0.0_JPRB
    ZW(JL,4)=0.0_JPRB
    ZW(JL,5)=0.0_JPRB
    ZW(JL,6)=0.0_JPRB
    PFD(JL,KLEV+1)=((1.0_JPRB-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
     & + PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)  
    PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
  ENDDO
  DO JL = KIDIA,KFDIA
    ZTA1(JL)=0.0_JPRB
    ZTO1(JL)=0.0_JPRB
  ENDDO
  DO JK = 1 , KLEV
    IKL = KLEV+1-JK
    DO JL = KIDIA,KFDIA
      ZRE = 1.0_JPRB/ZRMUE(JL,IKL)
      ZR0 = 1.0_JPRB/ZRMU0(JL,IKL)
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)*ZRE
      ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKL)*ZRE
      ZW(JL,3)=ZW(JL,3)+POZ(JL,  IKL)*ZRE
      ZW(JL,4)=ZW(JL,4)+PUD(JL,1,IKL)*ZR0
      ZW(JL,5)=ZW(JL,5)+PUD(JL,2,IKL)*ZR0
      ZW(JL,6)=ZW(JL,6)+POZ(JL,  IKL)*ZR0
    ENDDO
    
    CALL SWTT1 ( YDML_PHY_RAD%YRESWRT, KIDIA, KFDIA, KLON, KNU, 6,&
     & IIND,&
     & ZW,&
     & ZR                          )  

    DO JL = KIDIA,KFDIA
      ZRTMP(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZRJ(JL,JAJ,IKL)
      ZCLRTMP(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZRJ0(JL,JAJ,IKL)
      PFD(JL,IKL) = ((1.0_JPRB-PCLEAR(JL)) * ZRTMP(JL)&
       & +PCLEAR(JL)  * ZCLRTMP(JL)) * RSUN(KNU)  
      PCD(JL,IKL) = ZCLRTMP(JL) * RSUN(KNU)
      ZTA1(JL) = ZTA1(JL) + ZTAUAZ(JL,IKL)
      ZTO1(JL) = PTAU(JL,KNU,IKL)*(1.-(POMEGA(JL,KNU,IKL)*&
     &           PCG(JL,KNU,IKL)*PCG(JL,KNU,IKL))) + ZTO1(JL)
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    ZCLDIR = ZCLRTMP(JL)/ZRJ0(JL,JAJ,1)*EXP(-ZTA1(JL)/PRMU(JL))
    PDIRFS(JL) = ((1.0_JPRB-PCLEAR(JL))*ZCLDIR*EXP(-ZTO1(JL)/PRMU(JL))+&
   &              PCLEAR(JL)*ZCLDIR) * RSUN(KNU)
    PDIRFS(JL) = MIN(PFD(JL,1),PDIRFS(JL))
    PDIFFS(JL) = PFD(JL,1) - PDIRFS(JL)
    ZDIFT(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZTRCLD(JL)
    ZDIRT(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZTRCLR(JL)
    PSUDU1(JL) = ((1.0_JPRB-PCLEAR(JL)) * ZDIFT(JL)&
     & +PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)  
  ENDDO

!*         3.1.2  UPWARD FLUXES
!                 -------------

  DO JL = KIDIA,KFDIA
    PFU(JL,1) = ((1.0_JPRB-PCLEAR(JL))*ZRTMP(JL)*PALBD(JL,KNU)&
     & + PCLEAR(JL) *ZCLRTMP(JL)*PALBP(JL,KNU))&
     & * RSUN(KNU)  
    PCU(JL,1) = ZCLRTMP(JL) * PALBP(JL,KNU) * RSUN(KNU)
  ENDDO

  DO JK = 2 , KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW(JL,3)=ZW(JL,3)+POZ(JL,  IKM1)*1.66_JPRB
      ZW(JL,4)=ZW(JL,4)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW(JL,5)=ZW(JL,5)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW(JL,6)=ZW(JL,6)+POZ(JL,  IKM1)*1.66_JPRB
    ENDDO
    
    CALL SWTT1 ( YDML_PHY_RAD%YRESWRT, KIDIA, KFDIA, KLON, KNU, 6,&
     & IIND,&
     & ZW,&
     & ZR                          )  
  
    DO JL = KIDIA,KFDIA
      ZRTMP(JL) = ZR(JL,1)*ZR(JL,2)*ZR(JL,3)*ZRK(JL,JAJ,JK)
      ZCLRTMP(JL) = ZR(JL,4)*ZR(JL,5)*ZR(JL,6)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((1.0_JPRB-PCLEAR(JL)) * ZRTMP(JL)&
       & +PCLEAR(JL)  * ZCLRTMP(JL)) * RSUN(KNU)  
      PCU(JL,JK) = ZCLRTMP(JL) * RSUN(KNU)
    ENDDO
  ENDDO

ELSEIF (NSW == 6) THEN

!*         3.2   SIX SPECTRAL INTERVALS
!                ----------------------

  IIND(1)=1
  IIND(2)=2
  IIND(3)=1
  IIND(4)=2

!*         3.2,1  DOWNWARD FLUXES
!                 ---------------

  JAJ = 2

  DO JL = KIDIA,KFDIA
    ZW(JL,1)=0.0_JPRB
    ZW(JL,2)=0.0_JPRB
    ZW(JL,3)=0.0_JPRB
    ZW(JL,4)=0.0_JPRB
  
    ZO(JL,1)=0.0_JPRB
    ZO(JL,2)=0.0_JPRB
    PFD(JL,KLEV+1)=((1.0_JPRB-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
     & + PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)  
    PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
  ENDDO
  DO JL = KIDIA,KFDIA
    ZTA1(JL)=0.0_JPRB
    ZTO1(JL)=0.0_JPRB
  ENDDO
  DO JK = 1 , KLEV
    IKL = KLEV+1-JK
    DO JL = KIDIA,KFDIA
      ZRE = 1.0_JPRB/ZRMUE(JL,IKL)
      ZR0 = 1.0_JPRB/ZRMU0(JL,IKL)

      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)*ZRE
      ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKL)*ZRE
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKL)*ZR0
      ZW(JL,4)=ZW(JL,4)+PUD(JL,2,IKL)*ZR0
    
      ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKL)*ZRE
      ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKL)*ZR0
    ENDDO
 
    CALL SWTT1 ( YDML_PHY_RAD%YRESWRT, KIDIA, KFDIA, KLON, KNU, 4,&
     & IIND,&
     & ZW,&
     & ZR&
     & )  

    CALL SWUVO3 ( KIDIA, KFDIA, KLON, KNU, 2,&
     & YDML_PHY_RAD%YRESWRT%NEXPO3, YDML_PHY_RAD%YRESWRT%REXPO3,&
     & ZO,&
     & ZT&
     & )  

    DO JL = KIDIA,KFDIA
      ZRTMP(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZRJ(JL,JAJ,IKL)
      ZCLRTMP(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZRJ0(JL,JAJ,IKL)
      PFD(JL,IKL) = ((1.0_JPRB-PCLEAR(JL)) * ZRTMP(JL)&
       & +PCLEAR(JL)  * ZCLRTMP(JL)) * RSUN(KNU)  
      PCD(JL,IKL) = ZCLRTMP(JL) * RSUN(KNU)
      ZTA1(JL) = ZTA1(JL) + ZTAUAZ(JL,IKL)
      ZTO1(JL) = PTAU(JL,KNU,IKL)*(1.-(POMEGA(JL,KNU,IKL)*&
     &           PCG(JL,KNU,IKL)*PCG(JL,KNU,IKL))) + ZTO1(JL)
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    ZCLDIR = ZCLRTMP(JL)/ZRJ0(JL,JAJ,1)*EXP(-ZTA1(JL)/PRMU(JL))
    PDIRFS(JL) = ((1.0_JPRB-PCLEAR(JL))*ZCLDIR*EXP(-ZTO1(JL)/PRMU(JL))+&
   &              PCLEAR(JL)*ZCLDIR) * RSUN(KNU)
    PDIRFS(JL) = MIN(PFD(JL,1),PDIRFS(JL))
    PDIFFS(JL) = PFD(JL,1) - PDIRFS(JL)
    ZDIFT(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZTRCLD(JL)
    ZDIRT(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZTRCLR(JL)
    PSUDU1(JL) = ((1.0_JPRB-PCLEAR(JL)) * ZDIFT(JL)&
     & +PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)  
  ENDDO

!*         3.2.2  UPWARD FLUXES
!                 -------------

  DO JL = KIDIA,KFDIA
    PFU(JL,1) = ((1.0_JPRB-PCLEAR(JL))*ZRTMP(JL)*PALBD(JL,KNU)&
     & + PCLEAR(JL) *ZCLRTMP(JL)*PALBP(JL,KNU))&
     & * RSUN(KNU)  
    PCU(JL,1) = ZCLRTMP(JL) * PALBP(JL,KNU) * RSUN(KNU)
  ENDDO

  DO JK = 2 , KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW(JL,4)=ZW(JL,4)+PUD(JL,2,IKM1)*1.66_JPRB
      
      ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKM1)*1.66_JPRB
      ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKM1)*1.66_JPRB
    ENDDO

    CALL SWTT1 ( YDML_PHY_RAD%YRESWRT, KIDIA, KFDIA, KLON, KNU, 4,&
     & IIND,&
     & ZW,&
     & ZR&
     & )  

    CALL SWUVO3 ( KIDIA, KFDIA, KLON, KNU, 2,&
     & YDML_PHY_RAD%YRESWRT%NEXPO3, YDML_PHY_RAD%YRESWRT%REXPO3,&
     & ZO,&
     & ZT&
     & )  

    DO JL = KIDIA,KFDIA
      ZRTMP(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZRK(JL,JAJ,JK)
      ZCLRTMP(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((1.0_JPRB-PCLEAR(JL)) * ZRTMP(JL)&
       & +PCLEAR(JL)  * ZCLRTMP(JL)) * RSUN(KNU)  
      PCU(JL,JK) = ZCLRTMP(JL) * RSUN(KNU)
    ENDDO
  ENDDO
  
ENDIF  

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SW1S',1,ZHOOK_HANDLE)
END SUBROUTINE SW1S

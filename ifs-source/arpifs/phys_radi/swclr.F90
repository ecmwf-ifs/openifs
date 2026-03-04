! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SWCLR &
 & ( YDML_PHY_RAD,KIDIA , KFDIA , KLON  , KLEV  , KAER  , KNU,&
 & PAER  , PALBP , PDSIG , PRAYL , PSEC,&
 & PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ,&
 & PRK   , PRMU0 , PTAUAZ, PTRA1 , PTRA2 , PTRCLR, &
!++MODIFCODE
  & LDDUST,PPIZA_DST, PCGA_DST, PTAUREL_DST )
!--MODIFCODE

!**** *SWCLR* - CLEAR-SKY COLUMN COMPUTATIONS

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
!     CLEAR-SKY COLUMN

!**   INTERFACE.
!     ----------

!          *SWCLR* IS CALLED EITHER FROM *SW1S*
!                                OR FROM *SWNI*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!      JEAN-JACQUES MORCRETTE  *ECMWF*
!      ORIGINAL : 94-11-15

!     MODIFICATIONS.
!     --------------
!      03-10-10 Deborah Salmond and Marta Janiskova Optimisation
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      A.Grini (Meteo-France: 2005-11-10) 
!      Y.Seity 05-10-10 : add add 3 optional arg. for dust SW properties
!      Y.Seity 06-09-09 : add modset from O.Thouron (MesoNH) under NOVLP tests
!      M.Janiskova 08-03-18 : modifications for optimized TL/AD code + cleaning
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MODEL_PHYSICS_RADIATION_MOD, ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE YOERDU   , ONLY : REPSCT

IMPLICIT NONE

TYPE(MODEL_PHYSICS_RADIATION_TYPE), INTENT(IN)    :: YDML_PHY_RAD
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAER 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBP(KLON,YDML_PHY_RAD%YRERAD%NSW) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSIG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAYL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEC(KLON) 
!++MODIFCODE
LOGICAL           ,INTENT(IN)    :: LDDUST                   ! flag for DUST
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPIZA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCGA_DST(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUREL_DST(KLON,KLEV)
!--MODIFCODE
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCGAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPIZAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAY1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAY2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREFZ(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRJ(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRK(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRMU0(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAUAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRA1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRA2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTRCLR(KLON) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZC0I(KLON,KLEV+1)&
 & ,  ZCLEAR(KLON) &
 & ,  ZR21(KLON)&
 & ,  ZSS0(KLON) , ZSCAT(KLON)&
 & ,  ZTR(KLON,2,KLEV+1)  

INTEGER(KIND=JPIM) :: IKL, JA, JAE, JAJ, JK, JKL, JKLP1, JKM1, JL, INU1

REAL(KIND=JPRB) :: ZBMU0, ZBMU1, ZCORAE, ZDEN, ZDEN1, ZFACOA,&
 & ZFF, ZMUE, ZRATIO, ZRE11, &
 & ZTRAY, ZDENB
REAL(KIND=JPRB) :: ZIPIZAZ, ZIPTAUZ, ZIDENB, ZICLEAR, ZDIV1, ZDIV2
REAL(KIND=JPRB) :: ZRR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!++MODIFCODE
REAL(KIND=JPRB) ::ZFACOA_NEW(KLON,KLEV)
!--MODIFCODE

!     ------------------------------------------------------------------

!*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
!                --------------------------------------------

IF (LHOOK) CALL DR_HOOK('SWCLR',0,ZHOOK_HANDLE)
ASSOCIATE(NOVLP=>YDML_PHY_RAD%YRERAD%NOVLP, NSW=>YDML_PHY_RAD%YRERAD%NSW, &
 & REPCLC=>YDML_PHY_RAD%YRERDI%REPCLC, YDESWRT=>YDML_PHY_RAD%YRESWRT)
DO JK = 1 , KLEV+1
  DO JA = 1 , 6
    DO JL = KIDIA,KFDIA
      PRJ(JL,JA,JK) = 0.0_JPRB
      PRK(JL,JA,JK) = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO

! ------   NB: 'PAER' AEROSOLS ARE ENTERED FROM TOP TO BOTTOM

DO JK = 1 , KLEV
  IKL=KLEV+1-JK
  DO JL = KIDIA,KFDIA
    PCGAZ(JL,JK) = 0.0_JPRB
    PPIZAZ(JL,JK) =  0.0_JPRB
    PTAUAZ(JL,JK) = 0.0_JPRB
    ZFACOA_NEW(JL,JK) = 0.0_JPRB
  ENDDO

!++MODIFCODE  
  IF(NOVLP < 5)THEN !ECMWF VERSION
    DO JAE=1,6
      DO JL = KIDIA,KFDIA
        PTAUAZ(JL,JK)=PTAUAZ(JL,JK)+PAER(JL,JAE,IKL)*YDESWRT%RTAUA(KNU,JAE)
        PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL,JAE,IKL)&
         & * YDESWRT%RTAUA(KNU,JAE)*YDESWRT%RPIZA(KNU,JAE)  
        PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL,JAE,IKL)&
         & * YDESWRT%RTAUA(KNU,JAE)*YDESWRT%RPIZA(KNU,JAE)*YDESWRT%RCGA(KNU,JAE)  
      ENDDO
    ENDDO
  ELSE ! MESONH VERSION
    DO JAE=1,6
      DO JL = KIDIA,KFDIA
          !Special optical properties for dust
        IF (LDDUST.AND.(JAE==3)) THEN
          !Ponderation of aerosol optical properties:first step 
          !ti
          PTAUAZ(JL,JK)=PTAUAZ(JL,JK) + PAER(JL,JAE,IKL) * PTAUREL_DST(JL,IKL)
          !wi*ti
          PPIZAZ(JL,JK)=PPIZAZ(JL,JK) + PAER(JL,JAE,IKL)&
                   & *PTAUREL_DST(JL,IKL)*PPIZA_DST(JL,IKL)
          !wi*ti*gi
          PCGAZ(JL,JK) = PCGAZ(JL,JK) + PAER(JL,JAE,IKL)&
                &  *PTAUREL_DST(JL,IKL)*PPIZA_DST(JL,IKL)*PCGA_DST(JL,IKL)
          !wi*ti*(gi**2)
          ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)+PAER(JL, JAE, IKL)&
                & *PTAUREL_DST(JL,IKL) *PPIZA_DST(JL,IKL)*PCGA_DST(JL,IKL)*&
                & PCGA_DST(JL,IKL)
        ELSE
          !Ponderation of aerosol optical properties:first step 
          !ti
          PTAUAZ(JL,JK)=PTAUAZ(JL,JK)+PAER(JL, JAE, IKL)*YDESWRT%RTAUA(KNU,JAE)
          !wi*ti
          PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL, JAE, IKL)&
                &* YDESWRT%RTAUA(KNU,JAE)*YDESWRT%RPIZA(KNU,JAE)
          !wi*ti*gi
          PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL, JAE, IKL)&
                &* YDESWRT%RTAUA(KNU,JAE)*YDESWRT%RPIZA(KNU,JAE)*YDESWRT%RCGA(KNU,JAE)
          !wi*ti*(gi**2)
          ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)+PAER(JL, JAE, IKL)&
                &* YDESWRT%RTAUA(KNU,JAE)*YDESWRT%RPIZA(KNU,JAE)*YDESWRT%RCGA(KNU,JAE)*YDESWRT%RCGA(KNU,JAE)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
!--MODIFCODE  

!++MODIFCODE  
  IF (NOVLP < 5) THEN !ECMWF VERSION
    DO JL = KIDIA,KFDIA
      IF (KAER /= 0) THEN
        ZIPIZAZ = 1.0_JPRB/PPIZAZ(JL,JK)
        PCGAZ(JL,JK)=PCGAZ(JL,JK)*ZIPIZAZ
        ZIPTAUZ = 1.0_JPRB/PTAUAZ(JL,JK)
        PPIZAZ(JL,JK)=PPIZAZ(JL,JK)*ZIPTAUZ
        ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
!!!! wrong  ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))
!--     

        ZFF = PCGAZ(JL,JK) * PCGAZ(JL,JK)

!-- bug-fix: ZRATIO must be defined from the transformed value of 
!            optical thickness
        ZDENB = ZTRAY + PTAUAZ(JL,JK)*(1.0_JPRB-PPIZAZ(JL,JK)*ZFF)
        ZIDENB = 1.0_JPRB/ZDENB
        ZRATIO=ZTRAY*ZIDENB
 !--     
        PTAUAZ(JL,JK)=ZTRAY+PTAUAZ(JL,JK)*(1.0_JPRB-PPIZAZ(JL,JK)*ZFF)
        ZDIV1 = 1.0_JPRB/(1.0_JPRB + PCGAZ(JL,JK))
        PCGAZ(JL,JK) = PCGAZ(JL,JK) * (1.0_JPRB - ZRATIO) * ZDIV1
        ZDIV2 = 1.0_JPRB/(1.0_JPRB - PPIZAZ(JL,JK) * ZFF) 
        PPIZAZ(JL,JK) =ZRATIO+(1.0_JPRB-ZRATIO)*PPIZAZ(JL,JK)*(1.0_JPRB-ZFF)&
         & * ZDIV2  
      ELSE
        ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
        PTAUAZ(JL,JK) = ZTRAY
        PCGAZ(JL,JK) = 0.0_JPRB
        PPIZAZ(JL,JK) = 1.0_JPRB-REPSCT
      ENDIF
    ENDDO
  ELSE !MESONH VERSION
    DO JL = KIDIA,KFDIA
      IF (KAER /= 0) THEN
        ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
        ZRATIO =PPIZAZ(JL,JK)+ZTRAY
        !Ponderation G**2
        ZFACOA_NEW(JL,JK)= ZFACOA_NEW(JL,JK)/ZRATIO
        !Ponderation w
        PPIZAZ(JL,JK)=ZRATIO/(PTAUAZ(JL,JK)+ZTRAY)
        !Ponderation g
        PCGAZ(JL,JK)=PCGAZ(JL,JK)/ZRATIO
        !Ponderation+delta-modified parameters tau
        PTAUAZ(JL,JK)=(ZTRAY+PTAUAZ(JL,JK))*&
         &  (1.0_JPRB-PPIZAZ(JL,JK)*ZFACOA_NEW(JL,JK))
        !delta-modified parameters w
        PPIZAZ(JL,JK)=PPIZAZ(JL,JK)*(1.0_JPRB-ZFACOA_NEW(JL,JK))/&
         & (1.0_JPRB-ZFACOA_NEW(JL,JK)*PPIZAZ(JL,JK))     
        !delta-modified parameters g
        PCGAZ(JL,JK)=PCGAZ(JL,JK)/(1.0_JPRB+PCGAZ(JL,JK))
      
      ELSE
        ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
        ZFACOA_NEW(JL,JK)= 0.0_JPRB
        PTAUAZ(JL,JK) = ZTRAY
        PCGAZ(JL,JK) = 0.0_JPRB
        PPIZAZ(JL,JK) = 1.0_JPRB-REPSCT
      ENDIF
    ENDDO    
  ENDIF
!--MODIFCODE  
  
ENDDO

!     ------------------------------------------------------------------

!*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
!                ----------------------------------------------

DO JL = KIDIA,KFDIA
  ZC0I(JL,KLEV+1) = 0.0_JPRB
  ZCLEAR(JL) = 1.0_JPRB
  ZSCAT(JL) = 0.0_JPRB
ENDDO

JK = 1
JKL = KLEV+1 - JK
JKLP1 = JKL + 1
DO JL = KIDIA,KFDIA
!++MODIFCODE
  IF (NOVLP >= 5) THEN
    ZFACOA = PTAUAZ(JL,JK)
    ZCORAE = ZFACOA *  PSEC(JL)
  ELSE
    ZFACOA = 1.0_JPRB - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
    ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
  ENDIF
!--MODIFCODE
  ZR21(JL) = EXP(-ZCORAE   )
  ZSS0(JL) = 1.0_JPRB-ZR21(JL)

  IF (NOVLP == 1 .OR. NOVLP == 4) THEN
!* maximum-random
    ZICLEAR = 1.0_JPRB/(1.0_JPRB-MIN(ZSCAT(JL),1.0_JPRB-REPCLC)) 
    ZCLEAR(JL) = ZCLEAR(JL)&
     & *(1.0_JPRB-MAX(ZSS0(JL),ZSCAT(JL)))*ZICLEAR 
    ZC0I(JL,JKL) = 1.0_JPRB - ZCLEAR(JL)
    ZSCAT(JL) = ZSS0(JL)
  ELSEIF (NOVLP == 2) THEN
!* maximum
    ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
    ZC0I(JL,JKL) = ZSCAT(JL)
!++MODIFCODE
  ELSEIF ((NOVLP == 3).OR.(NOVLP  >=  5)) THEN
!--MODIFCODE
!* random
    ZCLEAR(JL)=ZCLEAR(JL)*(1.0_JPRB-ZSS0(JL))
    ZSCAT(JL) = 1.0_JPRB - ZCLEAR(JL)
    ZC0I(JL,JKL) = ZSCAT(JL)
  ENDIF
ENDDO

DO JK = 2 , KLEV
  JKL = KLEV+1 - JK
  JKLP1 = JKL + 1
  DO JL = KIDIA,KFDIA
!++MODIFCODE
    IF (NOVLP >= 5) THEN
      ZFACOA = PTAUAZ(JL,JK)
      ZCORAE = ZFACOA *  PSEC(JL)
    ELSE
      ZFACOA = 1.0_JPRB - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
    ENDIF
!--MODIFCODE
    ZR21(JL) = EXP(-ZCORAE   )
    ZSS0(JL) = 1.0_JPRB-ZR21(JL)

    IF (NOVLP == 1 .OR. NOVLP == 4) THEN
!* maximum-random
      ZICLEAR = 1.0_JPRB/(1.0_JPRB-MIN(ZSCAT(JL),1.0_JPRB-REPCLC)) 
      ZCLEAR(JL) = ZCLEAR(JL)&
       & *(1.0_JPRB-MAX(ZSS0(JL),ZSCAT(JL)))*ZICLEAR
      ZC0I(JL,JKL) = 1.0_JPRB - ZCLEAR(JL)
      ZSCAT(JL) = ZSS0(JL)
    ELSEIF (NOVLP == 2) THEN
!* maximum
      ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
      ZC0I(JL,JKL) = ZSCAT(JL)
!++MODIFCODE
    ELSEIF ((NOVLP == 3).OR.(NOVLP >= 5)) THEN
!--MODIFCODE
!* random
      ZCLEAR(JL)=ZCLEAR(JL)*(1.0_JPRB-ZSS0(JL))
      ZSCAT(JL) = 1.0_JPRB - ZCLEAR(JL)
      ZC0I(JL,JKL) = ZSCAT(JL)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
!                -----------------------------------------------

DO JL = KIDIA,KFDIA
  PRAY1(JL,KLEV+1) = 0.0_JPRB
  PRAY2(JL,KLEV+1) = 0.0_JPRB
  PREFZ(JL,2,1) = PALBP(JL,KNU)
  PREFZ(JL,1,1) = PALBP(JL,KNU)
  PTRA1(JL,KLEV+1) = 1.0_JPRB
  PTRA2(JL,KLEV+1) = 1.0_JPRB
ENDDO

DO JK = 2 , KLEV+1
  JKM1 = JK-1
  DO JL = KIDIA,KFDIA

!     ------------------------------------------------------------------

!*         3.1  EQUIVALENT ZENITH ANGLE
!               -----------------------

    ZMUE = (1.0_JPRB-ZC0I(JL,JK)) * PSEC(JL)+ ZC0I(JL,JK) * 1.66_JPRB
    PRMU0(JL,JK) = 1.0_JPRB/ZMUE

!     ------------------------------------------------------------------

!*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
!               ----------------------------------------------------

    ZBMU0 = 0.5_JPRB - 0.75_JPRB * PCGAZ(JL,JKM1) *PRMU0(JL,JK)
    ZDEN = 1.0_JPRB + (1.0_JPRB - PPIZAZ(JL,JKM1)&
     & + ZBMU0 * PPIZAZ(JL,JKM1)) * PTAUAZ(JL,JKM1) * ZMUE&
     & + (1-PPIZAZ(JL,JKM1)) * (1.0_JPRB - PPIZAZ(JL,JKM1)&
     & + 2.0_JPRB*ZBMU0*PPIZAZ(JL,JKM1))*PTAUAZ(JL,JKM1)*PTAUAZ(JL,JKM1)&
     & * ZMUE*ZMUE
    PTRA1(JL,JKM1) = 1.0_JPRB / ZDEN
    PRAY1(JL,JKM1) = ZBMU0 * PPIZAZ(JL,JKM1) * PTAUAZ(JL,JKM1)&
     & * ZMUE * PTRA1(JL,JKM1)

    ZBMU1 = 0.5_JPRB - 0.75_JPRB * PCGAZ(JL,JKM1) * 0.5_JPRB
    ZDEN1= 1.0_JPRB + (1.0_JPRB - PPIZAZ(JL,JKM1)&
     & + ZBMU1 * PPIZAZ(JL,JKM1)) * PTAUAZ(JL,JKM1) * 2.0_JPRB&
     & + (1-PPIZAZ(JL,JKM1)) * (1.0_JPRB - PPIZAZ(JL,JKM1)&
     & + 2.0_JPRB*ZBMU1*PPIZAZ(JL,JKM1))*PTAUAZ(JL,JKM1)&
     & * PTAUAZ(JL,JKM1)*4.0_JPRB  
    PTRA2(JL,JKM1) = 1.0_JPRB / ZDEN1
    PRAY2(JL,JKM1) = ZBMU1 * PPIZAZ(JL,JKM1) * PTAUAZ(JL,JKM1)&
     & * 2.0_JPRB *PTRA2(JL,JKM1)

    ZRR=1.0_JPRB/(1.0_JPRB-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1))
    PREFZ(JL,1,JK) = PRAY1(JL,JKM1)&
     & + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)&
     & * PTRA2(JL,JKM1)&
     & *ZRR  

    ZTR(JL,1,JKM1) = PTRA1(JL,JKM1)&
     & *ZRR  

    PREFZ(JL,2,JK) = PRAY1(JL,JKM1)&
     & + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
     & * PTRA2(JL,JKM1)   

    ZTR(JL,2,JKM1) = PTRA1(JL,JKM1)

  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZMUE = (1.0_JPRB-ZC0I(JL,1))*PSEC(JL)+ZC0I(JL,1)*1.66_JPRB
  PRMU0(JL,1)=1.0_JPRB/ZMUE
  PTRCLR(JL)=1.0_JPRB-ZC0I(JL,1)
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
    JKL = KLEV+1 - JK
    JKLP1 = JKL + 1
    DO JL = KIDIA,KFDIA
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
    ENDDO
  ENDDO

ELSE

  DO JAJ = 1 , 2
    DO JL = KIDIA,KFDIA
      PRJ(JL,JAJ,KLEV+1) = 1.0_JPRB
      PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
    ENDDO

    DO JK = 1 , KLEV
      JKL = KLEV+1 - JK
      JKLP1 = JKL + 1
      DO JL = KIDIA,KFDIA
        ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
        PRJ(JL,JAJ,JKL) = ZRE11
        PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
      ENDDO
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SWCLR',1,ZHOOK_HANDLE)
END SUBROUTINE SWCLR

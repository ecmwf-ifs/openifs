! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE UVCLR &
 &( YDERAD,YDERDI,KIDIA , KFDIA , KLON  , KLEV, &
 &  PALBP , PDSIG , PRAYL , PSEC, &
 &  PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ,&
 &  PRK   , PRMU0 , PTAUAZ, PTRA1 , PTRA2 , PTRCLR &
 &)  

!**** *UVCLR* - CLEAR-SKY COLUMN COMPUTATIONS

!     PURPOSE.
!     --------
!           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
!     CLEAR-SKY COLUMN

!**   INTERFACE.
!     ----------

!          *UVCLR* IS CALLED FROM *UVFLX*

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
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2005-10-04 based on SWCLR (but top to bottom)
   
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERAD   , ONLY : TERAD
USE YOERDI   , ONLY : TERDI

IMPLICIT NONE

TYPE(TERAD)       ,INTENT(INOUT):: YDERAD
TYPE(TERDI)       ,INTENT(INOUT):: YDERDI
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)  :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFDIA 

REAL(KIND=JPRB)   ,INTENT(IN)  :: PALBP(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDSIG(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PRAYL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PSEC(KLON) 

REAL(KIND=JPRB)   ,INTENT(IN)  :: PCGAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PPIZAZ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)  :: PTAUAZ(KLON,KLEV) 

REAL(KIND=JPRB)   ,INTENT(OUT) :: PRAY1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PRAY2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PREFZ(KLON,2,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PRJ(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PRK(KLON,6,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PRMU0(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PTRA1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PTRA2(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PTRCLR(KLON) 
!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!     ------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZC0I(KLON,KLEV+1)&
 & ,  ZCLE0(KLON,KLEV), ZCLEAR(KLON) &
 & ,  ZR21(KLON)&
 & ,  ZR23(KLON) , ZSS0(KLON) , ZSCAT(KLON)&
 & ,  ZTR(KLON,2,KLEV+1)  

INTEGER(KIND=JPIM) :: JA, JAJ, JK, JKP1, JL

REAL(KIND=JPRB) :: ZBMU0, ZBMU1, ZCORAE, ZDEN, ZDEN1, ZFACOA,&
 & ZGAP, ZMU1, ZMUE, ZRE11, ZTO, ZWW  
REAL(KIND=JPRB) :: ZRR,ZMU0,ZI2MU1,ZIMU1,ZIDEN,ZIDEN1

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UVCLR',0,ZHOOK_HANDLE)
ASSOCIATE(NOVLP=>YDERAD%NOVLP, &
 & REPCLC=>YDERDI%REPCLC)
!*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
!                --------------------------------------------

! N.B.: COMPUTED IN *UVRADI*

!     ------------------------------------------------------------------

!*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
!                ----------------------------------------------

DO JK = 1 , KLEV+1
  DO JA = 1 , 6
    DO JL = KIDIA,KFDIA
      PRJ(JL,JA,JK) = 0.0_JPRB
      PRK(JL,JA,JK) = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO

DO JL = KIDIA,KFDIA
  ZR23(JL) = 0.0_JPRB
  ZC0I(JL,1) = 0.0_JPRB
  ZCLEAR(JL) = 1.0_JPRB
  ZSCAT(JL) = 0.0_JPRB
ENDDO

DO JK = 1 , KLEV
  DO JL = KIDIA,KFDIA
    ZFACOA = 1.0_JPRB - PPIZAZ(JL,JK)*PCGAZ(JL,JK)*PCGAZ(JL,JK)
    ZCORAE = ZFACOA * PTAUAZ(JL,JK) * PSEC(JL)
    ZR21(JL) = EXP(-ZCORAE   )
    ZSS0(JL) = 1.0_JPRB-ZR21(JL)
    ZCLE0(JL,JK) = ZSS0(JL)

    IF (NOVLP == 1 .OR. NOVLP == 4) THEN
!* maximum-random      
      ZCLEAR(JL) = ZCLEAR(JL)&
       & *(1.0_JPRB-MAX(ZSS0(JL),ZSCAT(JL)))&
       & /(1.0_JPRB-MIN(ZSCAT(JL),1.0_JPRB-REPCLC))  
      ZC0I(JL,JK+1) = 1.0_JPRB - ZCLEAR(JL)
      ZSCAT(JL) = ZSS0(JL)
    ELSEIF (NOVLP == 2) THEN
!* maximum
      ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
      ZC0I(JL,JK+1) = ZSCAT(JL)
    ELSEIF (NOVLP == 3) THEN
!* random
      ZCLEAR(JL)=ZCLEAR(JL)*(1.0_JPRB-ZSS0(JL))
      ZSCAT(JL) = 1.0_JPRB - ZCLEAR(JL)
      ZC0I(JL,JK+1) = ZSCAT(JL)
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
!                -----------------------------------------------

DO JL = KIDIA,KFDIA
  PRAY1(JL,1) = 0.0_JPRB
  PRAY2(JL,1) = 0.0_JPRB
  PREFZ(JL,2,KLEV+1) = PALBP(JL)
  PREFZ(JL,1,KLEV+1) = PALBP(JL)
  PTRA1(JL,1) = 1.0_JPRB
  PTRA2(JL,1) = 1.0_JPRB
  PRMU0(JL,1) = 1.0_JPRB/PSEC(JL)
ENDDO
     
DO JK = KLEV,1,-1
  JKP1 = JK+1
  DO JL = KIDIA,KFDIA

!     ------------------------------------------------------------------

!*         3.1  EQUIVALENT ZENITH ANGLE
!               -----------------------

    ZMUE = (1.0_JPRB-ZC0I(JL,JK)) * PSEC(JL)+ ZC0I(JL,JK) * 1.66_JPRB
!-- no fiddling
!    ZMUE = PSEC(JL)
    PRMU0(JL,JK) = 1.0_JPRB/ZMUE
    ZMU0=PRMU0(JL,JK)

!     ------------------------------------------------------------------

!*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
!               ----------------------------------------------------

    ZGAP = PCGAZ(JL,JK)
    ZBMU0 = 0.5_JPRB - 0.75_JPRB * ZGAP *ZMU0
    ZWW = PPIZAZ(JL,JK)
    ZTO = PTAUAZ(JL,JK)
    ZDEN = 1.0_JPRB + (1.0_JPRB - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE&
     & + (1-ZWW) * (1.0_JPRB - ZWW +2.0_JPRB*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE  
    ZIDEN=1.0_JPRB / ZDEN
    PRAY1(JL,JKP1) = ZBMU0 * ZWW * ZTO * ZMUE * ZIDEN
    PTRA1(JL,JKP1) = ZIDEN

    ZMU1 =  0.5_JPRB
    ZIMU1 = 2.0_JPRB
    ZI2MU1= 4.0_JPRB
    ZBMU1 = 0.5_JPRB - 0.75_JPRB * ZGAP * ZMU1
    ZDEN1= 1.0_JPRB + (1.0_JPRB - ZWW + ZBMU1 * ZWW) * ZTO * ZIMU1&
     & + (1-ZWW) * (1.0_JPRB - ZWW +2.0_JPRB*ZBMU1*ZWW)*ZTO*ZTO*ZI2MU1  
    ZIDEN1=1.0_JPRB / ZDEN1
    PRAY2(JL,JKP1) = ZBMU1 * ZWW * ZTO * ZIMU1 *ZIDEN1
    PTRA2(JL,JKP1) = ZIDEN1

    ZRR=1.0_JPRB/(1.0_JPRB-PRAY2(JL,JKP1)*PREFZ(JL,1,JKP1))
    PREFZ(JL,1,JK) = PRAY1(JL,JKP1) + PREFZ(JL,1,JKP1) * PTRA1(JL,JKP1) * PTRA2(JL,JKP1) * ZRR  

    ZTR(JL,1,JKP1) = PTRA1(JL,JKP1) *ZRR  

    PREFZ(JL,2,JK) = PRAY1(JL,JKP1) + PREFZ(JL,2,JKP1) * PTRA1(JL,JKP1) * PTRA2(JL,JKP1)   

    ZTR(JL,2,JKP1) = PTRA1(JL,JKP1)

  ENDDO
ENDDO

DO JL = KIDIA,KFDIA
  ZMUE = (1.0_JPRB-ZC0I(JL,KLEV+1))*PSEC(JL)+ZC0I(JL,KLEV+1)*1.66_JPRB
!-- no fiddling
!  ZMUE = PSEC(JL)
  PRMU0(JL,KLEV+1)=1.0_JPRB/ZMUE
  PTRCLR(JL)=1.0_JPRB-ZC0I(JL,KLEV+1)
ENDDO

!     ------------------------------------------------------------------

!*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
!                 -------------------------------------------------

JAJ = 2
DO JL = KIDIA,KFDIA
  PRJ(JL,JAJ,1) = 1.0_JPRB
  PRK(JL,JAJ,1) = PREFZ(JL, 1,1)
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
IF (LHOOK) CALL DR_HOOK('UVCLR',1,ZHOOK_HANDLE)
END SUBROUTINE UVCLR


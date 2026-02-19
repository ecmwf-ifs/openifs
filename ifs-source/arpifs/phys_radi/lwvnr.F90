! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LWVNR &
 & ( KIDIA, KFDIA, KLON  , KLEV , KUAER,&
 & PABCU, PDBSL, PGA   , PGB,&
 & PADJD, PADJU, PCNTRB, PDBDT, PDWFSU &
 & )  

!**** *LWVNR*   - L.W., VERTICAL INTEGRATION, NEARBY LAYERS

!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION ON NEARBY LAYERS
!           TO GIVE LONGWAVE FLUXES OR RADIANCES

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU : (KLON,NUA,3*KLEV+1)  ; ABSORBER AMOUNTS
! PDBSL  : (KLON,KLEV*2)       ; SUB-LAYER PLANCK FUNCTION GRADIENT
! PGA, PGB                     ; PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PADJ.. : (KLON,KLEV+1)       ; CONTRIBUTION OF ADJACENT LAYERS
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
! PDBDT  : (KLON,NUA,KLEV)     ; LAYER PLANCK FUNCTION GRADIENT
! PDWFSU : (KLON,NSIL)         ; SPECTRAL DOWNWARD FLUX AT SURFACE

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
!     CONTRIBUTIONS OF THE ADJACENT LAYERS USING A GAUSSIAN QUADRATURE

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      JEAN-JACQUES MORCRETTE  *ECMWF*
!      ORIGINAL : 89-07-14

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Janiskova   22-Nov-2006 reduced version for H2O and CO2 based on lwvn
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOELW    , ONLY : NSIL     ,NIPD     ,NTRA     ,NUA      ,&
 &                    NG1      ,NG1P1    ,WG1  

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KUAER 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDBSL(KLON,NSIL,KLEV*2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGA(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGB(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PADJD(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PADJU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCNTRB(KLON,KLEV+1,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDBDT(KLON,NSIL,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDWFSU(KLON,NSIL) 

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZTT(KLON,NTRA), ZTT1(KLON,NTRA), ZTT2(KLON,NTRA)
REAL(KIND=JPRB) :: ZUU(KLON,NUA)

INTEGER(KIND=JPIM) :: IBS, IDD, IM12, IMU, IND, INU, IXD, IXU,&
 & JA, JG, JK, JK1, JK2, JL, JNU  

REAL(KIND=JPRB) :: ZWTR, ZWTR1, ZWTR2, ZWTR3, ZWTR4, ZWTR5, ZWTR6
REAL(KIND=JPRB) :: ZXN, ZXD, ZZ, ZXDIV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!#include "lwtt.intfb.h"

!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------

IF (LHOOK) CALL DR_HOOK('LWVNR',0,ZHOOK_HANDLE)
DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PADJD(JL,JK) = 0.0_JPRB
    PADJU(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO

!*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
!                  ---------------------------------

!DO JA = 1 , NTRA
DO JA = 1 , KUAER
  DO JL = KIDIA,KFDIA
    ZTT (JL,JA) = 1.0_JPRB
    ZTT1(JL,JA) = 1.0_JPRB
    ZTT2(JL,JA) = 1.0_JPRB
  ENDDO
ENDDO

!DO JA = 1 , NUA
DO JA = 1 , KUAER
  DO JL = KIDIA,KFDIA
    ZUU(JL,JA) = 0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.      VERTICAL INTEGRATION
!                  --------------------

!*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
!                  ---------------------------------

DO JK = 1 , KLEV

!*         2.1.1   DOWNWARD LAYERS
!                  ---------------

  IM12 = 2 * (JK - 1)
  IND = (JK - 1) * NG1P1 + 1
  IXD = IND
  INU = JK * NG1P1 + 1
  IXU = IND

  DO JG = 1 , NG1
    IBS = IM12 + JG
    IDD = IXD + JG

    DO JA = 1 , KUAER
      DO JL = KIDIA,KFDIA
        ZUU(JL,JA) = PABCU(JL,JA,IND) - PABCU(JL,JA,IDD)
      ENDDO
    ENDDO

    DO JA = 1 , 8
      DO JL = KIDIA,KFDIA
        ZZ = SQRT(ZUU(JL,JA))
        ZXD = PGB( JL,JA,1,JK) + ZZ* (PGB( JL,JA,2,JK) + ZZ )
        ZXN = PGA( JL,JA,1,JK) + ZZ*PGA( JL,JA,2,JK)
        ZXDIV = 1.0_JPRB/ZXD
        ZTT(JL,JA) = ZXN*ZXDIV
      ENDDO
    ENDDO
    DO JL = KIDIA,KFDIA
      IF (ZTT(JL,3) < 0.0_JPRB) THEN
        ZTT (JL,3) = 0.0_JPRB
      ENDIF
      ZTT (JL, 9) = ZTT (JL, 8)
    ENDDO

    DO JL = KIDIA,KFDIA
      ZWTR1=PDBSL(JL,1,IBS)*ZTT(JL,1)
      ZWTR2=PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)
      ZWTR3=PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)
      ZWTR4=PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)
      ZWTR5=PDBSL(JL,5,IBS)*ZTT(JL,3)
      ZWTR6=PDBSL(JL,6,IBS)*ZTT(JL,6)
      ZWTR=ZWTR1+ZWTR2+ZWTR3+ZWTR4+ZWTR5+ZWTR6
      PADJD(JL,JK) = PADJD(JL,JK) + ZWTR * WG1(JG)
      IF (JK == 1) THEN
        PDWFSU(JL,1)=PDWFSU(JL,1)+WG1(JG)*ZWTR1
        PDWFSU(JL,2)=PDWFSU(JL,2)+WG1(JG)*ZWTR2
        PDWFSU(JL,3)=PDWFSU(JL,3)+WG1(JG)*ZWTR3
        PDWFSU(JL,4)=PDWFSU(JL,4)+WG1(JG)*ZWTR4
        PDWFSU(JL,5)=PDWFSU(JL,5)+WG1(JG)*ZWTR5
        PDWFSU(JL,6)=PDWFSU(JL,6)+WG1(JG)*ZWTR6
      ENDIF
    ENDDO

!*         2.1.2   UPWARD LAYERS
!                  -------------

    IMU = IXU + JG
    DO JA = 1 , KUAER
      DO JL = KIDIA,KFDIA
        ZUU(JL,JA) = PABCU(JL,JA,IMU) - PABCU(JL,JA,INU)
      ENDDO
    ENDDO

    DO JA = 1 , 8
      DO JL = KIDIA,KFDIA
        ZZ = SQRT(ZUU(JL,JA))
        ZXD = PGB( JL,JA,1,JK) + ZZ* (PGB( JL,JA,2,JK) + ZZ )
        ZXN = PGA( JL,JA,1,JK) + ZZ*PGA( JL,JA,2,JK)
        ZXDIV = 1.0_JPRB/ZXD
        ZTT(JL,JA) = ZXN*ZXDIV
      ENDDO
    ENDDO
    DO JL = KIDIA,KFDIA
      IF (ZTT(JL,3) < 0.0_JPRB) THEN
        ZTT (JL,3) = 0.0_JPRB
      ENDIF
      ZTT (JL, 9) = ZTT (JL, 8)
    ENDDO

    DO JL = KIDIA,KFDIA
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1) &
       & +PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7) &
       & +PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8) &
       & +PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9) &
       & +PDBSL(JL,5,IBS)*ZTT(JL,3) &
       & +PDBSL(JL,6,IBS)*ZTT(JL,6)
      PADJU(JL,JK+1) = PADJU(JL,JK+1) + ZWTR * WG1(JG)
    ENDDO

  ENDDO

  DO JL = KIDIA,KFDIA
    PCNTRB(JL,JK+1,JK) = PADJD(JL,JK)
    PCNTRB(JL,JK,JK+1) = PADJU(JL,JK+1)
    PCNTRB(JL,JK  ,JK) = 0.0_JPRB
  ENDDO

ENDDO

DO JK = 1 , KLEV
  JK2 = 2 * JK
  JK1 = JK2 - 1

  DO JNU = 1 , NSIL
    DO JL = KIDIA,KFDIA
      PDBDT(JL,JNU,JK) = PDBSL(JL,JNU,JK1) + PDBSL(JL,JNU,JK2)
    ENDDO
  ENDDO
ENDDO

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LWVNR',1,ZHOOK_HANDLE)
END SUBROUTINE LWVNR

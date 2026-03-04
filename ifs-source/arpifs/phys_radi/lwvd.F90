! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LWVD &
 & ( YDPHNC,KIDIA,  KFDIA, KLON , KLEV  , KTRAER,&
 & PABCU,  PDBDT,&
 & PGA  ,  PGB,&
 & PCNTRB, PDISD, PDISU, PDWFSU,&
! for adjoint computation
 & PTTA , PTTB)

!**** *LWVD*   - L.W., VERTICAL INTEGRATION, DISTANT LAYERS

!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION ON THE DISTANT LAYERS

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU  : (KLON,NUA,3*KLEV+1) ; ABSORBER AMOUNTS
! PDBDT  : (KLON,KLEV)         ; LAYER PLANCK FUNCTION GRADIENT
! PGA, PGB                     ; PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PCNTRB : (KLON,KLEV+1,KLEV+1); ENERGY EXCHANGE MATRIX
! PDIS.. : (KLON,KLEV+1)       ; CONTRIBUTION BY DISTANT LAYERS
! PDWFSU : (KLON,NSIL)         ; SPECTRAL DOWNWARD FLUX AT SURFACE

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
!     CONTRIBUTIONS OF THE DISTANT LAYERS USING TRAPEZOIDAL RULE

!     EXTERNALS.
!     ----------

!          *LWTT*

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
!      M. Janiskova  16-Jan-2006 Cleaning
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOELW    , ONLY : NSIL     ,NIPD     ,NTRA     ,NUA      ,NG1P1
USE YOPHNC   , ONLY : TPHNC

IMPLICIT NONE

TYPE(TPHNC)       ,INTENT(IN)    :: YDPHNC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAER 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDBDT(KLON,NSIL,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGA(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGB(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCNTRB(KLON,KLEV+1,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISD(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDWFSU(KLON,NSIL)
 ! for adjoint computation
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTA(KLON,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTB(KLON,KTRAER,KLEV,KLEV+1)

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZTT(KLON,NTRA)
REAL(KIND=JPRB) :: ZTTP(KLON,NTRA,KLEV+1)

INTEGER(KIND=JPIM) :: IJKL, IKD1, IKD2, IKJ, IKJP1, IKM1, IKN,&
 & IKP1, IKU1, IKU2, JA, JK, JKJ, JL, JLK, IJA

REAL(KIND=JPRB) :: ZWW1, ZWW2, ZWW3, ZWW4, ZWW5, ZWW6
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "lwttm.intfb.h"

!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------

IF (LHOOK) CALL DR_HOOK('LWVD',0,ZHOOK_HANDLE)
ASSOCIATE(LWLOPT=>YDPHNC%LWLOPT)
DO JK = 1, KLEV+1
  DO JL = KIDIA,KFDIA
    PDISD(JL,JK) = 0.0_JPRB
    PDISU(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO

!*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
!                  ---------------------------------

DO JA = 1, NTRA
  DO JL = KIDIA,KFDIA
    ZTT (JL,JA) = 1.0_JPRB
    ZTTP(JL,JA,:) = 1.0_JPRB
  ENDDO
ENDDO

IF (LWLOPT) THEN
!for adjoint computation
  DO JA = 1, NTRA
    DO JL = KIDIA,KFDIA 
      PTTA(JL,JA,1,:) = 1.0_JPRB
      PTTB(JL,JA,1,:) = 1.0_JPRB
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*         2.      VERTICAL INTEGRATION
!                  --------------------

!*         2.2     CONTRIBUTION FROM DISTANT LAYERS
!                  ---------------------------------

!*         2.2.1   DISTANT AND ABOVE LAYERS
!                  ------------------------

!*         2.2.2   FIRST UPPER LEVEL
!                  -----------------

DO JK = 1 , KLEV-1
  IKP1=JK+1
  IKN=(JK-1)*NG1P1+1
  IKD1= JK  *NG1P1+1

  CALL LWTTM&
   & ( KIDIA         , KFDIA          , KLON,&
   & PGA(1,1,1,JK) , PGB(1,1,1,JK),&
   & PABCU(1,1,IKN), PABCU(1,1,IKD1), ZTTP(1,1,JK) )  

!*         2.2.3   HIGHER UP
!                  ---------

  DO JKJ=IKP1,KLEV
    IKJP1=JKJ+1
    IKD2= JKJ  *NG1P1+1

    CALL LWTTM&
     & ( KIDIA         , KFDIA          , KLON,&
     & PGA(1,1,1,JKJ), PGB(1,1,1,JKJ),&
     & PABCU(1,1,IKN), PABCU(1,1,IKD2), ZTTP(1,1,JKJ) )  

    DO JA = 1, KTRAER
      DO JL = KIDIA,KFDIA
        ZTT(JL,JA) = (ZTTP(JL,JA,JKJ-1)+ZTTP(JL,JA,JKJ))*0.5_JPRB
      ENDDO
    ENDDO

    IF (LWLOPT) THEN
! for adjoint computation
      DO JA = 1, KTRAER
        IJA = (JA-1)*KLEV+JKJ
        DO JL = KIDIA,KFDIA
          PTTA(JL,JA,JKJ,JK) = ZTT(JL,JA)
        ENDDO
      ENDDO
    ENDIF

    DO JL = KIDIA,KFDIA
      ZWW1=PDBDT(JL,1,JKJ)*ZTT(JL,1)          *ZTT(JL,10)
      ZWW2=PDBDT(JL,2,JKJ)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
      ZWW3=PDBDT(JL,3,JKJ)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
      ZWW4=PDBDT(JL,4,JKJ)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
      ZWW5=PDBDT(JL,5,JKJ)*ZTT(JL,3)          *ZTT(JL,14)
      ZWW6=PDBDT(JL,6,JKJ)*ZTT(JL,6)          *ZTT(JL,15)
      PCNTRB(JL,IKJP1,JK)=ZWW1+ZWW2+ZWW3+ZWW4+ZWW5+ZWW6
      PDISD(JL,JK)=PDISD(JL,JK)+PCNTRB(JL,IKJP1,JK)
      IF (JK == 1) THEN
        PDWFSU(JL,1)=PDWFSU(JL,1)+ZWW1
        PDWFSU(JL,2)=PDWFSU(JL,2)+ZWW2
        PDWFSU(JL,3)=PDWFSU(JL,3)+ZWW3
        PDWFSU(JL,4)=PDWFSU(JL,4)+ZWW4
        PDWFSU(JL,5)=PDWFSU(JL,5)+ZWW5
        PDWFSU(JL,6)=PDWFSU(JL,6)+ZWW6
      ENDIF
    ENDDO

  ENDDO
ENDDO

!*         2.2.4   DISTANT AND BELOW LAYERS
!                  ------------------------

!*         2.2.5   FIRST LOWER LEVEL
!                  -----------------

DO JK=3,KLEV+1
  IKN=(JK-1)*NG1P1+1
  IKM1=JK-1
  IKJ=JK-2
  IKU1= IKJ  *NG1P1+1

  CALL LWTTM&
   & ( KIDIA          , KFDIA         , KLON,&
   & PGA(1,1,1,IKJ) , PGB(1,1,1,IKJ),&
   & PABCU(1,1,IKU1), PABCU(1,1,IKN), ZTTP(1,1,1) )  

!*         2.2.6   DOWN BELOW
!                  ----------

  DO JLK=1,IKJ
    IJKL=IKM1-JLK
    IKU2=(IJKL-1)*NG1P1+1

    CALL LWTTM&
     & ( KIDIA          , KFDIA          , KLON,&
     & PGA(1,1,1,IJKL), PGB(1,1,1,IJKL),&
     & PABCU(1,1,IKU2), PABCU(1,1,IKN) , ZTTP(1,1,JLK+1) )  

    DO JA = 1, KTRAER
      DO JL = KIDIA,KFDIA
        ZTT(JL,JA) = (ZTTP(JL,JA,JLK)+ZTTP(JL,JA,JLK+1))*0.5_JPRB
      ENDDO
    ENDDO

    IF (LWLOPT) THEN
! for adjoint computation
      DO JA = 1, KTRAER
        IJA = (JA-1)*KLEV+JLK
        DO JL = KIDIA,KFDIA
          PTTB(JL,JA,JLK,JK) = ZTT(JL,JA)
        ENDDO
      ENDDO
    ENDIF

    DO JL = KIDIA,KFDIA
      PCNTRB(JL,IJKL,JK)=PDBDT(JL,1,IJKL)*ZTT(JL,1)*ZTT(JL,10)&
       & +PDBDT(JL,2,IJKL)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)&
       & +PDBDT(JL,3,IJKL)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)&
       & +PDBDT(JL,4,IJKL)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)&
       & +PDBDT(JL,5,IJKL)*ZTT(JL,3)          *ZTT(JL,14)&
       & +PDBDT(JL,6,IJKL)*ZTT(JL,6)          *ZTT(JL,15)  
      PDISU(JL,JK)=PDISU(JL,JK)+PCNTRB(JL,IJKL,JK)
    ENDDO

  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LWVD',1,ZHOOK_HANDLE)
END SUBROUTINE LWVD

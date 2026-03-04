! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUDUDV2 &
 & (YDECUMF2,KIDIA,    KFDIA,    KLON,     KLEV,&
 & KTOPM2,   KTYPE,    KCBOT,    KCTOP,    LDCUM,    PTSPHY,&
 & PAPH,     PUEN,     PVEN,     PMFU,     PMFD,&
 & PUU,      PUD,      PVU,      PVD,&
 & PTENU,    PTENV  )  

!**** *CUDUDV2* - UPDATES U AND V TENDENCIES,
!                 DOES GLOBAL DIAGNOSTIC OF DISSIPATION 
!                 VERSION FOR SIMPLIFIED MASS FLUX SCHEME

!          P. LOPEZ          E.C.M.W.F.     12/2003 

!          Originally duplicated from 
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89

!**   INTERFACE.
!     ----------

!          *CUDUDV2* IS CALLED FROM *CUMASTRN2*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
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
!    *KCTOP*        CLOUD TOP LEVEL

!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                      S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
!    *PUD*          U-VELOCITY IN DOWNDRAFTS                       M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S
!    *PVD*          V-VELOCITY IN DOWNDRAFTS                       M/S

!    UPDATED PARAMETERS (REAL):

!    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2

!          MODIFICATIONS
!          -------------
!        M.Hamrud      01-Oct-2003  CY28 Cleaning
!        P. Lopez      22-Oct-2007  Added implicit version
!----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG
USE YOECUMF2 , ONLY : TECUMF2

IMPLICIT NONE

TYPE(TECUMF2)     ,INTENT(INOUT) :: YDECUMF2
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOPM2 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENV(KLON,KLEV) 

REAL(KIND=JPRB) :: ZUEN(KLON,KLEV),   ZVEN(KLON,KLEV),&
                 & ZMFUU(KLON,KLEV),  ZMFDU(KLON,KLEV),&
                 & ZMFUV(KLON,KLEV),  ZMFDV(KLON,KLEV)  

INTEGER(KIND=JPIM) :: IK, IKB, JK, JL

REAL(KIND=JPRB) :: ZZP, ZDELP, ZP1, ZP2, ZIMP, ZTSPHY

! ALLOCATABLE ARAYS
REAL(KIND=JPRB),   DIMENSION(:,:), ALLOCATABLE :: ZDUDT, ZDVDT, ZGDP
REAL(KIND=JPRB),   DIMENSION(:,:), ALLOCATABLE :: ZB,  ZR1,  ZR2
LOGICAL,DIMENSION(:,:),   ALLOCATABLE :: LLCUMBAS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cubidiag.intfb.h"
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUDUDV2',0,ZHOOK_HANDLE)
ASSOCIATE(RMFSOLUV2=>YDECUMF2%RMFSOLUV2, RMFSOLRHS2=>YDECUMF2%RMFSOLRHS2)
ZIMP=1.0_JPRB-RMFSOLUV2
ZTSPHY=1.0_JPRB/PTSPHY

ALLOCATE(ZDUDT(KLON,KLEV))
ALLOCATE(ZDVDT(KLON,KLEV))
ALLOCATE(ZGDP(KLON,KLEV))

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF (LDCUM(JL)) THEN
      ZUEN(JL,JK)=PUEN(JL,JK)
      ZVEN(JL,JK)=PVEN(JL,JK)
      ZDELP=1.0_JPRB/(PAPH(JL,JK+1)-PAPH(JL,JK))
      ZGDP(JL,JK)=RG*ZDELP
    ENDIF
  ENDDO
ENDDO

!----------------------------------------------------------------------

!*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
!                  ----------------------------------------------

DO JK=KTOPM2,KLEV
  IK=JK-1
  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL)) THEN
      ZMFUU(JL,JK)=PMFU(JL,JK)*(PUU(JL,JK)-ZIMP*ZUEN(JL,IK))
      ZMFUV(JL,JK)=PMFU(JL,JK)*(PVU(JL,JK)-ZIMP*ZVEN(JL,IK))
      ZMFDU(JL,JK)=PMFD(JL,JK)*(PUD(JL,JK)-ZIMP*ZUEN(JL,IK))
      ZMFDV(JL,JK)=PMFD(JL,JK)*(PVD(JL,JK)-ZIMP*ZVEN(JL,IK))
    ENDIF
  ENDDO
ENDDO

! linear fluxes below cloud
IF (RMFSOLUV2 == 0.0_JPRB) THEN
  DO JK=KTOPM2,KLEV
!DIR$ IVDEP
!OCL NOVREC
    DO JL=KIDIA,KFDIA
      IF (LDCUM(JL) .AND. JK > KCBOT(JL)) THEN
        IKB=KCBOT(JL)
        ZP1=PAPH(JL,KLEV+1)-PAPH(JL,JK)
        ZP2=PAPH(JL,KLEV+1)-PAPH(JL,IKB)
        ZZP=ZP1/ZP2
        IF (KTYPE(JL) == 3) THEN
          ZZP=ZZP*ZZP
        ENDIF
        ZMFUU(JL,JK)=ZMFUU(JL,IKB)*ZZP
        ZMFUV(JL,JK)=ZMFUV(JL,IKB)*ZZP
        ZMFDU(JL,JK)=ZMFDU(JL,IKB)*ZZP
        ZMFDV(JL,JK)=ZMFDV(JL,IKB)*ZZP
      ENDIF
    ENDDO
  ENDDO
ENDIF

!*    1.2          COMPUTE TENDENCIES
!                  ------------------

DO JK=KTOPM2,KLEV

  IF (JK < KLEV) THEN
    DO JL=KIDIA,KFDIA
      IF (LDCUM(JL)) THEN
        ZDUDT(JL,JK)=ZGDP(JL,JK)*&
         & (ZMFUU(JL,JK+1)-ZMFUU(JL,JK)+&
         &  ZMFDU(JL,JK+1)-ZMFDU(JL,JK))  
        ZDVDT(JL,JK)=ZGDP(JL,JK)*&
         & (ZMFUV(JL,JK+1)-ZMFUV(JL,JK)+&
         &  ZMFDV(JL,JK+1)-ZMFDV(JL,JK))  
      ENDIF
    ENDDO

  ELSE
    DO JL=KIDIA,KFDIA
      IF (LDCUM(JL)) THEN
        ZDUDT(JL,JK)=-ZGDP(JL,JK)*(ZMFUU(JL,JK)+ZMFDU(JL,JK))
        ZDVDT(JL,JK)=-ZGDP(JL,JK)*(ZMFUV(JL,JK)+ZMFDV(JL,JK))
      ENDIF
    ENDDO
  ENDIF

ENDDO

IF (RMFSOLUV2 == 0.0_JPRB) THEN

!*    1.3          UPDATE TENDENCIES
!                  -----------------

  DO JK=KTOPM2,KLEV
    DO JL=KIDIA,KFDIA
      IF (LDCUM(JL)) THEN
        PTENU(JL,JK)=PTENU(JL,JK)+ZDUDT(JL,JK)
        PTENV(JL,JK)=PTENV(JL,JK)+ZDVDT(JL,JK)
      ENDIF
    ENDDO
  ENDDO

ELSE
!----------------------------------------------------------------------
   
!*      1.6          IMPLICIT SOLUTION
!                    -----------------

! Fill bi-diagonal Matrix vectors A=k-1, B=k;
! reuse ZMFUU=A and ZB=B; 
! ZDUDT and ZDVDT correspond to the RHS ("constants") of the equation
! The solution is in ZR1 and ZR2
  
   ALLOCATE(ZB(KLON,KLEV))
   ALLOCATE(ZR1(KLON,KLEV))
   ALLOCATE(ZR2(KLON,KLEV))
   ALLOCATE(LLCUMBAS(KLON,KLEV))

   LLCUMBAS(:,:)=.FALSE.
   ZB(:,:)=1.0_JPRB
   ZMFUU(:,:)=0.0_JPRB
  
! Fill vectors A, B and RHS 
  
   DO JK=KTOPM2,KLEV
      IK=JK+1
      DO JL=KIDIA,KFDIA
        LLCUMBAS(JL,JK)=LDCUM(JL) .AND. JK >= KCTOP(JL)-1 
        IF (LLCUMBAS(JL,JK)) THEN
          ZZP=RMFSOLUV2*PTSPHY*ZGDP(JL,JK)
          ZMFUU(JL,JK)=-ZZP*(PMFU(JL,JK)+PMFD(JL,JK))
          ZDUDT(JL,JK)=(ZDUDT(JL,JK)+PTENU(JL,JK)*RMFSOLRHS2)*PTSPHY+ZUEN(JL,JK)
          ZDVDT(JL,JK)=(ZDVDT(JL,JK)+PTENV(JL,JK)*RMFSOLRHS2)*PTSPHY+ZVEN(JL,JK)
          IF (JK < KLEV) THEN
            ZB(JL,JK)=1.0_JPRB+ZZP*(PMFU(JL,IK)+PMFD(JL,IK))
          ELSE
            ZB(JL,JK)=1.0_JPRB
          ENDIF
        ENDIF
      ENDDO
   ENDDO
  
   CALL CUBIDIAG &
      &( KIDIA, KFDIA, KLON, KLEV, &
      &  KCTOP, LLCUMBAS, &
      &  ZMFUU,    ZB,    ZDUDT,   ZR1 )
  
   CALL CUBIDIAG &
      &( KIDIA, KFDIA, KLON, KLEV, &
      &  KCTOP, LLCUMBAS, &
      &  ZMFUU,    ZB,    ZDVDT,   ZR2 )
  
! Compute tendencies
  
   DO JK=KTOPM2,KLEV
      DO JL=KIDIA,KFDIA
        IF (LLCUMBAS(JL,JK)) THEN
          PTENU(JL,JK)=PTENU(JL,JK)*(1.0_JPRB-RMFSOLRHS2)+(ZR1(JL,JK)-ZUEN(JL,JK))*ZTSPHY
          PTENV(JL,JK)=PTENV(JL,JK)*(1.0_JPRB-RMFSOLRHS2)+(ZR2(JL,JK)-ZVEN(JL,JK))*ZTSPHY
        ENDIF
      ENDDO
   ENDDO
  
   DEALLOCATE(LLCUMBAS)
   DEALLOCATE(ZR2)
   DEALLOCATE(ZR1)
   DEALLOCATE(ZB)
  
ENDIF
!----------------------------------------------------------------------

DEALLOCATE(ZGDP)
DEALLOCATE(ZDVDT)
DEALLOCATE(ZDUDT)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUDUDV2',1,ZHOOK_HANDLE)
END SUBROUTINE CUDUDV2

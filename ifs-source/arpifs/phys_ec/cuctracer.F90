! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUCTRACER &
 & (YDCUMFS, YDECUMF2, YDECUMF,&
 & KIDIA,    KFDIA,    KLON,     KLEV,  KTRAC,&
 & KCTOP,    KDTOP,    KTYPE,&
 & LDCUM,    LDDRAF,   PTSPHY,   PAPH,  PAP,&
 & PMFU,     PMFD,     PMFUO,    PMFDO, PUDRATE,  PDDRATE,&
 & PDMFUP,   PDMFDP,&
 & PCEN,     PTENC,    PSCAV  )  

!**** *CUCTRACER* - COMPUTE CONVECTIVE TRANSPORT OF CHEM. TRACERS
!                   IMPORTANT: ROUTINE IS FOR POSITIVE DEFINIT QUANTITIES
!                              SCAVENGING IS DONE IF SCAVENGING COEFFICIENT>0

!          P.BECHTOLD        E.C.M.W.F.              11/02/2004

!**   INTERFACE.
!     ----------

!          *CUTRACER* IS CALLED FROM *CUMASTR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF CHEMICAL TRACERS

!    *KCTOP*        CLOUD TOP  LEVEL
!    *KDTOP*        DOWNDRAFT TOP LEVEL
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION


!    INPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDDRAF*       FLAG: .TRUE. IF DOWNDRAFTS EXIST

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       PHYSICS TIME-STEP                              S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PAPH          PROVISIONAL PRESSURE ON FULL LEVELS            PA
!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATION
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUDRATE*      UPDRAFT DETRAINMENT                           KG/(M2*S)
!    *PDDRATE*      DOWNDRAFT DETRAINMENT                         KG/(M2*S)
!    *PDMFUP*       UPDRAUGT PRECIP (PRODUCTION) IN LAYER         KG/(M2*S)
!    *PDMFDP*       DDRAUGHT PRECIP (PRODUCTION) IN LAYER         KG/(M2*S)
!    *PSCAV*        WET SCAVENGING COEFFICIENTS                   UNITLESS 

!    UPDATED PARAMETERS (REAL):

!    *PTENC*        UPDATED TENDENCY OF CHEM. TRACERS              1/S

!          METHOD
!          -------
!     EXPLICIT UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
!     DEPENDING ON VALUE OF RMFSOLCT: 0=EXPLICIT 0-1 SEMI-IMPLICIT >=1 IMPLICIT

!     FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
!     FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
!                            TO CORRECT TENDENCIES BELOW CLOUD BASE


!------------------------------------------------------------------------------------
!     COMMENTS FOR OFFLINE USERS IN CHEMICAL TRANSPORT MODELS
!     (i.e. reading mass fluxes and detrainment rates from ECMWF archive:
!      ------------------------------------------------------------------
!     KCTOP IS FIRST LEVEL FROM TOP WHERE PMFU>0      
!     KDTOP IS FIRST LEVEL FROM TOP WHERE PMFD<0      
!     ATTENTION: ON ARCHIVE DETRAINMENT RATES HAVE UNITS KG/(M3*S), SO FOR USE
!                IN CURRENT ROUTINE YOU HAVE TO MULTIPLY ARCHIVED VALUES BY DZ
!     LDCUM  IS TRUE IF CONVECTION EXISTS, i.e. IF PMFU>0 IN COLUMN OR IF
!                       KCTOP>0 AND KCTOP<KLEV
!     LDDRAF IS TRUE IF DOWNDRAUGHTS EXIST IF PMFD<0 IN COLUMN OR IF
!                       KDTOP>0 AND KDTOP<KLEV
!     IF MASSFLUX SATISFIES CFL CRITERIUM M<=DP/Dt IT IS SUFFICIENT TO 
!     ONLY CONSIDER EXPLICIT SOLUTION (RMFSOLCT=0), IN THIS CASE
!     YOU CAN IGNORE IMPLICIT PART 7.0 OF CURRENT ROUTINE
!------------------------------------------------------------------------------------

!          EXTERNALS
!          ---------
!          CUBIDIAG

!          MODIFICATIONS
!          -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P.Bechtold    07-Oct-2009 add scavenging

!----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG
USE YOECUMF  , ONLY : TECUMF
USE YOECUMF2 , ONLY : TECUMF2
USE YOMCUMFS , ONLY : TCUMFS

IMPLICIT NONE

TYPE(TCUMFS)      ,INTENT(INOUT) :: YDCUMFS
TYPE(TECUMF)      ,INTENT(INOUT) :: YDECUMF
TYPE(TECUMF2)     ,INTENT(INOUT) :: YDECUMF2
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON) 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
LOGICAL           ,INTENT(IN)    :: LDDRAF(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFUO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFDO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUDRATE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDDRATE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDMFUP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDMFDP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCAV(KTRAC) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
INTEGER(KIND=JPIM) :: IK, IM, JK, JL, JN

REAL(KIND=JPRB) :: ZZP, ZMFA, ZIMP, ZERATE, ZPOSI, ZTSPHY, ZC, ZADVW(KLON)
REAL(KIND=JPRB) :: ZRMFCMIN, ZRMFSOLCT

!     ALLOCATABLE ARAYS
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZCEN, ZCU, ZCD, ZTENC, ZMFC
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZDP, ZB,  ZR1
LOGICAL, DIMENSION(:,:),  ALLOCATABLE :: LLCUMASK, LLCUMBAS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cubidiag.intfb.h"
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUCTRACER',0,ZHOOK_HANDLE)
ASSOCIATE(LECUMFS=>YDCUMFS%LECUMFS, &
 & RMFCMIN=>YDECUMF%RMFCMIN, RMFSOLCT=>YDECUMF%RMFSOLCT, &
 & RMFADVW=>YDECUMF%RMFADVW, RMFADVWDD=>YDECUMF%RMFADVWDD,&
 & RMFCMIN2=>YDECUMF2%RMFCMIN2, RMFSOLCT2=>YDECUMF2%RMFSOLCT2)
IF (LECUMFS) THEN
  ZRMFSOLCT=RMFSOLCT2
  ZRMFCMIN=RMFCMIN2
ELSE
  ZRMFSOLCT=RMFSOLCT
  ZRMFCMIN=RMFCMIN
ENDIF

ZIMP=1.0_JPRB-ZRMFSOLCT
ZTSPHY=1.0_JPRB/PTSPHY

DO JL=KIDIA,KFDIA
  ZADVW(JL)=0.0_JPRB
  IF(KTYPE(JL)==1) ZADVW(JL)=RMFADVW
ENDDO

ALLOCATE(ZCEN(KLON,KLEV,KTRAC)) !Half-level environmental values
ALLOCATE(ZCU(KLON,KLEV,KTRAC))  !Updraft values
ALLOCATE(ZCD(KLON,KLEV,KTRAC))  !Downdraft values
ALLOCATE(ZTENC(KLON,KLEV,KTRAC))!Tendency
ALLOCATE(ZMFC(KLON,KLEV,KTRAC)) !Fluxes
ALLOCATE(ZDP(KLON,KLEV))        !Pressure difference
ALLOCATE(LLCUMASK(KLON,KLEV))   !Mask for convection

! Initialize Cumulus mask + some setups

DO JK=2,KLEV
   DO JL=KIDIA,KFDIA
     LLCUMASK(JL,JK)=.FALSE.
     IF(LDCUM(JL)) THEN
       ZDP(JL,JK)=RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
       IF(JK>=KCTOP(JL)-1) THEN
          LLCUMASK(JL,JK)=.TRUE.
       ENDIF
     ENDIF
   ENDDO
ENDDO
!----------------------------------------------------------------------

DO JN=1,KTRAC 

!*    1.0          DEFINE TRACERS AT HALF LEVELS
!                  -----------------------------

  DO JK=2,KLEV
    IK=JK-1
    DO JL=KIDIA,KFDIA
      ZCEN(JL,JK,JN)=PCEN(JL,JK,JN)
      ZCD(JL,JK,JN) =PCEN(JL,IK,JN)
      ZCU(JL,JK,JN) =PCEN(JL,IK,JN)
      ZMFC(JL,JK,JN)=0.0_JPRB
      ZTENC(JL,JK,JN)=0.0_JPRB
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
      ZCU(JL,KLEV,JN) =PCEN(JL,KLEV,JN)
  ENDDO

!*    2.0          COMPUTE UPDRAFT VALUES
!                  ----------------------

  DO JK=KLEV-1,3,-1
    IK=JK+1
    DO JL=KIDIA,KFDIA
      IF ( LLCUMASK(JL,JK) ) THEN
        ZERATE=PMFU(JL,JK)-PMFU(JL,IK)+PUDRATE(JL,JK)
        ZMFA=1.0_JPRB/MAX(ZRMFCMIN,PMFU(JL,JK))
        IF (JK >=KCTOP(JL) )  THEN
          ZCU(JL,JK,JN)=( PMFU(JL,IK)*ZCU(JL,IK,JN)+ZERATE*PCEN(JL,JK,JN) &
           & -(PUDRATE(JL,JK)+PDMFUP(JL,JK)*PSCAV(JN))*ZCU(JL,IK,JN) )*ZMFA   
! if you have a source term dc/dt=dcdt write 
!             ZCU(JL,JK,JN)=( PMFU(JL,IK)*ZCU(JL,IK,JN)+ZERATE*PCEN(JL,JK,JN) &
!                           -PUDRATE(JL,JK)*ZCU(JL,IK,JN) )*ZMFA 
!                           +dcdt(jl,ik,jn)*ptsphy
        ENDIF
      ENDIF
    ENDDO
  ENDDO


!*    3.0          COMPUTE DOWNDRAFT VALUES
!                  ------------------------

  DO JK=3,KLEV
    IK=JK-1
    DO JL=KIDIA,KFDIA
        IF ( LDDRAF(JL).AND.JK==KDTOP(JL) ) THEN
         !Nota: in order to avoid final negative Tracer values at LFS the allowed value of ZCD
         !      depends on the jump in mass flux at the LFS
            !ZCD(JL,JK,JN)=0.5_JPRB*ZCU(JL,JK,JN)+0.5_JPRB*PCEN(JL,IK,JN)
             ZCD(JL,JK,JN)=0.1_JPRB*ZCU(JL,JK,JN)+0.9_JPRB*PCEN(JL,IK,JN)
        ELSEIF ( LDDRAF(JL).AND.JK>KDTOP(JL) ) THEN
             ZERATE=-PMFD(JL,JK)+PMFD(JL,IK)+PDDRATE(JL,JK)
             ZMFA=1._JPRB/MIN(-ZRMFCMIN,PMFD(JL,JK))
             ZCD(JL,JK,JN)=( PMFD(JL,IK)*ZCD(JL,IK,JN)-ZERATE*PCEN(JL,IK,JN) &
                           &+(PDDRATE(JL,JK)+PDMFDP(JL,JK)*PSCAV(JN))*ZCD(JL,IK,JN) )*ZMFA 
! if you have a source term dc/dt=dcdt write 
!             ZCD(JL,JK,JN)=( PMFD(JL,IK)*ZCD(JL,IK,JN)-ZERATE*PCEN(JL,IK,JN) &
!                           &+PDDRATE(JL,JK)*ZCD(JL,IK,JN) &
!                           &+dcdt(jl,ik,jn)*ptsphy
        ENDIF
    ENDDO
  ENDDO

! In order to avoid negative Tracer at KLEV adjust ZCD
  JK=KLEV
  IK=JK-1
  DO JL=KIDIA,KFDIA
    IF (LDDRAF(JL)) THEN
     ZPOSI=-ZDP(JL,JK)*(PMFU(JL,JK)*ZCU(JL,JK,JN)+PMFD(JL,JK)*ZCD(JL,JK,JN)&
                       &-(PMFU(JL,JK)+PMFD(JL,JK))*PCEN(JL,IK,JN) )
     IF( PCEN(JL,JK,JN)+ZPOSI*PTSPHY<0.0_JPRB ) THEN
        ZMFA=1._JPRB/MIN(-ZRMFCMIN,PMFD(JL,JK))
        ZCD(JL,JK,JN)=( (PMFU(JL,JK)+PMFD(JL,JK))*PCEN(JL,IK,JN)-PMFU(JL,JK)*ZCU(JL,JK,JN)&
                    &+PCEN(JL,JK,JN)/(PTSPHY*ZDP(JL,JK)) )*ZMFA
     ENDIF
    ENDIF
  ENDDO

ENDDO

!----------------------------------------------------------------------
 
  DO JN=1,KTRAC

!*    4.0          COMPUTE FLUXES
!                  --------------

    DO JK=2,KLEV
      IK=JK-1
      DO JL=KIDIA,KFDIA
        IF(LLCUMASK(JL,JK)) THEN
          ZMFA=PMFU(JL,JK)+PMFD(JL,JK)
          ZMFC(JL,JK,JN)=PMFU(JL,JK)*ZCU(JL,JK,JN)+PMFD(JL,JK)*ZCD(JL,JK,JN)&
           & -ZIMP*ZMFA*ZCEN(JL,IK,JN)   
        ENDIF
      ENDDO
    ENDDO

!*    5.0          COMPUTE TENDENCIES = RHS
!                  ------------------------

    DO JK=2,KLEV-1
      IK=JK+1
      DO JL=KIDIA,KFDIA
        IF(LLCUMASK(JL,JK)) THEN
          ZTENC(JL,JK,JN)=ZDP(JL,JK)*(ZMFC(JL,IK,JN)-ZMFC(JL,JK,JN))
         !ZTENC(JL,JK,JN)=ZTENC(JL,JK,JN)-ZDP(JL,JK)*PDMFUP(JL,JK)*PSCAV(JN)*ZCU(JL,JK,JN)
        ENDIF
      ENDDO
    ENDDO

    JK=KLEV
       DO JL=KIDIA,KFDIA
         IF(LDCUM(JL)) THEN
            ZTENC(JL,JK,JN)=-ZDP(JL,JK)*ZMFC(JL,JK,JN)
           !ZTENC(JL,JK,JN)=ZTENC(JL,JK,JN)-ZDP(JL,JK)*PDMFUP(JL,JK)*PSCAV(JN)*ZCU(JL,JK,JN)
         ENDIF
       ENDDO

  ENDDO

  IF ( ZRMFSOLCT==0.0_JPRB ) THEN


!*    6.0          UPDATE TENDENCIES
!                  -----------------

    DO JN=1,KTRAC
      DO JK=2,KLEV
        DO JL=KIDIA,KFDIA
          IF(LLCUMASK(JL,JK)) THEN
            PTENC(JL,JK,JN)=PTENC(JL,JK,JN)+ZTENC(JL,JK,JN)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  ELSE

!---------------------------------------------------------------------------

!*    7.0          IMPLICIT SOLUTION
!                  -----------------
 
   ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
   ! reuse ZMFC=A and ZB=B;
   ! ZTENC corresponds to the RHS ("constants") of the equation
   ! The solution is in ZR1 
 
    ALLOCATE(ZB(KLON,KLEV))
    ALLOCATE(ZR1(KLON,KLEV))   
    ALLOCATE(LLCUMBAS(KLON,KLEV))
    LLCUMBAS(:,:)=.FALSE.
    ZB(:,:)=1._JPRB
 
    DO JN=1,KTRAC

   ! Fill vectors A, B and RHS
 
      DO JK=2,KLEV
        IK=JK+1
        IM=JK-1
        DO JL=KIDIA,KFDIA
        ! LLCUMBAS(JL,JK)=LLCUMASK(JL,JK).AND.JK<=KCBOT(JL)
          LLCUMBAS(JL,JK)=LLCUMASK(JL,JK)
          IF(LLCUMBAS(JL,JK)) THEN
            ZZP=ZRMFSOLCT*ZDP(JL,JK)*PTSPHY
            ZMFC(JL,JK,JN)=-ZZP*(PMFU(JL,JK)+PMFD(JL,JK))
            IF(JK<KLEV) THEN
              ZB(JL,JK)=1.0_JPRB+ZZP*(PMFU(JL,IK)+PMFD(JL,IK))
            ELSE
              ZB(JL,JK)=1.0_JPRB
            ENDIF
            ZZP=RG*(PMFUO(JL,JK)+RMFADVWDD*PMFDO(JL,JK))/(PAP(JL,JK)-PAP(JL,IM))*PTSPHY*ZADVW(JL)
            ZC=ZZP*(PCEN(JL,IM,JN)-PCEN(JL,JK,JN))
            ZTENC(JL,JK,JN) = ZTENC(JL,JK,JN)*PTSPHY+PCEN(JL,JK,JN)-ZC
          ENDIF
        ENDDO
      ENDDO
 
      CALL CUBIDIAG&
       & ( KIDIA, KFDIA, KLON, KLEV,&
       & KCTOP, LLCUMBAS,&
       & ZMFC(:,:,JN),  ZB,   ZTENC(:,:,JN),   ZR1 )  
 
     ! Compute tendencies
 
      DO JK=2,KLEV
        DO JL=KIDIA,KFDIA
          IF(LLCUMBAS(JL,JK)) THEN
            PTENC(JL,JK,JN)=PTENC(JL,JK,JN)+(ZR1(JL,JK)-PCEN(JL,JK,JN))*ZTSPHY
         !  for implicit solution including tendency source term
         !  PTENC(JL,JK,JN)=(ZR1(JL,JK)-PCEN(JL,JK,JN))*ZTSPHY
          ENDIF
        ENDDO
      ENDDO

    ENDDO
 
    DEALLOCATE(LLCUMBAS)
    DEALLOCATE(ZB)
    DEALLOCATE(ZR1)

  ENDIF
!---------------------------------------------------------------------------

DEALLOCATE(LLCUMASK)  
DEALLOCATE(ZDP)
DEALLOCATE(ZMFC)
DEALLOCATE(ZTENC)
DEALLOCATE(ZCD)
DEALLOCATE(ZCU)
DEALLOCATE(ZCEN)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUCTRACER',1,ZHOOK_HANDLE)
END SUBROUTINE CUCTRACER

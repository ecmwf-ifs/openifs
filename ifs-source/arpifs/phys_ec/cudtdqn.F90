! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUDTDQN &
 & (YDEPHLI, YDPHNC,   YDECUMF,  YDEPHY,&
 & KIDIA,    KFDIA,    KLON,     KLEV,&
 & KTOPM2,   KTYPE,    KCTOP,    KCBOT,    KDTOP,    LDCUM,    LDDRAF,   LSCVFLAG, PTSPHY,&
 & PAPH,     PAP  ,    PGEOH,    PGEO,&
 & PTEN,     PTENH,    PQEN,     PQENH,    PQSEN,&
 & PLGLAC,   PLUDE,    PLUDELI,  PSNDE,    PMFU,     PMFD,&
 & PMFUS,    PMFDS,    PMFUQ,    PMFDQ,&
 & PMFUL,    PDMFUP,   PDPMEL,   PMFLXR,   PMFLXS,&
 & PTENT,    PTENQ,    PENTH )

!**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                DOES GLOBAL DIAGNOSTICS

!**   INTERFACE.
!     ----------

!          *CUDTDQ* IS CALLED FROM *CUMASTR*

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
!    *KCTOP*        CLOUD TOP LEVEL
!    *KCBOT*        CLOUD BASE LEVEL
!    *KDTOP*        TOP LEVEL OF DOWNDRAFTS

!    INPUT PARAMETERS (LOGICAL): 

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDDRAF*       FLAG: .TRUE. FOR DOWNDRAFT LEVEL

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PAP *         PROVISIONAL PRESSURE ON FULL LEVELS            PA
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PGEO*         GEOPOTENTIAL ON FULL LEVELS                   M2/S2
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PQSEN*        SATURATION ENV. SPEC. HUMIDITY (T+1)          KG/KG
!    *PLGLAC*       FLUX OF FROZEN CLOUDWATER IN UPDRAFTS         KG/(M2*S) 
!    *PLUDE*        DETRAINED CONDENSATE                          KG/(M2*S)
!    *PLUDELI*      DETRAINED LIQUID, ICE                         KG/(M2*S)
!    *PSNDE*        DETRAINED SNOW/RAIN                           KG/(M2*S)
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS        J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS          KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP.                    KG/(M2*S)
!    *PDPMEL*       CHANGE IN PRECIP.-FLUXES DUE TO MELTING       KG/(M2*S)

!    UPDATED PARAMETERS (REAL):

!    *PTENT*        TEMPERATURE TENDENCY                           K/S
!    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)

!    OUTPUT PARAMETERS (REAL):

!    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S) 

!----------------------------------------------------------------------

!     AUTHOR.
!     -------
!      M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      03-08-28       : Clean up LINPHYS  P.BECHTOLD
!      05-10-13       : implicit solution P.BECHTOLD
!      10-11-17       : allow to substract subbsidence (dyn) P. Bechtold
!      28-12-21       : correction for enthalpy conservation ice phase P. Bechtold
!---------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG       ,RCPD     ,RLVTT    ,RLSTT    ,RLMLT    ,RTT
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                    R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &                    RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
 &                    RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
USE YOECUMF  , ONLY : TECUMF
USE YOEPHY   , ONLY : TEPHY
USE YOEPHLI  , ONLY : TEPHLI
USE YOPHNC   , ONLY : TPHNC

IMPLICIT NONE

TYPE(TECUMF)      ,INTENT(INOUT) :: YDECUMF
TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
TYPE(TPHNC)       ,INTENT(INOUT) :: YDPHNC
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOPM2 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDTOP(KLON) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON)
LOGICAL           ,INTENT(IN)    :: LDDRAF(KLON) 
LOGICAL           ,INTENT(IN)    :: LSCVFLAG(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQEN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQENH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSEN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLGLAC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLUDE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLUDELI(KLON,KLEV,2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSNDE(KLON,KLEV,2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFUS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFDS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFUQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFDQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFUL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDMFUP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDPMEL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFLXR(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFLXS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PENTH(KLON,KLEV) 

LOGICAL :: LLTEST
INTEGER(KIND=JPIM) :: JK, IK, IM, JL
REAL(KIND=JPRB)    :: ZTSPHY, ZIMP, ZORCPD, ZOEALFA, ZTARG,&
                     &ZZP, ZGQ, ZGS, ZGH, ZS, ZQ
REAL(KIND=JPRB)    :: ZALV(KLON,KLEV), ZALV2(KLON,KLEV), ZINT(KLON), ZINT2(KLON)

REAL(KIND=JPRB) :: ZMFUS(KLON,KLEV), ZMFUQ(KLON,KLEV),&
                  &ZMFDS(KLON,KLEV), ZMFDQ(KLON,KLEV), ZADVW(KLON)

REAL(KIND=JPRB),   DIMENSION(:,:), ALLOCATABLE :: ZDTDT, ZDQDT, ZDP
REAL(KIND=JPRB),   DIMENSION(:,:), ALLOCATABLE :: ZB,    ZR1,   ZR2
LOGICAL,           DIMENSION(:,:), ALLOCATABLE :: LLCUMBAS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cubidiag.intfb.h"
#include "fcttre.func.h"

!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUDTDQN',0,ZHOOK_HANDLE)
ASSOCIATE(RMFSOLTQ=>YDECUMF%RMFSOLTQ, RMFSOLRHS=>YDECUMF%RMFSOLRHS, &
 & LPHYLIN=>YDEPHLI%LPHYLIN, RLPTRC=>YDEPHLI%RLPTRC, &
 & LEPCLD=>YDEPHY%LEPCLD, LENCLD2=>YDPHNC%LENCLD2, RMFADVW=>YDECUMF%RMFADVW, &
 & LEPCLD2=>YDPHNC%LEPCLD2, RMFADVWDD=>YDECUMF%RMFADVWDD, &
 & LMFENTHCONS=>YDECUMF%LMFENTHCONS)

!*    1.0          SETUP AND INITIALIZATIONS
!                  -------------------------

ZIMP=1.0_JPRB-RMFSOLTQ
ZTSPHY=1.0_JPRB/PTSPHY
ZORCPD=1.0_JPRB/RCPD
DO JL=KIDIA,KFDIA
  ZADVW(JL)=0.0_JPRB
  IF(KTYPE(JL)==1) ZADVW(JL)=RMFADVW
ENDDO

ALLOCATE(ZDTDT(KLON,KLEV))
ALLOCATE(ZDQDT(KLON,KLEV))
ALLOCATE(ZDP(KLON,KLEV))


DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PENTH(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!         MASS-FLUX APPROACH SWITCHED ON FOR DEEP CONVECTION ONLY
!         IN THE TANGENT-LINEAR AND ADJOINT VERSIONS

DO JL=KIDIA,KFDIA
  IF (KTYPE(JL) /= 1.AND.LPHYLIN) LDCUM(JL)=.FALSE.
ENDDO

! zero detrained liquid water if diagnostic cloud scheme to be used

! this means that detrained liquid water will be evaporated in the 
! cloud environment and not fed directly into a cloud liquid water 
! variable

LLTEST = (.NOT.LEPCLD.AND..NOT.LENCLD2.AND..NOT.LEPCLD2).OR. &
 &       (LPHYLIN.AND..NOT.LENCLD2.AND..NOT.LEPCLD2)

IF (LLTEST) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PLUDE(JL,JK)=0.0_JPRB
      PLUDELI(JL,JK,1)=0.0_JPRB
      PLUDELI(JL,JK,2)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL)) THEN
      ZDP(JL,JK)=RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
      ZMFUS(JL,JK)=PMFUS(JL,JK)
      ZMFDS(JL,JK)=PMFDS(JL,JK)
      ZMFUQ(JL,JK)=PMFUQ(JL,JK)
      ZMFDQ(JL,JK)=PMFDQ(JL,JK)
    ENDIF
  ENDDO
ENDDO

!------------------------------------------------------------------------------

IF ( RMFSOLTQ>0.0_JPRB ) THEN

!*    2.0          RECOMPUTE CONVECTIVE FLUXES IF IMPLICIT

    DO JK=KTOPM2,KLEV
      IK=JK-1
    !DIR$ IVDEP
    !OCL NOVREC
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL).AND.JK>=KCTOP(JL)-1) THEN
        ! compute interpolating coefficients ZGS and ZGQ for half-level values
          ZGQ=(PQENH(JL,JK)-PQEN(JL,IK))/PQSEN(JL,JK)
          ZGH =RCPD*PTEN(JL,JK)+PGEO(JL,JK)
          ZGS=(RCPD*(PTENH(JL,JK)-PTEN(JL,IK))+PGEOH(JL,JK)-PGEO(JL,IK))/ZGH

        !half-level environmental values for S and Q
          ZS =RCPD*(ZIMP*PTEN(JL,IK)+ZGS*PTEN(JL,JK))+PGEO(JL,IK)+ZGS*PGEO(JL,JK)
          ZQ =ZIMP*PQEN(JL,IK)+ZGQ*PQSEN(JL,JK)
          ZMFUS(JL,JK)=PMFUS(JL,JK)-PMFU(JL,JK)*ZS
          ZMFUQ(JL,JK)=PMFUQ(JL,JK)-PMFU(JL,JK)*ZQ
          IF(LDDRAF(JL).AND.JK >= KDTOP(JL)) THEN
            ZMFDS(JL,JK)=PMFDS(JL,JK)-PMFD(JL,JK)*ZS
            ZMFDQ(JL,JK)=PMFDQ(JL,JK)-PMFD(JL,JK)*ZQ
          ENDIF
        ENDIF
      ENDDO
    ENDDO

ENDIF

!*    3.0          COMPUTE TENDENCIES
!                  ------------------

ZINT=0._JPRB
DO JK=KTOPM2,KLEV
 DO JL=KIDIA,KFDIA
   ZALV(JL,JK)=FOELHMCU(PTEN(JL,JK))
   ZALV2(JL,JK)=ZALV(JL,JK)
 ENDDO
ENDDO
IF(LMFENTHCONS) THEN
! Integrals for conservation correction ice-phase
  DO JK=KTOPM2,KLEV-1
   DO JL=KIDIA,KFDIA
     IF(LSCVFLAG(JL)) ZALV(JL,JK)=RLVTT
     ZINT(JL)=ZINT(JL)+RLMLT*(PLGLAC(JL,JK)-PDPMEL(JL,JK))-ZALV(JL,JK)*&
             &(PMFUL(JL,JK+1)-PMFUL(JL,JK))+ZALV2(JL,JK)*PDMFUP(JL,JK)
   ENDDO
  ENDDO
  JK=KLEV
  DO JL=KIDIA,KFDIA
   IF(LSCVFLAG(JL)) ZALV(JL,JK)=RLVTT
   ZINT(JL)=ZINT(JL)+RLMLT*(PLGLAC(JL,JK)-PDPMEL(JL,JK))+&
           ZALV(JL,JK)*PMFUL(JL,JK)+ZALV2(JL,JK)*PDMFUP(JL,JK)
   ZINT2(JL)=ZINT(JL)-RLVTT*PMFLXR(JL,KLEV+1)-RLSTT*PMFLXS(JL,KLEV+1)
   IF(ABS(ZINT(JL))>0.1_JPRB) THEN
     ZINT(JL)=ZINT2(JL)/(ZINT(JL)+1.E-3_JPRB)
   ELSE
     ZINT(JL)=0.0_JPRB
   ENDIF
   IF(ABS(ZINT(JL))>2.0_JPRB) ZINT(JL)=0.0_JPRB
  ENDDO
ENDIF

DO JK=KTOPM2,KLEV

  IF(JK < KLEV) THEN
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        IF (LPHYLIN) THEN
          ZTARG=PTEN(JL,JK)
          ZOEALFA=MIN(1.0_JPRB,0.545_JPRB*(TANH(0.17_JPRB*(ZTARG-RLPTRC))+1.0_JPRB))
          ZALV(JL,JK)=ZOEALFA*RLVTT+(1.0_JPRB-ZOEALFA)*RLSTT
          ZALV2(JL,JK)=ZALV(JL,JK)
        ENDIF
        ZDTDT(JL,JK)=ZDP(JL,JK)*ZORCPD*&
         & (ZMFUS(JL,JK+1)-ZMFUS(JL,JK)+&
         & ZMFDS(JL,JK+1)-ZMFDS(JL,JK)&
         & +(RLMLT*(PLGLAC(JL,JK)-PDPMEL(JL,JK))&
         & -ZALV(JL,JK)*(PMFUL(JL,JK+1)-PMFUL(JL,JK))+ZALV2(JL,JK)*PDMFUP(JL,JK))*(1.0_JPRB-ZINT(JL))&
         & +RLVTT*(PLUDELI(JL,JK,1)+PSNDE(JL,JK,1)) &
         & +RLSTT*(PLUDELI(JL,JK,2)+PSNDE(JL,JK,2)) )

        ZDQDT(JL,JK)=ZDP(JL,JK)*&
         & (ZMFUQ(JL,JK+1)-ZMFUQ(JL,JK)+&
         & ZMFDQ(JL,JK+1)-ZMFDQ(JL,JK)+&
         & PMFUL(JL,JK+1)-PMFUL(JL,JK)-&
         & PLUDE(JL,JK)-PDMFUP(JL,JK)-PSNDE(JL,JK,1)-PSNDE(JL,JK,2))  
      ENDIF
    ENDDO

  ELSE
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        IF (LPHYLIN) THEN
          ZTARG=PTEN(JL,JK)
          ZOEALFA=MIN(1.0_JPRB,0.545_JPRB*(TANH(0.17_JPRB*(ZTARG-RLPTRC))+1.0_JPRB))
          ZALV(JL,JK)=ZOEALFA*RLVTT+(1.0_JPRB-ZOEALFA)*RLSTT
          ZALV2(JL,JK)=ZALV(JL,JK)
        ENDIF
        ZDTDT(JL,JK)=-ZDP(JL,JK)*ZORCPD*&
         & (ZMFUS(JL,JK)+ZMFDS(JL,JK)+(RLMLT*(PDPMEL(JL,JK)-PLGLAC(JL,JK))&
         &-ZALV(JL,JK)*PMFUL(JL,JK)-ZALV2(JL,JK)*PDMFUP(JL,JK))*(1.0_JPRB-ZINT(JL)) )

        ZDQDT(JL,JK)=-ZDP(JL,JK)*&
         & (ZMFUQ(JL,JK)+ZMFDQ(JL,JK)+&
         & (PMFUL(JL,JK)+PDMFUP(JL,JK)))  
      ENDIF
    ENDDO
  ENDIF

ENDDO

IF ( RMFSOLTQ==0.0_JPRB ) THEN

!*    3.1          UPDATE TENDENCIES
!                  -----------------

    DO JK=KTOPM2,KLEV
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL)) THEN
          PTENT(JL,JK)=PTENT(JL,JK)+ZDTDT(JL,JK)
          PTENQ(JL,JK)=PTENQ(JL,JK)+ZDQDT(JL,JK)
          PENTH(JL,JK)=ZDTDT(JL,JK)*RCPD
        ENDIF
      ENDDO
    ENDDO

ELSE
!----------------------------------------------------------------------

!*    3.2          IMPLICIT SOLUTION
!                  -----------------

     ! Fill bi-diagonal Matrix vectors A=k-1, B=k, C=k+1;
     ! reuse ZMFUS=A
     ! ZDTDT and ZDQDT correspond to the RHS ("constants") of the equation
     ! The solution is in ZR1 and ZR2

      ALLOCATE(ZB(KLON,KLEV))
      ALLOCATE(ZR1(KLON,KLEV))
      ALLOCATE(ZR2(KLON,KLEV))
      ALLOCATE(LLCUMBAS(KLON,KLEV))
      LLCUMBAS(:,:)=.FALSE.
      ZB(:,:)=1.0_JPRB
      ZMFUS(:,:)=0.0_JPRB

     ! Fill vectors A, B and RHS

      DO JK=KTOPM2,KLEV
         IK=JK+1
         IM=JK-1
         DO JL=KIDIA,KFDIA
           LLCUMBAS(JL,JK)=LDCUM(JL).AND.JK>=KCTOP(JL)-1
           IF(LLCUMBAS(JL,JK)) THEN
             ZZP=RMFSOLTQ*ZDP(JL,JK)*PTSPHY
             ZMFUS(JL,JK)=-ZZP*(PMFU(JL,JK)+PMFD(JL,JK))
             IF(JK<KLEV) THEN
               ZB(JL,JK)=1.0_JPRB+ZZP*(PMFU(JL,IK)+PMFD(JL,IK))
             ELSE
               ZB(JL,JK)=1.0_JPRB
             ENDIF
             ZZP=RG*(PMFU(JL,JK)+RMFADVWDD*PMFD(JL,JK))/(PAP(JL,JK)-PAP(JL,IM))*PTSPHY*ZADVW(JL)
             ZS=ZZP*(PTEN(JL,IM)-PTEN(JL,JK)+ZORCPD*(PGEO(JL,IM)-PGEO(JL,JK)) )
             ZQ=ZZP*(PQEN(JL,IM)-PQEN(JL,JK))
             ZDTDT(JL,JK) = (ZDTDT(JL,JK)+PTENT(JL,JK)*RMFSOLRHS)*PTSPHY+PTEN(JL,JK)-ZS
             ZDQDT(JL,JK) = (ZDQDT(JL,JK)+PTENQ(JL,JK)*RMFSOLRHS)*PTSPHY+PQEN(JL,JK)-ZQ
           ENDIF
         ENDDO
      ENDDO

      CALL CUBIDIAG&
         &( KIDIA,    KFDIA,   KLON,   KLEV, &
         &  KCTOP,    LLCUMBAS, &
         &  ZMFUS,    ZB,     ZDTDT,   ZR1 )

      CALL CUBIDIAG&
         &( KIDIA,    KFDIA,   KLON,   KLEV, &
         &  KCTOP,    LLCUMBAS, &
         &  ZMFUS,    ZB,     ZDQDT,   ZR2 )

     ! Compute tendencies

      DO JK=KTOPM2,KLEV
         DO JL=KIDIA,KFDIA
           IF(LLCUMBAS(JL,JK)) THEN
             PTENT(JL,JK)=PTENT(JL,JK)*(1.0_JPRB-RMFSOLRHS)+(ZR1(JL,JK)-PTEN(JL,JK))*ZTSPHY
             PTENQ(JL,JK)=PTENQ(JL,JK)*(1.0_JPRB-RMFSOLRHS)+(ZR2(JL,JK)-PQEN(JL,JK))*ZTSPHY
             PENTH(JL,JK)=(ZR1(JL,JK)-PTEN(JL,JK))*ZTSPHY
           ENDIF
         ENDDO
      ENDDO

      DEALLOCATE(LLCUMBAS)
      DEALLOCATE(ZR2)
      DEALLOCATE(ZR1)
      DEALLOCATE(ZB)

!----------------------------------------------------------------------
ENDIF

DEALLOCATE(ZDP)
DEALLOCATE(ZDQDT)
DEALLOCATE(ZDTDT)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUDTDQN',1,ZHOOK_HANDLE)
END SUBROUTINE CUDTDQN

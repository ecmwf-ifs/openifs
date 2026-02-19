! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUDTDQN2 &
 & (YDML_PHY_SLIN,     YDEPHY,&
 & KIDIA,    KFDIA,    KLON,     KLEV,&
 & KTOPM2,   KCTOP,    KDTOP,&
 & LDRAIN1D, LDCUM,    LDDRAF,   PTSPHY,&
 & PAPH,     PGEOH,    PGEO,&
 & PTEN,     PTENH,    PQEN,     PQENH,    PQSEN,&
 & PLGLAC,   PLUDE,    PMFU,     PMFD,&
 & PMFUS,    PMFDS,    PMFUQ,    PMFDQ,&
 & PMFUL,    PDMFUP,   PDMFDP,   PDPMEL,&
 & PTENT,    PTENQ,    PENTH )  

!**** *CUDTDQ2* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                 DOES GLOBAL DIAGNOSTICS (1dvar)

!**   INTERFACE.
!     ----------

!          *CUDTDQ2* IS CALLED FROM *CUMASTR*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KCTOP*        CLOUD TOP LEVEL
!    *KDTOP*        TOP LEVEL OF DOWNDRAFTS

!    INPUT PARAMETERS (LOGICAL): 

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDDRAF*       FLAG: .TRUE. FOR DOWNDRAFT LEVEL

!    INPUT PARAMETERS (REAL):

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PGEO*         GEOPOTENTIAL ON FULL LEVELS                   M2/S2
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PQSEN*        SATURATION ENV. SPEC. HUMIDITY (T+1)          KG/KG
!    *PLGLAC*       FLUX OF FROZEN CLOUDWATER IN UPDRAFTS         KG/(M2*S) 
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M3*S)
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS        J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS          KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS        KG/(M2*S)
!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS      KG/(M2*S)
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
!      03-10-21    : Duplication of CUDTDQN for simplified convection scheme    P. LOPEZ
!      11-01-2007  : Modified use of LDRAIN switch (1D-Var rain)   P. LOPEZ
!      22-10-2007  : Implicit version added   P. LOPEZ
!      01-10-2008  : LDRAIN1D name change to reflect usage  A. GEER
!      P. Lopez     10-05-2016  Revised tendencies to match full scheme
!----------------------------------------------------------------------

USE MODEL_PHYSICS_SIMPLINEAR_MOD , ONLY : MODEL_PHYSICS_SIMPLINEAR_TYPE
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG       ,RCPD     ,RLVTT    ,RLSTT    ,RLMLT    ,RTT
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &                    R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &                    RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU  ,&
 &                    RTWAT_RTICE_R      ,RTWAT_RTICECU_R  
USE YOEPHY   , ONLY : TEPHY

IMPLICIT NONE

TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_PHYSICS_SIMPLINEAR_TYPE),INTENT(INOUT):: YDML_PHY_SLIN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOPM2 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDTOP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDRAIN1D
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON) 
LOGICAL           ,INTENT(IN)    :: LDDRAF(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQEN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQENH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSEN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLGLAC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLUDE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFUS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFDS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFUQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFDQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFUL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDMFUP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDMFDP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDPMEL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PENTH(KLON,KLEV) 

INTEGER(KIND=JPIM) :: JK, JL, IK

REAL(KIND=JPRB) :: ZMFUS(KLON,KLEV), ZMFUQ(KLON,KLEV),&
                 & ZMFDS(KLON,KLEV), ZMFDQ(KLON,KLEV)

REAL(KIND=JPRB) :: ZALV, ZOEALFA, ZTARG
REAL(KIND=JPRB) :: ZGLMDP, ZGLVDP, ZRCPD, ZGDPCP
REAL(KIND=JPRB) :: ZFACT2, ZFACT3
REAL(KIND=JPRB) :: ZTSPHY, ZIMP, ZZP, ZGQ, ZGS, ZGH, ZS, ZQ, ZDELP

LOGICAL :: LLTEST

REAL(KIND=JPRB),   DIMENSION(:,:), ALLOCATABLE :: ZDTDT, ZDQDT, ZGDP
REAL(KIND=JPRB),   DIMENSION(:,:), ALLOCATABLE :: ZB,    ZR1,   ZR2
LOGICAL,           DIMENSION(:,:), ALLOCATABLE :: LLCUMBAS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "cubidiag.intfb.h"
#include "fcttre.func.h"

!----------------------------------------------------------------------

!*    1.0          INCREMENTATION OF T AND Q TENDENCIES
!                  ------------------------------------

IF (LHOOK) CALL DR_HOOK('CUDTDQN2',0,ZHOOK_HANDLE)
ASSOCIATE(RMFSOLTQ2=>YDML_PHY_SLIN%YRECUMF2%RMFSOLTQ2, RMFSOLRHS2=>YDML_PHY_SLIN%YRECUMF2%RMFSOLRHS2, &
 & LPHYLIN=>YDML_PHY_SLIN%YREPHLI%LPHYLIN, RLPTRC=>YDML_PHY_SLIN%YREPHLI%RLPTRC, &
 & LEPCLD=>YDEPHY%LEPCLD, &
 & LENCLD2=>YDML_PHY_SLIN%YRPHNC%LENCLD2, &
 & LEPCLD2=>YDML_PHY_SLIN%YRPHNC%LEPCLD2)
ZIMP=1.0_JPRB-RMFSOLTQ2
ZTSPHY=1.0_JPRB/PTSPHY
ZRCPD=1.0_JPRB/RCPD

ALLOCATE(ZDTDT(KLON,KLEV))
ALLOCATE(ZDQDT(KLON,KLEV))
ALLOCATE(ZGDP(KLON,KLEV))

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PENTH(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

! zero detrained liquid water if diagnostic cloud scheme to be used

! this means that detrained liquid water will be evaporated in the 
! cloud environment and not fed directly into a cloud liquid water 
! variable

LLTEST = ((.NOT.LEPCLD.AND..NOT.LENCLD2.AND..NOT.LEPCLD2).OR. &
       &  (LPHYLIN.AND..NOT.LENCLD2.AND..NOT.LEPCLD2)).AND..NOT.LDRAIN1D

IF (LLTEST) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PLUDE(JL,JK)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF(LDCUM(JL)) THEN
      ZDELP=1.0_JPRB/(PAPH(JL,JK+1)-PAPH(JL,JK))
      ZGDP(JL,JK)= RG*ZDELP
      ZMFUS(JL,JK)=PMFUS(JL,JK)
      ZMFDS(JL,JK)=PMFDS(JL,JK)
      ZMFUQ(JL,JK)=PMFUQ(JL,JK)
      ZMFDQ(JL,JK)=PMFDQ(JL,JK)
    ENDIF
  ENDDO
ENDDO

!------------------------------------------------------------------------------

IF (RMFSOLTQ2 > 0.0_JPRB) THEN

!*    2.0          RECOMPUTE CONVECTIVE FLUXES IF IMPLICIT

    DO JK=KTOPM2,KLEV
      IK=JK-1
!DIR$ IVDEP
!OCL NOVREC
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL).AND.JK>=KCTOP(JL)-1) THEN
! Compute interpolating coefficients ZGS and ZGQ for half-level values
          ZFACT2=1.0_JPRB/PQSEN(JL,JK)
          ZGQ=(PQENH(JL,JK)-PQEN(JL,IK))*ZFACT2
          ZGH =RCPD*PTEN(JL,JK)+PGEO(JL,JK)
          ZFACT3=1.0_JPRB/ZGH
          ZGS=(RCPD*(PTENH(JL,JK)-PTEN(JL,IK))+PGEOH(JL,JK)-PGEO(JL,IK))*ZFACT3

! Half-level environmental values for S and Q
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

DO JK=KTOPM2,KLEV

  IF(JK < KLEV) THEN
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        IF (LPHYLIN .OR. LDRAIN1D) THEN
          ZTARG=PTEN(JL,JK)
          ZOEALFA=MIN(1.0_JPRB,0.545_JPRB*(TANH(0.17_JPRB*(ZTARG-RLPTRC))+1.0_JPRB))
          ZALV=ZOEALFA*RLVTT+(1.0_JPRB-ZOEALFA)*RLSTT
        ELSE
          ZALV=FOELHMCU(PTEN(JL,JK))
        ENDIF
        ZGDPCP=ZGDP(JL,JK)*ZRCPD
        ZGLMDP=ZGDPCP*RLMLT
        ZGLVDP=ZGDPCP*ZALV 

        ZDTDT(JL,JK)=ZGDPCP*&
         & (ZMFUS(JL,JK+1)-ZMFUS(JL,JK)&
         & +ZMFDS(JL,JK+1)-ZMFDS(JL,JK))&
         & +ZGLMDP*(PLGLAC(JL,JK)-PDPMEL(JL,JK))&
         & -ZGLVDP*(PMFUL(JL,JK+1)-PMFUL(JL,JK)&
         & -PLUDE(JL,JK)&
         & -(PDMFUP(JL,JK)+PDMFDP(JL,JK)))  

        ZDQDT(JL,JK)=ZGDP(JL,JK)*&
         & (ZMFUQ(JL,JK+1)-ZMFUQ(JL,JK)&
         & +ZMFDQ(JL,JK+1)-ZMFDQ(JL,JK)&
         & +PMFUL(JL,JK+1)-PMFUL(JL,JK)&
         & -PLUDE(JL,JK)&
         & -(PDMFUP(JL,JK)+PDMFDP(JL,JK))) 
      ENDIF
    ENDDO

  ELSE
    DO JL=KIDIA,KFDIA
      IF(LDCUM(JL)) THEN
        IF (LPHYLIN .OR. LDRAIN1D) THEN
          ZTARG=PTEN(JL,JK)
          ZOEALFA=MIN(1.0_JPRB,0.545_JPRB*(TANH(0.17_JPRB*(ZTARG-RLPTRC))+1.0_JPRB))
          ZALV=ZOEALFA*RLVTT+(1.0_JPRB-ZOEALFA)*RLSTT
        ELSE
          ZALV=FOELHMCU(PTEN(JL,JK))
        ENDIF
        ZGDPCP=ZGDP(JL,JK)*ZRCPD
        ZGLMDP=ZGDPCP*RLMLT
        ZGLVDP=ZGDPCP*ZALV 

        ZDTDT(JL,JK)=-ZGDPCP*(ZMFUS(JL,JK)+ZMFDS(JL,JK))&
           &  -ZGLMDP*PDPMEL(JL,JK)&
           &  +ZGLVDP*(PMFUL(JL,JK)+PDMFUP(JL,JK)+PDMFDP(JL,JK))  

        ZDQDT(JL,JK)=-ZGDP(JL,JK)&
         & *(ZMFUQ(JL,JK)+ZMFDQ(JL,JK)&
         & +(PMFUL(JL,JK)+PDMFUP(JL,JK)+PDMFDP(JL,JK)))  
      ENDIF
    ENDDO
  ENDIF

ENDDO

IF (RMFSOLTQ2 == 0.0_JPRB) THEN

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

! Fill bi-diagonal Matrix vectors A=k-1, B=k, C=k+1; reuse ZMFUS=A
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
      DO JL=KIDIA,KFDIA
        LLCUMBAS(JL,JK)=LDCUM(JL) .AND. JK >= KCTOP(JL)-1
        IF(LLCUMBAS(JL,JK)) THEN
          ZZP=RMFSOLTQ2*ZGDP(JL,JK)*PTSPHY
          ZMFUS(JL,JK)=-ZZP*(PMFU(JL,JK)+PMFD(JL,JK))
          ZDTDT(JL,JK) = (ZDTDT(JL,JK)+PTENT(JL,JK)*RMFSOLRHS2)*PTSPHY+PTEN(JL,JK)
          ZDQDT(JL,JK) = (ZDQDT(JL,JK)+PTENQ(JL,JK)*RMFSOLRHS2)*PTSPHY+PQEN(JL,JK)
          IF(JK < KLEV) THEN
            ZB(JL,JK)=1.0_JPRB+ZZP*(PMFU(JL,IK)+PMFD(JL,IK))
          ELSE
            ZB(JL,JK)=1.0_JPRB
          ENDIF
        ENDIF
      ENDDO
   ENDDO

   CALL CUBIDIAG &
      &( KIDIA,    KFDIA,   KLON,   KLEV, &
      &  KCTOP,    LLCUMBAS, &
      &  ZMFUS,    ZB,     ZDTDT,   ZR1 )

   CALL CUBIDIAG &
      &( KIDIA,    KFDIA,   KLON,   KLEV, &
      &  KCTOP,    LLCUMBAS, &
      &  ZMFUS,    ZB,     ZDQDT,   ZR2 )

! Compute final tendencies

   DO JK=KTOPM2,KLEV
      DO JL=KIDIA,KFDIA
        IF(LLCUMBAS(JL,JK)) THEN
          PTENT(JL,JK)=PTENT(JL,JK)*(1.0_JPRB-RMFSOLRHS2)+(ZR1(JL,JK)-PTEN(JL,JK))*ZTSPHY
          PTENQ(JL,JK)=PTENQ(JL,JK)*(1.0_JPRB-RMFSOLRHS2)+(ZR2(JL,JK)-PQEN(JL,JK))*ZTSPHY
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

DEALLOCATE(ZGDP)
DEALLOCATE(ZDQDT)
DEALLOCATE(ZDTDT)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUDTDQN2',1,ZHOOK_HANDLE)
END SUBROUTINE CUDTDQN2

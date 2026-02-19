! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CUANCAPE2(YDECUMF,YDEPHLI, KIDIA,  KFDIA,  KLON,  KLEV,&
                   &PAP,    PAPH,   PT,    PQ,   PCAPE,   PCIN,  PPDEPL)

!***** CUANCAPE2 - COMPUTE APPROXIMATE CAPE,CIN  USING THETAE AND THETAES

!     E. HOLM + P. BECHTOLD     E.C.M.W.F.     13/10/2005

!     PURPOSE 
!     -------
!                 ESTIMATE CAPE FIRST FOR A MIXED-LAYER PARCEL, THEN
!                 LOOP OVER SUBSEQUENT DEPARTURE LAYERS IN LOWEST 350 hPa
!                 Theta_e =Theta*exp[L q_v/(C_p T)] 
!                         = T*(P0/P)**(R_d/C_p) * exp[L q_v/(C_p T)]
!                 -> THIS WILL BE THE UPDRAUGHT PARCEL (CONSERVING ITS
!                 PROPERTIES)  (no entrainment)
!                 CAPE    = Int ( g dTheta_v/Theta_v dz ) = 
!                   aprox = Int ( g (Theta_e_up-Theta_e_sat)/Theta_e_sat ) dz
!                 WITH THIS FORMULATION THE ACTUAL CAPE IS OVERESTIMATED  BY
!                 ROUGHLY 20%. DEEP CONVECTION CAN BE CONSIDERED FOR CAPE
!                 VALUES ABOVE 200-500 J/KG            


!     PARAMETER     DESCRIPTION                              UNITS
!     ---------     -----------                              -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (REAL):

!    *PAP*          PRESSURE ON FULL LEVELS                    PA
!    *PAPH*         PRESSURE ON HALF LEVELS                    PA
!    *PT*           TEMPERATURE ON FULL LEVELS                 K   
!    *PQ*           SPECIFIC HUMIDITY ON FULL LEVELS          KG/KG

!    OUTPUT PARAMETERS (REAL):

!    *PCAPE*        Different definitions of CAPE             J/KG
!         0=    max unstable using theta_e,theta_es approximationa s in ERA5
!         using theta_v:
!         1=for max unstable parcel
!         2=for 50  hPa low-level mixed layer parcel 
!         3=for 100 hPa low-level mixed layer parcel 
!    *PCIN *        CONVECTIVE INHIBITION  AS FOR PCAPE 1-3    J/KG
!    *PPDEPL*       DEPARTURE LEVEL (PA) OF MOST UNSTABLE PARCEL


!          MODIFICATIONS
!          -------------
!     24-06-2011 P. Bechtold :  Add CIN 
!     01-03-2019 P. Bechtold and I. Tsonevski:  Do computations correctly above LCL
!                                    enable accurate computations using Theta_v
!     28-12-2020 P. Bechtold and I. Tsonevski: Introduce max unstable and
!                                   50 and 100 hPa mixed layer CAPE/CIN using Theta_v

!-------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                    R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
 &                    RTWAT_RTICE_R, RTWAT_RTICECU_R
USE YOMCST   , ONLY : RETV, RLVTT, RLSTT, RTT, RD, RKAPPA, RATM, RESTT, RCPD, RV
USE YOECUMF  , ONLY : TECUMF
USE YOEPHLI  , ONLY : TEPHLI

IMPLICIT NONE

TYPE(TECUMF)      ,INTENT(IN)    :: YDECUMF
TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(KLON,0:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCIN(KLON,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPDEPL(KLON)

INTEGER(KIND=JPIM), PARAMETER    :: NL=2 ! additional types of mixed layers (50 and 100 hPa)
INTEGER(KIND=JPIM) :: JL, JK, JKK, JDEP(KLON,0:NL), JN, JNP

REAL(KIND=JPRB), DIMENSION(KLON) :: ZTMIX, ZTHVMIX, ZTHETEU
REAL(KIND=JPRB)                  :: ZCAPE(KLON,KLEV), ZCIN(KLON,KLEV), ZCIN2(KLON),&
                                   &ZTHMIX(KLON,KLEV), ZQMIX(KLON,KLEV), ZPMIXA(KLON),&
                                   &ZPLCL(KLON,KLEV), ZEXN(KLON,KLEV), ZPMIX(KLON,KLEV), ZPM(NL) 
REAL(KIND=JPRB) :: ZDP, ZTHETES, ZTH, ZDZ, ZTEMP, ZTVEMP, ZRPAP, ZORKAPPA,&
                  &ZTD, ZTLCL,  ZCST, ZP

REAL(KIND=JPRB) :: ZTHETAD(KLON,KLEV), ZTHV(KLON,KLEV), ZZDZ(KLON,KLEV), ZTU(KLON,KLEV), ZQS(KLON,KLEV)
LOGICAL         :: LLCAPE_REC=.TRUE.,LLZ(KLON), LLPAP(KLON,KLEV), LLDEP(KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcttre.func.h"
#include "cuadjtq.intfb.h"
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUANCAPE2',0,ZHOOK_HANDLE)
ASSOCIATE(NJKT1=>YDECUMF%NJKT1, NJKT2=>YDECUMF%NJKT2, NJKT4=>YDECUMF%NJKT4, RMINCIN=>YDECUMF%RMINCIN)

ZORKAPPA=1.0_JPRB/RKAPPA
ZCST=1.0_JPRB/(RESTT*RD/RV)
ZPM(1)=50.E2_JPRB
ZPM(2)=100.E2_JPRB

DO JL=KIDIA,KFDIA
  PCAPE(JL,:)=0.0_JPRB
  PCIN(JL,:)=RMINCIN
  PPDEPL(JL)=0.0_JPRB
  ZEXN(JL,KLEV)=(RATM/PAP(JL,KLEV))**RKAPPA  
ENDDO
  DO JK=KLEV-1,NJKT2,-1
     DO JL=KIDIA,KFDIA
       LLPAP(JL,JK)=( PAP(JL,JK)>80.E2_JPRB )
       IF(LLPAP(JL,JK)) THEN
          ZRPAP=1.0_JPRB/PAP(JL,JK)
          ZQS(JL,JK) = FOEEWM(PT(JL,JK))*ZRPAP
          ZQS(JL,JK) = MAX(1.E-8_JPRB,ZQS(JL,JK))
          ZQS(JL,JK) = ZQS(JL,JK)/(1.0_JPRB-RETV*ZQS(JL,JK)) ! small correction
          ZEXN(JL,JK)=(RATM*ZRPAP)**RKAPPA
          ZTH = PT(JL,JK)*ZEXN(JL,JK)
          ZTHETES=ZTH*EXP( FOELDCP(PT(JL,JK))*ZQS(JL,JK)/PT(JL,JK) ) 
          ZTHETAD(JL,JK)=1.0_JPRB/ZTHETES
          ZTHV(JL,JK)=1.0_JPRB/(ZTH*(1.0_JPRB+RETV*PQ(JL,JK)))
       ENDIF
     ENDDO
  ENDDO

  DO JK=NJKT2,KLEV-1
    DO JL=KIDIA,KFDIA
      ZRPAP=1.0_JPRB/PAP(JL,JK)
      ZZDZ(JL,JK)=(PAPH(JL,JK+1)-PAPH(JL,JK))*ZRPAP*RD*PT(JL,JK)*&
             &(1.0_JPRB+RETV*PQ(JL,JK))
    ENDDO
  ENDDO

! double vertical loop to compute departure level that produces maximum CAPE
! using simplified Theta_e algorithm , ie no stauration adjustment

DO JKK=KLEV-1,NJKT1,-1

  DO JL=KIDIA,KFDIA
    ZCAPE(JL,JKK)=0.0_JPRB
    ZCIN(JL,JKK) =0.0_JPRB
    IF (PAPH(JL,KLEV+1)-PAPH(JL,JKK-1)<60.E2_JPRB) THEN
      ZTMIX(JL)=0.0_JPRB
      ZTHMIX(JL,JKK)=0.0_JPRB
      ZQMIX(JL,JKK)=0.0_JPRB
      ZPMIX(JL,JKK)=0.0_JPRB
      DO JK=JKK+1,JKK-1,-1
        IF(ZPMIX(JL,JKK)<=30.E2_JPRB) THEN
          ZDP=PAPH(JL,JK+1)-PAPH(JL,JK)
          ZPMIX(JL,JKK)=ZPMIX(JL,JKK)+ZDP
          ZTHMIX(JL,JKK)=ZTHMIX(JL,JKK)+PT(JL,JK)*ZDP*ZEXN(JL,JK)
          ZQMIX(JL,JKK)=ZQMIX(JL,JKK)+PQ(JL,JK)*ZDP
        ENDIF
      ENDDO
      ZP=1.0_JPRB/ZPMIX(JL,JKK)
      ZQMIX(JL,JKK)=ZQMIX(JL,JKK)*ZP
      ZPMIX(JL,JKK)=PAPH(JL,JKK+2)-0.5_JPRB*ZPMIX(JL,JKK)
      ZTHMIX(JL,JKK)=ZTHMIX(JL,JKK)*ZP
      ZTMIX(JL)=ZTHMIX(JL,JKK)*(ZPMIX(JL,JKK)/RATM)**RKAPPA
    ELSE
      ZQMIX(JL,JKK)=PQ(JL,JKK)
      ZPMIX(JL,JKK)=PAP(JL,JKK)
      ZTMIX(JL)=PT(JL,JKK)
      ZTHMIX(JL,JKK)=PT(JL,JKK)*(RATM/ZPMIX(JL,JKK))**RKAPPA
    ENDIF
    ZTHETEU(JL)=ZTHMIX(JL,JKK)*&
               &EXP( FOELDCP(ZTMIX(JL))*ZQMIX(JL,JKK)/ZTMIX(JL) )
    ZTHVMIX(JL)=ZTHMIX(JL,JKK)*(1.0_JPRB+RETV*ZQMIX(JL,JKK))
    LLZ(JL)=(PAPH(JL,KLEV+1)-PAPH(JL,JKK))<350.E2_JPRB

    ! dewpoint temperature
    ZTEMP = LOG(MAX(1.E-4_JPRB,ZQMIX(JL,JKK))*ZPMIX(JL,JKK)*ZCST)
    ZTD = (R4LES*ZTEMP-R3LES*RTT)/(ZTEMP-R3LES)
    ! adiabatic saturation temperature
    ZTLCL  = ZTD - ( .212_JPRB + 1.571E-3_JPRB * ( ZTD - RTT )      &
               - 4.36E-4_JPRB* ( ZTMIX(JL) - RTT ) ) * ( ZTMIX(JL) - ZTD )
    ZTLCL  = MAX(160.0_JPRB,MIN( ZTLCL, ZTMIX(JL) ))
    ZPLCL(JL,JKK)  = RATM * ( ZTLCL / ZTHMIX(JL,JKK) ) ** ZORKAPPA
  ENDDO

  DO JK=JKK,NJKT2,-1
     DO JL=KIDIA,KFDIA
       IF(LLPAP(JL,JK) .AND. LLZ(JL)) THEN
          ZTEMP=ZTHETEU(JL)*ZTHETAD(JL,JK)-1.0_JPRB
          ZTVEMP=ZTHVMIX(JL)*ZTHV(JL,JK)-1.0_JPRB
          ZDZ=ZZDZ(JL,JK)
          IF(PAP(JL,JK)<=ZPLCL(JL,JKK).AND.ZTEMP>0.0_JPRB) THEN
            ZCAPE(JL,JKK)=ZCAPE(JL,JKK)+ZTEMP*ZDZ
          ENDIF
       ENDIF
     ENDDO
  ENDDO
  
ENDDO

      ! chose maximum CAPE and CIN values
JDEP(:,0)=KLEV-1
DO JK=KLEV-1,NJKT1,-1
  DO JL=KIDIA,KFDIA
    IF(ZCAPE(JL,JK)>PCAPE(JL,0)) THEN
      PCAPE(JL,0)=ZCAPE(JL,JK)
      JDEP(JL,0)=JK
    END IF
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  JKK=JDEP(JL,0)
  PPDEPL(JL)=PAP(JL,JKK)
ENDDO

!----------------------------------------------------------------------------------------------------------
! recompute most unstable CAPE/CIN using Tv (saturation adjustment) and mixed-layer 50 and 100 hPa  CAPE/CIN 
! as a small technical trick store mixed-layer 50 (100) on model level 1 (2), noting that the MUCAPE values are
! valid=stored on model level jk=JDEP(jl,0)
IF(LLCAPE_REC) THEN
  IF (NL>=1) THEN
      JK=KLEV
   ! compute mixed layer values for lowest 50 and 100 hPa layer
      DO JL=KIDIA,KFDIA
        ZDP=PAPH(JL,JK+1)-PAPH(JL,JK) 
        ZPMIX(JL,1)=ZDP
        ZPMIXA(JL)=ZDP
        ZTHMIX(JL,1)=PT(JL,JK)*ZDP*ZEXN(JL,JK)
        ZQMIX(JL,1)=PQ(JL,JK)*ZDP
      ENDDO
      DO JN=2,NL
        DO JL=KIDIA,KFDIA
          ZPMIX(JL,JN)=ZPMIX(JL,1)
          ZTHMIX(JL,JN)=ZTHMIX(JL,1)
          ZQMIX(JL,JN)=ZQMIX(JL,1)
        ENDDO
      ENDDO
      JNP=NL
      JN=1
      DO JK=KLEV-1,NJKT4,-1
        DO JL=KIDIA,KFDIA
            ZDP=PAPH(JL,JK+1)-PAPH(JL,JK)
            ZPMIXA(JL)=ZPMIXA(JL)+ZDP
            IF( ZPMIXA(JL)<=ZPM(JN)) THEN
              ZPMIX(JL,JN)=ZPMIX(JL,JN)+ZDP
              ZTHMIX(JL,JN)=ZTHMIX(JL,JN)+PT(JL,JK)*ZDP*ZEXN(JL,JK)
              ZQMIX(JL,JN)=ZQMIX(JL,JN)+PQ(JL,JK)*ZDP
            ENDIF
            IF( ZPMIXA(JL)<=ZPM(JNP)) THEN
              ZPMIX(JL,JNP)=ZPMIX(JL,JNP)+ZDP
              ZTHMIX(JL,JNP)=ZTHMIX(JL,JNP)+PT(JL,JK)*ZDP*ZEXN(JL,JK)
              ZQMIX(JL,JNP)=ZQMIX(JL,JNP)+PQ(JL,JK)*ZDP
            ENDIF
        ENDDO
      ENDDO
      DO JN=1,NL
        DO JL=KIDIA,KFDIA 
         ZP=1.0_JPRB/ZPMIX(JL,JN)
         ZTHMIX(JL,JN)=ZTHMIX(JL,JN)*ZP
         ZQMIX(JL,JN)=ZQMIX(JL,JN)*ZP
       ENDDO
      ENDDO

   !compute LCL for these mixed layer parcels
      DO JN=1,NL
        DO JL=KIDIA,KFDIA
          ZPMIX(JL,JN)=PAPH(JL,KLEV+1)-0.5_JPRB*ZPMIX(JL,JN)
        ! dewpoint temperature
          ZTMIX(JL)=ZTHMIX(JL,JN)*(ZPMIX(JL,JN)/RATM)**RKAPPA
          ZTEMP = LOG(MAX(1.E-4_JPRB,ZQMIX(JL,JN))*ZPMIX(JL,JN)*ZCST)
          ZTD = (R4LES*ZTEMP-R3LES*RTT)/(ZTEMP-R3LES)
        ! adiabatic saturation temperature
          ZTLCL  = ZTD - ( .212_JPRB + 1.571E-3_JPRB * ( ZTD - RTT )      &
                 - 4.36E-4_JPRB* ( ZTMIX(JL) - RTT ) ) * ( ZTMIX(JL) - ZTD )
          ZTLCL  = MAX(160.0_JPRB,MIN( ZTLCL, ZTMIX(JL) ))
          ZPLCL(JL,JN)  = RATM * ( ZTLCL / ZTHMIX(JL,JN) ) ** ZORKAPPA
          JDEP(JL,JN)=JN
        ENDDO
      ENDDO
  ENDIF

  ! loop over 0=maximum unstable and NL mixed layer parcels

  DO JN=0,NL
    JK=KLEV-1
    DO JL=KIDIA,KFDIA
       ZCIN2(JL)=0.0_JPRB
       JKK=JDEP(JL,JN)
       ZTU(JL,JK)=ZTHMIX(JL,JKK)/ZEXN(JL,JK)
       ZQS(JL,JK)=ZQMIX(JL,JKK)
       ZCAPE(JL,JKK)=0.0_JPRB
    ENDDO
    DO JK=KLEV-2,NJKT2,-1
       DO JL=KIDIA,KFDIA
          JKK=JDEP(JL,JN)
          LLZ(JL)=LLPAP(JL,JK) .AND. PAP(JL,JK)<=ZPLCL(JL,JKK)
       ENDDO
       DO JL=KIDIA,KFDIA
         IF(LLPAP(JL,JK)) THEN
           ZTU(JL,JK)=ZTU(JL,JK+1)*ZEXN(JL,JK+1)/ZEXN(JL,JK)
           ZQS(JL,JK)=ZQS(JL,JK+1)
         ENDIF
       ENDDO
  
       CALL CUADJTQ &
       & ( YDEPHLI, KIDIA,    KFDIA,    KLON,    KLEV, JK,&
       &   PAP(:,JK),   ZTU,   ZQS,     LLZ,   1)
  
       DO JL=KIDIA,KFDIA
         JKK=JDEP(JL,JN)
         LLDEP(JL)=LLPAP(JL,JK).AND.PAP(JL,JK)<=ZPMIX(JL,JKK)
         IF(LLDEP(JL)) THEN
           ZTVEMP=ZTU(JL,JK)*ZEXN(JL,JK)*(1.0_JPRB+RETV*ZQS(JL,JK))*ZTHV(JL,JK)-1.0_JPRB
           IF(PAP(JL,JK)<=ZPLCL(JL,JKK)) THEN
             ZCAPE(JL,JKK)=ZCAPE(JL,JKK)+MAX(0.0_JPRB,ZTVEMP)*ZZDZ(JL,JK)
           ENDIF
           IF(ZTVEMP < 0.0_JPRB.AND.ZCAPE(JL,JKK)<20.0_JPRB) THEN
             ZCIN2(JL)=ZCIN2(JL)+ZTVEMP*ZZDZ(JL,JK)
           ENDIF
         ENDIF
       ENDDO
    ENDDO
  
    DO JL=KIDIA,KFDIA
      JKK=JDEP(JL,JN)
      PCIN(JL,JN+1) =-MAX(RMINCIN,ZCIN2(JL))
      PCAPE(JL,JN+1)=ZCAPE(JL,JKK)
    ENDDO

 ENDDO

ENDIF


!-------------------------------------------------------------------------------
  
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUANCAPE2',1,ZHOOK_HANDLE)
END SUBROUTINE CUANCAPE2

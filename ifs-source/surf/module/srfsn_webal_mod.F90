! (C) Copyright 2015- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
MODULE SRFSN_WEBAL_MOD
CONTAINS
SUBROUTINE SRFSN_WEBAL(KIDIA,KFDIA,KLON,KLEVSN,&
 & PTMST,LLNOSNOW,PFRSN,&
 & PSSNM1M,PWSNM1M,PRSNM1M,PTSNM1M,&
 & PTSURF,PHFLUX,ZSNOTRS,PSNOWF,PRAINF,PEVAPSN,PSURFCOND,&
 & PAPRS,&
 & YDSOIL,YDCST,&
 & PSSN,PWSN,PTSN,&
 & PGSN,PRUSN,PMELTSN,PFREZSN,&
 & PDHTSS,PDHSSS)

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_SOIL , ONLY : TSOIL 
USE YOS_CST  , ONLY : TCST
USE ABORT_SURF_MOD

!**** *SRFSN_WEBAL* - Snow water & energy balance 
!     PURPOSE.
!     --------
!     THIS ROUTINE COMPUTES ENERGY AND WATER BALANCE IN THE SNOWPACK

!**   INTERFACE.
!     ----------
!          *SRFSN_WEBAL* IS CALLED FROM *SRFSN_DRIVER*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET
!    *KLEVSN*     VERTICAL SNOW LAYERS

!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)
!    *LLNOSNOW*   NO-SNOW/SNOW MASK (TRUE IF NO-SNOW)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PWLM1M*     SKIN RESERVOIR WATER CONTENT                 kg/m**2
!    *PFRSN*      total snow fraction tile 5 + 7
!    *PSSNM1M*    TOTAL SNOW MASS IN EACH LAYER (per unit area) kg/m**2
!    *PWSNM1M*    LIQUID WATER CONTENT IN SNOW                 kg/m**2
!    *PRSNM1M*    SNOW DENSITY in each layer                   kg/m**3
!    *PTSNM1M*    TEMPERATURE OF SNOW LAYER                    K
!    *PTSURF*     TEMPERATURE OF TOP SOIL LAYER                K
!    *PHFLUX*     CONDUCTIVE HEAT FLUX INTO THE SNOWPACK       W/m**2
!    *PSNOWF*     TOTAL SNOW FLUX AT THE SURFACE         KG/M**2/S
!    *PRAINF*     TOTAL RAIN FLUX AT THE SURFACE         KG/M**2/S
!    *PEVAPSN*    EVAPORATION FROM SNOW UNDER FOREST           KG/M2/S
!    *PSURFCOND*  THERMAL CONDUCTIVITY OF TOP SOIL LAYER
!    *ZSNOTRS*    SOLAR RADIATION FLUX INTO THE SNOWPACK      W/m**2
!    *PAPRS*      ATMOSPHERIC PRESSURE ON BOTTOM HALF LEVEL   Pa
!    

!     OUTPUT PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PSSN*        SNOW MASS each layer (per unit area)        kg/m**2
!    *PWSN*        LIQUID WATER CONTENT IN SNOW                 kg/m**2
!    *PTSN*        TEMPERATURE OF SNOW LAYER                    K

!    FLUXES FROM SNOW SCHEME:
!    *PGSN*       GROUND HEAT FLUX FROM SNOW DECK TO SOIL     W/M**2   (#)
!    *PRUSN*      FLUX OF MELT WATER FROM SNOW TO SOIL       KG/M**2/S (#)
!    *PMELTSN*    LATENT HEAT OF MELTED WATER                 J/m**2
!    *PFREZSN*    LATENT HEAT OF REFREEZE                     J/m**2



!     OUTPUT PARAMETERS (DIAGNOSTIC):
!    *PDHIIS*     Diagnostic array for interception layer (see module yomcdh)

!     METHOD.
!     -------
!          

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          

!     Modifications:
!     Original   E. Dutra      ECMWF     04/12/2015
!                G. Arduini              01/09/2021

!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments 
INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVSN
REAL(KIND=JPRB)   , INTENT(IN)   :: PTMST
LOGICAL           , INTENT(IN)   :: LLNOSNOW(:) 

REAL(KIND=JPRB)   , INTENT(IN)   :: PFRSN(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PSSNM1M(:,:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PWSNM1M(:,:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PRSNM1M(:,:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PTSNM1M(:,:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PTSURF(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PHFLUX(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: ZSNOTRS(:,:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PSNOWF(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PRAINF(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PEVAPSN(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PSURFCOND(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PAPRS(:)

TYPE(TSOIL)       , INTENT(IN)   :: YDSOIL
TYPE(TCST)        , INTENT(IN)   :: YDCST

REAL(KIND=JPRB)   , INTENT(OUT)  :: PSSN(:,:)
REAL(KIND=JPRB)   , INTENT(OUT)  :: PWSN(:,:)
REAL(KIND=JPRB)   , INTENT(OUT)  :: PTSN(:,:)
REAL(KIND=JPRB)   , INTENT(OUT)  :: PGSN(:)
REAL(KIND=JPRB)   , INTENT(OUT)  :: PRUSN(:)
REAL(KIND=JPRB)   , INTENT(OUT)  :: PMELTSN(:,:)
REAL(KIND=JPRB)   , INTENT(OUT)  :: PFREZSN(:,:)

REAL(KIND=JPRB)   , INTENT(OUT)  :: PDHTSS(:,:,:)
REAL(KIND=JPRB)   , INTENT(OUT)  :: PDHSSS(:,:,:)

! Local variables 
REAL(KIND=JPRB) :: ZDSN(KLEVSN)    ! actual snow depth
REAL(KIND=JPRB) :: ZSNHC(KLEVSN)   ! snow heat capacity 
REAL(KIND=JPRB) :: ZICE(KLEVSN)    ! snow ice content (PSSN-PWSN)
REAL(KIND=JPRB) :: ZSNCOND(KLEVSN) ! Snow thermal conductivity 
REAL(KIND=JPRB) :: ZSNCONDH(KLEVSN+1) ! THERMAL CONDUCTIVITY IN HALF LEVEL TERM 
REAL(KIND=JPRB) :: ZTA(KLEVSN),ZTB(KLEVSN),ZTC(KLEVSN),ZTR(KLEVSN)  ! TERMS TO TRI-DIAG
REAL(KIND=JPRB) :: ZTSTAR(KLEVSN)  ! New snow temperature 
REAL(KIND=JPRB) :: ZISTAR(KLEVSN)  ! New ice content 
REAL(KIND=JPRB) :: ZWSTAR(KLEVSN)  ! New liquid water content 
REAL(KIND=JPRB) :: ZLIQF(0:KLEVSN) ! LIQUID WATER FLUX 

REAL(KIND=JPRB) :: ZGSN,ZWCAP,ZQ,ZGSNRES

REAL(KIND=JPRB) :: ZSNVCOND
REAL(KIND=JPRB) :: ZTMP0,ZTMP1
REAL(KIND=JPRB)  :: ZTMST,ZIHCAP,ZSOILDEPTH1,ZEPSILON
REAL(KIND=JPRB) :: ZSOILRES, ZHOICE, ZSNRES, ZWHCAP
REAL(KIND=JPRB) :: ZDTM, ZDTH
INTEGER(KIND=JPIM) :: JL,JK,KLACT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!! INCLUDE FUNCTIONS

#include "fcsurf.h"

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_WEBAL_MOD:SRFSN_WEBAL',0,ZHOOK_HANDLE)

!    -----------------------------------------------------------------
ASSOCIATE(RHOCI=>YDSOIL%RHOCI,RHOICE=>YDSOIL%RHOICE,&
 & RTT=>YDCST%RTT,RLMLT=>YDCST%RLMLT,&
 & RALAMSN=>YDSOIL%RALAMSN, RLAMICE=>YDSOIL%RLAMICE,&
 & RDSNMAX=>YDSOIL%RDSNMAX, SNHCONDAV=>YDSOIL%SNHCONDAV, &
 & SNHCONDBV=>YDSOIL%SNHCONDBV, SNHCONDCV=>YDSOIL%SNHCONDCV, &
 & SNHCONDPOV=>YDSOIL%SNHCONDPOV )

! RLMLT (latent heat of fusion J Kg -1)
ZTMST = 1.0_JPRB/PTMST 
ZIHCAP = RHOCI/RHOICE  ! Ice heat capacity (J K-1 Kg-1)
ZSOILDEPTH1=YDSOIL%RDAW(1) ! 1sth Soil layer depth
ZEPSILON=10._JPRB*EPSILON(ZEPSILON)
ZWHCAP=4180._JPRB ! J K-1 Kg-1

DO JL=KIDIA,KFDIA
  IF (LLNOSNOW(JL)) THEN 
    PSSN(JL,1) = MAX(0._JPRB, PSSNM1M(JL,1) + PTMST*(PSNOWF(JL)+PEVAPSN(JL)))  ! update snow mass with snowfall
    PSSN(JL,2:KLEVSN) = 0._JPRB 
    PTSN(JL,:) = RTT
    PWSN(JL,:) = 0.0_JPRB
    PGSN(JL)   = 0.0_JPRB
    PRUSN(JL)  = 0.0_JPRB
    PMELTSN(JL,:) = 0.0_JPRB
    PFREZSN(JL,:) = 0.0_JPRB 
    
!! DDH DIAGNOSTICS 
    IF (SIZE(PDHTSS) > 0 .AND. SIZE(PDHSSS) > 0) THEN
      ! Snow heat capacity per unit surface
      !PDHTSS(JL,1:KLEVSN,1)=SUM(ZSNHC(1:KLEVSN)*PTMST)
      PDHTSS(JL,1:KLEVSN,1)=0._JPRB
      ! Snow temperature
      PDHTSS(JL,1:KLEVSN,2)=RTT
      ! Snow energy per unit surface
      PDHTSS(JL,1:KLEVSN,3)=0._JPRB
      ! Snow thermally active depth 
      PDHTSS(JL,1:KLEVSN,4)=0._JPRB
      ! Snow basal heat flux
      PDHTSS(JL,1:KLEVSN,13)=0._JPRB
      ! phase changes
      PDHTSS(JL,1:KLEVSN,14)=0._JPRB
      ! Snow heat content change
      PDHTSS(JL,1:KLEVSN,15)=0._JPRB
      
      ! SNOW MASS
      PDHSSS(JL,1:KLEVSN,1)=0._JPRB
      ! snow evap
      PDHSSS(JL,1:KLEVSN,4)=0._JPRB
      ! phase changes 
      PDHSSS(JL,1:KLEVSN,5)=0._JPRB
    ENDIF
  ELSE

!! Preparation
    DO JK=1,KLEVSN
      IF (PSSNM1M(JL,JK) > ZEPSILON ) KLACT=JK
    ENDDO
    DO JK=1,KLEVSN
      ZICE(JK)    = PSSNM1M(JL,JK) - PWSNM1M(JL,JK) ! Ice (Kg m-2)
      ! Limit thermally active snow depth to 1 meter of snow (for glaciers in particular)
      ZDSN(JK)    = MIN(RDSNMAX, (PSSNM1M(JL,JK) / PRSNM1M(JL,JK))) ! read snow depth (m)
      ZSNHC(JK)   = (ZIHCAP*PRSNM1M(JL,JK) * MIN(RDSNMAX, (ZICE(JK)/PRSNM1M(JL,JK))) + ZWHCAP*PWSNM1M(JL,JK) ) * ZTMST
      ! heat conductivity from water vapor transport into the snowpack
      ZSNVCOND=(SNHCONDPOV/PAPRS(JL))*MAX(0._JPRB,(SNHCONDAV-SNHCONDBV/(PTSNM1M(JL,JK)-SNHCONDCV)))
      ! snow heat conductivity 
      ZSNCOND(JK) = FSNTCOND(PRSNM1M(JL,JK))+ZSNVCOND  ! add the thermal cond from water vapor transfer
    ENDDO
    ! special case for 1 layer only
    IF (KLACT == 1 ) ZSNHC(1) = MAX(ZEPSILON, ZSNHC(1))
     
    ZSNCONDH(1)=PFRSN(JL)*ZSNCOND(1) / MAX(ZEPSILON,(0.5_JPRB*ZDSN(1))*(0.5_JPRB*ZDSN(1)))   ! ACTUALY NOT USED ! (W m-2 K-1)
    DO JK=2,KLACT
      ZSNCONDH(JK)=PFRSN(JL)*2._JPRB*(ZDSN(JK-1)*ZSNCOND(JK-1)+ZDSN(JK)*ZSNCOND(JK))/&
                  & MAX(ZEPSILON, (ZDSN(JK-1)+ZDSN(JK))*(ZDSN(JK-1)+ZDSN(JK)))
    ENDDO
    ! Divided by two to decrease heat flux between snow and soil (accounting for litter,organic etc). 
    ! This needs to be in double-precision to avoid some problematic computation between small numbers 
     ZSNCONDH(KLACT+1)=REAL(PFRSN(JL),KIND=JPRD)*&
                      &(REAL(ZDSN(KLACT),KIND=JPRD)*REAL(ZSNCOND(KLACT),KIND=JPRD)+ZSOILDEPTH1*PSURFCOND(JL))/&
                      &MAX(ZEPSILON, REAL(ZDSN(KLACT)+ZSOILDEPTH1,KIND=JPRD)*REAL(ZDSN(KLACT)+ZSOILDEPTH1,KIND=JPRD))
    ! ZSNCONDH(KLACT+1)=PFRSN(JL)*1._JPRB*(ZDSN(KLACT)*ZSNCOND(KLACT)+ZSOILDEPTH1*PSURFCOND(JL))/&
    !              &MAX(ZEPSILON, (ZDSN(KLACT)+ZSOILDEPTH1)*(ZDSN(KLACT)+ZSOILDEPTH1))

    IF (KLACT > 1 ) THEN 
      JK=1
      ZTA(JK)=0._JPRB
      ZTB(JK)=ZSNHC(JK)+ZSNCONDH(JK+1)
      ZTC(JK)=-ZSNCONDH(JK+1)
      ZTR(JK)=(PHFLUX(JL)-ZSNOTRS(JL,JK))+ZSNHC(JK)*PTSNM1M(JL,JK)
! Add radiative flux to second to bottom layer:
      JK=KLACT
      ZTA(JK)=-ZSNCONDH(JK)
      ZTB(JK)=ZSNHC(JK)+ZSNCONDH(JK)+ZSNCONDH(JK+1)
      ZTC(JK)=0._JPRB
      ZTR(JK)=ZSNHC(JK)*PTSNM1M(JL,JK)+ZSNCONDH(JK+1)*PTSURF(JL)+ZSNOTRS(JL,JK)
      DO JK=2,KLACT-1
        ZTA(JK)=-ZSNCONDH(JK)
        ZTB(JK)=ZSNHC(JK)+ZSNCONDH(JK)+ZSNCONDH(JK+1)
        ZTC(JK)=-ZSNCONDH(JK+1)
        ZTR(JK)=ZSNHC(JK)*PTSNM1M(JL,JK) + ZSNOTRS(JL,JK)
      ENDDO
    ENDIF
    
    ! SOLVE
    ZTSTAR(:) = PTSNM1M(JL,:)
    IF (KLACT == 1 ) THEN
      ! SINGLE LAYER
      ZTB(1)=ZSNHC(1)+ZSNCONDH(2)
      ZTR(1)=(PHFLUX(JL)-ZSNOTRS(JL,1))+ZSNHC(1)*PTSNM1M(JL,1)+ZSNCONDH(2)*PTSURF(JL)
      ZTSTAR(1)=ZTR(1)/ZTB(1)
      ZGSN=ZSNCONDH(2)*(ZTSTAR(1)-PTSURF(JL))
    ELSE
      ! MULTI LAYER
      CALL TRISOLVER(KLACT,ZTA,ZTB,ZTC,ZTR,ZTSTAR)
      ZGSN=ZSNCONDH(KLACT+1)*(ZTSTAR(KLACT)-PTSURF(JL))
    ENDIF 
    ! * SPECIAL CASE WHEN THERE IS MELTING ON THE 1ST LAYER
    IF ( ZTSTAR(1) > RTT ) THEN 
      IF (KLACT == 1 ) THEN
        ! SINGLE LAYER
        ZGSN=ZSNCONDH(2)*(RTT-PTSURF(JL))
        ZTSTAR(1)=(PHFLUX(JL)-ZSNOTRS(JL,1)+ZSNHC(1)*PTSNM1M(JL,1)-ZGSN)/MAX(ZEPSILON,ZSNHC(1))
      ELSE
        ! MULTI LAYER
        ! SOLVE SYSTEM SETTING T1 == RTT 
        ZTB(1) = 1._JPRB
        ZTC(1) = 0._JPRB 
        ZTR(1) = RTT 
        ZTSTAR(:) = PTSNM1M(JL,:)  
        CALL TRISOLVER(KLACT,ZTA,ZTB,ZTC,ZTR,ZTSTAR)
        ZGSN=ZSNCONDH(KLACT+1)*(ZTSTAR(KLACT)-PTSURF(JL))
        ! EXPLICIT 1ST LAYER TEMP CALCULATION
        ZTSTAR(1)=(PHFLUX(JL)-ZSNOTRS(JL,1)+ZSNHC(1)*PTSNM1M(JL,1)-ZSNCONDH(2)*(ZTSTAR(1)-ZTSTAR(2)))/MAX(ZEPSILON,ZSNHC(1))
      ENDIF
    ENDIF
    

!     !! UPDATE SOLID/LIQUID CONTENTS OF FIRST LAYER BEFORE FURTHER CALCULATIONS
    ZISTAR(:) = ZICE(:)
    ZWSTAR(:) = PWSNM1M(JL,:)
    ZLIQF(:)  = 0._JPRB
      !! INCLUDE SNOWFALL AND RAINFALL INTERCEPTION 
    ZISTAR(1) = ZICE(1)+PTMST*(PSNOWF(JL)+PEVAPSN(JL))
    ZLIQF(0) = PTMST*PRAINF(JL)
    
    !! PHASE CHANGES AND LIQUID WATER BALANCE STARTING IN LAYER 1
    PMELTSN(JL,:)=0._JPRB
    PFREZSN(JL,:)= 0._JPRB
    ZGSNRES=0._JPRB
    DO JK=1,KLACT
      ! UPDATE LIQUID WATER FROM FLUX ABOVER 
      ZWSTAR(JK) = PWSNM1M(JL,JK) + ZLIQF(JK-1)
      
      ! PHASE CHANGES 
      ZTMP0 = ZSNHC(JK)*(ZTSTAR(JK)-RTT)
      PMELTSN(JL,JK) = MAX(0._JPRB, MIN( ZTMP0 , RLMLT*ZTMST*ZISTAR(JK) ) )
      PFREZSN(JL,JK) = MIN(0._JPRB, MAX( ZTMP0 , -RLMLT*ZTMST*ZWSTAR(JK) ) )
      ZQ = PMELTSN(JL,JK) + PFREZSN(JL,JK)
      
      ! FINAL TEMP. UPDATE
      ZTSTAR(JK) = ZTSTAR(JK) - (ZQ+ZGSNRES)/MAX(ZEPSILON,ZSNHC(JK))
      IF ( ZTSTAR(JK) > RTT ) THEN
        ZGSNRES=ZSNHC(JK)*(RTT-ZTSTAR(JK))
        ZTSTAR(JK)=RTT
      ELSE
        ZGSNRES=0._JPRB
      ENDIF
      
      ! UPDATE SOLID / LIQUID MASS
      ZWSTAR(JK) = ZWSTAR(JK) + ZQ/RLMLT*PTMST 
      ZISTAR(JK) = ZISTAR(JK) - ZQ/RLMLT*PTMST 
      
      ! LIQUD WATER BUDGET
      ZWCAP = FLWC(ZISTAR(JK),PRSNM1M(JL,JK))
      ZLIQF(JK) = MAX( 0._JPRB, (ZWSTAR(JK)-ZWCAP))
      ZWSTAR(JK) = ZWSTAR(JK)-ZLIQF(JK)
    ENDDO
    DO JK=KLACT+1,KLEVSN
      ZTSTAR(JK) = ZTSTAR(KLACT)
    ENDDO  
      
    !! FINAL VALUES
    DO JK=1,KLEVSN
      PTSN(JL,JK) = MIN( RTT, ZTSTAR(JK) )
      PSSN(JL,JK) = MAX( 0._JPRB, ZWSTAR(JK) + ZISTAR(JK) )
      PWSN(JL,JK) = MIN( PSSN(JL,JK), MAX(0._JPRB, ZWSTAR(JK)) )
    ENDDO

    PRUSN(JL) = ZLIQF(KLACT)*ZTMST
    PGSN(JL) = ZGSN - ZGSNRES + ZSNOTRS(JL,KLACT+1)
    

    
    IF (ANY(PTSN(JL,:)>RTT+ZEPSILON)) THEN
      write(*,*) 'All snow melted...'
      CALL ABORT_SURF('ALL SNOW MELTED')
    ENDIF
    IF (ANY(PTSN(JL,:)<100._JPRB)) THEN
      write(*,*) 'Very snow cold temperature'
      write(*,*) 'Tsn-1',PTSNM1M(JL,:)
      write(*,*) 'Tsn',PTSN(JL,:)
      write(*,*) 'SWE-1',PSSNM1M(JL,:)
      write(*,*) 'SWE',PSSN(JL,:)
      write(*,*) 'Snow frac,heat,pg0',PFRSN(JL),PHFLUX(JL),PGSN(JL)
      
      CALL ABORT_SURF('Very snow cold temperature')
    ENDIF 
    
    
    !! DDH DIAGNOSTICS 
    IF (SIZE(PDHTSS) > 0 .AND. SIZE(PDHSSS) > 0) THEN
      DO JK=1,KLEVSN
      ! Snow heat capacity per unit surface
        PDHTSS(JL,JK,1)=ZSNHC(JK)*PTMST
      ! Snow temperature
        PDHTSS(JL,JK,2)=PTSNM1M(JL,JK)
      ! Snow energy per unit surface
        PDHTSS(JL,JK,3)=ZSNHC(JK)*PTMST*PTSNM1M(JL,JK)
      ! Snow thermally active depth 
        PDHTSS(JL,JK,4)=PSSNM1M(JL,JK) / PRSNM1M(JL,JK)
      ! Snow basal heat flux
        PDHTSS(JL,1,13)=PGSN(JL)
      ! phase changes
        PDHTSS(JL,JK,14)=PMELTSN(JL,JK)+PFREZSN(JL,JK)
      ! Snow heat content change
        PDHTSS(JL,JK,15)=ZSNHC(JK)*PTMST*(PTSN(JL,JK)-PTSNM1M(JL,JK))
      
      ! SNOW MASS
        PDHSSS(JL,JK,1)=PSSNM1M(JL,JK)
      ! snow evap
        PDHSSS(JL,JK,4)=PEVAPSN(JL)
      ! phase changes (snow melt/refreezing)
        PDHSSS(JL,JK,5)=PDHTSS(JL,JK,14)/RLMLT
!       ! snow liquid water content 
!       PDHSSS(JL,JK,6)=PWSNM1M(JL,JK)
!       ! snow runoff 
!       PDHSSS(JL,JK,8)=PRUSN(JL)
      ENDDO
    ENDIF 

  ENDIF 
ENDDO
                           
  
END ASSOCIATE
!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_WEBAL_MOD:SRFSN_WEBAL',1,ZHOOK_HANDLE)

END SUBROUTINE SRFSN_WEBAL

SUBROUTINE TRISOLVER(KSNLAC,ZA,ZB,ZC,ZF,ZTSTAR)
USE PARKIND1, ONLY : JPIM, JPRB
IMPLICIT NONE

!! SINGLE POINT TRIDIAGONAL SOLVER ! 

!* DECLARATION OF ARGUMENTS
INTEGER(KIND=JPIM)  ,INTENT(IN)  :: KSNLAC ! NUMBER OF ACTIVE SNOW LAYERS
REAL(KIND=JPRB),INTENT(IN),DIMENSION(:)  :: ZA,ZB,ZC,ZF
REAL(KIND=JPRB),INTENT(INOUT),DIMENSION(:) :: ZTSTAR

!* LOCAL VARIABLES  
INTEGER(KIND=JPIM)  :: JK
REAL(KIND=JPRB)     :: BET
REAL(KIND=JPRB),DIMENSION(SIZE(ZB)) :: GAM

BET=ZB(1)

ZTSTAR(1)=ZF(1)/BET
DO JK=2,KSNLAC
  GAM(JK)=ZC(JK-1)/BET
  BET=ZB(JK)-ZA(JK)*GAM(JK)
  ZTSTAR(JK)=(ZF(JK)-ZA(JK)*ZTSTAR(JK-1))/BET
ENDDO
DO JK=KSNLAC-1,1,-1
  ZTSTAR(JK)=ZTSTAR(JK)-GAM(JK+1)*ZTSTAR(JK+1)
ENDDO

END SUBROUTINE TRISOLVER


END MODULE SRFSN_WEBAL_MOD

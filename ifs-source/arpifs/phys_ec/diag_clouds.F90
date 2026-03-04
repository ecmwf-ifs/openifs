! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE DIAG_CLOUDS&
 & ( YDECUMF,KIDIA,  KFDIA,  KLON,   KLEV,   LDIAG_CL,&
 &   LDCUM,  KCBOT,  KCTOP,  PAP,    PAPH,   PGEO,   PGEOH,&   
 &   PT,     PQ,     PL,     PI,     PA,&
 &   PFPLCL, PFPLCN, PFPLSL, PFPLSN,&
 &   PRAINFRAC_TOPRFZ, P2T, P2D,&
 &   PCBASE, PCBASEA,PCTOPC, P0DEGL, PM10DEGL, PCONVIND,&
 &   PPRECTYPE, PFZRA, PZTWETB, PTROPOTP)

!    P.BECHTOLD AND RICHARD FORBES    E.C.M.W.F.     02/2010

!    PURPOSE
!    -------

!    TO PROVIDE DIAGNOSTICS FOR CLOUD BASE
!    AND ZERO DEGREE LEVEL AND PRECIPITATION TYPE

!    INTERFACE
!    ---------
!    THIS ROUTINE IS CALLED FROM *CALLPAR*.

!    METHOD.
!    --------
!    CLOUD BASE COMPUTATIONS TAKE INTO ACCOUNT "STRATIFORM" AND
!    CONVECTIVE CLOUD BASE
!    ZERO DEGREE LEVEL DETERMINATIONSEARCHES FOR FIRST SWAP FROM POSITIVE
!    TO NEGATIVE TEMPERATURES, CORRESPONDING HEIGHT IS SET TO MID-LEVEL VALUE

!    CLOUD BASE AND ZERO DEGREE LEVEL ARE GIVEN IN (M) ABOVE GROUND

!    PARAMETER     DESCRIPTION                                   UNITS 
!    ---------     -----------                                   ----- 
!    INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (LOGICAL):

!    *LDIAG_CL*     DO MOST COMPUTATIONS (APART FROM PRECIP DIAGNOSTICS) ONLY AT 
!                   TIME INTERVALS AS DEFINED IN CALLPAR
!    *LDCUM*        CONVECTION FLAG

 
!     INPUT PARAMETERS (REAL)

!    *PAP*          PRESSURE                                       PA
!    *PAPH*         PRESSURE ON HALF LEVELS                        PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PT*           TEMPERATURE                                    K
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG
!    *PL*           LIQUID WATER CONTENT                          KG/KG
!    *PI*           ICE WATER CONTENT                             KG/KG
!    *PA*           CLOUD FRACTION                                (0-1)
!    *PFPLCL*       CONVECTIVE LIQUID PRECIPITATION FLUX          KG/M2/S
!    *PFPLCN*       CONVECTIVE SNOW FLUX                          KG/M2/S
!    *PFPLSL*       RESOLVED LIQUID PRECIPITATION FLUX            KG/M2/S
!    *PFPLSN*       RESOLVED SNOW FLUX                            KG/M2/S
!    *PRAINFRAC_TOPRFZ* RAIN FRACTION AT TOP OF FREEZING LAYER
!    *P2T*          2M TEMPERATURE                                  K
!    *P2D*          2M DEWPOINT                                     K

!    OUTPUT PARAMETERS (REAL) 

!    *PCBASE*      CLOUD BASE HEIGHT                                M
!    *PCBASEA*     CLOUD BASE FOR AVIATION (5 OCTA CC REQUIRED)     M
!    *PTOPC*       CLOUD TOP HEIGHT CONVECTIVE                      M


!    *P0DEGL*      0 DEGREE CELSIUS LEVEL                           M 
!    *PM10DEGL*  -10 DEGREE CELSIUS LEVEL                           M 
!    *PCONVIND*    CONVECTIVE INDICES                             deg C 
!    *PPRECTYPE*   PRECIPITATION TYPE
!    *PFZRA*       ACCUMULATED FREEZING RAIN FLUX                 KG/M2
!    *PZTWETB*     HEIGHT OF 0 and 1C WET BULB TEMPERATURE          M

!          MODIFICATIONS
!          -------------
!    06/2011: Add convective Indices (Total Totals, K-Index) P.BECHTOLD
!    R.Forbes 01-Mar-2014 Added precipitation type
!    R.Forbes 10-Jan-2015 Added freezing rain FZRA  and Precipitation type
!    P.Bechtold 10-Nov-2015 Added convective cloud top, cloud base aviation
!                           height of 0,1C wet bulb T as in cloudsc.F90
!    R.Forbes 15-Dec-2019 Passed in dew point and 2m wet-bulb calculation
!                         Rewrite of T<0 precip type
!    P.Bechtold 26-Mar-2021 Added tropopause pressure (thermal)
!    R.Forbes   15-Sep-2021 Added freezing drizzle precip type
!
!----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RTT, RETV, RLVTT, RLSTT
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                    R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
 &                    RTWAT_RTICE_R, RTWAT_RTICECU_R
USE YOECUMF  , ONLY : TECUMF

IMPLICIT NONE

TYPE(TECUMF)      ,INTENT(INOUT) :: YDECUMF
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 

INTEGER(KIND=JPIM),INTENT(IN)    :: KCBOT(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCTOP(KLON)

LOGICAL,           INTENT(IN)    :: LDIAG_CL
LOGICAL,           INTENT(IN)    :: LDCUM(KLON)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLSL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLSN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAINFRAC_TOPRFZ(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P2T(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P2D(KLON) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCBASE(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCBASEA(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTOPC(KLON)

REAL(KIND=JPRB)   ,INTENT(OUT)   :: P0DEGL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PM10DEGL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCONVIND(KLON,2) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPRECTYPE(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFZRA(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZTWETB(KLON,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTROPOTP(KLON) 

LOGICAL :: LLSTB(KLON,2)
INTEGER(KIND=JPIM) :: JL, JK, JKS, IC0L(KLON), ICM10L(KLON), ILTROP(KLON)
REAL(KIND=JPRB) ::     ZRG, ZBASEC(KLON), ZBASES(KLON),&
                     & ZA, ZTD8, ZTD7, ZP, ZDZ
! Rain fraction and freezing rain variables
REAL(KIND=JPRB)    :: ZRAIN, ZSNOW, ZTOT, ZRAINFRAC_SFC

! Numerical fit to wet bulb temperature
REAL(KIND=JPRB),PARAMETER :: ZTW1 = 1329.31_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW2 = 0.0074615_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW3 = 0.85E5_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW4 = 40.637_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW5 = 275.0_JPRB
REAL(KIND=JPRB)           :: ZTWETB(KLON,KLEV+1), ZALFA, ZQSICE, ZSUBSAT, ZTM10
!REAL(KIND=JPRB)           :: Z2Q, Z2QS, Z2TWETB

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcttre.func.h"
#include "troplev.intfb.h"

!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIAG_CLOUDS',0,ZHOOK_HANDLE)
ASSOCIATE(NJKT2=>YDECUMF%NJKT2, NJKT4=>YDECUMF%NJKT4, NJKT5=>YDECUMF%NJKT5, &
 & NJKT6=>YDECUMF%NJKT6, RBASE0=>YDECUMF%RBASE0)

ZRG=1.0_JPRB/RG
ZTM10=RTT-10.0_JPRB

! Initialisations
DO JL=KIDIA,KFDIA
  PCBASE(JL)=RBASE0
  PCBASEA(JL)=RBASE0
  PCTOPC(JL)=RBASE0

  ZBASEC(JL)=RBASE0 ! convective cloud base
  ZBASES(JL)=RBASE0 ! "stratiform" or cloud scheme based cloud base
  P0DEGL(JL)=0.0_JPRB
  PM10DEGL(JL)=0.0_JPRB
  PCONVIND(JL,:)=0.0_JPRB
  ! Default precipitation type this timestep
  PPRECTYPE(JL) = 0.0_JPRB
  ! Default freezing rain accumulation this timestep 
  PFZRA(JL) = 0.0_JPRB
  PZTWETB(JL,:) = 0.0_JPRB   
  PTROPOTP(JL) = 0.0_JPRB   
ENDDO

!apart from precip diagnostics do computation only every hour

IF(LDIAG_CL) THEN

  ! Convective Cloud Base + Initialisations
  DO JL=KIDIA,KFDIA
    LLSTB(JL,:)=.TRUE.
    IC0L(JL)=0
    ICM10L(JL)=0
    IF (LDCUM(JL)) THEN
      JK=KCBOT(JL)
      ZBASEC(JL)=PGEOH(JL,JK)*ZRG
      JK=KCTOP(JL)
      PCTOPC(JL)=PGEOH(JL,JK)*ZRG
    ENDIF
  ENDDO
  
  ! Vertical Search loop 
  DO JK=KLEV-1,NJKT2,-1
  ! cloud base defined by liquid water/cloud fraction criterion
    DO JL=KIDIA,KFDIA
      IF(LLSTB(JL,1)) THEN
         IF(PA(JL,JK)>1.E-2_JPRB.AND.(PL(JL,JK)+PI(JL,JK))>1.E-6_JPRB) THEN
           ZBASES(JL)=PGEO(JL,JK)*ZRG
           LLSTB(JL,1)=.FALSE.
         ENDIF
      ENDIF
    ENDDO
  ! cloud base aviation, requiring 5 Octa cloud cover
    DO JL=KIDIA,KFDIA
      IF(LLSTB(JL,2)) THEN
         IF(PA(JL,JK)>5.E-1_JPRB) THEN
           PCBASEA(JL)=PGEO(JL,JK)*ZRG
           LLSTB(JL,2)=.FALSE.
         ENDIF
      ENDIF
    ENDDO
  
  ! zero and -10 degree level, exact zero crossing
    DO JL=KIDIA,KFDIA
      IF(IC0L(JL)<2) THEN
         IF(PT(JL,JK+1)>RTT.AND.PT(JL,JK)<=RTT) THEN
           ZDZ=(RTT-PT(JL,JK+1))*(PGEO(JL,JK)-PGEO(JL,JK+1))/MIN(-1.E-8_JPRB,PT(JL,JK)-PT(JL,JK+1))
           P0DEGL(JL)=(PGEO(JL,JK+1)+ZDZ)*ZRG
           IC0L(JL)=IC0L(JL)+1
         ENDIF
      ENDIF
      IF(ICM10L(JL)<2) THEN
         IF(PT(JL,JK+1)>ZTM10.AND.PT(JL,JK)<=ZTM10) THEN
           ZDZ=(ZTM10-PT(JL,JK+1))*(PGEO(JL,JK)-PGEO(JL,JK+1))/MIN(-1.E-8_JPRB,PT(JL,JK)-PT(JL,JK+1))
           PM10DEGL(JL)=(PGEO(JL,JK+1)+ZDZ)*ZRG
           ICM10L(JL)=ICM10L(JL)+1
         ENDIF
      ENDIF
    ENDDO
  ENDDO
  
  
  DO JL=KIDIA,KFDIA
    PCBASE(JL)=ZBASES(JL)
    IF(ZBASEC(JL)<RBASE0.AND.ZBASES(JL)<RBASE0) PCBASE(JL)=MIN(ZBASEC(JL),ZBASES(JL))
  ENDDO
  
  ! Convective Indices
  ZP=1.0_JPRB/R2ES
  DO JL=KIDIA,KFDIA
  ! Total Totals Index =T850+Td850-2*T500
    ZA=LOG(MAX(1.E-9_JPRB,PQ(JL,NJKT4))*PAP(JL,NJKT4)*ZP)
    ZTD8=(R4LES*ZA-R3LES*RTT)/(ZA-R3LES)
    PCONVIND(JL,1)=PT(JL,NJKT4)+ZTD8-2*PT(JL,NJKT5)
  ! K Index =T850-T500+Td850-(T700-Td700)
    ZA=LOG(MAX(1.E-9_JPRB,PQ(JL,NJKT6))*PAP(JL,NJKT6)*ZP)
    ZTD7=(R4LES*ZA-R3LES*RTT)/(ZA-R3LES)
    PCONVIND(JL,2)=PT(JL,NJKT4)-PT(JL,NJKT5)+ZTD8-PT(JL,NJKT6)+ZTD7-RTT
    IC0L(JL)=0
  ENDDO

  !Wet bulb temperature C exact crossing
  DO JK=KLEV,1,-1
    JKS=MIN(JK+1,KLEV)
    DO JL=KIDIA,KFDIA
      IF(PAP(JL,JK)>60.E2_JPRB.AND.IC0L(JL)<2) THEN
        !---------------------------------------------
        ! ice saturation T<273K
        ! liquid water saturation for T>273K 
        !---------------------------------------------
        ZALFA=FOEDELTA(PT(JL,JK))
        ZQSICE=MIN((ZALFA*FOEELIQ(PT(JL,JK))+ &
          &  (1.0_JPRB-ZALFA)*FOEEICE(PT(JL,JK)))/PAP(JL,JK),0.5_JPRB)
        ZQSICE=ZQSICE/(1.0_JPRB-RETV*ZQSICE)
        ZSUBSAT = MAX(ZQSICE-PQ(JL,JK),0.0_JPRB)       
        ZTWETB(JL,JK)=PT(JL,JK)-ZSUBSAT*(ZTW1+ZTW2*(PAP(JL,JK)-ZTW3)-ZTW4*(PT(JL,JK)-ZTW5))-RTT
        IF(JK<KLEV.AND.ZTWETB(JL,JKS)>0.0_JPRB.AND.ZTWETB(JL,JK)<=0.0_JPRB) THEN
          ZDZ=(0.0_JPRB-ZTWETB(JL,JKS))*(PGEO(JL,JK)-PGEO(JL,JKS))/&
                                        &MIN(-1.E-8_JPRB,ZTWETB(JL,JK)-ZTWETB(JL,JKS))
          PZTWETB(JL,1)=(PGEO(JL,JKS)+ZDZ)*ZRG
         !PZTWETB(JL,1)=PGEOH(JL,JK)*ZRG
          IC0L(JL)=IC0L(JL)+1
        ENDIF
        IF(JK<KLEV.AND.ZTWETB(JL,JKS)>1.0_JPRB.AND.ZTWETB(JL,JK)<=1.0_JPRB) THEN
          ZDZ=(1.0_JPRB-ZTWETB(JL,JKS))*(PGEO(JL,JK)-PGEO(JL,JKS))/&
                                        &MIN(-1.E-8_JPRB,ZTWETB(JL,JK)-ZTWETB(JL,JKS))
          PZTWETB(JL,2)=(PGEO(JL,JKS)+ZDZ)*ZRG
         !PZTWETB(JL,2)=PGEOH(JL,JK)*ZRG
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  !Tropopause pressure
  CALL TROPLEV(KLON,KIDIA,KFDIA,KLEV,.FALSE.,PT,PQ,PAP,ILTROP)
  DO JL=KIDIA,KFDIA
     PTROPOTP(JL)=PAP(JL,ILTROP(JL))
  ENDDO

ENDIF
!-------------------------------------------------------------------------------
!
!*   2.  Diagnose surface precipitation type and accumulate freezing rain
!        Uses WMO code table (4.201) for precipitation type
!        See phys_ec/sucldp.F90 for definitions
!
!-------------------------------------------------------------------------------

DO JL=KIDIA,KFDIA
  
  ! Add together stratiform and convective rain and snow categories
  ZRAIN = PFPLCL(JL,KLEV+1) + PFPLSL(JL,KLEV+1)
  ZSNOW = PFPLCN(JL,KLEV+1) + PFPLSN(JL,KLEV+1)
  ZTOT  = ZRAIN + ZSNOW

  IF (ZTOT > 1.E-10_JPRB) THEN   ! If there is precipitation

    ! Fraction of the total precipitation that is rain at the surface
    ZRAINFRAC_SFC = ZRAIN/ZTOT
    
    ! Calculate wet bulb temperature at 2m
    ! saturation wrt ice for T<273K, wrt water for T>273K 
    ! Commented out at the moment but keep for future 
    !Z2Q     = MIN(FOEEW(P2D(JL))/PAPH(JL,KLEV+1),0.5_JPRB)
    !Z2Q     = Z2Q/(1.0_JPRB-RETV*Z2Q)
    !Z2QS    = MIN(FOEEW(P2T(JL))/PAPH(JL,KLEV+1),0.5_JPRB)
    !Z2QS    = Z2QS/(1.0_JPRB-RETV*Z2QS)
    !ZSUBSAT = MAX(Z2QS-Z2Q,0.0_JPRB)       
    !Z2TWETB = P2T(JL)-ZSUBSAT*(ZTW1+ZTW2*(PAPH(JL,KLEV+1)-ZTW3)-ZTW4*(P2T(JL)-ZTW5))

    ! If 2m T above freezing (dry bulb at the moment although could use wet bulb)
    IF (P2T(JL) > RTT) THEN
    
      IF (ZRAINFRAC_SFC > 0.8_JPRB) THEN
        ! More than 80% of precipitation is rain

        PPRECTYPE(JL) = 1._JPRB  ! Rain
      
      ELSEIF (ZRAINFRAC_SFC > 0.2_JPRB .AND. ZRAINFRAC_SFC  <=0.8_JPRB) THEN
        ! Rain fraction is between 20% and 80%

        PPRECTYPE(JL) = 7._JPRB  ! Melting snow (sleet)
      
      ELSEIF (ZRAINFRAC_SFC <= 0.2_JPRB) THEN
        ! Rain fraction is less than 20% (mainly snow with some melted water)

        PPRECTYPE(JL) = 6._JPRB  ! Wet snow
      
      ELSEIF (ZRAINFRAC_SFC <= 0.01_JPRB) THEN

        ! Rain fraction is essentially zero
        PPRECTYPE(JL) = 5._JPRB  ! Dry snow
      
      ENDIF
    
    ! T below freezing
    ELSE

      !-----------------------------------------------------------------------
      ! No warm layer aloft or very shallow (so <20% rain fraction)
      !-----------------------------------------------------------------------
      IF (PRAINFRAC_TOPRFZ(JL) < 0.01_JPRB) THEN

        IF (ZRAINFRAC_SFC < 0.01_JPRB) THEN
          ! If rain fraction close to zero then dry snow (T<0)

          PPRECTYPE(JL) = 5._JPRB  ! Dry snow

        ELSEIF (ZRAINFRAC_SFC < 0.5_JPRB) THEN 
          ! Rain fraction is less than 50% (mainly snow with some melted water)
          ! but temperatuere is less than zero, so will freeze and not wet snow

          PPRECTYPE(JL) = 5._JPRB  ! Dry snow
      
        ELSE ! IF (ZRAINFRAC_SFC >= 0.5_JPRB) THEN
          ! If rain fraction at the surface is >50% then there must be 
          ! significant supercooled rain production in the cloud 
          ! -> freezing drizzle

          PPRECTYPE(JL) = 12._JPRB ! Freezing drizzle
          PFZRA(JL)     = ZRAIN    ! Mark as freezing rain rate 
         
        ENDIF

      ELSEIF (PRAINFRAC_TOPRFZ(JL) < 0.2_JPRB) THEN
     
        IF (ZRAINFRAC_SFC < 0.01_JPRB) THEN
          ! If rain fraction close to zero then dry snow (T<0)

          PPRECTYPE(JL) = 5._JPRB  ! Dry snow

        ELSEIF (ZRAINFRAC_SFC < 0.5_JPRB) THEN 
          ! Rain fraction is less than 50% (mainly snow with some melted water)
          ! but temperatuere is less than zero, so will freeze and not wet snow

          PPRECTYPE(JL) = 5._JPRB  ! Dry snow
      
        ELSE ! IF (ZRAINFRAC_SFC >= 0.5_JPRB) THEN
          ! Could be freezing drizzle production below warm layer,
          ! but class as freezing rain
          
          PPRECTYPE(JL) = 3._JPRB ! Freezing rain
          PFZRA(JL)     = ZRAIN   ! Freezing rain rate 
         
        ENDIF

      !-----------------------------------------------------------------------
      ! There is a warm layer aloft, but not so deep that all snow has melted
      ! Potential for rapid refreezing and formation of ice pellets 
      !-----------------------------------------------------------------------
      ELSEIF (PRAINFRAC_TOPRFZ(JL) >= 0.2_JPRB .AND. &
            & PRAINFRAC_TOPRFZ(JL) < 0.8_JPRB) THEN 
      
        IF (ZRAINFRAC_SFC <= 0.5_JPRB) THEN
          ! If rain fraction at the surface is <50% then ice pellets
          
          PPRECTYPE(JL) = 8._JPRB  ! Ice pellets
         
        ELSE ! IF (ZRAINFRAC_SFC > 0.5_JPRB) THEN
          ! If rain fraction at the surface is >50% then freezing rain

          PPRECTYPE(JL) = 3._JPRB  ! Freezing rain
          PFZRA(JL)     = ZRAIN    ! Freezing rain rate 
         
        ENDIF
      
      !-----------------------------------------------------------------------
      ! There is a warm layer aloft and it is deep so most snow has melted 
      !-----------------------------------------------------------------------
      ELSE ! IF (PRAINFRAC_TOPRFZ(JL) >= 0.8_JPRB) THEN

        IF (ZRAINFRAC_SFC <= 0.5_JPRB) THEN
          ! If rain fraction at the surface is <50% then there has been 
          ! significant refreezing (ice pellets) or snow production
          
          PPRECTYPE(JL) = 8._JPRB  ! Ice pellets
         
        ELSE ! IF (ZRAINFRAC_SFC > 0.5_JPRB) THEN
          ! If rain fraction at the surface is >50% then most of the 
          ! supercooled rain from the warm layer remains supercooled

          PPRECTYPE(JL) = 3._JPRB  ! Freezing rain
          PFZRA(JL)     = ZRAIN    ! Freezing rain rate 
         
        ENDIF    

      ENDIF ! on PRAINFRAC_TOPRFZ 

    ENDIF ! on 2m temperature

  ENDIF ! on precipitation present

ENDDO

!-----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DIAG_CLOUDS',1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_CLOUDS

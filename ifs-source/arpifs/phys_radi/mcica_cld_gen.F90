! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE MCICA_CLD_GEN &
 &( YDEMCICA,YDRIP,PNUAGE    , PLIQWCIN, PICEWCIN, PSNOWCIN, PRLC_CF, PRLC_CW, PGLAT , PGLON ,&
 &  PSIGMA_QCW, PT, PAP, KLON, KIDIA, KFDIA, KLEV, KNX_LOC, &
 &  LDPPH     , LDMAXRAN  , &
 &  KNCLDY    , PLIQWCIN_S, PICEWCIN_S, PSNOWCIN_S)

! from Barker, Cole and Raisanen

! JJMorcrette 20060115  adapted to ECMWF system by 
! PBechtold+NSemane 13082012 add RG/RD for small planet
! N. Semane+P.Bechtold 04-10-2012 add RPLRG for small planet
! ABozzo Apr 2014 fixed bug in leyer height computation

! --------------------------------------------------------------------
! Driver to call stochastic cloud generator and produce subcolumns of
! cloud liquid and ice water contents that can be used by the McICA 
! radiative transfer routine.
! --------------------------------------------------------------------   

USE YOE_MCICA , ONLY : TEMCICA
USE YOMRIP    , ONLY : TRIP
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST    , ONLY : RG, RD
USE YOMDYNCORE, ONLY : RPLRG

IMPLICIT NONE

!EXTERNAL McICA_CLD_GENERATOR

! Note: KLEV    => Number of layers
! Note: KNX_LOC => Number of subcolumns to generate


! PARAMETER


TYPE(TEMCICA)      ,INTENT(INOUT):: YDEMCICA
TYPE(TRIP)         ,INTENT(INOUT):: YDRIP
REAL(KIND=JPRB), PARAMETER :: &
 & M2KM = 1.0_JPRB/1000.0_JPRB, & ! Convert meters to kilometers
 & CUT  = 0.001_JPRB              ! Cutoff for cloud amount    

! KNX_LOC, LDPPH, LDMAXRAN should be passed in somehow not defined as PARAMETERS

!-- all the following are now passed as arguments
!================================================
!INTEGER(KIND=JPIM), PARAMETER :: &
! & KNX_LOC = 150       ! Number of subcolumns to generate

!LOGICAL, PARAMETER :: &
! & LDPPH    = .TRUE., &
! & LDMAXRAN = .TRUE.

!-----------------------------------------------------------

! INPUT DATA

INTEGER(KIND=JPIM), INTENT(IN) :: KNX_LOC
LOGICAL, INTENT(IN) :: LDPPH, LDMAXRAN


REAL(KIND=JPRB), INTENT(IN) :: &
  & PNUAGE(KLON,KLEV),     & ! Column cloud fraction
  & PICEWCIN(KLON,KLEV),   & ! Column in-cloud ice water mixing ratio profile    (kg/kg)
  & PSNOWCIN(KLON,KLEV),   & ! Column in-cloud snow water mixing ratio profile    (kg/kg)
  & PLIQWCIN(KLON,KLEV),   & ! Column in-cloud liquid water mixing ratio profile (kg/kg)
  & PGLAT(KLON),          & ! Latitude
  & PGLON(KLON)             ! Longitude

INTEGER(KIND=JPIM), INTENT(IN) :: & ! Counter and array sizes
  & KLON,                 & ! Length of latitude band   
  & KIDIA,                & ! Starting index for latitude
  & KFDIA,                & ! Ending point for latitude
  & KLEV                    ! Number of layers

REAL(KIND=JPRB), INTENT(IN) :: &
  & PRLC_CF(KLON,KLEV),     & ! Cloud fraction decorrelation length               (km)
  & PRLC_CW(KLON,KLEV),     & ! Cloud condensate decorrelation length             (km)
  & PSIGMA_QCW(KLON,KLEV),  & ! Normalized standard deviation of cloud condensate (unitless)      
  & PT(KLON,KLEV),          & ! Column temperature at layer midpoint              (K)
  & PAP(KLON,KLEV)            ! FULL LEVEL PRESSURE                               (Pa)  

! OUTPUT DATA


REAL(KIND=JPRB), INTENT(OUT) :: & ! Subcolumns
  & PICEWCIN_S(KLON,KLEV,KNX_LOC), & ! Column ice water mixing ratio profile          (g/m^3)
  & PSNOWCIN_S(KLON,KLEV,KNX_LOC), & ! Column snow water mixing ratio profile          (g/m^3)
  & PLIQWCIN_S(KLON,KLEV,KNX_LOC)    ! Column liquid water mixing ratio profile       (g/m^3) 

INTEGER(KIND=JPIM),INTENT(OUT) :: &
  & KNCLDY(KLON)                 ! Number of cloudy subcolumns


! LOCAL DATA


INTEGER(KIND=JPIM) :: &
  & JL,     & ! Counter over GCM columns
  & II,     & ! Counter over NX subcolumns
  & JK        ! Counter over lay vertical layers

REAL(KIND=JPRB) :: &
  & ZRHO      ! Density of air                                 (g/m^3)

REAL(KIND=JPRB) :: &
  & ZPLAY,  & ! Pressure ratio                                 
  & ZROG      ! Dry air gas constant/gravity (RD/RG)

REAL(KIND=JPRB) :: &
  & ZQI_PROF(KLON,KLEV), & ! Temporary array to convert ice water content     (kg/kg) -> (g/m^3)
  & ZQS_PROF(KLON,KLEV), & ! Temporary array to convert snow water content    (kg/kg) -> (g/m^3)
  & ZQC_PROF(KLON,KLEV), & ! Temporary array to convert liquid water content  (kg/kg) -> (g/m^3)
  & ZM(KLON,KLEV)     , &  ! layer depth                                      (km)
  & ZDMULT(KLON,KLEV)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "mcica_cld_generator.intfb.h"
!---------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MCICA_CLD_GEN',0,ZHOOK_HANDLE)

! Zero out fields
DO JL = KIDIA,KFDIA
  KNCLDY(JL) = 0
ENDDO ! JL
      
DO JK = 1 , KLEV 
  DO JL = KIDIA,KFDIA
    ZQC_PROF(JL,JK) = 0.0_JPRB
    ZQI_PROF(JL,JK) = 0.0_JPRB
    ZQS_PROF(JL,JK) = 0.0_JPRB
  ENDDO ! JL
ENDDO ! JK

DO II = 1, KNX_LOC
  DO JK = 1, KLEV
    DO JL = KIDIA,KFDIA
      PICEWCIN_S(JL,JK,II) = 0.0_JPRB
      PSNOWCIN_S(JL,JK,II) = 0.0_JPRB
      PLIQWCIN_S(JL,JK,II) = 0.0_JPRB
    ENDDO
  ENDDO
ENDDO

! Compute the heights of mid-layers

ZROG=RD/(RG/RPLRG)

DO JK = 1, KLEV-1
  DO JL = KIDIA, KFDIA
    ZPLAY = PAP(JL,JK+1)/PAP(JL,JK)
    ZM(JL,JK) = ZROG*0.5*(PT(JL,JK)+PT(JL,JK+1))*LOG(ZPLAY)*M2KM 
  ENDDO
ENDDO

! Convert the cloud condensate from kg/kg to g/m^3 (is the cloud condensate cloud or domain mean)?
DO JK = 1, KLEV
  DO JL = KIDIA, KFDIA
! Compute layer height               
    IF (PNUAGE(JL,JK) > CUT) THEN
! ZRHO in kg/m^3                  
      ZRHO            = 1000.0_JPRB * PAP(JL,JK)/(RD*PT(JL,JK))
      ZDMULT(JL,JK)   = ZRHO !/PNUAGE(JL,JK) ! Divide by cloud amount if not provided in-cloud water contents
      ZQI_PROF(JL,JK) = PICEWCIN(JL,JK)*ZDMULT(JL,JK)
      ZQS_PROF(JL,JK) = PSNOWCIN(JL,JK)*ZDMULT(JL,JK)
      ZQC_PROF(JL,JK) = PLIQWCIN(JL,JK)*ZDMULT(JL,JK)
    ENDIF
  ENDDO  ! JL
ENDDO ! JK
         
! Call cloud generator

CALL MCICA_CLD_GENERATOR &
  &( YDEMCICA,YDRIP, ZM, PNUAGE, ZQC_PROF, ZQI_PROF, ZQS_PROF, PGLAT, PGLON, &
  &  PRLC_CF, PRLC_CW, PSIGMA_QCW,   &
  &  KLON, KIDIA, KFDIA, KLEV, KNX_LOC,  &
  &  LDPPH, LDMAXRAN, &
  &  PICEWCIN_S, PSNOWCIN_S, PLIQWCIN_S, KNCLDY)

! Note that all of the cloudy subcolumns in the arrays PICEWCIN_S, PSNOWCIN_S 
! and PLIQWCIN_S should be at the "front" of the arrays and the array 
! KNCLDY contains the number of cloud subcolumns


!-- put back the cloud water in same unit as inputs!
DO JK= 1, KLEV
  DO JL= KIDIA, KFDIA
    IF (PNUAGE(JL,JK) > CUT) THEN
      DO II= 1, KNX_LOC
        PICEWCIN_S(JL,JK,II)=PICEWCIN_S(JL,JK,II)/ZDMULT(JL,JK)
        PSNOWCIN_S(JL,JK,II)=PSNOWCIN_S(JL,JK,II)/ZDMULT(JL,JK)
        PLIQWCIN_S(JL,JK,II)=PLIQWCIN_S(JL,JK,II)/ZDMULT(JL,JK)
      ENDDO
    ENDIF
  ENDDO
ENDDO 

!---------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MCICA_CLD_GEN',1,ZHOOK_HANDLE)
END SUBROUTINE MCICA_CLD_GEN

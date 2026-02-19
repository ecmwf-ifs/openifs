! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE MCICA_CLD_GENERATOR &
  &( YDEMCICA,YDRIP,PZF    , PAVG_CF, PAVG_QLW, PAVG_QIW, PAVG_QSW, PGLAT, PGLON, &
  &  PRLC_CF, PRLC_CW, PSIGMA_QCW, &
  &  KLON   , KIDIA   , KFDIA      , KLEV    , KNX_LOC , &
  &  LDPPH     , LDMAXRAN       , &
  &  PICEWCIN_S, PSNOWCIN_S, PLIQWCIN_S     , KNCLDY)

! --------------------------------------------------------------------
! --------------------------------------------------------------------
! Generate a NX column by KLEV layer subcolumns.
! Method and model evaluation is described in:
! Raisanen et al, 2004, Stochastic generation of subgrid-scale cloudy 
! columns for large-scale models,
! Quart. J. Roy. Meteorol. Soc., 2004 , 130 , 2047-2068
! Input profiles must be from top of the model to the bottom as is the output.

!  Original:
!  ---------
!   ???

!  Modifications:
!  --------------
!      K. Yessad (July 2014): Move some variables.
!      R. Hogan  (Dec 2014):  Changed the ISEED definition to avoid 
!                             diagonal stripes in instantaneous fluxes
! --------------------------------------------------------------------   

USE PARKIND1, ONLY : JPIM, JPRB, JPRD
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOE_MCICA,ONLY : TEMCICA

!-- Mike Fisher's random number generator
USE YOMCT3  , ONLY : NSTEP
USE YOMRIP0 , ONLY : NINDAT
USE YOMRIP  , ONLY : TRIP

USE RANDOM_NUMBERS_MIX  , ONLY : RANDOMNUMBERSTREAM, UNIFORM_DISTRIBUTION

! --------------------------------------------------------------------   

IMPLICIT NONE

! This subroutine needs the array xcw, and dimensions of the array NMcI1 and NMcI2.
! These are provided within YOE_McICA.

! Note: KLEV    => Number of layers in GCM
!       KNX_LOC => Number of columns to generate


! PARAMETERS


TYPE(TEMCICA)      ,INTENT(INOUT):: YDEMCICA
TYPE(TRIP)         ,INTENT(INOUT):: YDRIP
REAL(KIND=JPRB), PARAMETER ::&
  & CUT = 0.001 ! Cutoff for cloud amount

! KNX_LOC, LDPPH, LDMAXRAN should be passed in somehow not defined as PARAMETERS

!INTEGER(KIND=JPIM), PARAMETER :: &
!  & KNX_LOC = 150       ! Number of subcolumns to generate

!LOGICAL, PARAMETER :: &
!  & LDPPH = .TRUE.,     & ! Switch for assuming cloud to be PPH (.TRUE.) or not (.FALSE.)
!  & LDMAXRAN = .TRUE.     ! Switch for assuming cloud to be maximum-random overlap (.TRUE.) or not (.FALSE.)

!----------------------------------------------------------------------------------


! INPUT DATA

INTEGER(KIND=JPIM), INTENT(IN) :: KNX_LOC
LOGICAL, INTENT(IN) :: LDPPH, LDMAXRAN

INTEGER(KIND=JPIM), INTENT(IN) ::&! Counters and array sizes
  & KLON,&! Length of latitude band
  & KIDIA,&! Starting index for latitude
  & KFDIA,&! Ending point for latitude
  & KLEV                     ! Number of layers

REAL(KIND=JPRB), INTENT(IN) ::&
  & PZF(KLON,KLEV),&! Layer depth                           (km)
  & PAVG_CF(KLON,KLEV),&! Cloud fraction for each layer         (unitless)
  & PAVG_QLW(KLON,KLEV),&! Cloud mean liquid water mixing ratio  (g/m^3)
  & PAVG_QIW(KLON,KLEV),&! Cloud mean ice mixing ratio           (g/m^3)
  & PAVG_QSW(KLON,KLEV),&! Cloud mean snow mixing ratio          (g/m^3)
  & PSIGMA_QCW(KLON,KLEV),&! Normalized cloud condensate std. dev. (unitless) 
  & PRLC_CF(KLON,KLEV),&! Cloud fraction decorrelation length   (km)
  & PRLC_CW(KLON,KLEV),&! Cloud condensate decorrelation length (km)
  & PGLAT(KLON),&! Latitude in degrees
  & PGLON(KLON)              ! Longitude in degrees

! OUTPUT DATA


REAL(KIND=JPRB),INTENT(OUT) ::&
  & PICEWCIN_S(KLON,KLEV,KNX_LOC),&! Column ice water mixing ratio profile    (g/m^3)
  & PSNOWCIN_S(KLON,KLEV,KNX_LOC),&! Column snow water mixing ratio profile    (g/m^3)
  & PLIQWCIN_S(KLON,KLEV,KNX_LOC)    ! Column liquid water mixing ratio profile (g/m^3)

INTEGER(KIND=JPIM), INTENT(OUT) ::&
  & KNCLDY(KLON)                   ! Number of cloudy subcolumns


! LOCAL DATA


INTEGER(KIND=JPIM) ::&
  & IL,&! Counter over KLON GCM columns
  & I,&! Counter
  & JK,&! Counter over KLEV vertical layers
  & I_TOP(KLON),&! Index of top most cloud layer
  & I_BASE(KLON),&! Index of lowest cloud layer
  & IND1,&! Index in variability calculation 
  & IND2,&! Index in variability calculation
  & I_LOC(KLON)   ! Counter to place the new subcolumns into the arrays starting from the front

REAL(KIND=JPRB) ::&
  & ZALPHA(KLON,KLEV),&! Fraction of maximum/random cloud overlap
  & ZRCORR(KLON,KLEV),&! Fraction of maximum/random cloud condensate overlap
  & ZRIND1,&! Real index in variability calculation
  & ZRIND2,&! Real index in variability calculation
  & ZCW                  ! Ratio of cloud condensate mixing ratio for this cell to its layer cloud-mean value

REAL(KIND=JPRB) :: & ! Random number vectors
  & ZX(KLON),  &
  & ZY(KLON),  &
  & ZX1(KLON), &
  & ZY1(KLON), &
  & ZX2(KLON), &
  & ZY2(KLON)  

REAL(KIND=JPRB) :: & ! Random number vectors
  & ZZX (KLEV,KLON), &
  & ZZY (KLEV,KLON), &
  & ZZX1(KLEV,KLON), &
  & ZZY1(KLEV,KLON), &
  & ZZX2(KLEV,KLON), &
  & ZZY2(KLEV,KLON)

LOGICAL :: LL_CLD(KLON)  ! Flag if cloudy subcolumn was generated

TYPE(RANDOMNUMBERSTREAM) :: YL_RANDOM_STREAM(KLON)
INTEGER(KIND=JPIM) :: ISEED, ITIM, IDAY

REAL(KIND=JPRB) :: ZTMP1, ZTMP2
INTEGER(KIND=JPIM) :: ISTOP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! --------------------------------------------------------------------   

#include "fcttim.func.h"
#include "setran.intfb.h"

! --------------------------------------------------------------------   
IF (LHOOK) CALL DR_HOOK('MCICA_CLD_GENERATOR',0,ZHOOK_HANDLE)
ASSOCIATE(NMCI1=>YDEMCICA%NMCI1, NMCI2=>YDEMCICA%NMCI2, XCW=>YDEMCICA%XCW, &
 & TSTEP=>YDRIP%TSTEP)
! --------------------------------------------------------------------   

! Create the relevant seed from date and time
! get the starting day and number of minutes since start
IDAY=NDD(NINDAT)
ITIM=NINT(NSTEP*TSTEP/60._JPRB)


! Initialize the arrays

DO IL = KIDIA, KFDIA
  I_LOC(IL) = 1
  LL_CLD(IL) = .FALSE.
ENDDO

! Find uppermost cloudy layer

DO IL = KIDIA, KFDIA
  DO JK=1,KLEV
    I_TOP(IL) = JK
    IF (PAVG_CF(IL,JK) > CUT) EXIT
  ENDDO ! JK
ENDDO ! IL


! Find lowermost cloudy layer

DO IL = KIDIA, KFDIA
  DO JK=KLEV,1,-1
    I_BASE(IL) = JK
    IF (PAVG_CF(IL,JK) > CUT) THEN
      EXIT
    ENDIF
  ENDDO ! JK
ENDDO ! IL

! Calculate overlap factors ZALPHA for cloud fraction and ZRCORR for cloud
! condensate based on layer midpoint distances and decorrelation depths

ISTOP=0

ZALPHA=0.0_JPRB
ZRCORR=0.0_JPRB

DO IL = KIDIA, KFDIA
  DO JK=I_TOP(IL),I_BASE(IL)-1
    IF (PRLC_CF(IL,JK) > 0.0) THEN
      ZTMP1 = PZF(IL,JK) / PRLC_CF(IL,JK)
      ZALPHA(IL,JK) = EXP(-ZTMP1)
    ELSE
      ZALPHA(IL,JK) = 0.0
    ENDIF
    IF (PRLC_CW(IL,JK) > 0.0) THEN
      ZTMP2 = PZF(IL,JK) / PRLC_CW(IL,JK)
      ZRCORR(IL,JK) = EXP(-ZTMP2)
    ELSE
      ZRCORR(IL,JK) = 0.0
    ENDIF
!    IF (ISTOP > 0) THEN
!      PRINT 9002,ISTOP,IL,JK,ZTMP1,ZTMP2,PZF(IL,JK),PZF(IL,JK+1),PRLC_CF(IL,JK),PRLC_CW(IL,JK),ZALPHA(IL,JK),ZRCORR(IL,K)
9002  FORMAT(1X,'CldGen:',3I4,10E12.5)
!    ENDIF
  ENDDO ! JK
ENDDO ! IL

DO IL = KIDIA, KFDIA
  ! Original method produced diagonal stripes across the globe
  !    ISEED=NINT( PGLON(IL)+PGLAT(IL)+ITIM+IDAY )
  ! New method gives unique value for roughly every 1-km square on the globe
  ISEED= NINT(PGLON(IL)*9000.0_JPRD + PGLAT(IL)*100.0_JPRD) +ITIM+IDAY
  CALL SETRAN ( ISEED, YL_RANDOM_STREAM(IL) )
ENDDO

DO I=1,KNX_LOC

  DO IL = KIDIA, KFDIA
    CALL UNIFORM_DISTRIBUTION (ZZX(1:KLEV,IL) , YL_RANDOM_STREAM(IL) )
    CALL UNIFORM_DISTRIBUTION (ZZY(1:KLEV,IL) , YL_RANDOM_STREAM(IL) )
    CALL UNIFORM_DISTRIBUTION (ZZX1(1:KLEV,IL), YL_RANDOM_STREAM(IL) )
    CALL UNIFORM_DISTRIBUTION (ZZY1(1:KLEV,IL), YL_RANDOM_STREAM(IL) )
    CALL UNIFORM_DISTRIBUTION (ZZX2(1:KLEV,IL), YL_RANDOM_STREAM(IL) )
    CALL UNIFORM_DISTRIBUTION (ZZY2(1:KLEV,IL), YL_RANDOM_STREAM(IL) )
  ENDDO

!! Need to check if a cloudy subcolumn was generated
!  DO IL = KIDIA, KFDIA
!    IF (LL_CLD(IL)) THEN
!      I_LOC(IL) = I_LOC(IL) + 1
!      LL_CLD(IL) = .FALSE.
!    ENDIF
!  ENDDO

! Generate all subcolumns for latitude chain

DO IL=KIDIA,KFDIA
  DO JK = I_TOP(IL), I_BASE(IL)

! Generate all of the needed random numbers

    IF (JK == I_TOP(IL)) THEN
!      CALL RANDOM_NUMBER(ZX)
!      CALL RANDOM_NUMBER(ZY)
!     CALL UNIFORM_DISTRIBUTION ( ZX, YL_RANDOM_STREAM )
!     CALL UNIFORM_DISTRIBUTION ( ZY, YL_RANDOM_STREAM )
      ZX(IL)=ZZX(JK,IL)
      ZY(IL)=ZZY(JK,IL)
    ENDIF
    
!    CALL RANDOM_NUMBER(ZX1)
!    CALL RANDOM_NUMBER(ZY1)
!   CALL UNIFORM_DISTRIBUTION ( ZX1, YL_RANDOM_STREAM )
!   CALL UNIFORM_DISTRIBUTION ( ZY1, YL_RANDOM_STREAM )
    ZX1(IL)=ZZX1(JK,IL)
    ZY1(IL)=ZZY1(JK,IL)

!    CALL RANDOM_NUMBER(ZX2)
!    CALL RANDOM_NUMBER(ZY2)
!   CALL UNIFORM_DISTRIBUTION ( ZX2, YL_RANDOM_STREAM )
!   CALL UNIFORM_DISTRIBUTION ( ZY2, YL_RANDOM_STREAM )
    ZX2(IL)=ZZX2(JK,IL)
    ZY2(IL)=ZZY2(JK,IL)
            

! Maximum-random overlap
      IF (LDMAXRAN) THEN
        IF (ZX(IL) < 1.0_JPRB - PAVG_CF(IL,JK-1)) THEN !It is clear above
          ZX(IL) = ZX1(IL) * (1.0_JPRB - PAVG_CF(IL,JK-1))
          ZY(IL) = ZY1(IL)
        ENDIF
! Generalized overlap
      ELSE
        IF (ZX1(IL) > ZALPHA(IL,JK-1)) ZX(IL) = ZX2(IL)
        IF (ZY1(IL) > ZRCORR(IL,JK-1)) ZY(IL) = ZY2(IL)
      ENDIF

! Treatment of cloudy cells
      IF (ZX(IL) > 1.0_JPRB - PAVG_CF(IL,JK)) THEN ! Generate cloud in this layer

        IF (LDPPH) THEN ! Homogeneous clouds
          ZCW = 1.0_JPRB
        ELSE
! Horizontally variable clouds:
! Determine ZCW = ratio of cloud condensate miximg ratio QC for this cell to
! its mean value for all cloudy cells in this layer.
! Use bilinear interpolation of ZCW tabulated in array XCW as a function of
!    * cumulative probability Y
!    * relative standard deviation PSIGMA
! Take care that the definition of ZRIND2 is consistent with subroutine
! TABULATE_XCW

          ZRIND1 = ZY(IL) * (NMCI1 - 1) + 1.0_JPRB
          IND1  = MAX(1, MIN(INT(ZRIND1), NMCI1-1))
          ZRIND1 = ZRIND1 - IND1
          ZRIND2 = 40.0_JPRB * PSIGMA_QCW(IL,JK) - 3.0_JPRB
          IND2  = MAX(1, MIN(INT(ZRIND2), NMCI2-1))
          ZRIND2 = ZRIND2 - IND2

          ZCW = (1.0-ZRIND1) * (1.0_JPRB - ZRIND2) * XCW(IND1,IND2)&
            &               + (1.0_JPRB - ZRIND1) * ZRIND2       * XCW(IND1,IND2+1)&
            &               + ZRIND1 * (1.0_JPRB - ZRIND2)       * XCW(IND1+1,IND2)&
            &               + ZRIND1 * ZRIND2             * XCW(IND1+1,IND2+1) 
        ENDIF

! A horizontally constant IWC/LWC ratio is assumed for each layer so far
! -- BUG all cloudy bits on same spectral side
!        PLIQWCIN_S(IL,JK,I_LOC(IL)) = ZCW * PAVG_QLW(IL,JK)
!        PICEWCIN_S(IL,JK,I_LOC(IL)) = ZCW * PAVG_QIW(IL,JK)
        PLIQWCIN_S(IL,JK, I ) = ZCW * PAVG_QLW(IL,JK)
        PICEWCIN_S(IL,JK, I ) = ZCW * PAVG_QIW(IL,JK)
        PSNOWCIN_S(IL,JK, I ) = ZCW * PAVG_QSW(IL,JK)

! Note that a cloud subcolumn was generated
        LL_CLD(IL)             = .TRUE.
      ENDIF
               
    ENDDO                ! JK
  ENDDO              ! IL

!-- JasonC 20060110 (moved from start to end of loop over I)
! Need to check if a cloudy subcolumn was generated
  DO IL = KIDIA, KFDIA
    IF (LL_CLD(IL)) THEN
      I_LOC(IL) = I_LOC(IL) + 1
      LL_CLD(IL) = .FALSE.
    ENDIF
  ENDDO

ENDDO                  ! I

! Record the number of cloudy subcolumns generated
! The variable KNCLDY can be used later to randomly select 
! only cloudy columns
DO IL = KIDIA, KFDIA
  KNCLDY(IL) = I_LOC(IL)-1
ENDDO ! IL

!print 9701,KNX_LOC,(KNCLDY(IL),IL=KIDIA,KFDIA)
9701 FORMAT(1X,'mcgen1 ',I6,1X,50I4)
!IL=KIDIA
!do JK=I_TOP(IL),I_BASE(IL)
!  print 9702,(PLIQWCIN_S(IL,JK,I),I=1,KNX_LOC,8)
!  print 9703,(PICEWCIN_S(IL,JK,I),I=1,KNX_LOC,8)
9702 FORMAT(1X,'mcgen2 ',20E12.4)
9703 FORMAT(1X,'mcgen3 ',20E12.4)
!enddo


!---------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MCICA_CLD_GENERATOR',1,ZHOOK_HANDLE)
END SUBROUTINE MCICA_CLD_GENERATOR

! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_SO2SO4 &
  &( YDRIP,YDEAERSNK,KIDIA, KFDIA , KLON  , KLEV  , &
  &  PSO2 , PITSO2, PGELAT, PGELAM, PTSPHY, PRHCL, PT, &
  &  PTSO2, PTSO4, PFSO2, PFSO4,PDP)

!*** * AER_SO2SO4* - GAS-TO-PARTICLE (SULPHATE AEROSOLS)

!**   INTERFACE.
!     ----------
!          *AER_SO2SO4* IS CALLED FROM *AER_PHY3*.

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        FROM O.BOUCHER and N.HUNNEUS's gastoparticle 

!     SOURCE.
!     -------

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 20080508
!        SR 03/2016 add in-cloud conv rate increase
!        SR 03/2016 add diurnal cycle and dependence on temperature
!        SR 05/2017 externalize diurnal cycle
!-----------------------------------------------------------------------

USE YOMRIP   , ONLY : TRIP
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RDAY , RMSO2, RMSO4, RG
USE YOEAERSNK, ONLY : TEAERSNK

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERSNK)    ,INTENT(INOUT) :: YDEAERSNK
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSO2(KLON,KLEV), PITSO2(KLON,KLEV), PGELAT(KLON), PGELAM(KLON), PDP(KLON, KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHCL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSO2(KLON,KLEV), PTSO4(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSO2(KLON), PFSO4(KLON)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL, JK
REAL(KIND=JPRB) :: ZSO4SO2       , ZFACT(KLON), ZINCR,ZFACT_T
REAL(KIND=JPRB) :: ZAERCONV(KLON), ZTAUCHEM0(KLON),ZSCALE(KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------

#include "compo_diurnal.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_SO2SO4',0,ZHOOK_HANDLE)
ASSOCIATE(RSO2CV1=>YDEAERSNK%RSO2CV1, RSO2CV2=>YDEAERSNK%RSO2CV2)
ZSO4SO2=RMSO4/RMSO2
! In cloud conversion rate is increased by a factor 2
ZINCR=2.0_JPRB

!-- Huneeus et al., 2009, $3.1, p.215, eq.2

CALL COMPO_DIURNAL(YDRIP,KIDIA,KFDIA,KLON,'Sine',PGELAM,PGELAT,ZSCALE,PAMPLITUDE=1.0_JPRB,PHOURPEAK=12.0_JPRB)


DO JL=KIDIA,KFDIA
  ZTAUCHEM0(JL) = RDAY * (RSO2CV1 - RSO2CV2 * COS(PGELAT(JL)) )
! naj clipping and squaring of ZSCALE may need review 
  ZSCALE(JL)=MAX(0.5_JPRB, ZSCALE(JL))
  ZFACT(JL)=(1._JPRB-EXP(-PTSPHY/ZTAUCHEM0(JL))) / PTSPHY
  ZFACT(JL)=ZFACT(JL)*ZSCALE(JL)**2
ENDDO

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ! dependency on temperature following Eatough et al. (1994)
    ZFACT_T=EXP(32.37_JPRB - 9000_JPRB/PT(JL,JK))
    ZFACT_T=MAX(0.3_JPRB, ZFACT_T)
    ZFACT_T=MIN(1.5_JPRB, ZFACT_T)
    ZAERCONV(JL) =  (PSO2(JL,JK) + PTSPHY * PITSO2(JL,JK)) *ZFACT(JL)*ZFACT_T
    IF (PRHCL(JL,JK) > 0.98_JPRB ) THEN
      PTSO4(JL,JK) =  ZAERCONV(JL) * ZSO4SO2 * ZINCR
      PTSO2(JL,JK) = -ZAERCONV(JL) * ZINCR
    ELSE
      PTSO4(JL,JK) =  ZAERCONV(JL) * ZSO4SO2
      PTSO2(JL,JK) = -ZAERCONV(JL)    
    ENDIF
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  PFSO4(JL)=0.0_JPRB
  PFSO2(JL)=0.0_JPRB
ENDDO

DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
      PFSO4(JL) = PFSO4(JL) + PTSO4(JL,JK)*(PDP(JL,JK))/RG
      PFSO2(JL) = PFSO2(JL) + PTSO2(JL,JK)*(PDP(JL,JK))/RG
   ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_SO2SO4',1,ZHOOK_HANDLE)
END SUBROUTINE AER_SO2SO4

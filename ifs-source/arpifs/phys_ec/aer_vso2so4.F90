! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_VSO2SO4 &
  &( YDEAERSNK,KIDIA, KFDIA , KLON  , KLEV  , &
  &  PSO2 , PITSO2, PGELAT, PTSPHY, &
  &  PTSO2, PTSO4 )

!*** * AER_SO2SO4* - GAS-TO-PARTICLE (VOLCANIC SULPHATE AEROSOLS)

!**   INTERFACE.
!     ----------
!          *AER_VSO2SO4* IS CALLED FROM *AER_PHY3*.

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        FROM O.BOUCHER and N.HUNNEUS's gastoparticle 

!     SOURCE.
!     -------

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 20130729  COPY OF *AER_SO2SO4*

!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST    ,ONLY : RDAY , RMSO2, RMSO4
USE YOEAERSNK ,ONLY : TEAERSNK

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERSNK)    ,INTENT(INOUT) :: YDEAERSNK
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSO2(KLON,KLEV), PITSO2(KLON,KLEV), PGELAT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSO2(KLON,KLEV), PTSO4(KLON,KLEV)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL, JK
REAL(KIND=JPRB) :: ZSO4SO2       , ZFACT(KLON)
REAL(KIND=JPRB) :: ZAERCONV(KLON), ZTAUCHEM0(KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_VSO2SO4',0,ZHOOK_HANDLE)
ASSOCIATE(RVSO2CV1=>YDEAERSNK%RVSO2CV1, RVSO2CV2=>YDEAERSNK%RVSO2CV2)
ZSO4SO2=RMSO4/RMSO2

!-- Huneeus et al., 2009, $3.1, p.215, eq.2

DO JL=KIDIA,KFDIA
  ZTAUCHEM0(JL) = RDAY * (RVSO2CV1 - RVSO2CV2 * COS(PGELAT(JL)) )
  ZFACT(JL)=(1._JPRB-EXP(-PTSPHY/ZTAUCHEM0(JL))) / PTSPHY
ENDDO

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZAERCONV(JL) =  (PSO2(JL,JK) + PTSPHY * PITSO2(JL,JK)) * ZFACT(JL)

    PTSO4(JL,JK) =  ZAERCONV(JL) * ZSO4SO2
    PTSO2(JL,JK) = -ZAERCONV(JL)    
  ENDDO
ENDDO
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_VSO2SO4',1,ZHOOK_HANDLE)
END SUBROUTINE AER_VSO2SO4

! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_CGROWTH &
  &( YDEAERATM,KIDIA   , KFDIA    , KLON     , KLEV  , &
  &  PAERPHO , PITAERPHO, PTSPHY  ,&
  &  PTAERPHI, PTAERPHO & 
  &)

!*** * AER_CGROWTH* - GROWTH OF O.M. AND B.C. AEROSOLS

!**   INTERFACE.
!     ----------
!          *AER_CGROWTH* IS CALLED FROM *AER_PHY3*.

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        FROM O.BOUCHER's cgrowth 

!     SOURCE.
!     -------

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2006-11-21
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEAERATM ,ONLY : TEAERATM

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERATM)    ,INTENT(INOUT) :: YDEAERATM
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV

REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERPHO(KLON,KLEV) , PITAERPHO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAERPHI(KLON,KLEV), PTAERPHO(KLON,KLEV)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL, JK
REAL(KIND=JPRB) :: ZCOEF, ZAERCONV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_CGROWTH',0,ZHOOK_HANDLE)
ASSOCIATE(RGRATE=>YDEAERATM%RGRATE)
ZCOEF= RGRATE

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZAERCONV = (PAERPHO(JL,JK) + PTSPHY * PITAERPHO(JL,JK)) * ZCOEF
    PTAERPHI(JL,JK) =  ZAERCONV
    PTAERPHO(JL,JK) = -ZAERCONV    
  ENDDO
ENDDO
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_CGROWTH',1,ZHOOK_HANDLE)
END SUBROUTINE AER_CGROWTH

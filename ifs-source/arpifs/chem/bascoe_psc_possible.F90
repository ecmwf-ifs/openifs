! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_PSC_POSSIBLE( KIDIA, KFDIA, KLON, KLEV, KJULDAY, PTP,PRSF1,PLAT,KTROPOP, &
    &  LD_PSC_POSSIBLE,KTOP_PSC,KBOT_PSC)

!**   DESCRIPTION 
!   Set llpsc_possible, jtop_psc and jbot_psc (and initialize psc_diags)
!   i.e. determine when and where PSC can appear:
!    - levels below 10 hPa and above tropopause
!    - S.H. from May 1 to Oct 31, N.H. from Dec 1 to April 30
!    - Temperature < 200 K !!
!-----------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
INTEGER(KIND=JPIM),INTENT(IN)      :: KIDIA, KFDIA, KLON, KLEV
REAL(KIND=JPRB),INTENT(IN)         :: PTP(KLON,KLEV),PRSF1(KLON,KLEV),PLAT(KLON)
INTEGER(KIND=JPIM), INTENT(IN)     :: KJULDAY
INTEGER(KIND=JPIM), INTENT(IN)     :: KTROPOP(KLON)
LOGICAL, INTENT(OUT)               :: LD_PSC_POSSIBLE(KLON)
INTEGER(KIND=JPIM)   , INTENT(OUT) :: KTOP_PSC(KLON),KBOT_PSC(KLON)

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------

!*       0.5   LOCAL VARIABLES
!              ---------------
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JL,JK,J1,J2
INTEGER(KIND=JPIM), DIMENSION(KLON)      ::  JTOP

IF (LHOOK) CALL DR_HOOK('BASCOE_PSC_POSSIBLE',0,ZHOOK_HANDLE )


LD_PSC_POSSIBLE(KIDIA:KFDIA) = .FALSE.
KTOP_PSC(KIDIA:KFDIA) = -998
KBOT_PSC(KIDIA:KFDIA) = -999

JTOP(KIDIA:KFDIA) = 1
! Find 10 hPa level

DO JL=KIDIA,KFDIA
  DO JK=1,KLEV/2
    IF (PRSF1(JL,JK) <= 1000_JPRB)  JTOP(JL) =  JK
  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  ! Use tropopause level
  J1 = JTOP(JL)
  J2 = KTROPOP(JL) - 2
  J2 = MAX(J1+1,J2)
  IF( ANY( PTP(JL,J1:J2) < 220. ) .AND. &
&       ( PLAT(JL)<-50. .and. (KJULDAY>120 .and. KJULDAY<305) ) &
&  .OR. ( PLAT(JL)> 50. .and. (KJULDAY>334 .or.  KJULDAY<121) ) &
&           ) THEN
    LD_PSC_POSSIBLE(JL) = .TRUE.
    KTOP_PSC(JL) = J1
    KBOT_PSC(JL) = J2
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('BASCOE_PSC_POSSIBLE',1,ZHOOK_HANDLE )
END SUBROUTINE BASCOE_PSC_POSSIBLE

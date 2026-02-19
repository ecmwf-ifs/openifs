! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE METHOX(KIDIA,  KFDIA,  KLON,  KLEV,PQ,     PTENQ,  PAP )

!**** *METHOX*   - Calculate humidity tendencies from methane
!                  oxidation and photolysis

!**   INTERFACE.
!     ----------
!        CALL *METHOX* FROM *CALLPAR*
!              ------        -------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS

!     INPUT PARAMETERS (REAL):

!    *PAP*          PRESSURE                                      PA
!    *PQ*           SPECIFIC HUMIDITY                             KG/KG

!     UPDATED PARAMETERS (REAL):

!    *PTENQ*        TENDENCY OF SPECIFIC HUMIDITY                 KG/(KG*S)

!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        MODULE YOEMETH
!        MODULE YOMCST

!     METHOD.
!     -------
!        SEE RD-MEMO R60.1/AJS/31

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        SEE RD-MEMO R60.1/AJS/31

!     AUTHOR.
!     -------
!        C.JAKOB   *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 98-04-07
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P Bechtold 18/05/2012   Use RDAYI for 86400
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOEMETH   , ONLY : RALPHA1 ,RALPHA2  ,RQLIM   ,&
 & RPBOTOX,  RPBOTPH ,RPTOPOX  ,RPTOPPH ,&
 & RALPHA3,  RLOGPPH  
USE YOMCST   , ONLY :  RPI, RDAYI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
LOGICAL :: LLOXID,         LLPHOTO

INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPRB) :: ZARG, ZPRATIO, ZTAU1, ZTAU2, ZTDAYS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('METHOX',0,ZHOOK_HANDLE)
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

    PTENQ(JL,JK)= 0.0_JPRB

    LLOXID=PAP(JL,JK) < RPBOTOX.AND.PQ(JL,JK) < RQLIM
    LLPHOTO=PAP(JL,JK) < RPBOTPH

!     METHANE OXIDATION

    IF(LLOXID) THEN
      IF(PAP(JL,JK) <= RPTOPOX) THEN
        ZTDAYS=100._JPRB
      ELSE
        ZPRATIO=(LOG(PAP(JL,JK)/RPTOPOX))**4._JPRB/LOG(RPBOTOX/PAP(JL,JK))
        ZTDAYS=100._JPRB*(1+RALPHA1*ZPRATIO)
      ENDIF
      ZTAU1=RDAYI*ZTDAYS
      PTENQ(JL,JK)=PTENQ(JL,JK)+(RQLIM-PQ(JL,JK))/ZTAU1
    ENDIF

!     PHOTOLYSIS

    IF(LLPHOTO) THEN
      IF(PAP(JL,JK) <= RPTOPPH) THEN
        ZTDAYS=3._JPRB
      ELSE
        ZARG=RALPHA2-RALPHA3*(1+COS((RPI*LOG(PAP(JL,JK)/RPBOTPH))/RLOGPPH))
        ZTDAYS=1.0_JPRB/(EXP(ZARG)-0.01_JPRB)
      ENDIF
      ZTAU2=RDAYI*ZTDAYS
      PTENQ(JL,JK)=PTENQ(JL,JK)-PQ(JL,JK)/ZTAU2
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('METHOX',1,ZHOOK_HANDLE)
END SUBROUTINE METHOX
